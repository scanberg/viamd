#include <viamd_event.h>
#include <viamd.h>
#include <representation.h>
#include <color_utils.h>
#include <event.h>
#include <serialization_utils.h>
#include <display_units.h>
#include <plot_series.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_log.h>
#include <md_system.h>
#include <md_util.h>
#include <md_nonbonded.h>

#include <imgui.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>

#include <string>
#include <algorithm>

namespace dataset {

struct AtomElementMapping {
    char lbl[31] = "";
    md_element_t elem = 0;
};

// What an atom type looked like straight out of the loader, before the user got to it. Loading the same
// file again reproduces exactly this, so it is the baseline a workspace stores modifications against:
// whatever differs from it is something the user changed, and nothing else needs to be written.
struct AtomTypeLoadState {
    md_atomic_number_t z = 0;
    float      radius = 0;
    float      mass   = 0;
    uint32_t   color  = 0;
    md_atom_type_flags_t flags = MD_ATOM_TYPE_FLAG_NONE;
};

// The user may turn an atom type into a coarse grained bead and back (an atom). A virtual site stays one.
static inline bool type_is_bead(md_atom_type_flags_t flags) {
    return md_atom_type_flags_particle_kind(flags) == MD_PARTICLE_BEAD;
}

static inline md_atom_type_flags_t type_with_bead(md_atom_type_flags_t flags, bool bead) {
    if (bead) return md_atom_type_flags_set_particle_kind(flags, MD_PARTICLE_BEAD);
    return type_is_bead(flags) ? md_atom_type_flags_set_particle_kind(flags, MD_PARTICLE_ATOM) : flags;
}

// What a particle is, is topology: whatever was derived from it is stale when it changes
static inline void set_type_bead(md_system_t& sys, size_t i, bool bead) {
    const md_atom_type_flags_t flags = type_with_bead(sys.atom.type.flags[i], bead);
    if (flags != sys.atom.type.flags[i]) {
        sys.atom.type.flags[i] = flags;
        md_system_topology_changed(&sys);
    }
}

// A single deserialized [AtomType] section. The workspace is parsed before the system is loaded, so these
// are buffered here and applied once the atom types actually exist (see apply_pending_atom_type_overrides).
struct AtomTypeOverride {
    char name[32] = "";             // Name of the atom type
    char ff_type[64] = "";          // Force field type, empty when the type has none
    md_atomic_number_t load_z = 0;  // Element the type was assigned upon load. Name, force field type and element identify the type

    // Only the fields which were present in the section are applied
    bool has_z              = false;
    bool has_radius         = false;
    bool has_mass           = false;
    bool has_color          = false;
    bool has_coarse_grained = false;
    bool has_use_defaults   = false;

    md_atomic_number_t z = 0;
    float    radius = 0;
    float    mass   = 0;
    uint32_t color  = 0;
    bool     coarse_grained = false;
    bool     use_defaults   = true;
};

// We use this to represent a single entity within the loaded system, e.g. a residue type
// This is used to represent multiple types, so all fields are not used in all cases
struct DatasetItem {
    char label[32] = "";
    uint32_t count = 0;
    float fraction = 0;
    
    // Extended metadata for popups
    uint64_t key = 0;            // Unique key of the type
    md_array(int) indices = 0;   // Indices into the corresponding structures which are represented by this item: i.e. chain or residue indices (for highlighting)

    // Atom type only
    bool use_defaults = true; // Flag if this particular atom type should be linked to the default values stemming for the element (only applicable to atom types with an element, i.e. not coarse grained)

    AtomTypeLoadState load = {}; // What the loader assigned, the baseline user modifications are stored against
};

struct ElementDefault {
    vec4_t color;
    float radius;
    float mass;
};

// Which properties of an element default the user has changed away from the values built into mdlib.
// Only these are stored in a workspace, and only these are pushed onto the atom types linked to the element:
// a property the user never touched must not overwrite whatever the loader supplied for a type.
struct ElementDefaultDelta {
    bool color  = false;
    bool radius = false;
    bool mass   = false;

    explicit operator bool() const { return color || radius || mass; }
};

// ## Component types
//
// Two tiers, as for the atom types: the component types of the system (sys.component.type) are linked to the defaults of
// their standard monomer (md_monomer_t, as an atom type is to its element) and follow them, or carry a color of their own.
// A type which is no form of a standard monomer (a ligand, water, an ion) always carries its own.

// The defaults of a standard monomer, mdlib's built in ones until the user edits them
struct MonomerDefault {
    vec4_t color;
};

struct ComponentTypeItem {
    char label[32] = "";
    uint32_t count = 0;             // Components of the type
    uint32_t atom_count = 0;
    float fraction = 0;             // Of all atoms
    md_array(int) comp_indices = 0; // The components of the type (for highlighting)
    bool use_defaults = false;      // Linked to the defaults of its monomer, only possible for a type with one
    uint32_t load_color = 0;        // The color it was given when it was inferred, the baseline modifications are stored against
};

// A single deserialized [ComponentType] section, buffered until the system is loaded (as AtomTypeOverride)
struct ComponentTypeOverride {
    char name[16] = "";
    md_component_kind_t kind = MD_COMPONENT_KIND_OTHER;  // Name and kind identify the type
    bool has_color = false;
    bool has_use_defaults = false;
    uint32_t color = 0;
    bool use_defaults = true;
};

struct ComponentTypeDelta {
    bool color    = false;
    bool defaults = false;

    explicit operator bool() const { return color || defaults; }
};

static inline bool comp_type_has_monomer(const md_component_type_data_t& type, size_t i) {
    return md_component_type_monomer(&type, i) != MD_MONOMER_NONE;
}

struct Dataset : viamd::EventHandler {
    bool show_window = false;
    char series_filter[64] = "";

    // Cached at ViamdInitialize: the serialize event only carries the serialization state,
    // so we need our own handle on the system in order to diff the atom types.
    ApplicationState* app_state = nullptr;
    
    // Dataset data (moved from ApplicationState.dataset)
    md_array(AtomElementMapping) atom_element_remappings = 0;
    md_array(DatasetItem) atom_types = 0;
    md_array(ComponentTypeItem) comp_types = 0;     // One per type of sys.component.type
    md_allocator_i* arena = 0;

    char atom_type_filter[64] = "";
    md_array(int) atom_type_order = 0;      // Display order of the atom types as sorted by the table, allocated in `arena`
    bool atom_type_resort     = false;      // The order is to be sorted again, e.g. after the atom types changed
    int  atom_type_selected   = -1;         // Atom type selected in the table, edited below it
    int  element_popup_z = -1;              // Element whose defaults are being edited

    ElementDefault element_defaults[MD_Z_Count] = {};

    char comp_type_filter[64] = "";
    md_array(int) comp_type_order = 0;      // Display order of the component types as sorted by the table, allocated in `arena`
    bool comp_type_resort   = false;
    int  comp_type_selected = -1;           // Component type selected in the table, edited below it
    int  monomer_popup      = -1;           // Monomer whose defaults are being edited

    MonomerDefault monomer_defaults[MD_MONOMER_COUNT] = {};

    // What the last Optimize colors of the atom and residue types did, shown under the button
    struct ColorOptimizePanel {
        bool has_result = false;
        size_t num_changed = 0;     // Types whose color it changed
        float min_distance = 0.0f;  // optimize_colors
    };
    ColorOptimizePanel atom_type_optimize;
    ColorOptimizePanel comp_type_optimize;

    // Atom type overrides read from a workspace, waiting for a system to be applied to.
    // Deliberately heap allocated: `arena` is reset by init_dataset_items, which runs after these are parsed.
    md_array(AtomTypeOverride) pending_overrides = 0;
    md_array(ComponentTypeOverride) pending_comp_type_overrides = 0;
    char pending_workspace[1024] = "";  // Workspace the pending overrides stem from, used to discard stale ones

    Dataset() { 
        viamd::event_system_register_handler(*this); 
    }
    
    ~Dataset() {
        // The arena will be cleaned up when the persistent allocator is destroyed
        // No need for explicit cleanup here
    }

    void init_element_defaults() {
        for (md_atomic_number_t z = 0; z < MD_Z_Count; ++z) {
            element_defaults[z].color   = vec4_from_u32(md_atomic_number_cpk_color(z));
            element_defaults[z].mass    = md_atomic_number_mass(z);
            element_defaults[z].radius  = md_atomic_number_vdw_radius(z);
        }
    }

    void init_monomer_defaults() {
        for (int m = 0; m < MD_MONOMER_COUNT; ++m) {
            monomer_defaults[m].color = vec4_from_u32(md_monomer_color((md_monomer_t)m));
        }
    }

    void clear_dataset_items() {
        comp_types = 0;
        atom_types = 0;
        atom_type_optimize = {};
        comp_type_optimize = {};
        atom_type_order = 0;
        atom_type_selected = -1;
        comp_type_order = 0;
        comp_type_selected = -1;
        if (arena) {
            md_arena_allocator_reset(arena);
        }
    }

    void init_dataset_items(ApplicationState& data) {
        if (!arena) {
            arena = md_arena_allocator_create(data.allocator.persistent, MEGABYTES(1));
        }
        clear_dataset_items();

        const md_system_t& sys = data.mold.sys;

        size_t type_count = md_system_atom_type_count(&sys);
        size_t atom_count = md_system_atom_count(&sys);

        if (atom_count == 0) return;

        md_array_resize(atom_types, type_count, arena);

        // Map atom types into dataset items
        for (size_t i = 0; i < type_count; ++i) {
            str_t atom_type_name = md_atom_type_name(&sys.atom.type, i);
            str_t ff_type = md_atom_type_ff_type(&sys.atom.type, i);
            DatasetItem item = { .key = i };
            if (!str_empty(ff_type) && !str_eq(ff_type, atom_type_name)) {
                // The force field type is what tells apart types of the same name (a Martini SC1 is a
                // different bead in every residue), so it is shown whenever it adds something
                snprintf(item.label, sizeof(item.label), STR_FMT " (" STR_FMT ")", STR_ARG(atom_type_name), STR_ARG(ff_type));
            } else {
                snprintf(item.label, sizeof(item.label), STR_FMT, STR_ARG(atom_type_name));
            }

            // Snapshot what the loader assigned, so that we can tell later on what the user has changed
            item.load.z      = sys.atom.type.z[i];
            item.load.radius = sys.atom.type.radius[i];
            item.load.mass   = sys.atom.type.mass[i];
            item.load.color  = sys.atom.type.color[i];
            item.load.flags  = sys.atom.type.flags[i];

            // Coarse grained types have no element to inherit from, so they always carry custom properties
            item.use_defaults = !type_is_bead(item.load.flags);

            atom_types[i] = item;
        }

        // Count and set indices for each atom type
        for (size_t i = 0; i < atom_count; ++i) {
            md_atom_type_idx_t type_idx = sys.atom.type_idx[i]; 
            atom_types[type_idx].count += 1;
            md_array_push(atom_types[type_idx].indices, (int)i, arena);
        }
        
        // Calculate fractions
        for (size_t i = 0; i < type_count; ++i) {
            atom_types[i].fraction = atom_types[i].count / (float)atom_count;
        }

        init_component_type_items(data);
    }

    void init_component_type_items(ApplicationState& data) {
        const md_system_t& sys = data.mold.sys;
        const md_component_type_data_t& type = sys.component.type;
        const size_t type_count = md_component_type_count(&type);
        const size_t comp_count = md_system_component_count(&sys);
        const size_t atom_count = md_system_atom_count(&sys);
        if (type_count == 0 || !sys.component.type_idx || atom_count == 0) return;

        md_array_resize(comp_types, type_count, arena);
        for (size_t i = 0; i < type_count; ++i) {
            ComponentTypeItem item = {};
            snprintf(item.label, sizeof(item.label), STR_FMT, STR_ARG(md_component_type_name(&type, i)));
            item.load_color   = md_component_type_color(&type, i);
            item.use_defaults = comp_type_has_monomer(type, i);
            comp_types[i] = item;
        }
        for (size_t c = 0; c < comp_count; ++c) {
            const md_component_type_idx_t t = md_component_type_idx(&sys.component, c);
            if (t < 0 || (size_t)t >= type_count) continue;
            ComponentTypeItem& item = comp_types[t];
            item.count += 1;
            item.atom_count += (uint32_t)md_component_atom_count(&sys.component, c);
            md_array_push(item.comp_indices, (int)c, arena);
        }
        for (size_t i = 0; i < type_count; ++i) {
            comp_types[i].fraction = comp_types[i].atom_count / (float)atom_count;
        }
    }

    ElementDefaultDelta compute_element_default_delta(md_atomic_number_t z) const {
        const ElementDefault& def = element_defaults[z];
        ElementDefaultDelta delta;
        delta.color  = u32_from_vec4(def.color) != md_atomic_number_cpk_color(z);
        delta.radius = def.radius != md_atomic_number_vdw_radius(z);
        delta.mass   = def.mass   != md_atomic_number_mass(z);
        return delta;
    }

    // Write one [ElementDefault] section per element the user has customized. The table is not tied to the
    // loaded system, so every customized element is written, not just the ones the current system happens to use.
    void serialize_element_defaults(viamd::serialization_state_t& state) {
        for (int z = 0; z < MD_Z_Count; ++z) {
            ElementDefaultDelta delta = compute_element_default_delta((md_atomic_number_t)z);
            if (!delta) continue;

            const ElementDefault& def = element_defaults[z];
            viamd::write_section_header(state, STR_LIT("ElementDefault"));
            viamd::write_int(state, STR_LIT("Element"), z);

            if (delta.color)  viamd::write_vec4(state, STR_LIT("Color"),  def.color);
            if (delta.radius) viamd::write_flt (state, STR_LIT("Radius"), def.radius);
            if (delta.mass)   viamd::write_flt (state, STR_LIT("Mass"),   def.mass);
        }
    }

    // Unlike the atom types these need no buffering: the element defaults are application state with no
    // dependency on the loaded system, so a section can be applied the moment it is parsed.
    void deserialize_element_default(viamd::deserialization_state_t& state) {
        int  z = -1;
        bool has_color = false, has_radius = false, has_mass = false;
        vec4_t color = {};
        float radius = 0, mass = 0;

        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "Element")) {
                viamd::extract_int(z, arg);
            } else if (str_eq_cstr(ident, "Color")) {
                has_color = viamd::extract_flt_vec(color.elem, 4, arg);
            } else if (str_eq_cstr(ident, "Radius")) {
                has_radius = viamd::extract_flt(radius, arg);
            } else if (str_eq_cstr(ident, "Mass")) {
                has_mass = viamd::extract_flt(mass, arg);
            }
        }

        if (z < 0 || z >= MD_Z_Count) {
            MD_LOG_INFO("Dataset: skipping [ElementDefault] section with a missing or out of range Element entry");
            return;
        }

        ElementDefault& def = element_defaults[z];
        if (has_color)  def.color  = color;
        if (has_radius) def.radius = radius;
        if (has_mass)   def.mass   = mass;
    }

    // Push the customized element defaults onto every atom type which is linked to them, mirroring what editing
    // the table in the UI does. Property by property, so a mass the loader took from a force field is not
    // clobbered by the element's mass just because the user recolored that element.
    void apply_element_defaults_to_atom_types(ApplicationState& data, bool& radius_changed, bool& color_changed) {
        md_atom_type_data_t& type = data.mold.sys.atom.type;

        size_t type_count = md_system_atom_type_count(&data.mold.sys);
        if (type_count == 0 || md_array_size(atom_types) != type_count) return;

        for (size_t i = 0; i < type_count; ++i) {
            if (!atom_types[i].use_defaults) continue;

            const md_atomic_number_t z = type.z[i];
            const ElementDefaultDelta delta = compute_element_default_delta(z);
            if (!delta) continue;

            const ElementDefault& def = element_defaults[z];
            if (delta.radius) {
                radius_changed |= type.radius[i] != def.radius;
                type.radius[i] = def.radius;
            }
            if (delta.mass) {
                type.mass[i] = def.mass;
            }
            if (delta.color) {
                const uint32_t color = u32_from_vec4(def.color);
                color_changed |= type.color[i] != color;
                type.color[i] = color;
            }
        }
    }

    // Which properties of an atom type the user has changed since it was loaded.
    // Converts to true if any of them have, i.e. if the type needs to be stored in the workspace at all.
    struct AtomTypeDelta {
        bool elem     = false;
        bool radius   = false;
        bool mass     = false;
        bool color    = false;
        bool coarse   = false;
        bool defaults = false;

        explicit operator bool() const { return elem || radius || mass || color || coarse || defaults; }
    };

    AtomTypeDelta compute_atom_type_delta(const md_atom_type_data_t& type, const DatasetItem& item, size_t i) const {
        const AtomTypeLoadState& load = item.load;

        // The state this type is restored to before any [AtomType] entry is applied: what the loader assigned,
        // with the customized element defaults pushed on top if the type is linked to them. Diffing against that
        // is what keeps an element level edit in [ElementDefault] instead of duplicated across every type using it.
        float    base_radius = load.radius;
        float    base_mass   = load.mass;
        uint32_t base_color  = load.color;

        if (item.use_defaults) {
            const ElementDefaultDelta ed = compute_element_default_delta(type.z[i]);
            const ElementDefault& def = element_defaults[type.z[i]];
            if (ed.radius) base_radius = def.radius;
            if (ed.mass)   base_mass   = def.mass;
            if (ed.color)  base_color  = u32_from_vec4(def.color);
        }

        AtomTypeDelta delta;
        delta.elem   = type.z[i]      != load.z;
        delta.radius = type.radius[i] != base_radius;
        delta.mass   = type.mass[i]   != base_mass;
        delta.color  = type.color[i]  != base_color;
        delta.coarse = type_is_bead(type.flags[i]) != type_is_bead(load.flags);

        // use_defaults is not stored in the system, its implicit default follows whether the type is a bead
        delta.defaults = item.use_defaults != !type_is_bead(load.flags);
        return delta;
    }

    // ## Component types: their monomer defaults, what the user changed, and the workspace

    bool compute_monomer_default_delta(int m) const {
        return u32_from_vec4(monomer_defaults[m].color) != md_monomer_color((md_monomer_t)m);
    }

    uint32_t monomer_default_color(md_monomer_t m) const {
        return u32_from_vec4(monomer_defaults[m].color);
    }

    // The color a component type is restored to: that of the defaults of its monomer when it is linked to them, else
    // the one it was inferred with. Diffing against it keeps an edit of the defaults in [MonomerDefault], instead of
    // duplicated in every type which follows them.
    uint32_t comp_type_base_color(const md_component_type_data_t& type, size_t i) const {
        return comp_types[i].use_defaults ? monomer_default_color(md_component_type_monomer(&type, i)) : comp_types[i].load_color;
    }

    ComponentTypeDelta compute_comp_type_delta(const md_component_type_data_t& type, size_t i) const {
        ComponentTypeDelta delta;
        delta.color    = md_component_type_color(&type, i) != comp_type_base_color(type, i);
        // use_defaults is not stored in the system, its implicit default follows whether the type has a monomer
        delta.defaults = comp_types[i].use_defaults != comp_type_has_monomer(type, i);
        return delta;
    }

    // Push the defaults of its monomer onto a component type which is linked to them
    void link_comp_type_to_monomer(md_component_type_data_t& type, size_t i, bool& color_changed) {
        if (!comp_type_has_monomer(type, i)) return;
        const uint32_t color = monomer_default_color(md_component_type_monomer(&type, i));
        color_changed |= type.color[i] != color;
        type.color[i] = color;
    }

    // Restore a component type to its state right after loading: linked when it has a monomer
    void reset_comp_type(md_component_type_data_t& type, size_t i, bool& color_changed) {
        ComponentTypeItem& item = comp_types[i];
        item.use_defaults = comp_type_has_monomer(type, i);
        const uint32_t color = comp_type_base_color(type, i);
        color_changed |= type.color[i] != color;
        type.color[i] = color;
    }

    void apply_monomer_defaults_to_comp_types(ApplicationState& data, bool& color_changed) {
        md_component_type_data_t& type = data.mold.sys.component.type;
        const size_t type_count = md_component_type_count(&type);
        if (type_count == 0 || md_array_size(comp_types) != type_count) return;
        for (size_t i = 0; i < type_count; ++i) {
            if (comp_types[i].use_defaults) link_comp_type_to_monomer(type, i, color_changed);
        }
    }

    // Apply the buffered [ComponentType] sections onto the freshly inferred types, identified by name and kind
    void apply_pending_comp_type_overrides(ApplicationState& data, bool& color_changed) {
        const size_t num_overrides = md_array_size(pending_comp_type_overrides);
        if (num_overrides == 0) return;

        md_component_type_data_t& type = data.mold.sys.component.type;
        const size_t type_count = md_component_type_count(&type);
        if (type_count == 0 || md_array_size(comp_types) != type_count) return;

        for (size_t j = 0; j < num_overrides; ++j) {
            const ComponentTypeOverride& ovr = pending_comp_type_overrides[j];
            const md_component_type_idx_t t = md_component_type_find(&type, str_from_cstr(ovr.name), ovr.kind);
            if (t < 0) {
                MD_LOG_INFO("Dataset: workspace contains overrides for component type '%s' (%s) which is not present in the loaded system", ovr.name, md_component_kind_name(ovr.kind));
                continue;
            }
            if (ovr.has_use_defaults) {
                comp_types[t].use_defaults = ovr.use_defaults && comp_type_has_monomer(type, t);
            }
            if (ovr.has_color) {
                color_changed |= type.color[t] != ovr.color;
                type.color[t] = ovr.color;
            }
        }
    }

    // One [MonomerDefault] section per monomer whose defaults the user changed, identified by its id ("HIS")
    void serialize_monomer_defaults(viamd::serialization_state_t& state) {
        for (int m = MD_MONOMER_NONE + 1; m < MD_MONOMER_COUNT; ++m) {
            if (!compute_monomer_default_delta(m)) continue;
            viamd::write_section_header(state, STR_LIT("MonomerDefault"));
            viamd::write_str(state, STR_LIT("Monomer"), md_monomer_id((md_monomer_t)m));
            viamd::write_vec4(state, STR_LIT("Color"), monomer_defaults[m].color);
        }
    }

    // Applied the moment it is parsed: the defaults do not depend on the loaded system
    void deserialize_monomer_default(viamd::deserialization_state_t& state) {
        md_monomer_t m = MD_MONOMER_NONE;
        bool has_color = false;
        vec4_t color = {};

        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "Monomer")) {
                m = md_monomer_from_id(arg);
            } else if (str_eq_cstr(ident, "Color")) {
                has_color = viamd::extract_flt_vec(color.elem, 4, arg);
            }
        }
        if (m == MD_MONOMER_NONE) {
            MD_LOG_INFO("Dataset: skipping [MonomerDefault] section without a known Monomer entry");
            return;
        }
        if (has_color) monomer_defaults[m].color = color;
    }

    // One [ComponentType] section per component type the user modified: its name and kind, then what differs
    void serialize_component_types(viamd::serialization_state_t& state) {
        if (!app_state) return;
        const md_component_type_data_t& type = app_state->mold.sys.component.type;
        const size_t type_count = md_component_type_count(&type);
        if (type_count == 0 || md_array_size(comp_types) != type_count) return;

        for (size_t i = 0; i < type_count; ++i) {
            const ComponentTypeDelta delta = compute_comp_type_delta(type, i);
            if (!delta) continue;
            viamd::write_section_header(state, STR_LIT("ComponentType"));
            viamd::write_str(state, STR_LIT("Name"), md_component_type_name(&type, i));
            viamd::write_int(state, STR_LIT("Kind"), (int)md_component_type_kind(&type, i));
            if (delta.color)    viamd::write_vec4(state, STR_LIT("Color"), vec4_from_u32(type.color[i]));
            if (delta.defaults) viamd::write_bool(state, STR_LIT("UseDefaults"), comp_types[i].use_defaults);
        }
    }

    // Overrides are buffered per workspace. Sections of another workspace than the one buffered never made it onto a
    // system (it failed to load) and are stale.
    void begin_pending_overrides(const viamd::deserialization_state_t& state) {
        if (!str_eq_cstr(state.filename, pending_workspace)) {
            free_pending_overrides();
            str_copy_to_char_buf(pending_workspace, sizeof(pending_workspace), state.filename);
        }
    }

    void deserialize_component_type(viamd::deserialization_state_t& state) {
        begin_pending_overrides(state);
        ComponentTypeOverride ovr = {};

        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "Name")) {
                viamd::extract_to_char_buf(ovr.name, sizeof(ovr.name), arg);
            } else if (str_eq_cstr(ident, "Kind")) {
                int kind = 0;
                if (viamd::extract_int(kind, arg) && kind >= 0 && kind <= (int)MD_COMPONENT_FLAG_KIND_MASK) {
                    ovr.kind = (md_component_kind_t)kind;
                }
            } else if (str_eq_cstr(ident, "Color")) {
                vec4_t color = {};
                if (viamd::extract_flt_vec(color.elem, 4, arg)) {
                    ovr.color = u32_from_vec4(color);
                    ovr.has_color = true;
                }
            } else if (str_eq_cstr(ident, "UseDefaults")) {
                ovr.has_use_defaults = viamd::extract_bool(ovr.use_defaults, arg);
            }
        }
        if (ovr.name[0] == '\0') {
            MD_LOG_INFO("Dataset: skipping [ComponentType] section without a Name entry");
            return;
        }
        md_array_push(pending_comp_type_overrides, ovr, md_get_heap_allocator());
    }

    // Write one [AtomType] section per atom type the user has modified, containing the identifying name +
    // loaded element followed by only those properties which actually differ from what was loaded.
    void serialize_atom_types(viamd::serialization_state_t& state) {
        if (!app_state) return;

        const md_system_t& sys = app_state->mold.sys;
        const md_atom_type_data_t& type = sys.atom.type;

        size_t type_count = md_system_atom_type_count(&sys);
        if (type_count == 0 || md_array_size(atom_types) != type_count) return;

        for (size_t i = 0; i < type_count; ++i) {
            const DatasetItem& item = atom_types[i];
            if (i == 0 && item.count == 0) {
                // Skip sentinel "unknown" atom type if unused
                continue;
            }

            AtomTypeDelta delta = compute_atom_type_delta(type, item, i);
            if (!delta) continue;

            viamd::write_section_header(state, STR_LIT("AtomType"));
            viamd::write_str(state, STR_LIT("Name"), md_atom_type_name(&type, i));
            const str_t ff_type = md_atom_type_ff_type(&type, i);
            if (!str_empty(ff_type)) {
                viamd::write_str(state, STR_LIT("ForceFieldType"), ff_type);
            }
            viamd::write_int(state, STR_LIT("LoadedElement"), item.load.z);

            if (delta.elem)     viamd::write_int (state, STR_LIT("Element"),       type.z[i]);
            if (delta.radius)   viamd::write_flt (state, STR_LIT("Radius"),        type.radius[i]);
            if (delta.mass)     viamd::write_flt (state, STR_LIT("Mass"),          type.mass[i]);
            if (delta.color)    viamd::write_vec4(state, STR_LIT("Color"),         vec4_from_u32(type.color[i]));
            if (delta.coarse)   viamd::write_bool(state, STR_LIT("CoarseGrained"), type_is_bead(type.flags[i]));
            if (delta.defaults) viamd::write_bool(state, STR_LIT("UseDefaults"),   item.use_defaults);
        }
    }

    // Parse a single [AtomType] section. The system is not loaded at this point, so we only buffer it.
    void deserialize_atom_type(viamd::deserialization_state_t& state) {
        begin_pending_overrides(state);

        AtomTypeOverride ovr = {};

        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "Name")) {
                viamd::extract_to_char_buf(ovr.name, sizeof(ovr.name), arg);
            } else if (str_eq_cstr(ident, "ForceFieldType")) {
                viamd::extract_to_char_buf(ovr.ff_type, sizeof(ovr.ff_type), arg);
            } else if (str_eq_cstr(ident, "LoadedElement")) {
                int z;
                if (viamd::extract_int(z, arg)) {
                    ovr.load_z = (md_atomic_number_t)z;
                }
            } else if (str_eq_cstr(ident, "Element")) {
                int z;
                if (viamd::extract_int(z, arg)) {
                    ovr.z = (md_atomic_number_t)z;
                    ovr.has_z = true;
                }
            } else if (str_eq_cstr(ident, "Radius")) {
                ovr.has_radius = viamd::extract_flt(ovr.radius, arg);
            } else if (str_eq_cstr(ident, "Mass")) {
                ovr.has_mass = viamd::extract_flt(ovr.mass, arg);
            } else if (str_eq_cstr(ident, "Color")) {
                vec4_t color = {};
                if (viamd::extract_flt_vec(color.elem, 4, arg)) {
                    ovr.color = u32_from_vec4(color);
                    ovr.has_color = true;
                }
            } else if (str_eq_cstr(ident, "CoarseGrained")) {
                ovr.has_coarse_grained = viamd::extract_bool(ovr.coarse_grained, arg);
            } else if (str_eq_cstr(ident, "UseDefaults")) {
                ovr.has_use_defaults = viamd::extract_bool(ovr.use_defaults, arg);
            }
        }

        if (ovr.name[0] == '\0') {
            MD_LOG_INFO("Dataset: skipping [AtomType] section without a Name entry");
            return;
        }

        md_array_push(pending_overrides, ovr, md_get_heap_allocator());
    }

    void free_pending_overrides() {
        md_array_free(pending_overrides, md_get_heap_allocator());
        pending_overrides = 0;
        md_array_free(pending_comp_type_overrides, md_get_heap_allocator());
        pending_comp_type_overrides = 0;
        pending_workspace[0] = '\0';
    }

    // Apply the buffered overrides onto the freshly loaded atom types. Must run after init_dataset_items so
    // that the load state has been snapshotted: it stays the baseline, the user's modifications sit on top.
    // A property a section does not carry was never modified, so it keeps whatever the loader assigned.
    void apply_pending_atom_type_overrides(ApplicationState& data, bool& radius_changed, bool& color_changed) {
        size_t num_overrides = md_array_size(pending_overrides);
        if (num_overrides == 0) return;

        md_system_t& sys = data.mold.sys;
        md_atom_type_data_t& type = sys.atom.type;

        size_t type_count = md_system_atom_type_count(&sys);
        if (type_count == 0 || md_array_size(atom_types) != type_count) return;

        for (size_t j = 0; j < num_overrides; ++j) {
            const AtomTypeOverride& ovr = pending_overrides[j];
            str_t name = str_from_cstr(ovr.name);
            str_t ff_type = str_from_cstr(ovr.ff_type);

            // Identify the type by its name, its force field type (empty for sources without one, which
            // is also what a workspace written before force field types existed holds) and the element
            // it had upon load
            bool matched = false;
            for (size_t i = 0; i < type_count; ++i) {
                DatasetItem& item = atom_types[i];
                if (item.load.z != ovr.load_z) continue;
                if (!str_eq(md_atom_type_name(&type, i), name)) continue;
                if (!str_eq(md_atom_type_ff_type(&type, i), ff_type)) continue;

                if (ovr.has_z) {
                    type.z[i] = ovr.z;
                }
                if (ovr.has_radius) {
                    radius_changed |= type.radius[i] != ovr.radius;
                    type.radius[i] = ovr.radius;
                }
                if (ovr.has_mass) {
                    type.mass[i] = ovr.mass;
                }
                if (ovr.has_color) {
                    color_changed |= type.color[i] != ovr.color;
                    type.color[i] = ovr.color;
                }
                if (ovr.has_coarse_grained) {
                    set_type_bead(sys, i, ovr.coarse_grained);
                }
                if (ovr.has_use_defaults) {
                    item.use_defaults = ovr.use_defaults;
                }

                matched = true;
                break;
            }

            if (!matched) {
                MD_LOG_INFO("Dataset: workspace contains overrides for atom type '%s' which is not present in the loaded system", ovr.name);
            }
        }
    }

    void on_system_init(ApplicationState& data) {
        init_dataset_items(data);

        bool radius_changed = false;
        bool color_changed  = false;

        // Order matters. The per type overrides carry the UseDefaults flags, and those decide which types the
        // element defaults are allowed to touch, so they have to be in place before the defaults are pushed.
        apply_pending_atom_type_overrides(data, radius_changed, color_changed);
        apply_element_defaults_to_atom_types(data, radius_changed, color_changed);

        // The same for the component types and the defaults of their monomers
        apply_pending_comp_type_overrides(data, color_changed);
        apply_monomer_defaults_to_comp_types(data, color_changed);

        // Applied, or for a system they do not belong to
        free_pending_overrides();

        if (radius_changed) {
            data.mold.dirty_gpu_buffers |= MolBit_DirtyRadius;
        }
        if (color_changed) {
            flag_all_representations_as_dirty(&data);
        }
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                // Initialize component
                app_state = (ApplicationState*)e.payload;
                init_element_defaults();
                init_monomer_defaults();
                workspace_register_window("System", &show_window);
                break;
            }
            case viamd::EventType_ViamdDeserializeBegin:
                // A workspace without [ElementDefault], [AtomType], [MonomerDefault] or [ComponentType] sections has the built in ones
                init_element_defaults();
                init_monomer_defaults();
                free_pending_overrides();
                break;
            case viamd::EventType_ViamdShutdown:
                // Cleanup
                clear_dataset_items();
                free_pending_overrides();
                app_state = nullptr;
                break;
            case viamd::EventType_ViamdSystemInit: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                on_system_init(state);
                break;
            }
            case viamd::EventType_ViamdFrameTick: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                draw(state);
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("System", &show_window);
                break;
            case viamd::EventType_ViamdSerialize: {
                viamd::serialization_state_t& state = *(viamd::serialization_state_t*)e.payload;
                serialize_element_defaults(state);
                serialize_atom_types(state);
                serialize_monomer_defaults(state);
                serialize_component_types(state);
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t &state = *(viamd::deserialization_state_t*)e.payload;
                str_t section = viamd::section_header(state);
                if (str_eq_cstr(section, "ElementDefault")) {
                    deserialize_element_default(state);
                } else if (str_eq_cstr(section, "AtomType")) {
                    deserialize_atom_type(state);
                } else if (str_eq_cstr(section, "MonomerDefault")) {
                    deserialize_monomer_default(state);
                } else if (str_eq_cstr(section, "ComponentType")) {
                    deserialize_component_type(state);
                }
                break;
            }
            default:
                break;
            }
        }
    }

    static inline void handle_item_click(ApplicationState& state) {
        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
            // Shift + click on element header to select all atom types that use this element
            // Since the hover should already contain the correct selection we can just copy it to the filter mask to achieve this
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                md_bitmask_or(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
            } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                md_bitmask_andnot(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
            }
        }
    }

    // Return the row in the periodic table for a given atomic number
    static inline int element_row(int z) {
        if (z == 0) return 9; // Placeholder for unknown elements

        if (z <= 2) return 0;
        if (z <= 10) return 1;
        if (z <= 18) return 2;
        if (z <= 36) return 3;
        if (z <= 54) return 4;

        if (57 <= z && z <= 71) return 8; // Lanthanides
        if (89 <= z && z <= 103) return 9; // Actinides

        if (z <= 86) return 5;

        return 6; // For actinides and beyond
    }

    // Return the column in the periodic table for a given atomic number
    static inline int element_col(int z) {
        if (z == 0) return 0; // Placeholder for unknown elements
        if (z == 1) return 0;
        if (z == 2) return 17;
        if (z <= 4)  return z - 3;
        if (z <= 10) return z + 7;

        if (z <= 12) return z - 11;
        if (z <= 18) return z - 1;
        if (z <= 36) return (z - 1) % 18;
        if (z <= 54) return (z - 1) % 18;
        if (z <= 56) return (z - 1) % 18;
        if (z <= 71) return (z - 57) % 18 + 2; // Lanthanides
        if (z <= 86) return z - 69;
        if (z <= 88) return z - 87;
        if (z <= 103) return (z - 89) % 18 + 2; // Actinides
        if (z <= 118) return (z - 101) % 18;

        return 0;
    }

    struct PeriodicTableResult {
        bool hovered = false;
        bool clicked = false;
        int  z = -1;
    };

    static bool element_button(const char* lbl, vec4_t color) {
        const float button_size = ImGui::GetFontSize() * 1.5f;
        const vec4_t color_hover  = vec4_clamp(color + vec4_set1(0.2f), vec4_set1(0.0f), vec4_set1(1.0f));
        const vec4_t color_active = color * vec4_set(0.8f, 0.8f, 0.8f, 1.0f);

        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(color.x, color.y, color.z, color.w));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0, 0, 0, 1));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(color_hover.x,  color_hover.y,  color_hover.z,  1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(color_active.x, color_active.y, color_active.z, 1.0f));
		ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0, 0, 0, 1));

        bool clicked = ImGui::Button(lbl, ImVec2(button_size, button_size));

        ImGui::PopStyleColor(5);
		ImGui::PopStyleVar();

		return clicked;
    }

    // If `enabled_mask` is non-null: bit=1 means enabled.
    // Returns per-frame interaction info.
    static PeriodicTableResult periodic_table_widget(const ElementDefault* elem_defs, const uint64_t enabled_mask[2] = nullptr) {
        PeriodicTableResult res = {};

        const int width  = 18;
        const int height = 10;
        const float button_size = ImGui::GetFontSize() * 1.5f;
        const float spacing = ImGui::GetStyle().ItemSpacing.x * 0.25f;

        ImVec2 offset = ImGui::GetCursorPos();
        ImGui::Dummy(ImVec2(width * (button_size + spacing), height * (button_size + spacing)));

        for (int z = 0; z < MD_Z_Count; ++z) {
            int row = element_row(z);
            int col = element_col(z);

            ImGui::PushID(z);
            ImGui::SetCursorPos(offset + ImVec2(col * (button_size + spacing), row * (button_size + spacing)));

            bool disabled = enabled_mask ? !(enabled_mask[z / 64] & (1ULL << (z % 64))) : false;

            const ElementDefault& def = elem_defs[z];
            const vec4_t color = def.color;
            const str_t sym = md_atomic_number_symbol((md_atomic_number_t)z);

			ImGui::PushStyleVar(ImGuiStyleVar_DisabledAlpha, 0.25f);

            if (disabled) ImGui::BeginDisabled();
			if (element_button(str_ptr(sym), color)) {
                res.clicked = true;
                res.z = z;
            }
            if (disabled) ImGui::EndDisabled();

            ImGui::PopStyleVar();

            // Interaction is queried after the item
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                res.hovered = true;
                res.z = z;
            }

            str_t name = md_atomic_number_name((md_atomic_number_t)z);
            ImGui::SetItemTooltip("%d: " STR_FMT, z, STR_ARG(name));

            ImGui::PopID();
        }

        return res;
    }

    // ## Files: what the system is composed of

    struct FileRow {
        const char* role;
        str_t path;
        char details[128];
        bool missing = false;
    };

    static void file_row(const FileRow& row) {
        ImGui::TableNextRow();
        ImGui::PushID(row.path.ptr, row.path.ptr + row.path.len);

        ImGui::TableSetColumnIndex(0);
        ImGui::TextUnformatted(row.role);

        ImGui::TableSetColumnIndex(1);
        str_t name = row.path;
        extract_file(&name, row.path);
        char name_buf[256];
        str_copy_to_char_buf(name_buf, sizeof(name_buf), name);
        ImGui::Selectable(name_buf, false, ImGuiSelectableFlags_SpanAllColumns);
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
            ImGui::SetTooltip(STR_FMT, STR_ARG(row.path));
        }
        if (ImGui::BeginPopupContextItem("##file")) {
            if (ImGui::MenuItem("Copy path")) {
                char buf[1024];
                str_copy_to_char_buf(buf, sizeof(buf), row.path);
                ImGui::SetClipboardText(buf);
            }
            ImGui::EndPopup();
        }

        ImGui::TableSetColumnIndex(2);
        if (row.missing) {
            ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f), ICON_FA_TRIANGLE_EXCLAMATION " Not loaded");
            ImGui::SetItemTooltip("The workspace names this file, but it could not be loaded.\nIt is kept in the workspace when it is saved: put the file back and load the workspace again.");
        } else {
            ImGui::TextDisabled("%s", row.details);
        }
        ImGui::PopID();
    }

    static const char* file_role_label(const FileEntry& e, FileRole role) {
        switch (role) {
        case FileRole_System:     return "Structure";
        case FileRole_Trajectory: return "Trajectory";
        case FileRole_Supplement: return e.type == LoaderType_ITP ? "Topology" : "Supplement";
        case FileRole_Series:     return e.type == LoaderType_EDR ? "Energies" : "Series";
        default:                  return "";
        }
    }

    // The files the system was built from, in the order they were loaded: its own, which build and fill it
    // in, and those of its run. Loading them again in their places builds it again, which is what a workspace does.
    void draw_files(ApplicationState& data) {
        if (!ImGui::CollapsingHeader("Files", ImGuiTreeNodeFlags_DefaultOpen)) return;

        const FileList* list = &data.files.list;
        if (file_list_count(list) == 0) {
            ImGui::TextDisabled("Nothing loaded. Open files with File > Open File..., or drop them on the window.");
            return;
        }

        const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp;
        if (!ImGui::BeginTable("##files", 3, flags)) return;
        ImGui::TableSetupColumn("Role", ImGuiTableColumnFlags_WidthStretch, 1.2f);
        ImGui::TableSetupColumn("File", ImGuiTableColumnFlags_WidthStretch, 3.0f);
        ImGui::TableSetupColumn("Holds", ImGuiTableColumnFlags_WidthStretch, 2.5f);
        ImGui::TableHeadersRow();

        str_t groups[64];
        const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), &data), ARRAY_SIZE(groups));

        for (size_t i = 0; i < file_list_count(list); ++i) {
            const FileEntry& e = *file_list_get(list, i);
            const FileRole role = file_list_role(list, i);
            FileRow row = { file_role_label(e, role), str_from_cstr(e.path), "" };
            row.missing = e.missing;

            if (!e.missing) {
                char* d = row.details;
                const size_t cap = sizeof(row.details);
                switch (role) {
                case FileRole_System:
                    snprintf(d, cap, "%zu atoms%s", data.mold.sys.atom.count, md_system_is_coarse_grained(&data.mold.sys) ? ", coarse grained" : "");
                    break;
                case FileRole_Trajectory: {
                    const size_t num_frames = run_num_frames(&data);
                    const double* times = run_frame_times(&data);
                    if (num_frames > 0 && times && !md_unit_is_none(run_time_unit(&data))) {
                        char unit_buf[32];
                        const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), run_time_unit(&data));
                        snprintf(d, cap, "%zu frames, %.4g - %.4g %s", num_frames, times[0] * scl, times[num_frames - 1] * scl, unit_buf);
                    } else {
                        snprintf(d, cap, "%zu frames", num_frames);
                    }
                    break;
                }
                case FileRole_Series:
                    // The group the file published, which names its source
                    for (size_t g = 0; g < num_groups; ++g) {
                        if (str_eq(system_series_group_source(&data, groups[g]), row.path)) {
                            snprintf(d, cap, "%zu series", system_series_members(nullptr, 0, &data, groups[g]));
                            break;
                        }
                    }
                    break;
                case FileRole_Supplement:
                    snprintf(d, cap, STR_FMT, STR_ARG(loader::type_name(e.type)));
                    break;
                default:
                    break;
                }
            }
            file_row(row);
        }

        ImGui::EndTable();
    }

    // ## Composition: what the system is made of
    //
    // One row per entity - a kind of molecule - with how many molecules of it there are, rather than a
    // header per instance: a solvated protein is a handful of rows, not sixteen thousand waters.
    // Molecules are counted as what the bonds connect (sys.structure), since an instance is whatever
    // the file or the inference grouped: all the waters of a box are typically one. Computed when the
    // system changes, never per frame.

    struct EntityRow {
        int         entity = -1;        // -1: the atoms no instance holds
        const char* kind = "";
        char        name[96] = "";
        uint32_t    num_instances = 0;
        uint32_t    res_min = 0;        // components per instance
        uint32_t    res_max = 0;
        size_t      num_residues = 0;   // components, all instances
        size_t      num_molecules = 0;  // connected structures whose first atom is in the entity
        size_t      num_atoms = 0;
        double      mass = 0.0;         // amu
        double      charge = 0.0;       // e
        bool        polymer = false;
    };

    struct Composition {
        uint64_t key = 0;
        md_array(EntityRow) rows = nullptr;

        size_t num_atoms = 0;
        size_t num_components = 0;
        size_t num_molecules = 0;
        double mass = 0.0;
        bool   has_charge = false;
        double charge = 0.0;

        size_t bonds_topology = 0;
        size_t bonds_inferred = 0;
        size_t bonds_user = 0;
        size_t bonds_file = 0;
        size_t bonds_elastic = 0;   // Not bonds of the molecules, and so in none of the above

        size_t atoms_without_element = 0;   // neither an element nor a coarse grained bead, and with mass
        size_t unresolved_amino = 0;        // named like amino acids, backbone not found
        size_t unresolved_nucleic = 0;
    } composition;

    int  expanded_entity = -1;
    int  expanded_instance = -1;
    bool use_short_labels = true;

    // What an entity is: its kind, and for a small molecule what its residue is
    struct EntityCounts {
        size_t amino = 0;       // components which are amino acids
        size_t nucleic = 0;     // nucleotides
        size_t ion = 0;         // ions
        size_t beads = 0;       // coarse grained beads
    };

    static const char* entity_kind(md_entity_kind_t kind, const EntityCounts& n, size_t residues) {
        switch (kind) {
        case MD_ENTITY_KIND_WATER:      return "Water";
        case MD_ENTITY_KIND_PEPTIDE:    return n.beads ? "Protein (CG)" : "Protein";
        case MD_ENTITY_KIND_DNA:        return "DNA";
        case MD_ENTITY_KIND_RNA:        return "RNA";
        case MD_ENTITY_KIND_NUCLEIC:    return "Nucleic acid";
        case MD_ENTITY_KIND_BRANCHED:   return "Oligosaccharide";
        case MD_ENTITY_KIND_POLYMER:    return "Polymer";
        default: break;
        }
        if (n.ion && n.ion == residues)         return "Ion";
        if (n.amino && n.amino == residues)     return "Amino acid";
        if (n.nucleic && n.nucleic == residues) return "Nucleotide";
        if (n.beads)                            return "Coarse grained";
        return "Other";
    }

    // The inference marks the entities it named with " (*)"; the kind column says what they are
    static str_t entity_display_name(str_t desc) {
        if (str_ends_with(desc, STR_LIT(" (*)"))) desc = str_substr(desc, 0, desc.len - 4);
        return desc;
    }

    static uint64_t composition_key(const ApplicationState& data) {
        const md_system_t& sys = data.mold.sys;
        const size_t counts[] = { sys.atom.count, sys.entity.count, sys.instance.count, sys.component.count, sys.bond.count, sys.structure.count };
        uint64_t key = md_hash64(counts, sizeof(counts), 0);
        const md_attribute_t* charge = md_attributes_find(&sys.attributes, STR_LIT("atom/charge"));
        const uint64_t charge_version = charge ? charge->version : 0;
        key = md_hash64(&charge_version, sizeof(charge_version), key);
        return md_hash64(&data.mold.sys.atom.type_idx, sizeof(void*), key);
    }

    void compute_composition(ApplicationState& data) {
        const md_system_t& sys = data.mold.sys;
        Composition& c = composition;
        md_array_free(c.rows, data.allocator.persistent);
        c = Composition{};
        c.key = composition_key(data);

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        // Charges, when the file has them (a topology does)
        float* charge = nullptr;
        if (const md_attribute_t* attr = md_attributes_find(&sys.attributes, STR_LIT("atom/charge"))) {
            if (attr->format.rank == 1 && md_attribute_element_count(&attr->format) == sys.atom.count) {
                charge = md_temp_alloc_array(temp, float, sys.atom.count + 1);
                c.has_charge = md_attribute_extract_f32(charge, sys.atom.count, attr, md_attribute_slice_all(), md_unit_none()) == sys.atom.count;
                if (!c.has_charge) charge = nullptr;
            }
        }

        c.num_atoms = sys.atom.count;
        c.num_components = sys.component.count;
        c.num_molecules = sys.structure.count;

        // The row each atom counts toward; the last row is the atoms of no instance
        int32_t* atom_row = md_temp_alloc_array(temp, int32_t, sys.atom.count + 1);
        for (size_t a = 0; a < sys.atom.count; ++a) atom_row[a] = -1;
        EntityCounts* counts = md_temp_alloc_array(temp, EntityCounts, sys.entity.count + 1);
        for (size_t e = 0; e <= sys.entity.count; ++e) counts[e] = EntityCounts{};

        for (size_t e = 0; e < sys.entity.count; ++e) {
            EntityRow row = {};
            row.entity = (int)e;
            row.kind = "Other";
            row.polymer = md_entity_kind_is_polymer(md_entity_kind(&sys.entity, e));
            str_t desc = entity_display_name(md_entity_description(&sys.entity, e));
            if (str_empty(desc)) desc = md_entity_id(&sys.entity, e);
            str_copy_to_char_buf(row.name, sizeof(row.name), desc);
            row.res_min = UINT32_MAX;
            md_array_push(c.rows, row, data.allocator.persistent);
        }

        for (size_t i = 0; i < sys.instance.count; ++i) {
            const int e = md_instance_entity_idx(&sys.instance, i);
            if (e < 0 || (size_t)e >= md_array_size(c.rows)) continue;
            EntityRow& row = c.rows[e];
            const md_urange_t atoms = md_system_instance_atom_range(&sys, i);
            const md_urange_t comps = md_system_instance_comp_range(&sys, i);
            const uint32_t num_comps = comps.end - comps.beg;
            row.num_instances += 1;
            row.res_min = MIN(row.res_min, num_comps);
            row.res_max = MAX(row.res_max, num_comps);
            row.num_residues += num_comps;
            row.num_atoms += atoms.end - atoms.beg;
            for (uint32_t a = atoms.beg; a < atoms.end && a < sys.atom.count; ++a) {
                row.mass += md_atom_mass(&sys.atom, a);
                if (charge && !atom_attribute_value_absent(charge[a])) row.charge += charge[a];
                atom_row[a] = e;
                counts[e].beads += md_atom_particle_kind(&sys.atom, a) == MD_PARTICLE_BEAD;
            }
            for (uint32_t k = comps.beg; k < comps.end; ++k) {
                const md_component_kind_t kind = md_component_kind(&sys.component, k);
                counts[e].amino   += kind == MD_COMPONENT_KIND_AMINO_ACID;
                counts[e].nucleic += kind == MD_COMPONENT_KIND_NUCLEOTIDE;
                counts[e].ion     += kind == MD_COMPONENT_KIND_ION;
            }
        }

        // Atoms no instance holds, so the rows add up to the system
        EntityRow rest = {};
        rest.kind = "Unassigned";
        snprintf(rest.name, sizeof(rest.name), "%s", sys.entity.count ? "In no entity" : "All atoms");
        const int32_t rest_row = (int32_t)md_array_size(c.rows);
        for (size_t a = 0; a < sys.atom.count; ++a) {
            if (atom_row[a] >= 0) continue;
            atom_row[a] = rest_row;
            rest.num_atoms += 1;
            rest.mass += md_atom_mass(&sys.atom, a);
            // An atom without a charge (a QM atom beside an embedding's sites) adds none
            if (charge && !atom_attribute_value_absent(charge[a])) rest.charge += charge[a];
        }
        if (rest.num_atoms > 0) {
            md_array_push(c.rows, rest, data.allocator.persistent);
        }

        // Molecules: each connected structure counts toward the row of its first atom
        if (sys.structure.count > 0 && sys.structure.offset && sys.structure.atom_idx) {
            for (size_t k = 0; k < sys.structure.count; ++k) {
                const int32_t root = sys.structure.atom_idx[sys.structure.offset[k]];
                if (root < 0 || (size_t)root >= sys.atom.count) continue;
                const int32_t r = atom_row[root];
                if (r >= 0 && (size_t)r < md_array_size(c.rows)) c.rows[r].num_molecules += 1;
            }
        }

        for (size_t i = 0; i < md_array_size(c.rows); ++i) {
            EntityRow& row = c.rows[i];
            if (row.res_min == UINT32_MAX) row.res_min = 0;
            if (row.entity >= 0) {
                row.kind = entity_kind(md_entity_kind(&sys.entity, row.entity), counts[row.entity], row.num_residues);
            }
            c.mass += row.mass;
            c.charge += row.charge;
        }

        for (size_t b = 0; b < sys.bond.count; ++b) {
            const uint32_t f = sys.bond.flags ? (uint32_t)sys.bond.flags[b] : 0u;
            if (f & MD_BOND_FLAG_ELASTIC) {
                c.bonds_elastic += 1;
                continue;
            }
            switch (md_bond_origin((md_bond_flags_t)f)) {
            case MD_BOND_ORIGIN_USER:     c.bonds_user     += 1; break;
            case MD_BOND_ORIGIN_TOPOLOGY: c.bonds_topology += 1; break;
            case MD_BOND_ORIGIN_INFERRED: c.bonds_inferred += 1; break;
            default:                      c.bonds_file     += 1; break;
            }
        }

        for (size_t a = 0; a < sys.atom.count; ++a) {
            if (atom_without_element(sys, a)) c.atoms_without_element += 1;
        }
        // The backbones of coarse grained residues are beads, which atom names cannot resolve, so they are not counted
        for (size_t i = 0; i < sys.component.count && !md_system_is_coarse_grained(&sys); ++i) {
            const md_component_flags_t f = md_component_flags(&sys.component, i);
            if (f & MD_COMPONENT_FLAG_RESOLVED) continue;
            c.unresolved_amino   += md_component_flags_kind(f) == MD_COMPONENT_KIND_AMINO_ACID;
            c.unresolved_nucleic += md_component_flags_kind(f) == MD_COMPONENT_KIND_NUCLEOTIDE;
        }
    }

    // No element, not a bead, and not a massless virtual site (TIP4P's M carries no element by design)
    static bool atom_without_element(const md_system_t& sys, size_t a) {
        if (md_atom_atomic_number(&sys.atom, a) != 0) return false;
        if (md_atom_particle_kind(&sys.atom, a) == MD_PARTICLE_BEAD) return false;
        return md_atom_mass(&sys.atom, a) > 0.0f;
    }

    static void highlight_entity(ApplicationState& data, int entity) {
        const md_system_t& sys = data.mold.sys;
        md_bitmask_t* hb = &data.selection.highlight_mask;
        md_bitmask_clear(hb);
        if (entity < 0) {
            // The atoms of no instance
            md_bitmask_set_range(hb, 0, (uint32_t)sys.atom.count);
            for (size_t i = 0; i < sys.instance.count; ++i) {
                const md_urange_t r = md_system_instance_atom_range(&sys, i);
                md_bitmask_clear_range(hb, r.beg, r.end);
            }
        } else {
            for (size_t i = 0; i < sys.instance.count; ++i) {
                if (md_instance_entity_idx(&sys.instance, i) == entity) {
                    const md_urange_t r = md_system_instance_atom_range(&sys, i);
                    md_bitmask_set_range(hb, r.beg, r.end);
                }
            }
        }
    }

    // A polymer's sequence, one chip per component, wrapping to the window
    void draw_sequence(ApplicationState& data, size_t inst_idx, bool short_labels) {
        const md_system_t& sys = data.mold.sys;
        const ImGuiStyle& style = ImGui::GetStyle();
        const md_urange_t range = md_system_instance_comp_range(&sys, inst_idx);
        const uint32_t max_shown = 4000;
        const uint32_t end = MIN(range.end, range.beg + max_shown);

        for (uint32_t comp_idx = range.beg; comp_idx < end; ++comp_idx) {
            str_t comp_name = md_component_name(&sys.component, comp_idx);
            const md_component_kind_t comp_kind = md_component_kind(&sys.component, comp_idx);
            const bool residue = comp_kind == MD_COMPONENT_KIND_AMINO_ACID || comp_kind == MD_COMPONENT_KIND_NUCLEOTIDE;

            // The one letter code of the standard monomer of its type (M for MSE), the full name of one none is known of
            const md_component_type_idx_t type_idx = md_component_type_idx(&sys.component, comp_idx);
            const md_monomer_t monomer = type_idx >= 0 ? md_component_type_monomer(&sys.component.type, (size_t)type_idx) : MD_MONOMER_NONE;
            const char letter[2] = { md_monomer_letter(monomer), '\0' };
            const bool short_label = short_labels && residue && monomer != MD_MONOMER_NONE && monomer != MD_MONOMER_UNK && monomer != MD_MONOMER_N;
            if (short_label) {
                comp_name = str_from_cstr(letter);
                ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
            }
            // In the color of its type, as the Res Type coloring of the representations shows it
            const bool colored = residue && type_idx >= 0;
            if (colored) {
                const uint32_t rgb = md_component_type_color(&sys.component.type, (size_t)type_idx) & 0x00FFFFFFU;
                const uint32_t bg  = rgb | 0xAA000000U;
                ImGui::PushStyleColor(ImGuiCol_Header, bg);
                ImGui::PushStyleColor(ImGuiCol_HeaderActive, bg);
                ImGui::PushStyleColor(ImGuiCol_Text, color_is_light(rgb) ? IM_COL32(0, 0, 0, 255) : IM_COL32(255, 255, 255, 255));
            }
            const ImVec2 text_sz = ImGui::CalcTextSize(str_beg(comp_name), str_end(comp_name));
            const ImVec2 item_sz(text_sz.x + style.ItemSpacing.x * 2.0f, text_sz.y + style.ItemSpacing.x * 2.0f);
            // Flow: stay on the same line only if it fits to the right of the previous chip
            if (comp_idx != range.beg) {
                const float last_x = ImGui::GetItemRectMax().x;
                const float max_x  = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
                if (last_x + style.ItemSpacing.x + item_sz.x <= max_x) {
                    ImGui::SameLine();
                }
            }
            char label[32];
            str_copy_to_char_buf(label, sizeof(label), comp_name);
            ImGui::PushID((int)comp_idx);
            ImGui::Selectable(label, true, 0, text_sz);
            ImGui::PopID();
            if (short_label) {
                ImGui::PopStyleVar();
            }
            if (colored) {
                ImGui::PopStyleColor(3);
            }
            if (ImGui::IsItemHovered()) {
                const str_t full = md_component_name(&sys.component, comp_idx);
                ImGui::SetTooltip(STR_FMT " %d", STR_ARG(full), md_component_seq_id(&sys.component, comp_idx));
                const md_urange_t atoms = md_system_component_atom_range(&sys, comp_idx);
                md_bitmask_assign(&data.selection.highlight_mask, md_bitset_range(atoms.beg, atoms.end));
                handle_item_click(data);
            }
        }
        if (end < range.end) {
            ImGui::TextDisabled("... and %u more", range.end - end);
        }
    }

    void draw_composition(ApplicationState& data) {
        if (!ImGui::CollapsingHeader("Composition", ImGuiTreeNodeFlags_DefaultOpen)) return;
        const md_system_t& sys = data.mold.sys;
        if (sys.atom.count == 0) {
            ImGui::TextDisabled("No system loaded");
            return;
        }
        if (composition.key != composition_key(data)) {
            compute_composition(data);
            expanded_entity = -1;
            expanded_instance = -1;
        }
        const Composition& c = composition;

        ImGui::TextDisabled("Hover a row to highlight its atoms, shift click to select them, click to list its instances.");

        const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp;
        const int num_cols = c.has_charge ? 7 : 6;
        if (ImGui::BeginTable("##composition", num_cols, flags)) {
            ImGui::TableSetupColumn("Kind",     ImGuiTableColumnFlags_WidthStretch, 1.3f);
            ImGui::TableSetupColumn("Name",     ImGuiTableColumnFlags_WidthStretch, 2.0f);
            ImGui::TableSetupColumn("Molecules", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("Residues",  ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("Atoms",    ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("Mass",     ImGuiTableColumnFlags_WidthStretch, 1.2f);
            if (c.has_charge) ImGui::TableSetupColumn("Charge", ImGuiTableColumnFlags_WidthStretch, 0.9f);
            // The header, with what the less obvious columns count
            ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
            const char* header_tips[7] = { nullptr, nullptr, "Separate molecules: groups of atoms connected by bonds", "Residues (components), all molecules together",
                                           nullptr, "Of all molecules together", "Net charge of all molecules together, from the topology" };
            for (int col = 0; col < num_cols; ++col) {
                if (!ImGui::TableSetColumnIndex(col)) continue;
                ImGui::TableHeader(ImGui::TableGetColumnName(col));
                if (header_tips[col] && ImGui::IsItemHovered()) ImGui::SetTooltip("%s", header_tips[col]);
            }

            for (size_t r = 0; r < md_array_size(c.rows); ++r) {
                const EntityRow& row = c.rows[r];
                ImGui::TableNextRow();
                ImGui::PushID((int)r);
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Selectable(row.kind, expanded_entity == (int)r, ImGuiSelectableFlags_SpanAllColumns)) {
                    if (!ImGui::IsKeyDown(ImGuiMod_Shift)) {
                        expanded_entity = expanded_entity == (int)r ? -1 : (int)r;
                        expanded_instance = -1;
                    }
                }
                if (ImGui::IsItemHovered()) {
                    highlight_entity(data, row.entity);
                    handle_item_click(data);
                }
                ImGui::TableSetColumnIndex(1); ImGui::TextUnformatted(row.name);
                ImGui::TableSetColumnIndex(2);
                if (row.num_molecules > 0) ImGui::Text("%zu", row.num_molecules); else ImGui::TextDisabled("-");
                ImGui::TableSetColumnIndex(3);
                if (row.entity >= 0) ImGui::Text("%zu", row.num_residues); else ImGui::TextDisabled("-");
                ImGui::TableSetColumnIndex(4); ImGui::Text("%zu", row.num_atoms);
                ImGui::TableSetColumnIndex(5); ImGui::Text("%.4g Da", row.mass);
                if (c.has_charge) {
                    ImGui::TableSetColumnIndex(6);
                    ImGui::Text("%+.3g", fabs(row.charge) < 5e-4 ? 0.0 : row.charge);
                }
                ImGui::PopID();
            }
            ImGui::EndTable();
        }

        // The copies of the expanded entity, and a polymer's sequence
        if (expanded_entity >= 0 && (size_t)expanded_entity < md_array_size(c.rows) && c.rows[expanded_entity].entity >= 0) {
            const EntityRow& row = c.rows[expanded_entity];
            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };
            md_array(uint32_t) instances = 0;
            for (size_t i = 0; i < sys.instance.count; ++i) {
                if (md_instance_entity_idx(&sys.instance, i) == row.entity) md_array_push(instances, (uint32_t)i, md_temp_allocator(temp));
            }
            const int count = (int)md_array_size(instances);

            ImGui::Indent();
            ImGui::Text("%s: %d %s", row.name, count, count == 1 ? "instance" : "instances");
            ImGui::SetItemTooltip("An instance is a chain or a group of molecules as the file (or the inference) grouped them");
            if (row.polymer) {
                ImGui::SameLine();
                ImGui::Checkbox("One letter codes", &use_short_labels);
            }

            if (count == 1 && row.polymer) {
                draw_sequence(data, instances[0], use_short_labels);
            } else {
                const float line = ImGui::GetTextLineHeightWithSpacing();
                const float height = line * (float)MIN(count, 10) + ImGui::GetStyle().WindowPadding.y * 2;
                if (ImGui::BeginChild("##instances", ImVec2(0, height), ImGuiChildFlags_Borders)) {
                    ImGuiListClipper clipper;
                    clipper.Begin(count);
                    if (expanded_instance >= 0) clipper.IncludeItemByIndex(expanded_instance);
                    while (clipper.Step()) {
                        for (int k = clipper.DisplayStart; k < clipper.DisplayEnd; ++k) {
                            const uint32_t inst = instances[k];
                            const str_t id   = md_instance_id(&sys.instance, inst);
                            const str_t auth = md_instance_auth_id(&sys.instance, inst);
                            const md_urange_t comps = md_system_instance_comp_range(&sys, inst);
                            const md_urange_t atoms = md_system_instance_atom_range(&sys, inst);
                            char label[128];
                            if (!str_empty(auth)) {
                                snprintf(label, sizeof(label), STR_FMT " (" STR_FMT ")   %u residues, %u atoms", STR_ARG(id), STR_ARG(auth), comps.end - comps.beg, atoms.end - atoms.beg);
                            } else {
                                snprintf(label, sizeof(label), STR_FMT "   %u residues, %u atoms", STR_ARG(id), comps.end - comps.beg, atoms.end - atoms.beg);
                            }
                            ImGui::PushID(k);
                            if (ImGui::Selectable(label, expanded_instance == k) && !ImGui::IsKeyDown(ImGuiMod_Shift) && row.polymer) {
                                expanded_instance = expanded_instance == k ? -1 : k;
                            }
                            if (ImGui::IsItemHovered()) {
                                md_bitmask_assign(&data.selection.highlight_mask, md_bitset_range((uint32_t)(atoms.beg), (uint32_t)(atoms.end)));
                                handle_item_click(data);
                            }
                            if (expanded_instance == k && row.polymer) {
                                ImGui::Indent();
                                draw_sequence(data, inst, use_short_labels);
                                ImGui::Unindent();
                            }
                            ImGui::PopID();
                        }
                    }
                }
                ImGui::EndChild();
            }
            ImGui::Unindent();
        }

        draw_totals(data);
        draw_warnings(data);
        draw_forcefield(data);
        ImGui::Spacing();
    }

    void draw_totals(ApplicationState& data) {
        const Composition& c = composition;
        ImGui::Spacing();
        ImGui::Text("%zu atoms, %zu residues, %zu molecules", c.num_atoms, c.num_components, c.num_molecules);
        ImGui::SetItemTooltip("Molecules: groups of atoms connected by bonds");

        char mass_buf[64];
        snprintf(mass_buf, sizeof(mass_buf), "Mass %.6g Da", c.mass);
        ImGui::TextUnformatted(mass_buf);

        // Density, with a box: a quick check that the box and the contents belong together
        const md_unitcell_t& cell = data.mold.state.unitcell;
        mat3_t A = {0};
        md_unitcell_A_extract_float(A.elem, &cell);
        const double volume = fabs((double)mat3_determinant(A));   // Angstrom^3
        if ((cell.flags & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) && volume > 0.0) {
            const double density = c.mass / volume * 1.66053906660;    // amu/A^3 -> g/cm^3
            ImGui::SameLine();
            ImGui::Text(",  density %.4g g/cm\xc2\xb3", density);
            ImGui::SetItemTooltip("Mass over the volume of the box in the current frame. Water is close to 1.");
        }

        if (c.has_charge) {
            ImGui::Text("Net charge %+.3f e", fabs(c.charge) < 5e-4 ? 0.0 : c.charge);
        }

        ImGui::Text("%zu bonds", c.bonds_topology + c.bonds_inferred + c.bonds_user + c.bonds_file);
        if (c.bonds_topology + c.bonds_inferred + c.bonds_user + c.bonds_file > 0) {
            char buf[256];
            int len = 0;
            auto part = [&](size_t n, const char* what) {
                if (n == 0 || len >= (int)sizeof(buf)) return;
                len += snprintf(buf + len, sizeof(buf) - len, "%s%zu %s", len ? ", " : "", n, what);
            };
            part(c.bonds_topology, "from the topology");
            part(c.bonds_file,     "from the file");
            part(c.bonds_inferred, "guessed from distances");
            part(c.bonds_user,     "added by you");
            ImGui::SameLine();
            ImGui::TextDisabled("(%s)", buf);
        }
        if (c.bonds_elastic > 0) {
            ImGui::Text("%zu elastic bonds", c.bonds_elastic);
            ImGui::SameLine();
            ImGui::TextDisabled("(springs of the force field: shown by licorice and ball and stick when asked to)");
        }

        const md_unitcell_t& uc = cell;
        if (uc.flags & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) {
            char unit_buf[32];
            const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), md_unit_angstrom());
            const bool px = uc.flags & MD_UNITCELL_PBC_X, py = uc.flags & MD_UNITCELL_PBC_Y, pz = uc.flags & MD_UNITCELL_PBC_Z;
            char periodic[32] = "not periodic";
            if (px || py || pz) snprintf(periodic, sizeof(periodic), "periodic in %s%s%s", px ? "x" : "", py ? "y" : "", pz ? "z" : "");
            if (uc.flags & MD_UNITCELL_ORTHO) {
                ImGui::Text("Box %.4g x %.4g x %.4g %s, %s", uc.x * scl, uc.y * scl, uc.z * scl, unit_buf, periodic);
            } else {
                ImGui::Text("Triclinic box, %s", periodic);
                ImGui::SetItemTooltip("x %.4g, y %.4g, z %.4g\nxy %.4g, xz %.4g, yz %.4g %s", uc.x * scl, uc.y * scl, uc.z * scl, uc.xy * scl, uc.xz * scl, uc.yz * scl, unit_buf);
            }
        } else {
            ImGui::TextDisabled("No box");
        }
    }

    // Only what applies is shown; each highlights its atoms on hover
    void draw_warnings(ApplicationState& data) {
        const Composition& c = composition;
        const md_system_t& sys = data.mold.sys;
        const ImVec4 warn = ImVec4(1.0f, 0.75f, 0.2f, 1.0f);

        if (c.atoms_without_element > 0) {
            ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " %zu atoms have no element", c.atoms_without_element);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Their element could not be told from the file: they are drawn as unknown and\nhave no radius or colour of their own. Set it under Atom Types.");
                md_bitmask_t* hb = &data.selection.highlight_mask;
                md_bitmask_clear(hb);
                for (size_t a = 0; a < sys.atom.count; ++a) {
                    if (atom_without_element(sys, a)) md_bitmask_set_bit(hb, (uint32_t)a);
                }
                handle_item_click(data);
            }
        }
        auto residue_warning = [&](size_t count, md_component_kind_t kind, const char* what) {
            if (count == 0) return;
            ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " %zu residues are named like %s but their backbone was not found", count, what);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("They are not part of any chain, so cartoons, secondary structure and\nbackbone angles skip them. Their atom names may not follow the usual convention.");
                md_bitmask_t* hb = &data.selection.highlight_mask;
                md_bitmask_clear(hb);
                for (size_t i = 0; i < sys.component.count; ++i) {
                    const md_component_flags_t f = md_component_flags(&sys.component, i);
                    if (md_component_flags_kind(f) == kind && !(f & MD_COMPONENT_FLAG_RESOLVED)) {
                        const md_urange_t r = md_system_component_atom_range(&sys, i);
                        md_bitmask_set_range(hb, r.beg, r.end);
                    }
                }
                handle_item_click(data);
            }
        };
        residue_warning(c.unresolved_amino,   MD_COMPONENT_KIND_AMINO_ACID, "amino acids");
        residue_warning(c.unresolved_nucleic, MD_COMPONENT_KIND_NUCLEOTIDE, "nucleotides");

        if (c.has_charge && fabs(c.charge) >= 0.01) {
            ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " The system is not neutral (%+.3f e)", c.charge);
            ImGui::SetItemTooltip("With Ewald electrostatics (PME) a net charge is compensated by a uniform background,\nwhich is usually a sign of missing counter ions.");
        }
        if (c.bonds_user > 0) {
            ImGui::TextDisabled(ICON_FA_CIRCLE_INFO " %zu bonds were added by you", c.bonds_user);
            if (ImGui::IsItemHovered()) {
                md_bitmask_t* hb = &data.selection.highlight_mask;
                md_bitmask_clear(hb);
                for (size_t b = 0; b < sys.bond.count; ++b) {
                    if (sys.bond.flags && md_bond_origin(sys.bond.flags[b]) == MD_BOND_ORIGIN_USER) {
                        md_bitmask_set_bit(hb, (uint32_t)sys.bond.pairs[b].idx[0]);
                        md_bitmask_set_bit(hb, (uint32_t)sys.bond.pairs[b].idx[1]);
                    }
                }
            }
        }
    }

    // The interactions the simulation used, when the system came with them (a .tpr)
    void draw_forcefield(ApplicationState& data) {
        const md_nb_forcefield_t* ff = data.mold.sys.nonbonded;
        if (!ff) return;
        const md_nb_potential_t& pot = ff->potential;
        const char* coulomb = "none";
        switch (pot.coulomb) {
        case MD_NB_COULOMB_CUTOFF:          coulomb = "plain cut-off"; break;
        case MD_NB_COULOMB_REACTION_FIELD:  coulomb = "reaction field"; break;
        case MD_NB_COULOMB_EWALD:           coulomb = "Ewald (PME), real space"; break;
        default: break;
        }
        char unit_buf[32];
        const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), md_unit_nanometer());
        ImGui::TextDisabled("Force field: %zu non-bonded types, Lennard-Jones cut-off %.3g %s, electrostatics %s, cut-off %.3g %s",
            ff->num_types, sqrt(pot.lj_cutoff2) * scl, unit_buf, coulomb, sqrt(pot.coulomb_cutoff2) * scl, unit_buf);
        ImGui::SetItemTooltip("From the run input file. The Contacts window evaluates interaction energies with these.");
    }

    // ## Trajectory: how the frames were sampled and what they carry

    struct TrajectoryInfo {
        uint64_t key = 0;
        size_t   num_frames = 0;
        double   t0 = 0, t1 = 0;            // run's time unit
        double   dt_min = 0, dt_max = 0, dt_mean = 0;
        size_t   backwards_at = 0;          // first frame whose time is before the previous, 0 for none
        bool     has_time = false;
        bool     has_velocity = false;
        bool     has_force = false;
        char     sections[128] = "";        // written at their own interval (a .trr's velocities, forces)
        bool     has_box = false;
        double   vol_min = 0, vol_max = 0;  // Angstrom^3
    } traj;

    uint64_t trajectory_key(const ApplicationState& data) {
        const md_attribute_t* axis = run_time_axis(&data);
        const md_attribute_t* cell = run_attribute(&data, STR_LIT("unitcell"));
        const uint64_t v[3] = { axis ? axis->version : 0, cell ? cell->version : 0, (uint64_t)md_attributes_count(&data.mold.sys.attributes) };
        return md_hash64(v, sizeof(v), md_hash64_str(str_from_cstr(data.mold.run), 0));
    }

    void compute_trajectory(ApplicationState& data) {
        traj = TrajectoryInfo{};
        traj.key = trajectory_key(data);
        const md_attributes_t* attrs = &data.mold.sys.attributes;
        const size_t n = run_num_frames(&data);
        const double* t = run_frame_times(&data);
        traj.num_frames = n;
        traj.has_time = n > 0 && t && !md_unit_is_none(run_time_unit(&data));
        if (n > 0 && t) {
            traj.t0 = t[0];
            traj.t1 = t[n - 1];
            traj.dt_min = DBL_MAX;
            traj.dt_max = -DBL_MAX;
            for (size_t i = 1; i < n; ++i) {
                const double dt = t[i] - t[i - 1];
                if (dt < 0 && !traj.backwards_at) traj.backwards_at = i;
                traj.dt_min = MIN(traj.dt_min, dt);
                traj.dt_max = MAX(traj.dt_max, dt);
            }
            traj.dt_mean = n > 1 ? (t[n - 1] - t[0]) / (double)(n - 1) : 0.0;
            if (n < 2) traj.dt_min = traj.dt_max = 0;
        }

        traj.has_velocity = run_attribute(&data, STR_LIT("atom/velocity")) != nullptr;
        traj.has_force    = run_attribute(&data, STR_LIT("atom/force")) != nullptr;

        // Sections with an axis of their own: "<run>/trr/<name>/time"
        char buf[256];
        const str_t trr = run_attribute_path(buf, sizeof(buf), &data, STR_LIT("trr"));
        if (!str_empty(trr)) {
            int len = 0;
            for (md_attribute_iter_t it = md_attributes_iter_children(attrs, trr); md_attributes_next(&it);) {
                const md_attribute_t* axis = md_attributes_find_in(attrs, it.child_path, STR_LIT("time"));
                if (!axis) continue;
                const size_t m = axis->format.shape[0];
                len += snprintf(traj.sections + len, sizeof(traj.sections) - len, "%s" STR_FMT " in %zu frames", len ? ", " : "", STR_ARG(it.child), m);
                if (len >= (int)sizeof(traj.sections)) break;
            }
        }

        // The box, frame by frame
        const md_attribute_t* cell = run_attribute(&data, STR_LIT("unitcell"));
        if (cell && cell->data && md_attribute_element_count(&cell->format) == n * 9 && n > 0) {
            md_temp_scope_t temp = md_temp_begin();
            float* m = md_temp_alloc_array(temp, float, n * 9);
            if (md_attribute_extract_f32(m, n * 9, cell, md_attribute_slice_all(), md_unit_none()) == n * 9) {
                traj.vol_min = DBL_MAX;
                traj.vol_max = 0;
                for (size_t f = 0; f < n; ++f) {
                    const float* a = m + f * 9;
                    const double det = a[0] * ((double)a[4] * a[8] - (double)a[5] * a[7])
                                     - a[1] * ((double)a[3] * a[8] - (double)a[5] * a[6])
                                     + a[2] * ((double)a[3] * a[7] - (double)a[4] * a[6]);
                    const double v = fabs(det);
                    if (v <= 0) continue;
                    traj.has_box = true;
                    traj.vol_min = MIN(traj.vol_min, v);
                    traj.vol_max = MAX(traj.vol_max, v);
                }
            }
            md_temp_end(temp);
        }
    }

    void draw_trajectory(ApplicationState& data) {
        if (run_num_frames(&data) == 0) return;
        if (!ImGui::CollapsingHeader("Trajectory", ImGuiTreeNodeFlags_DefaultOpen)) return;
        if (traj.key != trajectory_key(data)) {
            compute_trajectory(data);
        }

        char tu[32] = "";
        const double ts = display_units::factor_print(tu, sizeof(tu), run_time_unit(&data));
        const ImVec4 warn = ImVec4(1.0f, 0.75f, 0.2f, 1.0f);

        if (traj.has_time) {
            ImGui::Text("%zu frames, %.6g - %.6g %s", traj.num_frames, traj.t0 * ts, traj.t1 * ts, tu);
        } else {
            ImGui::Text("%zu frames", traj.num_frames);
            ImGui::SameLine();
            ImGui::TextDisabled("(the file has no time)");
        }

        if (traj.num_frames > 1 && traj.has_time) {
            // Uniform when the steps agree to a part in a thousand; float times of an xtc wobble below that
            const bool uniform = (traj.dt_max - traj.dt_min) <= 1e-3 * fabs(traj.dt_mean);
            if (uniform) {
                ImGui::Text("A frame every %.4g %s", traj.dt_mean * ts, tu);
            } else {
                ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " Frames are not evenly spaced: %.4g - %.4g %s apart", traj.dt_min * ts, traj.dt_max * ts, tu);
                ImGui::SetItemTooltip("Typically trajectories from restarts joined together, or frames written at a changed interval.\nPlots over time are drawn at the actual times; playback steps frame by frame.");
            }
            if (traj.backwards_at) {
                ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " Time goes backwards at frame %zu", traj.backwards_at);
                ImGui::SetItemTooltip("Usually runs concatenated with overlap. Time based lookups (energy files, plots) will be unreliable.");
            }
        }

        char holds[256] = "positions";
        if (traj.has_velocity) strncat(holds, ", velocities", sizeof(holds) - strlen(holds) - 1);
        if (traj.has_force)    strncat(holds, ", forces",     sizeof(holds) - strlen(holds) - 1);
        ImGui::Text("Each frame holds %s", holds);
        if (traj.sections[0]) {
            ImGui::TextDisabled("Written at their own interval: %s", traj.sections);
        }

        if (traj.has_box) {
            char lu[32];
            const double ls = display_units::factor_print(lu, sizeof(lu), md_unit_angstrom());
            const double vs = ls * ls * ls;
            const double rel = traj.vol_max > 0 ? (traj.vol_max - traj.vol_min) / traj.vol_max : 0.0;
            if (rel < 1e-6) {
                ImGui::Text("Box: constant, volume %.5g %s\xc2\xb3", traj.vol_min * vs, lu);
            } else {
                ImGui::Text("Box: changes, volume %.5g - %.5g %s\xc2\xb3 (%.2g%%)", traj.vol_min * vs, traj.vol_max * vs, lu, rel * 100.0);
                ImGui::SetItemTooltip("A box that changes size: a constant pressure (NPT) run");
            }
        } else {
            ImGui::TextDisabled("Box: none in the frames");
        }

        const double frame = data.animation.frame;
        const size_t fi = (size_t)CLAMP(frame + 0.5, 0.0, (double)(traj.num_frames - 1));
        const double* t = run_frame_times(&data);
        if (traj.has_time && t) {
            ImGui::Text("Showing frame %zu, %.6g %s", fi, t[fi] * ts, tu);
        } else {
            ImGui::Text("Showing frame %zu", fi);
        }

        // What was loaded along the frames, and how it lines up with them
        str_t groups[64];
        const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), &data), ARRAY_SIZE(groups));
        for (size_t g = 0; g < num_groups; ++g) {
            const md_attribute_t* axis = md_attributes_find_in(&data.mold.sys.attributes, groups[g], STR_LIT("time"));
            char label[256];
            system_series_group_label(label, sizeof(label), &data, groups[g]);
            if (!axis || axis->format.shape[0] < 1) {
                ImGui::TextDisabled("%s: one value per frame", label);
                continue;
            }
            const size_t m = axis->format.shape[0];
            double ends[2] = {0, 0};
            md_attribute_extract_f64(&ends[0], 1, axis, md_attribute_slice_1(0),                  run_time_unit(&data));
            md_attribute_extract_f64(&ends[1], 1, axis, md_attribute_slice_1((uint32_t)(m - 1)), run_time_unit(&data));
            const double every = m > 1 ? (ends[1] - ends[0]) / (double)(m - 1) : 0.0;
            const bool covers = traj.has_time && ends[0] <= traj.t0 + 1e-6 * fabs(traj.t1 - traj.t0) && ends[1] >= traj.t1 - 1e-6 * fabs(traj.t1 - traj.t0);
            ImGui::TextDisabled("%s: every %.4g %s, %s", label, every * ts, tu,
                covers ? "covers every frame" : "covers part of the trajectory");
            if (!covers) {
                ImGui::SetItemTooltip("%.6g - %.6g %s, the frames span %.6g - %.6g %s", ends[0] * ts, ends[1] * ts, tu, traj.t0 * ts, traj.t1 * ts, tu);
            }
        }
        ImGui::Spacing();
    }

    // ASCII case insensitive 'needle occurs in haystack'
    static bool contains_ignore_case(str_t haystack, str_t needle) {
        if (needle.len == 0) return true;
        if (needle.len > haystack.len) return false;
        for (size_t i = 0; i + needle.len <= haystack.len; ++i) {
            if (str_eq_ignore_case(str_substr(haystack, i, needle.len), needle)) return true;
        }
        return false;
    }

    // Quantities loaded along the run - an energy file, the columns of an .xvg or a .csv - one node
    // per file. Each row is a series the Timelines and Distributions windows can plot: drag it into
    // a subplot of either, double click it to toggle it in the first timeline subplot, or pick one
    // from its context menu.
    void draw_series(ApplicationState& data) {
        str_t groups[64];
        const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), &data), ARRAY_SIZE(groups));
        if (num_groups == 0) return;

        if (!ImGui::CollapsingHeader("Series", ImGuiTreeNodeFlags_DefaultOpen)) return;

        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputTextWithHint("##series_filter", "Filter by name", series_filter, sizeof(series_filter));
        const str_t filter = str_from_cstr(series_filter);

        const md_attributes_t* table = &data.mold.sys.attributes;
        const str_t run = str_from_cstr(data.mold.run);

        for (size_t g = 0; g < num_groups; ++g) {
            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };

            const size_t num = system_series_members(nullptr, 0, &data, groups[g]);
            str_t* paths = md_temp_alloc_array(temp, str_t, num + 1);
            system_series_members(paths, num, &data, groups[g]);

            char group_label[256];
            system_series_group_label(group_label, sizeof(group_label), &data, groups[g]);

            ImGui::PushID((int)md_hash64_str(groups[g], 0));
            defer { ImGui::PopID(); };

            if (!ImGui::TreeNodeEx("##group", ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_SpanAvailWidth, "%s   %zu series", group_label, num)) {
                continue;
            }
            defer { ImGui::TreePop(); };

            const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp;
            if (!ImGui::BeginTable("##series", 4, table_flags)) {
                continue;
            }
            ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthStretch, 3.0f);
            ImGui::TableSetupColumn("Unit", ImGuiTableColumnFlags_WidthStretch, 1.5f);
            ImGui::TableSetupColumn("Samples", ImGuiTableColumnFlags_WidthStretch, 1.0f);
            ImGui::TableSetupColumn("Time", ImGuiTableColumnFlags_WidthStretch, 2.0f);
            ImGui::TableHeadersRow();

            for (size_t i = 0; i < num; ++i) {
                const md_attribute_t* attr = md_attributes_find(table, paths[i]);
                if (!attr) continue;

                const SeriesKey key = series_key(SeriesSource_System, paths[i]);
                char name[64];
                series_label(name, sizeof(name), &data, key);
                if (!contains_ignore_case(str_from_cstr(name), filter) && !contains_ignore_case(md_attribute_leaf(attr), filter)) {
                    continue;
                }

                // How attr() in the script reads it: the path below the run
                char expr[SERIES_PATH_CAP + 16] = "";
                if (str_begins_with(paths[i], run) && paths[i].len > run.len + 1) {
                    const str_t rel = str_substr(paths[i], run.len + 1);
                    snprintf(expr, sizeof(expr), "attr(\"" STR_FMT "\")", STR_ARG(rel));
                }

                bool in_subplot[PLOT_MAX_SUBPLOTS] = {};
                bool in_any = false;
                for (int s = 0; s < data.timeline.num_subplots; ++s) {
                    in_subplot[s] = plot_find_series(data.timeline.subplots[s], key) != -1;
                    in_any |= in_subplot[s];
                }

                const ImVec4 color = series_default_color(&data, key);

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::PushID(paths[i].ptr, paths[i].ptr + paths[i].len);

                ImPlot::ItemIcon(color);
                ImGui::SameLine();
                if (ImGui::Selectable(name, in_any, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowDoubleClick)) {
                    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                        const int idx = plot_find_series(data.timeline.subplots[0], key);
                        if (idx != -1) {
                            plot_remove_series(data.timeline.subplots[0], idx);
                        } else {
                            plot_add_series(&data, data.timeline.subplots[0], key);
                            data.timeline.show_window = true;
                        }
                    }
                }
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
                    ImGui::BeginTooltip();
                    ImGui::TextUnformatted(name);
                    if (!str_empty(attr->description)) {
                        ImGui::TextDisabled(STR_FMT, STR_ARG(attr->description));
                    }
                    ImGui::TextDisabled(STR_FMT, STR_ARG(attr->path));
                    if (expr[0]) ImGui::TextDisabled("In the script: %s", expr);
                    ImGui::Separator();
                    ImGui::TextUnformatted("Drag into a timeline or distribution subplot, double click to toggle it in the first timeline");
                    ImGui::EndTooltip();
                }
                if (ImGui::BeginDragDropSource()) {
                    series_set_drag_payload(TIMELINE_SERIES_DND, key, -1, name, color);
                    ImGui::EndDragDropSource();
                }
                if (ImGui::BeginPopupContextItem("##context")) {
                    ImGui::TextDisabled("Show in timeline");
                    for (int s = 0; s < data.timeline.num_subplots; ++s) {
                        char item[32];
                        snprintf(item, sizeof(item), "Subplot %d", s + 1);
                        if (ImGui::MenuItem(item, nullptr, in_subplot[s])) {
                            const int idx = plot_find_series(data.timeline.subplots[s], key);
                            if (idx != -1) {
                                plot_remove_series(data.timeline.subplots[s], idx);
                            } else {
                                plot_add_series(&data, data.timeline.subplots[s], key);
                                data.timeline.show_window = true;
                            }
                        }
                    }
                    if (data.timeline.num_subplots < PLOT_MAX_SUBPLOTS && ImGui::MenuItem("New subplot")) {
                        const int s = data.timeline.num_subplots++;
                        plot_add_series(&data, data.timeline.subplots[s], key);
                        data.timeline.show_window = true;
                    }
                    ImGui::Separator();
                    {
                        PlotSubplot& dist = data.distributions.subplots[0];
                        const int idx = plot_find_series(dist, key);
                        if (ImGui::MenuItem("Show distribution", nullptr, idx != -1)) {
                            if (idx != -1) {
                                plot_remove_series(dist, idx);
                            } else {
                                plot_add_series(&data, dist, key);
                                data.distributions.show_window = true;
                            }
                        }
                    }
                    if (expr[0]) {
                        ImGui::Separator();
                        if (ImGui::MenuItem("Copy script expression", expr)) {
                            ImGui::SetClipboardText(expr);
                        }
                    }
                    ImGui::EndPopup();
                }
                ImGui::PopID();

                ImGui::TableSetColumnIndex(1);
                {
                    char unit_buf[32] = "";
                    display_units::factor_print(unit_buf, sizeof(unit_buf), attr->unit);
                    ImGui::TextUnformatted(unit_buf[0] ? unit_buf : "-");
                }

                const uint32_t n = attr->format.shape[0];
                ImGui::TableSetColumnIndex(2);
                {
                    const size_t per_sample = md_attribute_element_count(&attr->format) / n;
                    if (per_sample > 1) {
                        ImGui::Text("%u x %zu", n, per_sample);
                    } else {
                        ImGui::Text("%u", n);
                    }
                }

                ImGui::TableSetColumnIndex(3);
                {
                    const md_attribute_t* axis = md_attributes_axis(table, attr);
                    double t[2] = {0, 0};
                    const md_attribute_slice_t first = md_attribute_slice_1(0);
                    const md_attribute_slice_t last  = md_attribute_slice_1(n - 1);
                    if (axis && axis->format.shape[0] == n &&
                        md_attribute_extract_f64(&t[0], 1, axis, first, md_unit_none()) == 1 &&
                        md_attribute_extract_f64(&t[1], 1, axis, last,  md_unit_none()) == 1)
                    {
                        char unit_buf[32] = "";
                        const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), axis->unit);
                        ImGui::Text("%.4g - %.4g %s", t[0] * scl, t[1] * scl, unit_buf);
                    } else {
                        ImGui::TextUnformatted("-");
                    }
                }
            }
            ImGui::EndTable();
        }
    }

    // ## Formatting helpers

    // Integer with thousands separators, e.g. 12,345
    static const char* fmt_count(char* buf, size_t cap, size_t n) {
        char tmp[32];
        const int len = snprintf(tmp, sizeof(tmp), "%zu", n);
        size_t out = 0;
        for (int i = 0; i < len && out + 2 < cap; ++i) {
            if (i > 0 && (len - i) % 3 == 0) buf[out++] = ',';
            buf[out++] = tmp[i];
        }
        buf[out] = '\0';
        return buf;
    }

    // Whether black text reads better than white on the color (RGBA, R in the lowest byte, sRGB)
    static bool color_is_light(uint32_t rgba) {
        const float r = (float)((rgba >>  0) & 0xFF) / 255.0f;
        const float g = (float)((rgba >>  8) & 0xFF) / 255.0f;
        const float b = (float)((rgba >> 16) & 0xFF) / 255.0f;
        return 0.2126f * r + 0.7152f * g + 0.0722f * b > 0.55f;
    }

    // Key value row for the tooltip tables
    static void info_row(const char* key) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::TextDisabled("%s", key);
        ImGui::TableSetColumnIndex(1);
    }

    // ## Optimized colors
    //
    // A button for the atom types and one for the residue types: the colors of the types the system uses made as
    // distinct as they can be (optimize_colors). Every type in use is recolored, starting from its default, so that the
    // result depends on the types the system uses only. A type with a color of meaning (an element, a standard residue,
    // the element of an ion) is pulled back to it and moves a little at most, the others move freely. What comes out
    // are ordinary colors of the types, unlinked from their defaults where they moved: stored in the workspace and
    // restored by Reset like any edit.

    static constexpr float color_optimize_pull = 0.5f;

    static ColorOptimizeDesc color_optimize_desc(const ApplicationState& data) {
        ColorOptimizeDesc desc;
        desc.linear_input = data.visuals.shading.colors == postprocess_pipeline::ColorSpace_Linear;
        return desc;
    }

    // The row of the table with the filter, Reset and Optimize colors, and under it what the last optimization did.
    // Returns true when the colors are to be optimized.
    static bool draw_type_table_toolbar(const char* id, bool show_filter, char* filter, size_t filter_cap, const char* filter_hint,
                                        size_t num_modified, bool& reset_clicked, const char* reset_tip,
                                        ColorOptimizePanel& panel, const char* optimize_tip) {
        const ImGuiStyle& style = ImGui::GetStyle();
        char reset_lbl[48] = "";
        float reset_w = 0;
        if (num_modified) {
            snprintf(reset_lbl, sizeof(reset_lbl), "Reset %zu modified", num_modified);
            reset_w = ImGui::CalcTextSize(reset_lbl).x + style.FramePadding.x * 2 + style.ItemSpacing.x;
        }
        const char* opt_lbl = "Optimize colors";
        const float opt_w = ImGui::CalcTextSize(opt_lbl).x + style.FramePadding.x * 2;

        ImGui::PushID(id);
        defer { ImGui::PopID(); };
        if (show_filter) {
            ImGui::SetNextItemWidth(-(reset_w + opt_w + style.ItemSpacing.x));
            ImGui::InputTextWithHint("##filter", filter_hint, filter, filter_cap);
            ImGui::SameLine();
        }
        reset_clicked = false;
        if (num_modified) {
            reset_clicked = ImGui::Button(reset_lbl);
            ImGui::SetItemTooltip("%s", reset_tip);
            ImGui::SameLine();
        }
        const bool run = ImGui::Button(opt_lbl);
        ImGui::SetItemTooltip("%s", optimize_tip);

        if (reset_clicked) panel.has_result = false;
        if (panel.has_result && !run) {
            char changed[48];
            if (panel.num_changed == 0) snprintf(changed, sizeof(changed), "No color changed");
            else snprintf(changed, sizeof(changed), "%zu type%s recolored", panel.num_changed, panel.num_changed == 1 ? "" : "s");
            if (panel.min_distance > 0.0f) {
                ImGui::TextDisabled("%s, closest pair %.3f apart (OKLab)", changed, panel.min_distance);
            } else {
                ImGui::TextDisabled("%s", changed);     // Fewer than two colors
            }
            if (panel.min_distance > 0.0f && panel.min_distance < 0.03f) {
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f), "Too many types to tell apart by color alone");
            }
        }
        return run;
    }

    // ## Atom types

    static constexpr float atom_type_min_mass   = 1.0f;
    static constexpr float atom_type_max_mass   = 500.0f;
    static constexpr float atom_type_min_radius = 0.1f;
    static constexpr float atom_type_max_radius = 20.0f;

    // Restore an atom type to its state right after loading: what the loader assigned, with the customized element
    // defaults pushed on top if it is linked to them (see on_system_init). This is the baseline compute_atom_type_delta
    // diffs against, so a reset type is no longer written to the workspace.
    void reset_atom_type(ApplicationState& data, size_t i, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_atom_type_data_t& type = data.mold.sys.atom.type;
        DatasetItem& item = atom_types[i];
        const AtomTypeLoadState& load = item.load;

        float    radius = load.radius;
        float    mass   = load.mass;
        uint32_t color  = load.color;

        item.use_defaults = !type_is_bead(load.flags);
        if (item.use_defaults) {
            const ElementDefaultDelta ed = compute_element_default_delta(load.z);
            const ElementDefault& def = element_defaults[load.z];
            if (ed.radius) radius = def.radius;
            if (ed.mass)   mass   = def.mass;
            if (ed.color)  color  = u32_from_vec4(def.color);
        }

        radius_changed |= type.radius[i] != radius;
        mass_changed   |= type.mass[i]   != mass;
        color_changed  |= type.color[i]  != color || type.z[i] != load.z;

        type.z[i]      = load.z;
        type.radius[i] = radius;
        type.mass[i]   = mass;
        type.color[i]  = color;
        set_type_bead(data.mold.sys, i, type_is_bead(load.flags));
    }

    // Push the element defaults onto an atom type which is linked to them
    void link_atom_type_to_element(ApplicationState& data, size_t i, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_atom_type_data_t& type = data.mold.sys.atom.type;
        const ElementDefault& def = element_defaults[type.z[i]];
        const uint32_t color = u32_from_vec4(def.color);

        radius_changed |= type.radius[i] != def.radius;
        mass_changed   |= type.mass[i]   != def.mass;
        color_changed  |= type.color[i]  != color;

        type.radius[i] = def.radius;
        type.mass[i]   = def.mass;
        type.color[i]  = color;
    }

    void atom_type_tooltip(const ApplicationState& data, size_t i) const {
        const md_atom_type_data_t& type = data.mold.sys.atom.type;
        const DatasetItem& item = atom_types[i];
        const bool cg = type_is_bead(type.flags[i]);
        const str_t name    = md_atom_type_name(&type, i);
        const str_t ff_type = md_atom_type_ff_type(&type, i);
        char num[32];

        ImGui::BeginTooltip();
        ImGui::Text(STR_FMT, STR_ARG(name));
        if (!str_empty(ff_type) && !str_eq(ff_type, name)) {
            ImGui::SameLine();
            ImGui::TextDisabled("force field type " STR_FMT, STR_ARG(ff_type));
        }
        ImGui::Separator();

        if (ImGui::BeginTable("##atom_type_info", 2, ImGuiTableFlags_SizingFixedFit)) {
            info_row("Element");
            if (cg) {
                ImGui::TextUnformatted("None, coarse grained bead");
            } else {
                const md_atomic_number_t z = type.z[i];
                ImGui::Text(STR_FMT " (" STR_FMT ", %d)", STR_ARG(md_atomic_number_name(z)), STR_ARG(md_atomic_number_symbol(z)), (int)z);
            }
            info_row("Radius");
            ImGui::Text("%.3f \xC3\x85", type.radius[i]);
            info_row("Mass");
            ImGui::Text("%.3f u", type.mass[i]);
            info_row("Atoms");
            ImGui::Text("%s (%.2f%%)", fmt_count(num, sizeof(num), item.count), item.fraction * 100.0f);
            info_row("Properties");
            ImGui::TextUnformatted(cg ? "Custom (coarse grained)" : item.use_defaults ? "From the element defaults" : "Custom");

            const AtomTypeDelta delta = compute_atom_type_delta(type, item, i);
            if (delta) {
                char buf[128] = "";
                size_t len = 0;
                const char* parts[] = {
                    delta.elem     ? "element"        : nullptr,
                    delta.radius   ? "radius"         : nullptr,
                    delta.mass     ? "mass"           : nullptr,
                    delta.color    ? "color"          : nullptr,
                    delta.coarse   ? "coarse grained" : nullptr,
                    delta.defaults ? "linking"        : nullptr,
                };
                for (const char* part : parts) {
                    if (!part) continue;
                    const int n = snprintf(buf + len, sizeof(buf) - len, "%s%s", len ? ", " : "", part);
                    if (n < 0 || (size_t)n >= sizeof(buf) - len) break;
                    len += (size_t)n;
                }
                info_row("Modified");
                ImGui::TextUnformatted(buf);
            }
            ImGui::EndTable();
        }
        ImGui::Separator();
        ImGui::TextDisabled("Shift + click to select, Shift + right click to deselect");
        ImGui::EndTooltip();
    }

    // Columns of the atom type table. Used as the column user IDs, so sorting keeps working when columns are hidden
    enum AtomTypeColumn {
        AtomTypeCol_Type,
        AtomTypeCol_Element,
        AtomTypeCol_Radius,
        AtomTypeCol_Mass,
        AtomTypeCol_Atoms,
        AtomTypeCol_Count
    };

    static int compare3(double a, double b) { return (a > b) - (a < b); }

    int compare_atom_types(const md_atom_type_data_t& type, int a, int b, ImGuiID column) const {
        const bool cg_a = type_is_bead(type.flags[a]);
        const bool cg_b = type_is_bead(type.flags[b]);
        switch (column) {
        case AtomTypeCol_Type:     return strcmp(atom_types[a].label, atom_types[b].label);
        case AtomTypeCol_Element:  return compare3(cg_a ? (int)MD_Z_Count : (int)type.z[a], cg_b ? (int)MD_Z_Count : (int)type.z[b]); // Beads have no element and go last
        case AtomTypeCol_Radius:   return compare3(type.radius[a], type.radius[b]);
        case AtomTypeCol_Mass:     return compare3(type.mass[a], type.mass[b]);
        case AtomTypeCol_Atoms:    return compare3(atom_types[a].count, atom_types[b].count);
        default: return 0;
        }
    }

    // Order the atom types by the sort specs of the table. Without any specs (tri-state sorting) the loader order is restored.
    void sort_atom_types(const md_atom_type_data_t& type, const ImGuiTableSortSpecs* specs) {
        const size_t n = md_array_size(atom_type_order);
        for (size_t k = 0; k < n; ++k) atom_type_order[k] = (int)k;
        if (!specs || specs->SpecsCount == 0) return;

        std::stable_sort(atom_type_order, atom_type_order + n, [&](int a, int b) {
            for (int s = 0; s < specs->SpecsCount; ++s) {
                const ImGuiTableColumnSortSpecs& spec = specs->Specs[s];
                const int c = compare_atom_types(type, a, b, spec.ColumnUserID);
                if (c) return spec.SortDirection == ImGuiSortDirection_Descending ? c > 0 : c < 0;
            }
            return false;
        });
    }

    static void cell_text_right(const char* text, bool dim) {
        const float w = ImGui::CalcTextSize(text).x;
        const float avail = ImGui::GetContentRegionAvail().x;
        if (avail > w) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + avail - w);
        ImGui::AlignTextToFramePadding();
        if (dim) ImGui::TextDisabled("%s", text);
        else     ImGui::TextUnformatted(text);
    }

    // Editor of the selected atom type, shown below the table
    void draw_atom_type_editor(ApplicationState& data, size_t i, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_atom_type_data_t& type = data.mold.sys.atom.type;
        DatasetItem& item = atom_types[i];

        ImGui::PushID("##atom_type_editor");
        ImGui::PushID((int)i);
        defer { ImGui::PopID(); ImGui::PopID(); };

        {
            char title[64];
            snprintf(title, sizeof(title), "%s%s", item.label, compute_atom_type_delta(type, item, i) ? " *" : "");
            ImGui::SeparatorText(title);
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                atom_type_tooltip(data, i);
            }
        }

        bool coarse_grained = type_is_bead(type.flags[i]);
        if (ImGui::Checkbox("Coarse grained", &coarse_grained)) {
            set_type_bead(data.mold.sys, i, coarse_grained);
        }

        if (!coarse_grained) {
            const md_atomic_number_t z = type.z[i];
            ImGui::SameLine();
            ImGui::TextUnformatted("Element");
            ImGui::SameLine();
            if (element_button(str_ptr(md_atomic_number_symbol(z)), element_defaults[z].color)) {
                ImGui::OpenPopup("##element_popup");
            }
            ImGui::SetItemTooltip(STR_FMT ", click to change", STR_ARG(md_atomic_number_name(z)));
            if (ImGui::BeginPopup("##element_popup")) {
                const PeriodicTableResult res = periodic_table_widget(element_defaults);
                if (res.clicked && res.z >= 0) {
                    color_changed |= type.z[i] != (md_atomic_number_t)res.z;
                    type.z[i] = (md_atomic_number_t)res.z;
                    if (item.use_defaults) {
                        link_atom_type_to_element(data, i, radius_changed, color_changed, mass_changed);
                    }
                    ImGui::CloseCurrentPopup();
                }
                ImGui::EndPopup();
            }

            ImGui::SameLine();
            if (ImGui::Checkbox("Use element defaults", &item.use_defaults) && item.use_defaults) {
                // Relinking takes the element defaults as they currently are, including the user's edits of them
                link_atom_type_to_element(data, i, radius_changed, color_changed, mass_changed);
            }
            ImGui::SetItemTooltip("Take radius, mass and color from the element defaults and follow them when those are edited");
        } else {
            item.use_defaults = false; // Coarse grained types always carry custom properties
        }

        // Editing a value of a linked type unlinks it, as the value no longer follows the element defaults
        const ImGuiSliderFlags slider_flags = ImGuiSliderFlags_Logarithmic | ImGuiSliderFlags_NoRoundToFormat | ImGuiSliderFlags_AlwaysClamp;
        bool edited = false;
        if (ImGui::SliderFloat("Radius", &type.radius[i], atom_type_min_radius, atom_type_max_radius, "%.3f \xC3\x85", slider_flags)) {
            radius_changed = edited = true;
        }
        ImGui::SetItemTooltip("Ctrl + click to enter a value");
        if (ImGui::SliderFloat("Mass", &type.mass[i], atom_type_min_mass, atom_type_max_mass, "%.3f u", slider_flags)) {
            mass_changed = edited = true;
        }
        ImGui::SetItemTooltip("Ctrl + click to enter a value");
        ImVec4 color = ImColor(type.color[i]);
        if (ImGui::ColorEdit4("Color", &color.x)) {
            type.color[i] = ImColor(color);
            color_changed = edited = true;
        }
        if (edited) item.use_defaults = false;

        if (item.use_defaults) {
            ImGui::TextDisabled("The values follow the element defaults, editing one unlinks the type");
        }

        if (compute_atom_type_delta(type, item, i)) {
            if (ImGui::Button("Reset")) {
                reset_atom_type(data, i, radius_changed, color_changed, mass_changed);
            }
            ImGui::SetItemTooltip("Restore what was loaded");
        }
    }

    void optimize_atom_type_colors(ApplicationState& data, bool& color_changed) {
        md_atom_type_data_t& type = data.mold.sys.atom.type;
        const size_t num_types = md_system_atom_type_count(&data.mold.sys);
        if (num_types == 0 || md_array_size(atom_types) != num_types) return;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_array(ColorOptimizeItem) items = 0;
        md_array(size_t) index = 0;
        for (size_t i = 0; i < num_types; ++i) {
            if (atom_types[i].count == 0) continue;     // Only the types in use
            // A type of an element is pulled to the color of its defaults, a bead or a type without an element starts
            // from the color it was loaded with and is free
            const bool has_element = !type_is_bead(type.flags[i]) && type.z[i] != 0;
            ColorOptimizeItem it;
            it.color  = has_element ? u32_from_vec4(element_defaults[type.z[i]].color) : atom_types[i].load.color;
            it.pull   = has_element ? color_optimize_pull : 0.0f;
            it.weight = (float)atom_types[i].count;
            md_array_push(items, it, md_temp_allocator(temp));
            md_array_push(index, i, md_temp_allocator(temp));
        }
        const ColorOptimizeResult res = optimize_colors(items, md_array_size(items), color_optimize_desc(data));
        atom_type_optimize.has_result   = true;
        atom_type_optimize.min_distance = res.min_distance;
        atom_type_optimize.num_changed  = 0;

        for (size_t k = 0; k < md_array_size(items); ++k) {
            const size_t i = index[k];
            DatasetItem& item = atom_types[i];
            if (items[k].pull > 0.0f) {
                const ElementDefault& def = element_defaults[type.z[i]];
                if (items[k].color == u32_from_vec4(def.color)) {
                    // It kept the color of its element: linked again, when that changes nothing else
                    if (!item.use_defaults && type.radius[i] == def.radius && type.mass[i] == def.mass) item.use_defaults = true;
                } else {
                    item.use_defaults = false;  // A color of its own now
                }
            }
            atom_type_optimize.num_changed += type.color[i] != items[k].color ? 1 : 0;
            color_changed |= type.color[i] != items[k].color;
            type.color[i] = items[k].color;
        }
    }

    void draw_atom_types(ApplicationState& data, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_system_t& sys = data.mold.sys;
        md_atom_type_data_t& type = sys.atom.type;
        const size_t num_types = md_system_atom_type_count(&sys);
        if (num_types == 0 || md_array_size(atom_types) != num_types) return;

        // The sentinel "unknown" type is only listed when used
        const size_t first = atom_types[0].count == 0 ? 1 : 0;

        char header[64];
        snprintf(header, sizeof(header), "Atom Types (%zu)###AtomTypes", num_types - first);
        if (!ImGui::CollapsingHeader(header, ImGuiTreeNodeFlags_DefaultOpen)) return;

        if (atom_type_selected >= (int)num_types) atom_type_selected = -1;

        size_t num_modified = 0;
        for (size_t i = first; i < num_types; ++i) {
            num_modified += compute_atom_type_delta(type, atom_types[i], i) ? 1 : 0;
        }

        const bool show_filter = num_types - first > 8;
        {
            bool reset_clicked = false;
            const bool optimize = draw_type_table_toolbar("##atom_type_toolbar", show_filter, atom_type_filter, sizeof(atom_type_filter), "Filter by name or element",
                num_modified, reset_clicked, "Restore the modified atom types (marked *) to what was loaded", atom_type_optimize,
                "Make the colors of the atom types in use as distinct as they can be. Every type starts from its default:\n"
                "colors you set are replaced, and Reset undoes it.\n\n"
                "Types of an element are pulled back to its color and move a little at most, where two are too close to tell apart.\n"
                "Coarse grained beads and types without an element move freely.");
            if (reset_clicked) {
                for (size_t i = first; i < num_types; ++i) {
                    if (compute_atom_type_delta(type, atom_types[i], i)) {
                        reset_atom_type(data, i, radius_changed, color_changed, mass_changed);
                    }
                }
            }
            if (optimize) {
                optimize_atom_type_colors(data, color_changed);
            }
        }
        const str_t filter = show_filter ? str_trim(str_from_cstr(atom_type_filter)) : STR_LIT("");

        auto passes_filter = [&](size_t i) {
            if (i < first) return false;
            if (str_empty(filter)) return true;
            const bool cg = type_is_bead(type.flags[i]);
            return contains_ignore_case(str_from_cstr(atom_types[i].label), filter) || (!cg && str_eq_ignore_case(md_atomic_number_symbol(type.z[i]), filter));
        };

        size_t num_rows = 0;
        for (size_t i = 0; i < num_types; ++i) {
            num_rows += passes_filter(i) ? 1 : 0;
        }

        if (num_rows == 0) {
            ImGui::TextDisabled("No atom type matches the filter");
        } else {
            draw_atom_type_table(data, num_rows, passes_filter, radius_changed, color_changed, mass_changed);
        }

        // The editor stays when the filter hides the selected type, its header tells which type it is
        if (atom_type_selected >= 0) {
            draw_atom_type_editor(data, (size_t)atom_type_selected, radius_changed, color_changed, mass_changed);
        } else {
            ImGui::TextDisabled("Select an atom type to edit it");
        }
        ImGui::Spacing();
    }

    template <typename Filter>
    void draw_atom_type_table(ApplicationState& data, size_t num_rows, const Filter& passes_filter, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_atom_type_data_t& type = data.mold.sys.atom.type;
        const size_t num_types = md_system_atom_type_count(&data.mold.sys);

        if (md_array_size(atom_type_order) != num_types) {
            md_array_resize(atom_type_order, num_types, arena);
            for (size_t k = 0; k < num_types; ++k) atom_type_order[k] = (int)k;
            atom_type_resort = true;
        }

        // Rows hold frame sized widgets, the table scrolls once there are more than fit in max_rows
        const ImGuiStyle& style = ImGui::GetStyle();
        const float  frame_h  = ImGui::GetFrameHeight();
        const float  row_h    = frame_h + style.CellPadding.y * 2.0f;
        const float  header_h = ImGui::GetFontSize() + style.CellPadding.y * 2.0f;
        const size_t max_rows = 12;
        const ImVec2 outer_size(0.0f, header_h + row_h * (float)MIN(num_rows, max_rows) + 2.0f);

        const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersInnerV |
            ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable | ImGuiTableFlags_Hideable | ImGuiTableFlags_Sortable | ImGuiTableFlags_SortTristate |
            ImGuiTableFlags_SizingFixedFit;
        if (!ImGui::BeginTable("##atom_types", AtomTypeCol_Count, table_flags, outer_size)) return;

        const float value_w = ImGui::CalcTextSize("000.000").x;
        const float swatch  = floorf(ImGui::GetTextLineHeight() * 0.75f);
        const float elem_w  = swatch + style.ItemSpacing.x + ImGui::CalcTextSize("Wm").x;

        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Type",    ImGuiTableColumnFlags_WidthStretch | ImGuiTableColumnFlags_NoHide, 0.0f, AtomTypeCol_Type);
        ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthFixed, elem_w,  AtomTypeCol_Element);
        ImGui::TableSetupColumn("Radius",  ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_PreferSortDescending, value_w, AtomTypeCol_Radius);
        ImGui::TableSetupColumn("Mass",    ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_PreferSortDescending, value_w, AtomTypeCol_Mass);
        ImGui::TableSetupColumn("Atoms",   ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_PreferSortDescending, 0.0f,    AtomTypeCol_Atoms);

        // Headers are submitted one by one to give them tooltips
        static const char* header_tips[AtomTypeCol_Count] = {
            "Name of the atom type (force field type), * marks modified types\nClick a row to edit the type, right click for more options",
            "Color and element, CG for coarse grained beads which have no element",
            "Radius in \xC3\x85ngstr\xC3\xB6m, dimmed when taken from the element defaults",
            "Mass in atomic mass units, dimmed when taken from the element defaults",
            "Number of atoms of the type (share of all atoms)",
        };
        ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
        for (int c = 0; c < AtomTypeCol_Count; ++c) {
            if (!ImGui::TableSetColumnIndex(c)) continue;
            ImGui::TableHeader(ImGui::TableGetColumnName(c));
            ImGui::SetItemTooltip("%s", header_tips[c]);
        }

        // Sorted when the specs change, not when values are edited, so a row does not move away while it is being edited
        if (ImGuiTableSortSpecs* specs = ImGui::TableGetSortSpecs()) {
            if (specs->SpecsDirty || atom_type_resort) {
                sort_atom_types(type, specs);
                specs->SpecsDirty = false;
                atom_type_resort  = false;
            }
        }

        if (ImGui::IsWindowHovered()) {
            md_bitmask_clear(&data.selection.highlight_mask);
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_array(int) rows = 0;
        md_array_ensure(rows, num_rows, md_temp_allocator(temp));
        for (size_t k = 0; k < num_types; ++k) {
            const int i = atom_type_order[k];
            if (passes_filter((size_t)i)) md_array_push(rows, i, md_temp_allocator(temp));
        }

        ImGuiListClipper clipper;
        clipper.Begin((int)md_array_size(rows), row_h);
        while (clipper.Step()) {
            for (int r = clipper.DisplayStart; r < clipper.DisplayEnd; ++r) {
                const size_t i = (size_t)rows[r];
                DatasetItem& item = atom_types[i];
                const bool cg       = type_is_bead(type.flags[i]);
                const bool linked   = item.use_defaults && !cg;
                const bool modified = (bool)compute_atom_type_delta(type, item, i);
                const bool selected = atom_type_selected == (int)i;

                ImGui::PushID((int)i);
                ImGui::TableNextRow(ImGuiTableRowFlags_None, row_h);

                // Type: spans the row for hovering and selection
                ImGui::TableSetColumnIndex(AtomTypeCol_Type);
                {
                    char lbl[48];
                    snprintf(lbl, sizeof(lbl), "%s%s", item.label, modified ? " *" : "");
                    ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.0f, 0.5f));
                    // Shift + click is taken by the selection of atoms (handle_item_click)
                    if (ImGui::Selectable(lbl, selected, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap, ImVec2(0, frame_h)) && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
                        atom_type_selected = selected ? -1 : (int)i;
                    }
                    ImGui::PopStyleVar();

                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenOverlappedByItem)) {
                        md_bitmask_clear(&data.selection.highlight_mask);
                        md_bitmask_set_indices(&data.selection.highlight_mask, (uint32_t*)item.indices, md_array_size(item.indices));
                        handle_item_click(data);
                        // Shift + right click is taken by deselect
                        if (!ImGui::IsKeyDown(ImGuiMod_Shift) && ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                            ImGui::OpenPopup("##atom_type_menu");
                        }
                    }
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                        atom_type_tooltip(data, i);
                    }
                    if (ImGui::BeginPopup("##atom_type_menu")) {
                        ImGui::TextUnformatted(item.label);
                        ImGui::Separator();
                        if (ImGui::MenuItem("Reset", nullptr, false, modified)) {
                            reset_atom_type(data, i, radius_changed, color_changed, mass_changed);
                        }
                        ImGui::SetItemTooltip("Restore what was loaded");
                        if (ImGui::MenuItem("Link to element defaults", nullptr, false, !cg && !item.use_defaults)) {
                            item.use_defaults = true;
                            link_atom_type_to_element(data, i, radius_changed, color_changed, mass_changed);
                        }
                        ImGui::EndPopup();
                    }
                }

                // Element: swatch in the color of the type and the symbol
                if (ImGui::TableSetColumnIndex(AtomTypeCol_Element)) {
                    const ImVec2 p  = ImGui::GetCursorScreenPos();
                    const float off = floorf((frame_h - swatch) * 0.5f);
                    ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(p.x, p.y + off), ImVec2(p.x + swatch, p.y + off + swatch), type.color[i] | IM_COL32_A_MASK, 2.0f);
                    ImGui::Dummy(ImVec2(swatch, frame_h));
                    ImGui::SameLine();
                    ImGui::AlignTextToFramePadding();
                    if (cg) {
                        ImGui::TextDisabled("CG");
                    } else {
                        ImGui::Text(STR_FMT, STR_ARG(md_atomic_number_symbol(type.z[i])));
                    }
                }

                // Values taken from the element defaults are dimmed
                char buf[48];
                if (ImGui::TableSetColumnIndex(AtomTypeCol_Radius)) {
                    snprintf(buf, sizeof(buf), "%.3f", type.radius[i]);
                    cell_text_right(buf, linked);
                }
                if (ImGui::TableSetColumnIndex(AtomTypeCol_Mass)) {
                    snprintf(buf, sizeof(buf), "%.3f", type.mass[i]);
                    cell_text_right(buf, linked);
                }
                if (ImGui::TableSetColumnIndex(AtomTypeCol_Atoms)) {
                    char num[32];
                    snprintf(buf, sizeof(buf), "%s (%.1f%%)", fmt_count(num, sizeof(num), item.count), item.fraction * 100.0f);
                    cell_text_right(buf, false);
                }

                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }

    // ## Element defaults

    void draw_element_defaults(ApplicationState& data, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_system_t& sys = data.mold.sys;
        md_atom_type_data_t& type = sys.atom.type;
        const size_t num_types = md_system_atom_type_count(&sys);
        if (num_types == 0 || md_array_size(atom_types) != num_types) return;

        // Elements used by atom types linked to the defaults, the others are shown disabled
        uint64_t elem_mask[2] = { 0 };
        for (size_t i = 0; i < num_types; ++i) {
            const int z = type.z[i];
            if (atom_types[i].use_defaults && atom_types[i].count > 0) {
                elem_mask[z / 64] |= (1ULL << (z % 64));
            }
        }
        if (!(elem_mask[0] || elem_mask[1])) return;
        if (!ImGui::CollapsingHeader("Element Defaults")) return;

        ImGui::Indent();
        defer { ImGui::Unindent(); };

        const PeriodicTableResult table_res = periodic_table_widget(element_defaults, elem_mask);
        if (table_res.hovered) {
            md_bitmask_t* hb = &data.selection.highlight_mask;
            md_bitmask_clear(hb);
            for (size_t i = 0; i < num_types; ++i) {
                const DatasetItem& item = atom_types[i];
                if (type.z[i] == table_res.z) {
                    md_bitmask_set_indices(hb, (uint32_t*)item.indices, md_array_size(item.indices));
                }
            }
            handle_item_click(data);

            if (table_res.clicked && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
                ImGui::OpenPopup("##element_default_popup");
                element_popup_z = table_res.z;
            }
        }
        ImGui::TextDisabled("Click an element to edit its defaults, which apply to all atom types linked to them");

        if (element_popup_z >= 0 && element_popup_z < MD_Z_Count && ImGui::BeginPopup("##element_default_popup")) {
            const int z = element_popup_z;
            const str_t sym  = md_atomic_number_symbol((md_atomic_number_t)z);
            const str_t name = md_atomic_number_name((md_atomic_number_t)z);
            ImGui::Text("%d: " STR_FMT " (" STR_FMT ")", z, STR_ARG(name), STR_ARG(sym));
            ImGui::Separator();

            ElementDefault& def = element_defaults[z];
            bool edited = false;
            edited |= ImGui::ColorEdit3("Color", def.color.elem);
            if (ImGui::InputFloat("Van der Waals Radius", &def.radius)) {
                def.radius = MAX(def.radius, 0.01f);
                edited = true;
            }
            if (ImGui::InputFloat("Atomic Mass", &def.mass)) {
                def.mass = MAX(def.mass, 0.001f);
                edited = true;
            }
            if (compute_element_default_delta((md_atomic_number_t)z)) {
                if (ImGui::Button("Reset")) {
                    def.color  = vec4_from_u32(md_atomic_number_cpk_color((md_atomic_number_t)z));
                    def.radius = md_atomic_number_vdw_radius((md_atomic_number_t)z);
                    def.mass   = md_atomic_number_mass((md_atomic_number_t)z);
                    edited = true;
                }
                ImGui::SetItemTooltip("Restore the built in values");
            }

            if (edited) {
                // Push onto every atom type of this element which is linked to the defaults
                for (size_t i = 0; i < num_types; ++i) {
                    if (type.z[i] == z && atom_types[i].use_defaults) {
                        link_atom_type_to_element(data, i, radius_changed, color_changed, mass_changed);
                    }
                }
            }
            ImGui::EndPopup();
        }
    }

    void draw_element_mappings() {
        const size_t num_mappings = md_array_size(atom_element_remappings);
        if (num_mappings && ImGui::CollapsingHeader("Atom Element Mappings")) {
            for (size_t i = 0; i < num_mappings; ++i) {
                const auto& mapping = atom_element_remappings[i];
                ImGui::Text("%s -> %s (%s)", mapping.lbl, md_util_element_name(mapping.elem).ptr, md_util_element_symbol(mapping.elem).ptr);
            }
        }
    }

    // The atom types and what they are drawn with: an editor, kept out of the overview
    void draw_atom_types_tab(ApplicationState& data) {
        bool radius_changed = false;
        bool color_changed  = false;
        bool mass_changed   = false;

        draw_atom_types(data, radius_changed, color_changed, mass_changed);
        draw_element_defaults(data, radius_changed, color_changed, mass_changed);

        if (radius_changed) {
            data.mold.dirty_gpu_buffers |= MolBit_DirtyRadius;
        }
        if (color_changed) {
            // @NOTE: Only the color within representations needs to be updated, not the filter.
            flag_all_representations_as_dirty(&data);
        }
        (void)mass_changed; // Mass is not used for rendering, but tracked in case it is used for other purposes in the future

        draw_element_mappings();
    }

    // ## Component types (residue types)

    void highlight_comp_type(ApplicationState& data, size_t i) {
        const md_system_t& sys = data.mold.sys;
        md_bitmask_t* hb = &data.selection.highlight_mask;
        md_bitmask_clear(hb);
        const ComponentTypeItem& item = comp_types[i];
        for (size_t k = 0; k < md_array_size(item.comp_indices); ++k) {
            const md_urange_t r = md_system_component_atom_range(&sys, (size_t)item.comp_indices[k]);
            md_bitmask_set_range(hb, r.beg, r.end);
        }
    }

    void comp_type_tooltip(const ApplicationState& data, size_t i) const {
        const md_component_type_data_t& type = data.mold.sys.component.type;
        const ComponentTypeItem& item = comp_types[i];
        const md_monomer_t m = md_component_type_monomer(&type, i);
        char num[32];

        ImGui::BeginTooltip();
        ImGui::Text(STR_FMT, STR_ARG(md_component_type_name(&type, i)));
        ImGui::Separator();
        if (ImGui::BeginTable("##comp_type_info", 2, ImGuiTableFlags_SizingFixedFit)) {
            info_row("Kind");
            ImGui::TextUnformatted(md_component_kind_name(md_component_type_kind(&type, i)));
            info_row("Standard");
            if (m != MD_MONOMER_NONE) {
                ImGui::Text(STR_FMT " (" STR_FMT ", %c)", STR_ARG(md_monomer_id(m)), STR_ARG(md_monomer_name(m)), md_monomer_letter(m));
            } else {
                ImGui::TextDisabled("None");
            }
            info_row("Residues");
            ImGui::TextUnformatted(fmt_count(num, sizeof(num), item.count));
            info_row("Atoms");
            ImGui::Text("%s (%.2f%%)", fmt_count(num, sizeof(num), item.atom_count), item.fraction * 100.0f);
            info_row("Color");
            ImGui::TextUnformatted(m == MD_MONOMER_NONE ? "Custom (no standard residue)" : item.use_defaults ? "From the residue defaults" : "Custom");
            const ComponentTypeDelta delta = compute_comp_type_delta(type, i);
            if (delta) {
                info_row("Modified");
                ImGui::TextUnformatted(delta.color && delta.defaults ? "color, linking" : delta.color ? "color" : "linking");
            }
            ImGui::EndTable();
        }
        ImGui::Separator();
        ImGui::TextDisabled("Shift + click to select, Shift + right click to deselect");
        ImGui::EndTooltip();
    }

    // Columns of the component type table, used as the column user IDs
    enum CompTypeColumn {
        CompTypeCol_Type,
        CompTypeCol_Kind,
        CompTypeCol_Standard,
        CompTypeCol_Residues,
        CompTypeCol_Count
    };

    int compare_comp_types(const md_component_type_data_t& type, int a, int b, ImGuiID column) const {
        switch (column) {
        case CompTypeCol_Type:     return strcmp(comp_types[a].label, comp_types[b].label);
        case CompTypeCol_Kind:     return compare3((int)md_component_type_kind(&type, a), (int)md_component_type_kind(&type, b));
        case CompTypeCol_Standard: {
            // Those without a standard residue go last
            const int ma = md_component_type_monomer(&type, a), mb = md_component_type_monomer(&type, b);
            return compare3(ma ? ma : (int)MD_MONOMER_COUNT, mb ? mb : (int)MD_MONOMER_COUNT);
        }
        case CompTypeCol_Residues: return compare3(comp_types[a].count, comp_types[b].count);
        default: return 0;
        }
    }

    void sort_comp_types(const md_component_type_data_t& type, const ImGuiTableSortSpecs* specs) {
        const size_t n = md_array_size(comp_type_order);
        for (size_t k = 0; k < n; ++k) comp_type_order[k] = (int)k;
        if (!specs || specs->SpecsCount == 0) return;

        std::stable_sort(comp_type_order, comp_type_order + n, [&](int a, int b) {
            for (int s = 0; s < specs->SpecsCount; ++s) {
                const ImGuiTableColumnSortSpecs& spec = specs->Specs[s];
                const int c = compare_comp_types(type, a, b, spec.ColumnUserID);
                if (c) return spec.SortDirection == ImGuiSortDirection_Descending ? c > 0 : c < 0;
            }
            return false;
        });
    }

    // Editor of the selected component type, shown below the table
    void draw_comp_type_editor(ApplicationState& data, size_t i, bool& color_changed) {
        md_component_type_data_t& type = data.mold.sys.component.type;
        ComponentTypeItem& item = comp_types[i];
        const md_monomer_t m = md_component_type_monomer(&type, i);

        ImGui::PushID("##comp_type_editor");
        ImGui::PushID((int)i);
        defer { ImGui::PopID(); ImGui::PopID(); };

        {
            char title[64];
            snprintf(title, sizeof(title), "%s%s", item.label, compute_comp_type_delta(type, i) ? " *" : "");
            ImGui::SeparatorText(title);
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                comp_type_tooltip(data, i);
            }
        }

        if (m != MD_MONOMER_NONE) {
            ImGui::AlignTextToFramePadding();
            ImGui::Text("Standard residue " STR_FMT " (" STR_FMT ")", STR_ARG(md_monomer_id(m)), STR_ARG(md_monomer_name(m)));
            ImGui::SameLine();
            if (ImGui::Checkbox("Use residue defaults", &item.use_defaults) && item.use_defaults) {
                // Relinking takes the defaults as they currently are, including the user's edits of them
                link_comp_type_to_monomer(type, i, color_changed);
            }
            ImGui::SetItemTooltip("Take the color from the defaults of the standard residue, and follow them when those are edited");
        } else {
            item.use_defaults = false;
            ImGui::TextDisabled("No standard residue: the type carries a color of its own");
        }

        // Editing the color of a linked type unlinks it, as it no longer follows the defaults
        ImVec4 color = ImColor(type.color[i]);
        if (ImGui::ColorEdit3("Color", &color.x)) {
            color.w = 1.0f;
            type.color[i] = ImColor(color);
            color_changed = true;
            item.use_defaults = false;
        }
        if (item.use_defaults) {
            ImGui::TextDisabled("The color follows the residue defaults, editing it unlinks the type");
        }

        if (compute_comp_type_delta(type, i)) {
            if (ImGui::Button("Reset")) {
                reset_comp_type(type, i, color_changed);
            }
            ImGui::SetItemTooltip("Restore what was loaded");
        }
    }

    template <typename Filter>
    void draw_comp_type_table(ApplicationState& data, size_t num_rows, const Filter& passes_filter, bool& color_changed) {
        md_component_type_data_t& type = data.mold.sys.component.type;
        const size_t num_types = md_component_type_count(&type);

        if (md_array_size(comp_type_order) != num_types) {
            md_array_resize(comp_type_order, num_types, arena);
            for (size_t k = 0; k < num_types; ++k) comp_type_order[k] = (int)k;
            comp_type_resort = true;
        }

        const ImGuiStyle& style = ImGui::GetStyle();
        const float  frame_h  = ImGui::GetFrameHeight();
        const float  row_h    = frame_h + style.CellPadding.y * 2.0f;
        const float  header_h = ImGui::GetFontSize() + style.CellPadding.y * 2.0f;
        const size_t max_rows = 12;
        const ImVec2 outer_size(0.0f, header_h + row_h * (float)MIN(num_rows, max_rows) + 2.0f);

        const ImGuiTableFlags table_flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersInnerV |
            ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable | ImGuiTableFlags_Hideable | ImGuiTableFlags_Sortable | ImGuiTableFlags_SortTristate |
            ImGuiTableFlags_SizingFixedFit;
        if (!ImGui::BeginTable("##comp_types", CompTypeCol_Count, table_flags, outer_size)) return;

        const float swatch  = floorf(ImGui::GetTextLineHeight() * 0.75f);
        const float std_w   = swatch + style.ItemSpacing.x + ImGui::CalcTextSize("WWW (W)").x;
        const float kind_w  = ImGui::CalcTextSize("amino acid").x;

        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Type",     ImGuiTableColumnFlags_WidthStretch | ImGuiTableColumnFlags_NoHide, 0.0f, CompTypeCol_Type);
        ImGui::TableSetupColumn("Kind",     ImGuiTableColumnFlags_WidthFixed, kind_w, CompTypeCol_Kind);
        ImGui::TableSetupColumn("Standard", ImGuiTableColumnFlags_WidthFixed, std_w,  CompTypeCol_Standard);
        ImGui::TableSetupColumn("Residues", ImGuiTableColumnFlags_WidthFixed | ImGuiTableColumnFlags_PreferSortDescending, 0.0f, CompTypeCol_Residues);

        static const char* header_tips[CompTypeCol_Count] = {
            "Name of the residue type, * marks modified types\nClick a row to edit the type, right click for more options",
            "What the residues of the type are",
            "Color, and the standard residue (one letter code) the type is a form of, whose defaults it can follow",
            "Number of residues of the type (share of all atoms)",
        };
        ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
        for (int c = 0; c < CompTypeCol_Count; ++c) {
            if (!ImGui::TableSetColumnIndex(c)) continue;
            ImGui::TableHeader(ImGui::TableGetColumnName(c));
            ImGui::SetItemTooltip("%s", header_tips[c]);
        }

        if (ImGuiTableSortSpecs* specs = ImGui::TableGetSortSpecs()) {
            if (specs->SpecsDirty || comp_type_resort) {
                sort_comp_types(type, specs);
                specs->SpecsDirty = false;
                comp_type_resort  = false;
            }
        }

        if (ImGui::IsWindowHovered()) {
            md_bitmask_clear(&data.selection.highlight_mask);
        }

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_array(int) rows = 0;
        md_array_ensure(rows, num_rows, md_temp_allocator(temp));
        for (size_t k = 0; k < num_types; ++k) {
            const int i = comp_type_order[k];
            if (passes_filter((size_t)i)) md_array_push(rows, i, md_temp_allocator(temp));
        }

        ImGuiListClipper clipper;
        clipper.Begin((int)md_array_size(rows), row_h);
        while (clipper.Step()) {
            for (int r = clipper.DisplayStart; r < clipper.DisplayEnd; ++r) {
                const size_t i = (size_t)rows[r];
                ComponentTypeItem& item = comp_types[i];
                const md_monomer_t m = md_component_type_monomer(&type, i);
                const bool modified = (bool)compute_comp_type_delta(type, i);
                const bool selected = comp_type_selected == (int)i;

                ImGui::PushID((int)i);
                ImGui::TableNextRow(ImGuiTableRowFlags_None, row_h);

                ImGui::TableSetColumnIndex(CompTypeCol_Type);
                {
                    char lbl[48];
                    snprintf(lbl, sizeof(lbl), "%s%s", item.label, modified ? " *" : "");
                    ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.0f, 0.5f));
                    if (ImGui::Selectable(lbl, selected, ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap, ImVec2(0, frame_h)) && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
                        comp_type_selected = selected ? -1 : (int)i;
                    }
                    ImGui::PopStyleVar();

                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenOverlappedByItem)) {
                        highlight_comp_type(data, i);
                        handle_item_click(data);
                        if (!ImGui::IsKeyDown(ImGuiMod_Shift) && ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                            ImGui::OpenPopup("##comp_type_menu");
                        }
                    }
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
                        comp_type_tooltip(data, i);
                    }
                    if (ImGui::BeginPopup("##comp_type_menu")) {
                        ImGui::TextUnformatted(item.label);
                        ImGui::Separator();
                        if (ImGui::MenuItem("Reset", nullptr, false, modified)) {
                            reset_comp_type(type, i, color_changed);
                        }
                        ImGui::SetItemTooltip("Restore what was loaded");
                        if (ImGui::MenuItem("Link to residue defaults", nullptr, false, m != MD_MONOMER_NONE && !item.use_defaults)) {
                            item.use_defaults = true;
                            link_comp_type_to_monomer(type, i, color_changed);
                        }
                        ImGui::EndPopup();
                    }
                }

                if (ImGui::TableSetColumnIndex(CompTypeCol_Kind)) {
                    ImGui::AlignTextToFramePadding();
                    ImGui::TextDisabled("%s", md_component_kind_name(md_component_type_kind(&type, i)));
                }

                // Standard: swatch in the color of the type, the monomer dimmed when the color follows its defaults
                if (ImGui::TableSetColumnIndex(CompTypeCol_Standard)) {
                    const ImVec2 p  = ImGui::GetCursorScreenPos();
                    const float off = floorf((frame_h - swatch) * 0.5f);
                    ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(p.x, p.y + off), ImVec2(p.x + swatch, p.y + off + swatch), type.color[i] | IM_COL32_A_MASK, 2.0f);
                    ImGui::Dummy(ImVec2(swatch, frame_h));
                    ImGui::SameLine();
                    ImGui::AlignTextToFramePadding();
                    if (m != MD_MONOMER_NONE) {
                        if (item.use_defaults) ImGui::TextDisabled(STR_FMT " (%c)", STR_ARG(md_monomer_id(m)), md_monomer_letter(m));
                        else                   ImGui::Text(STR_FMT " (%c)", STR_ARG(md_monomer_id(m)), md_monomer_letter(m));
                    } else {
                        ImGui::TextDisabled("-");
                    }
                }

                if (ImGui::TableSetColumnIndex(CompTypeCol_Residues)) {
                    char buf[48], num[32];
                    snprintf(buf, sizeof(buf), "%s (%.1f%%)", fmt_count(num, sizeof(num), item.count), item.fraction * 100.0f);
                    cell_text_right(buf, false);
                }

                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }

    void optimize_comp_type_colors(ApplicationState& data, bool& color_changed) {
        md_component_type_data_t& type = data.mold.sys.component.type;
        const size_t num_types = md_component_type_count(&type);
        if (num_types == 0 || md_array_size(comp_types) != num_types) return;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };
        md_array(ColorOptimizeItem) items = 0;
        md_array(size_t) index = 0;
        for (size_t i = 0; i < num_types; ++i) {
            if (comp_types[i].count == 0) continue;     // Only the types in use
            // A form of a standard residue is pulled to the color of its defaults, a monatomic ion to that of its element
            // (which it was inferred with), unless it is a bead without one. The others start from the color they were
            // inferred with and are free.
            const bool has_monomer = comp_type_has_monomer(type, i);
            bool element_ion = false;
            if (!has_monomer && md_component_type_kind(&type, i) == MD_COMPONENT_KIND_ION) {
                const md_urange_t atoms = md_system_component_atom_range(&data.mold.sys, (size_t)comp_types[i].comp_indices[0]);
                element_ion = atoms.end > atoms.beg && md_atom_atomic_number(&data.mold.sys.atom, atoms.beg) != 0;
            }
            ColorOptimizeItem it;
            it.color  = has_monomer ? monomer_default_color(md_component_type_monomer(&type, i)) : comp_types[i].load_color;
            it.pull   = has_monomer || element_ion ? color_optimize_pull : 0.0f;
            it.weight = (float)comp_types[i].atom_count;
            md_array_push(items, it, md_temp_allocator(temp));
            md_array_push(index, i, md_temp_allocator(temp));
        }
        const ColorOptimizeResult res = optimize_colors(items, md_array_size(items), color_optimize_desc(data));
        comp_type_optimize.has_result   = true;
        comp_type_optimize.min_distance = res.min_distance;
        comp_type_optimize.num_changed  = 0;

        for (size_t k = 0; k < md_array_size(items); ++k) {
            const size_t i = index[k];
            if (comp_type_has_monomer(type, i)) {
                // Linked to its defaults as long as it kept their color
                comp_types[i].use_defaults = items[k].color == monomer_default_color(md_component_type_monomer(&type, i));
            }
            comp_type_optimize.num_changed += type.color[i] != items[k].color ? 1 : 0;
            color_changed |= type.color[i] != items[k].color;
            type.color[i] = items[k].color;
        }
    }

    void draw_comp_types(ApplicationState& data, bool& color_changed) {
        md_system_t& sys = data.mold.sys;
        md_component_type_data_t& type = sys.component.type;
        const size_t num_types = md_component_type_count(&type);
        if (num_types == 0 || md_array_size(comp_types) != num_types) {
            ImGui::TextDisabled("No residues in this system");
            return;
        }

        // Types no residue is of any more (reclassified when the system was inferred again) are not listed
        auto in_use = [&](size_t i) { return comp_types[i].count > 0; };
        size_t num_used = 0;
        size_t num_modified = 0;
        for (size_t i = 0; i < num_types; ++i) {
            num_used     += in_use(i) ? 1 : 0;
            num_modified += compute_comp_type_delta(type, i) ? 1 : 0;
        }

        char header[64];
        snprintf(header, sizeof(header), "Residue Types (%zu)###ResidueTypes", num_used);
        if (!ImGui::CollapsingHeader(header, ImGuiTreeNodeFlags_DefaultOpen)) return;

        if (comp_type_selected >= (int)num_types) comp_type_selected = -1;

        const bool show_filter = num_used > 8;
        {
            bool reset_clicked = false;
            const bool optimize = draw_type_table_toolbar("##comp_type_toolbar", show_filter, comp_type_filter, sizeof(comp_type_filter), "Filter by name, standard residue or kind",
                num_modified, reset_clicked, "Restore the modified residue types (marked *) to what was loaded", comp_type_optimize,
                "Make the colors of the residue types in use as distinct as they can be. Every type starts from its default:\n"
                "colors you set are replaced, and Reset undoes it.\n\n"
                "Standard residues are pulled back to their palette color, and ions to the color of their element: they move a little\n"
                "at most, where two are too close to tell apart. The rest (ligands, lipids, water) move freely.");
            if (reset_clicked) {
                for (size_t i = 0; i < num_types; ++i) {
                    if (compute_comp_type_delta(type, i)) reset_comp_type(type, i, color_changed);
                }
            }
            if (optimize) {
                optimize_comp_type_colors(data, color_changed);
            }
        }
        const str_t filter = show_filter ? str_trim(str_from_cstr(comp_type_filter)) : STR_LIT("");

        auto passes_filter = [&](size_t i) {
            if (!in_use(i)) return false;
            if (str_empty(filter)) return true;
            const md_monomer_t m = md_component_type_monomer(&type, i);
            return contains_ignore_case(str_from_cstr(comp_types[i].label), filter) ||
                   (m != MD_MONOMER_NONE && str_eq_ignore_case(md_monomer_id(m), filter)) ||
                   contains_ignore_case(str_from_cstr(md_component_kind_name(md_component_type_kind(&type, i))), filter);
        };

        size_t num_rows = 0;
        for (size_t i = 0; i < num_types; ++i) {
            num_rows += passes_filter(i) ? 1 : 0;
        }

        if (num_rows == 0) {
            ImGui::TextDisabled("No residue type matches the filter");
        } else {
            draw_comp_type_table(data, num_rows, passes_filter, color_changed);
        }

        if (comp_type_selected >= 0) {
            draw_comp_type_editor(data, (size_t)comp_type_selected, color_changed);
        } else {
            ImGui::TextDisabled("Select a residue type to edit it");
        }
        ImGui::Spacing();
    }

    // A button in the color of a monomer, labeled with its id
    static bool monomer_button(const char* lbl, vec4_t color, float width) {
        const vec4_t color_hover  = vec4_clamp(color + vec4_set1(0.2f), vec4_set1(0.0f), vec4_set1(1.0f));
        const vec4_t color_active = color * vec4_set(0.8f, 0.8f, 0.8f, 1.0f);
        const bool light = color_is_light(u32_from_vec4(color));

        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
        ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(color.x, color.y, color.z, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_Border,        ImVec4(0, 0, 0, 1));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(color_hover.x,  color_hover.y,  color_hover.z,  1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(color_active.x, color_active.y, color_active.z, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_Text,          light ? ImVec4(0, 0, 0, 1) : ImVec4(1, 1, 1, 1));
        const bool clicked = ImGui::Button(lbl, ImVec2(width, ImGui::GetFrameHeight()));
        ImGui::PopStyleColor(5);
        ImGui::PopStyleVar();
        return clicked;
    }

    // The defaults of the standard monomers (tier 1), which the residue types linked to them follow
    void draw_monomer_defaults(ApplicationState& data, bool& color_changed) {
        md_component_type_data_t& type = data.mold.sys.component.type;
        const size_t num_types = md_component_type_count(&type);
        if (num_types == 0 || md_array_size(comp_types) != num_types) return;

        // Monomers of residue types linked to the defaults, the others are shown disabled
        bool used[MD_MONOMER_COUNT] = {};
        bool any = false;
        for (size_t i = 0; i < num_types; ++i) {
            if (comp_types[i].use_defaults && comp_types[i].count > 0) {
                used[md_component_type_monomer(&type, i)] = true;
                any = true;
            }
        }
        if (!any) return;
        if (!ImGui::CollapsingHeader("Residue Defaults")) return;

        ImGui::Indent();
        defer { ImGui::Unindent(); };

        const ImGuiStyle& style = ImGui::GetStyle();
        const float spacing = style.ItemSpacing.x * 0.25f;
        const float width   = ImGui::CalcTextSize("WWW").x + style.FramePadding.x * 2.0f;
        const float max_x   = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;

        auto draw_group = [&](const char* label, int beg, int end) {
            ImGui::TextDisabled("%s", label);
            for (int m = beg; m < end; ++m) {
                if (m != beg) {
                    if (ImGui::GetItemRectMax().x + spacing + width <= max_x) ImGui::SameLine(0, spacing);
                }
                char lbl[16];
                snprintf(lbl, sizeof(lbl), STR_FMT "##monomer", STR_ARG(md_monomer_id((md_monomer_t)m)));
                ImGui::BeginDisabled(!used[m]);
                if (monomer_button(lbl, monomer_defaults[m].color, width) && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
                    monomer_popup = m;
                    ImGui::OpenPopup("##monomer_default_popup");
                }
                ImGui::EndDisabled();
                if (ImGui::IsItemHovered()) {
                    // The atoms of the residue types which follow these defaults
                    md_bitmask_t* hb = &data.selection.highlight_mask;
                    md_bitmask_clear(hb);
                    for (size_t i = 0; i < num_types; ++i) {
                        if (!comp_types[i].use_defaults || md_component_type_monomer(&type, i) != m) continue;
                        for (size_t k = 0; k < md_array_size(comp_types[i].comp_indices); ++k) {
                            const md_urange_t r = md_system_component_atom_range(&data.mold.sys, (size_t)comp_types[i].comp_indices[k]);
                            md_bitmask_set_range(hb, r.beg, r.end);
                        }
                    }
                    handle_item_click(data);
                    ImGui::SetTooltip(STR_FMT " (%c)", STR_ARG(md_monomer_name((md_monomer_t)m)), md_monomer_letter((md_monomer_t)m));
                }
            }
        };
        draw_group("Amino acids", MD_MONOMER_ALA, MD_MONOMER_DA);
        draw_group("Nucleotides", MD_MONOMER_DA, MD_MONOMER_COUNT);
        ImGui::TextDisabled("Click a residue to edit its defaults, which apply to all residue types linked to them");

        if (monomer_popup > MD_MONOMER_NONE && monomer_popup < MD_MONOMER_COUNT && ImGui::BeginPopup("##monomer_default_popup")) {
            const md_monomer_t m = (md_monomer_t)monomer_popup;
            ImGui::Text(STR_FMT ": " STR_FMT " (%c)", STR_ARG(md_monomer_id(m)), STR_ARG(md_monomer_name(m)), md_monomer_letter(m));
            ImGui::Separator();

            MonomerDefault& def = monomer_defaults[m];
            bool edited = ImGui::ColorEdit3("Color", def.color.elem);
            if (compute_monomer_default_delta(m)) {
                if (ImGui::Button("Reset")) {
                    def.color = vec4_from_u32(md_monomer_color(m));
                    edited = true;
                }
                ImGui::SetItemTooltip("Restore the built in color");
            }
            if (edited) {
                def.color.w = 1.0f;
                for (size_t i = 0; i < num_types; ++i) {
                    if (comp_types[i].use_defaults && md_component_type_monomer(&type, i) == m) {
                        link_comp_type_to_monomer(type, i, color_changed);
                    }
                }
            }
            ImGui::EndPopup();
        }
    }

    void draw_comp_types_tab(ApplicationState& data) {
        bool color_changed = false;
        draw_comp_types(data, color_changed);
        draw_monomer_defaults(data, color_changed);
        if (color_changed) {
            flag_all_representations_as_dirty(&data);
        }
    }

    void draw(ApplicationState& data) {
        if (!show_window) return;

        ImGui::SetNextWindowSize(ImVec2(560, 700), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("System", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::IsWindowHovered()) {
                md_bitmask_clear(&data.selection.highlight_mask);
            }
            if (ImGui::BeginTabBar("##system_tabs")) {
                if (ImGui::BeginTabItem("Overview")) {
                    draw_files(data);
                    draw_composition(data);
                    draw_trajectory(data);
                    draw_series(data);
                    ImGui::EndTabItem();
                }
                if (ImGui::BeginTabItem("Atom Types")) {
                    draw_atom_types_tab(data);
                    ImGui::EndTabItem();
                }
                if (ImGui::BeginTabItem("Residue Types")) {
                    draw_comp_types_tab(data);
                    ImGui::EndTabItem();
                }
                ImGui::EndTabBar();
            }
        }
        ImGui::End();
    }
};

static Dataset instance;

}  // namespace dataset
