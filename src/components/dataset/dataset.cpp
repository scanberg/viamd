#include <viamd_event.h>
#include <viamd.h>
#include <event.h>
#include <serialization_utils.h>
#include <display_units.h>
#include <plot_series.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_log.h>
#include <md_system.h>
#include <md_util.h>
#include <md_nonbonded.h>

#include <imgui.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>

#include <string>

namespace dataset {

// Helper function to convert amino acids and nucleotides from three-letter to single-letter codes
static str_t convert_to_short(str_t in_str) {
    
    // Standard amino acids
    if (str_eq_cstr(in_str, "ALA")) return STR_LIT("A");
    if (str_eq_cstr(in_str, "ARG")) return STR_LIT("R");
    if (str_eq_cstr(in_str, "ASN")) return STR_LIT("N");
    if (str_eq_cstr(in_str, "ASP")) return STR_LIT("D");
    if (str_eq_cstr(in_str, "CYS")) return STR_LIT("C");
    if (str_eq_cstr(in_str, "GLU")) return STR_LIT("E");
    if (str_eq_cstr(in_str, "GLN")) return STR_LIT("Q");
    if (str_eq_cstr(in_str, "GLY")) return STR_LIT("G");
    if (str_eq_cstr(in_str, "HIS")) return STR_LIT("H");
    if (str_eq_cstr(in_str, "ILE")) return STR_LIT("I");
    if (str_eq_cstr(in_str, "LEU")) return STR_LIT("L");
    if (str_eq_cstr(in_str, "LYS")) return STR_LIT("K");
    if (str_eq_cstr(in_str, "MET")) return STR_LIT("M");
    if (str_eq_cstr(in_str, "PHE")) return STR_LIT("F");
    if (str_eq_cstr(in_str, "PRO")) return STR_LIT("P");
    if (str_eq_cstr(in_str, "SER")) return STR_LIT("S");
    if (str_eq_cstr(in_str, "THR")) return STR_LIT("T");
    if (str_eq_cstr(in_str, "TRP")) return STR_LIT("W");
    if (str_eq_cstr(in_str, "TYR")) return STR_LIT("Y");
    if (str_eq_cstr(in_str, "VAL")) return STR_LIT("V");

    // Standard nucleotides
    if (str_eq_cstr(in_str, "DA")) return STR_LIT("A");
    if (str_eq_cstr(in_str, "DC")) return STR_LIT("C");
    if (str_eq_cstr(in_str, "DG")) return STR_LIT("G");
    if (str_eq_cstr(in_str, "DT")) return STR_LIT("T");
    if (str_eq_cstr(in_str, "DU")) return STR_LIT("U");
    
    // Failed to map it, return itself
    return in_str;
}

// Essentially a bswap to go from RGBA to ABGR
#define RGBA_HEX(x) BSWAP32(x)

// Helper function to assign colors to components based on their type
static uint32_t component_color(str_t str) {
    if (str_eq_cstr(str, "DG")) return RGBA_HEX(0xD5B3EFAA); // Purple
    if (str_eq_cstr(str, "DT")) return RGBA_HEX(0x9FF3A0AA); // Green
    if (str_eq_cstr(str, "DC")) return RGBA_HEX(0xF8EE5CAA); // Yellow
    if (str_eq_cstr(str, "DU")) return RGBA_HEX(0xFE9D2DAA); // Orange
    if (str_eq_cstr(str, "DA")) return RGBA_HEX(0xFC697AAA); // Red

    return 0;
}

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
    md_flags_t flags  = MD_FLAG_NONE;
};

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
    md_array(int) sub_items = 0; // Indices into the items of the subcatagories: i.e. for chain -> unique residues types within that chain

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

struct Dataset : viamd::EventHandler {
    bool show_window = false;
    char series_filter[64] = "";

    // Cached at ViamdInitialize: the serialize event only carries the serialization state,
    // so we need our own handle on the system in order to diff the atom types.
    ApplicationState* app_state = nullptr;
    
    // Dataset data (moved from ApplicationState.dataset)
    md_array(AtomElementMapping) atom_element_remappings = 0;
    md_array(DatasetItem) inst_types = 0;
    md_array(DatasetItem) comp_types = 0;
    md_array(DatasetItem) atom_types = 0;
    md_allocator_i* arena = 0;

    ElementDefault element_defaults[MD_Z_Count] = {};

    // Atom type overrides read from a workspace, waiting for a system to be applied to.
    // Deliberately heap allocated: `arena` is reset by init_dataset_items, which runs after these are parsed.
    md_array(AtomTypeOverride) pending_overrides = 0;
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

    void clear_dataset_items() {
        inst_types = 0;
        comp_types = 0;
        atom_types = 0;
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
        size_t comp_count = md_system_component_count(&sys);
        size_t inst_count = md_system_instance_count(&sys);

        if (atom_count == 0) return;

        md_temp_scope_t temp_scope = md_temp_begin();
        md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
        defer { md_temp_end(temp_scope); };
        
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
            item.use_defaults = !(item.load.flags & MD_FLAG_COARSE_GRAINED);

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

        md_array(int) sequence = 0;
        md_array(int) comp_idx_type = 0;

        if (comp_count > 0) {
			md_array_resize(comp_idx_type, comp_count, temp_arena);
			MEMSET(comp_idx_type, -1, md_array_bytes(comp_idx_type));
        }

        // Process components - group by name + atom type sequence
        for (size_t i = 0; i < comp_count; ++i) {
            str_t comp_name = md_component_name(&sys.component, i);

            md_array_shrink(sequence, 0);
            md_urange_t range = md_component_atom_range(&sys.component, i);
            for (uint32_t j = range.beg; j < range.end; ++j) {
                int ai = sys.atom.type_idx[j];
                md_array_push(sequence, ai, temp_arena);
            }

            // Create combined string for hash (label + sequence of atom types)
            uint64_t hash = md_hash64_str(comp_name, md_hash64(sequence, md_array_bytes(sequence), 0));

            // Check if we already have this chain type
            DatasetItem* item = nullptr;
            for (size_t j = 0; j < md_array_size(comp_types); ++j) {
                if (comp_types[j].key == hash) {
                    item = &comp_types[j];
                    md_array_push(comp_idx_type, (int)j, temp_arena);
                    break;
                }
            }
            if (!item) {
                DatasetItem it = { .key = hash };
                snprintf(it.label, sizeof(it.label), STR_FMT, STR_ARG(comp_name));
                md_array_push(comp_idx_type, (int)md_array_size(comp_types), temp_arena);
                md_array_push(comp_types, it, arena);
                item = md_array_last(comp_types);
                md_array_push_array(item->sub_items, sequence, md_array_size(sequence), arena);
            }

            size_t comp_atom_count = md_component_atom_count(&sys.component, i);

            item->count += 1;
            item->fraction += (float)(comp_atom_count / (double)atom_count);
            md_array_push(item->indices, (int)i, arena);
        }

        // Process chains - key = residue type sequence
        if (comp_count > 0 && inst_count > 0) {
            for (size_t i = 0; i < inst_count; ++i) {
                md_array_shrink(sequence, 0);
                md_urange_t range = md_instance_component_range(&sys.instance, i);
                for (uint32_t j = range.beg; j < range.end; ++j) {
                    int res_type_idx = (int)comp_idx_type[j];
                    md_array_push(sequence, res_type_idx, temp_arena);
                }

                // Create combined string for hash (label + sequence of residue types)
                uint64_t hash = md_hash64(sequence, md_array_bytes(sequence), 0);

                // Check if we already have this chain type
                DatasetItem* item = nullptr;
                for (size_t j = 0; j < md_array_size(inst_types); ++j) {
                    if (inst_types[j].key == hash) {
                        item = &inst_types[j];
                        break;
                    }
                }
                if (!item) {
					int c_idx = (int)md_array_size(inst_types);
                    DatasetItem it = { .key = hash };
                    snprintf(it.label, sizeof(it.label), "Type %i", c_idx + 1);
                    md_array_push(inst_types, it, arena);
                    item = md_array_last(inst_types);
                    md_array_push_array(item->sub_items, sequence, md_array_size(sequence), arena);
                }

                size_t chain_atom_count = md_system_instance_atom_count(&sys, i);

                item->count += 1;
                item->fraction += (float)(chain_atom_count / (double)atom_count);
                md_array_push(item->indices, (int)i, arena);
            }
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
        delta.coarse = (type.flags[i] & MD_FLAG_COARSE_GRAINED) != (load.flags & MD_FLAG_COARSE_GRAINED);

        // use_defaults is not stored in the system, its implicit default follows the coarse grained flag
        delta.defaults = item.use_defaults != !(load.flags & MD_FLAG_COARSE_GRAINED);
        return delta;
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
            if (delta.coarse)   viamd::write_bool(state, STR_LIT("CoarseGrained"), (type.flags[i] & MD_FLAG_COARSE_GRAINED) != 0);
            if (delta.defaults) viamd::write_bool(state, STR_LIT("UseDefaults"),   item.use_defaults);
        }
    }

    // Parse a single [AtomType] section. The system is not loaded at this point, so we only buffer it.
    void deserialize_atom_type(viamd::deserialization_state_t& state) {
        // Overrides are buffered per workspace. If this section stems from another workspace than the one
        // currently buffered, then those never made it onto a system (the molecule failed to load) and are stale.
        if (!str_eq_cstr(state.filename, pending_workspace)) {
            free_pending_overrides();
            str_copy_to_char_buf(pending_workspace, sizeof(pending_workspace), state.filename);
        }

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
        pending_workspace[0] = '\0';
    }

    // Apply the buffered overrides onto the freshly loaded atom types. Must run after init_dataset_items so
    // that the load state has been snapshotted: it stays the baseline, the user's modifications sit on top.
    // A property a section does not carry was never modified, so it keeps whatever the loader assigned.
    void apply_pending_atom_type_overrides(ApplicationState& data, bool& radius_changed, bool& color_changed) {
        size_t num_overrides = md_array_size(pending_overrides);
        if (num_overrides == 0) return;
        defer { free_pending_overrides(); };

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
                    if (ovr.coarse_grained) {
                        type.flags[i] |=  MD_FLAG_COARSE_GRAINED;
                    } else {
                        type.flags[i] &= ~MD_FLAG_COARSE_GRAINED;
                    }
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
                workspace_register_window("System", &show_window);
                break;
            }
            case viamd::EventType_ViamdDeserializeBegin:
                // A workspace without [ElementDefault] or [AtomType] sections has the built in ones
                init_element_defaults();
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
                break;
            }
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t &state = *(viamd::deserialization_state_t*)e.payload;
                str_t section = viamd::section_header(state);
                if (str_eq_cstr(section, "ElementDefault")) {
                    deserialize_element_default(state);
                } else if (str_eq_cstr(section, "AtomType")) {
                    deserialize_atom_type(state);
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
                md_bitfield_or_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
            } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                md_bitfield_andnot_inplace(&state.selection.selection_mask, &state.selection.highlight_mask);
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
        ImGui::TextDisabled("%s", row.details);
        ImGui::PopID();
    }

    // The files the system was built from, in the order they compose it: the structure, the
    // trajectory it moves along, and what was loaded along that trajectory
    void draw_files(ApplicationState& data) {
        if (!ImGui::CollapsingHeader("Files", ImGuiTreeNodeFlags_DefaultOpen)) return;

        const bool has_structure = data.files.molecule[0] != '\0';
        str_t groups[64];
        const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), &data), ARRAY_SIZE(groups));
        if (!has_structure && num_groups == 0) {
            ImGui::TextDisabled("Nothing loaded. Open files with File > Open File..., or drop them on the window.");
            return;
        }

        const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp;
        if (!ImGui::BeginTable("##files", 3, flags)) return;
        ImGui::TableSetupColumn("Role", ImGuiTableColumnFlags_WidthStretch, 1.2f);
        ImGui::TableSetupColumn("File", ImGuiTableColumnFlags_WidthStretch, 3.0f);
        ImGui::TableSetupColumn("Holds", ImGuiTableColumnFlags_WidthStretch, 2.5f);
        ImGui::TableHeadersRow();

        if (has_structure) {
            FileRow row = { "Structure", str_from_cstr(data.files.molecule), "" };
            snprintf(row.details, sizeof(row.details), "%zu atoms%s", data.mold.sys.atom.count, data.files.coarse_grained ? ", coarse grained" : "");
            file_row(row);
        }

        if (data.files.trajectory[0] != '\0') {
            FileRow row = { "Trajectory", str_from_cstr(data.files.trajectory), "" };
            const size_t num_frames = run_num_frames(&data);
            const double* times = run_frame_times(&data);
            if (num_frames > 0 && times && !md_unit_is_none(run_time_unit(&data))) {
                char unit_buf[32];
                const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), run_time_unit(&data));
                snprintf(row.details, sizeof(row.details), "%zu frames, %.4g - %.4g %s", num_frames, times[0] * scl, times[num_frames - 1] * scl, unit_buf);
            } else {
                snprintf(row.details, sizeof(row.details), "%zu frames", num_frames);
            }
            file_row(row);
        }

        // Loaded along the trajectory: named by what they hold, which the group says
        const str_t run = str_from_cstr(data.mold.run);
        for (size_t g = 0; g < num_groups; ++g) {
            const str_t source = system_series_group_source(&data, groups[g]);
            if (str_empty(source)) continue;
            str_t rel = groups[g];
            if (!str_empty(run) && str_begins_with(rel, run) && rel.len > run.len + 1) {
                rel = str_substr(rel, run.len + 1);
            }
            const char* role = str_begins_with(rel, STR_LIT("edr")) ? "Energies" : "Series";
            FileRow row = { role, source, "" };
            snprintf(row.details, sizeof(row.details), "%zu series", system_series_members(nullptr, 0, &data, groups[g]));
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

        size_t atoms_without_element = 0;   // neither an element nor a coarse grained bead, and with mass
        size_t unresolved_amino = 0;        // named like amino acids, backbone not found
        size_t unresolved_nucleic = 0;
    } composition;

    int  expanded_entity = -1;
    int  expanded_instance = -1;
    bool use_short_labels = true;

    // What an entity is, from its flags and what its residues were recognised as. The flags alone
    // miss a peptide one of whose residues has an unusual backbone.
    static const char* entity_kind(md_flags_t flags, size_t amino, size_t nucleic, size_t residues) {
        if (flags & MD_FLAG_WATER)          return "Water";
        if (flags & MD_FLAG_ION)            return "Ion";
        if (flags & MD_FLAG_COARSE_GRAINED) return amino ? "Protein (CG)" : "Coarse grained";
        if ((flags & MD_FLAG_POLYPEPTIDE) || (amino && 2 * amino >= residues))     return residues > 1 ? "Protein" : "Amino acid";
        if ((flags & MD_FLAG_NUCLEOTIDE)  || (nucleic && 2 * nucleic >= residues)) return residues > 1 ? "Nucleic acid" : "Nucleotide";
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
                c.has_charge = md_attribute_extract_f32(charge, sys.atom.count, attr, md_unit_none()) == sys.atom.count;
                if (!c.has_charge) charge = nullptr;
            }
        }

        c.num_atoms = sys.atom.count;
        c.num_components = sys.component.count;
        c.num_molecules = sys.structure.count;

        // The row each atom counts toward; the last row is the atoms of no instance
        int32_t* atom_row = md_temp_alloc_array(temp, int32_t, sys.atom.count + 1);
        for (size_t a = 0; a < sys.atom.count; ++a) atom_row[a] = -1;
        size_t* amino   = md_temp_alloc_array(temp, size_t, sys.entity.count + 1);
        size_t* nucleic = md_temp_alloc_array(temp, size_t, sys.entity.count + 1);
        MEMSET(amino,   0, (sys.entity.count + 1) * sizeof(size_t));
        MEMSET(nucleic, 0, (sys.entity.count + 1) * sizeof(size_t));

        for (size_t e = 0; e < sys.entity.count; ++e) {
            EntityRow row = {};
            row.entity = (int)e;
            const md_flags_t flags = md_entity_flags(&sys.entity, e);
            row.kind = "Other";
            row.polymer = flags & (MD_FLAG_POLYMER | MD_FLAG_POLYPEPTIDE | MD_FLAG_NUCLEOTIDE);
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
                if (charge) row.charge += charge[a];
                atom_row[a] = e;
            }
            for (uint32_t k = comps.beg; k < comps.end; ++k) {
                const md_flags_t f = md_component_flags(&sys.component, k);
                if (f & MD_FLAG_AMINO_ACID) amino[e] += 1;
                if (f & MD_FLAG_NUCLEOTIDE) nucleic[e] += 1;
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
            if (charge) rest.charge += charge[a];
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
                row.kind = entity_kind(md_entity_flags(&sys.entity, row.entity), amino[row.entity], nucleic[row.entity], row.num_residues);
            }
            c.mass += row.mass;
            c.charge += row.charge;
        }

        for (size_t b = 0; b < sys.bond.count; ++b) {
            const uint32_t f = sys.bond.flags ? (uint32_t)sys.bond.flags[b] : 0u;
            if      (f & MD_BOND_FLAG_USER_DEFINED) c.bonds_user += 1;
            else if (f & MD_BOND_FLAG_TOPOLOGY)     c.bonds_topology += 1;
            else if (f & MD_BOND_FLAG_INFERRED)     c.bonds_inferred += 1;
            else                                    c.bonds_file += 1;
        }

        for (size_t a = 0; a < sys.atom.count; ++a) {
            if (atom_without_element(sys, a)) c.atoms_without_element += 1;
        }
        for (size_t i = 0; i < sys.component.count; ++i) {
            const md_flags_t f = md_component_flags(&sys.component, i);
            if ((f & MD_FLAG_AMINO_ACID) && !(f & MD_FLAG_POLYPEPTIDE)) c.unresolved_amino += 1;
            if ((f & MD_FLAG_NUCLEOTIDE) && !(f & MD_FLAG_NUCLEIC_ACID)) c.unresolved_nucleic += 1;
        }
    }

    // No element, not a bead, and not a massless virtual site (TIP4P's M carries no element by design)
    static bool atom_without_element(const md_system_t& sys, size_t a) {
        if (md_atom_atomic_number(&sys.atom, a) != 0) return false;
        if (sys.atom.flags && (sys.atom.flags[a] & MD_FLAG_COARSE_GRAINED)) return false;
        return md_atom_mass(&sys.atom, a) > 0.0f;
    }

    static void highlight_entity(ApplicationState& data, int entity) {
        const md_system_t& sys = data.mold.sys;
        md_bitfield_clear(&data.selection.highlight_mask);
        if (entity < 0) {
            // The atoms of no instance
            md_bitfield_set_range(&data.selection.highlight_mask, 0, sys.atom.count);
            for (size_t i = 0; i < sys.instance.count; ++i) {
                const md_urange_t r = md_system_instance_atom_range(&sys, i);
                md_bitfield_clear_range(&data.selection.highlight_mask, r.beg, r.end);
            }
            return;
        }
        for (size_t i = 0; i < sys.instance.count; ++i) {
            if (md_instance_entity_idx(&sys.instance, i) == entity) {
                const md_urange_t r = md_system_instance_atom_range(&sys, i);
                md_bitfield_set_range(&data.selection.highlight_mask, r.beg, r.end);
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
            const md_flags_t comp_flags = md_component_flags(&sys.component, comp_idx);
            const bool short_label = short_labels && (comp_flags & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE));
            if (short_label) {
                const uint32_t color = component_color(comp_name);
                comp_name = convert_to_short(comp_name);
                ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
                ImGui::PushStyleColor(ImGuiCol_Header, color);
                ImGui::PushStyleColor(ImGuiCol_HeaderActive, color);
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
                ImGui::PopStyleColor(2);
            }
            if (ImGui::IsItemHovered()) {
                const str_t full = md_component_name(&sys.component, comp_idx);
                ImGui::SetTooltip(STR_FMT " %d", STR_ARG(full), md_component_seq_id(&sys.component, comp_idx));
                md_bitfield_clear(&data.selection.highlight_mask);
                const md_urange_t atoms = md_system_component_atom_range(&sys, comp_idx);
                md_bitfield_set_range(&data.selection.highlight_mask, atoms.beg, atoms.end);
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
                                md_bitfield_clear(&data.selection.highlight_mask);
                                md_bitfield_set_range(&data.selection.highlight_mask, atoms.beg, atoms.end);
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
                md_bitfield_clear(&data.selection.highlight_mask);
                for (size_t a = 0; a < sys.atom.count; ++a) {
                    if (atom_without_element(sys, a)) md_bitfield_set_bit(&data.selection.highlight_mask, a);
                }
                handle_item_click(data);
            }
        }
        auto residue_warning = [&](size_t count, md_flags_t has, md_flags_t lacks, const char* what) {
            if (count == 0) return;
            ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " %zu residues are named like %s but their backbone was not found", count, what);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("They are not part of any chain, so cartoons, secondary structure and\nbackbone angles skip them. Their atom names may not follow the usual convention.");
                md_bitfield_clear(&data.selection.highlight_mask);
                for (size_t i = 0; i < sys.component.count; ++i) {
                    const md_flags_t f = md_component_flags(&sys.component, i);
                    if ((f & has) && !(f & lacks)) {
                        const md_urange_t r = md_system_component_atom_range(&sys, i);
                        md_bitfield_set_range(&data.selection.highlight_mask, r.beg, r.end);
                    }
                }
                handle_item_click(data);
            }
        };
        residue_warning(c.unresolved_amino,   MD_FLAG_AMINO_ACID, MD_FLAG_POLYPEPTIDE, "amino acids");
        residue_warning(c.unresolved_nucleic, MD_FLAG_NUCLEOTIDE, MD_FLAG_NUCLEIC_ACID, "nucleotides");

        if (c.has_charge && fabs(c.charge) >= 0.01) {
            ImGui::TextColored(warn, ICON_FA_TRIANGLE_EXCLAMATION " The system is not neutral (%+.3f e)", c.charge);
            ImGui::SetItemTooltip("With Ewald electrostatics (PME) a net charge is compensated by a uniform background,\nwhich is usually a sign of missing counter ions.");
        }
        if (c.bonds_user > 0) {
            ImGui::TextDisabled(ICON_FA_CIRCLE_INFO " %zu bonds were added by you", c.bonds_user);
            if (ImGui::IsItemHovered()) {
                md_bitfield_clear(&data.selection.highlight_mask);
                for (size_t b = 0; b < sys.bond.count; ++b) {
                    if (sys.bond.flags && (sys.bond.flags[b] & MD_BOND_FLAG_USER_DEFINED)) {
                        md_bitfield_set_bit(&data.selection.highlight_mask, sys.bond.pairs[b].idx[0]);
                        md_bitfield_set_bit(&data.selection.highlight_mask, sys.bond.pairs[b].idx[1]);
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
        char buf[256];
        const md_attribute_t* cell = md_attributes_find(&data.mold.sys.attributes, run_attribute_path(buf, sizeof(buf), &data, STR_LIT("unitcell")));
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

        char buf[256];
        traj.has_velocity = md_attributes_find(attrs, run_attribute_path(buf, sizeof(buf), &data, STR_LIT("atom/velocity"))) != nullptr;
        traj.has_force    = md_attributes_find(attrs, run_attribute_path(buf, sizeof(buf), &data, STR_LIT("atom/force"))) != nullptr;

        // Sections with an axis of their own: "<run>/trr/<name>/time"
        const str_t trr = run_attribute_path(buf, sizeof(buf), &data, STR_LIT("trr"));
        if (!str_empty(trr)) {
            str_t names[8];
            const size_t num = MIN(md_attributes_query_children(names, ARRAY_SIZE(names), attrs, trr), ARRAY_SIZE(names));
            int len = 0;
            for (size_t i = 0; i < num; ++i) {
                char tbuf[300];
                const int tl = snprintf(tbuf, sizeof(tbuf), STR_FMT "/" STR_FMT "/time", STR_ARG(trr), STR_ARG(names[i]));
                const md_attribute_t* axis = (tl > 0 && (size_t)tl < sizeof(tbuf)) ? md_attributes_find(attrs, str_t{tbuf, (size_t)tl}) : nullptr;
                if (!axis) continue;
                const size_t m = axis->format.shape[0];
                len += snprintf(traj.sections + len, sizeof(traj.sections) - len, "%s" STR_FMT " in %zu frames", len ? ", " : "", STR_ARG(names[i]), m);
                if (len >= (int)sizeof(traj.sections)) break;
            }
        }

        // The box, frame by frame
        const md_attribute_t* cell = md_attributes_find(attrs, run_attribute_path(buf, sizeof(buf), &data, STR_LIT("unitcell")));
        if (cell && cell->data && md_attribute_element_count(&cell->format) == n * 9 && n > 0) {
            md_temp_scope_t temp = md_temp_begin();
            float* m = md_temp_alloc_array(temp, float, n * 9);
            if (md_attribute_extract_f32(m, n * 9, cell, md_unit_none()) == n * 9) {
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
            char axis_buf[512];
            const int al = snprintf(axis_buf, sizeof(axis_buf), STR_FMT "/time", STR_ARG(groups[g]));
            const md_attribute_t* axis = (al > 0 && (size_t)al < sizeof(axis_buf)) ? md_attributes_find(&data.mold.sys.attributes, str_t{axis_buf, (size_t)al}) : nullptr;
            char label[256];
            system_series_group_label(label, sizeof(label), &data, groups[g]);
            if (!axis || axis->format.shape[0] < 1) {
                ImGui::TextDisabled("%s: one value per frame", label);
                continue;
            }
            const size_t m = axis->format.shape[0];
            double ends[2] = {0, 0};
            const md_attribute_slice_t s0 = md_attribute_slice_1(0), s1 = md_attribute_slice_1((uint32_t)(m - 1));
            md_attribute_extract_slice_f64(&ends[0], 1, axis, &s0, run_time_unit(&data));
            md_attribute_extract_slice_f64(&ends[1], 1, axis, &s1, run_time_unit(&data));
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
                        md_attribute_extract_slice_f64(&t[0], 1, axis, &first, md_unit_none()) == 1 &&
                        md_attribute_extract_slice_f64(&t[1], 1, axis, &last,  md_unit_none()) == 1)
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

    // The atom types and what they are drawn with: an editor, kept out of the overview
    void draw_atom_types(ApplicationState& data) {
			size_t num_atom_types = md_system_atom_type_count(&data.mold.sys);

            if (num_atom_types) {
                const float min_mass = 1.0f;
                const float max_mass = 500.0f;

                const float min_radius = 0.1f;
                const float max_radius = 20.0f;

                bool radius_changed = false;
                bool color_changed  = false;
                bool mass_changed   = false;

                if (ImGui::CollapsingHeader("Atom Types", ImGuiTreeNodeFlags_DefaultOpen)) {
                    ImGui::Indent();
                    for (size_t i = 0; i < num_atom_types; ++i) {
                        DatasetItem& item = atom_types[i];
                        if (i == 0 && item.count == 0) {
                            // Skip sentinel "unknown" atom type if unused
                            continue;
                        }

                        ImGui::PushID((int)i);
                        defer{ ImGui::PopID(); };
                        char buf_num[32];
						snprintf(buf_num, sizeof(buf_num), "%d (%.2f%%)", item.count, item.fraction * 100.0f);
                        char buf_tot[256];
                        snprintf(buf_tot, sizeof(buf_tot), "%-4s %10s", item.label, buf_num);
                        bool expand = ImGui::CollapsingHeader(buf_tot);
                        if (ImGui::IsItemHovered()) {
                            md_bitfield_clear(&data.selection.highlight_mask);
                            md_bitfield_set_indices_u32(&data.selection.highlight_mask, (uint32_t*)item.indices, md_array_size(item.indices));
                            if (ImGui::IsItemClicked()) {
                                handle_item_click(data);
                            }
                        }
                        if (expand) {
                            bool coarse_grained = data.mold.sys.atom.type.flags[i] & MD_FLAG_COARSE_GRAINED;
                            if (ImGui::Checkbox("Coarse Grained", &coarse_grained)) {
                                if (coarse_grained) {
                                    data.mold.sys.atom.type.flags[i] |=  MD_FLAG_COARSE_GRAINED;
                                } else {
                                    data.mold.sys.atom.type.flags[i] &= ~MD_FLAG_COARSE_GRAINED;
                                }
                            }
                            if (!(data.mold.sys.atom.type.flags[i] & MD_FLAG_COARSE_GRAINED)) {
								str_t symbol = md_atomic_number_symbol((md_atomic_number_t)data.mold.sys.atom.type.z[i]);
                                if (ImGui::Checkbox("Use element defaults", &item.use_defaults)) {
                                    if (item.use_defaults) {
										// If the user enables the use_defaults flag then we should set the values back to the element defaults
                                        md_element_t elem = data.mold.sys.atom.type.z[i];
										data.mold.sys.atom.type.radius[i] = md_util_element_vdw_radius(elem);
										data.mold.sys.atom.type.mass[i]   = md_util_element_atomic_mass(elem);
										data.mold.sys.atom.type.color[i]  = md_util_element_cpk_color(elem);
										radius_changed = true;
                                        color_changed = true;
										mass_changed = true;
                                    }
                                }
                                if (item.use_defaults) {
                                    ImGui::SameLine();
                                    int z = data.mold.sys.atom.type.z[i];
                                    vec4_t color = element_defaults[z].color;
                                    if (element_button(str_ptr(symbol), color)) {
                                        ImGui::OpenPopup("Element Popup");
                                    }
                                }
                                if (ImGui::BeginPopup("Element Popup")) {
                                    PeriodicTableResult table_res = periodic_table_widget(element_defaults);
                                    if (table_res.clicked) {
                                        int z = table_res.z;
                                        data.mold.sys.atom.type.z[i] = (md_atomic_number_t)table_res.z;
                                        data.mold.sys.atom.type.color[i] = u32_from_vec4(element_defaults[z].color);
                                        data.mold.sys.atom.type.radius[i] = element_defaults[z].radius;
                                        data.mold.sys.atom.type.mass[i] = element_defaults[z].mass;

										radius_changed = true;
										color_changed = true;
										mass_changed = true;
                                        ImGui::CloseCurrentPopup();
                                    }
									ImGui::EndPopup();
                                }
                            } else {
								item.use_defaults = false; // Coarse grained types always have custom properties
                            }

                            if (item.use_defaults) {
                                ImGui::PushDisabled();
                            }
                            float* radius = &data.mold.sys.atom.type.radius[i];
                            radius_changed |= ImGui::SliderFloat("Radius", radius, min_radius, max_radius);

                            float* mass = &data.mold.sys.atom.type.mass[i];
                            mass_changed |= ImGui::SliderFloat("Mass", mass, min_mass, max_mass);

                            ImVec4 color = ImColor(data.mold.sys.atom.type.color[i]);
                            if (ImGui::ColorEdit4("Color", &color.x)) {
                                data.mold.sys.atom.type.color[i] = ImColor(color);
                                color_changed = true;
                            }

                            if (item.use_defaults) {
                                ImGui::PopDisabled();
							}
                        }
                    }
                    ImGui::Unindent();
                }

                // Keep track of what elements are used in the system
				uint64_t elem_mask[2] = { 0 };

                for (size_t i = 0; i < num_atom_types; ++i) {
                    int z = data.mold.sys.atom.type.z[i];
                    if (atom_types[i].use_defaults && atom_types[i].count > 0) {
						elem_mask[z / 64] |= (1ULL << (z % 64));
                    }
                }

                if ((elem_mask[0] || elem_mask[1]) && ImGui::CollapsingHeader("Element Defaults")) {
                    ImGui::Indent();

                    static int z = -1;
					PeriodicTableResult table_res = periodic_table_widget(element_defaults, elem_mask);
					if (table_res.hovered) {
                        md_bitfield_clear(&data.selection.highlight_mask);
                        for (size_t i = 0; i < num_atom_types; ++i) {
							const DatasetItem& item = atom_types[i];
                            if (data.mold.sys.atom.type.z[i] == table_res.z) {
                                md_bitfield_set_indices_u32(&data.selection.highlight_mask, (uint32_t*)item.indices, md_array_size(item.indices));
                            }
                        }
                        handle_item_click(data);

                        if (table_res.clicked && !ImGui::IsKeyDown(ImGuiMod_Shift)) {
                            ImGui::OpenPopup("Element Popup");
							z = table_res.z;
                        }
                    }

                    if (ImGui::BeginPopup("Element Popup")) {
                        str_t sym = md_atomic_number_symbol((md_atomic_number_t)z);
                        str_t name = md_atomic_number_name((md_atomic_number_t)z);
                        char buf[64];
                        snprintf(buf, sizeof(buf), "%d: %s (%s)", z, str_ptr(name), str_ptr(sym));
                        ImGui::Text("%s", buf);
                        ImGui::Separator();

                        ElementDefault& elem_def = element_defaults[z];

                        if (ImGui::ColorEdit3("Color", elem_def.color.elem)) {
                            // Iterate and set color for all atom types that use this element and have use_defaults = true
                            for (size_t i = 0; i < num_atom_types; ++i) {
                                if (data.mold.sys.atom.type.z[i] == z && atom_types[i].use_defaults) {
                                    data.mold.sys.atom.type.color[i] = u32_from_vec4(elem_def.color);
                                }
                            }
                            color_changed = true;
                        }
                        if (ImGui::InputFloat("Van der Waals Radius", &elem_def.radius)) {
                            // Iterate and set radius for all atom types that use this element and have use_defaults = true
                            for (size_t i = 0; i < num_atom_types; ++i) {
                                if (data.mold.sys.atom.type.z[i] == z && atom_types[i].use_defaults) {
                                    data.mold.sys.atom.type.radius[i] = elem_def.radius;
                                }
                            }
                            radius_changed = true;
                        }
                        if (ImGui::InputFloat("Atomic Mass", &elem_def.mass)) {
                            // Iterate and set mass for all atom types that use this element and have use_defaults = true
                            for (size_t i = 0; i < num_atom_types; ++i) {
                                if (data.mold.sys.atom.type.z[i] == z && atom_types[i].use_defaults) {
                                    data.mold.sys.atom.type.mass[i] = elem_def.mass;
                                }
                            }
                            mass_changed = true;
                        }

                        ImGui::EndPopup();
                    }
                    ImGui::Unindent();
                }

                if (radius_changed) {
                    data.mold.dirty_gpu_buffers |= MolBit_DirtyRadius;
                }

                if (color_changed) {
                    // @NOTE: Only the color within representations needs to be updated, not the filter.
                    flag_all_representations_as_dirty(&data);
                }
                (void)mass_changed; // Currently mass is not used for rendering, but we track changes in case it's used for other purposes in the future
            }


            // Draw the three sections
            //draw_dataset_section("Chain Types",     inst_types,      md_array_size(inst_types),   0);
            //draw_dataset_section("Residue Types",   comp_types,    md_array_size(comp_types), 1);  
            //draw_dataset_section("Atom Types",      atom_types,       md_array_size(atom_types),    2);

            // Atom Element Mappings section (keep existing functionality)
            const size_t num_mappings = md_array_size(atom_element_remappings);
            if (num_mappings) {
                if (ImGui::CollapsingHeader("Atom Element Mappings")) {
                    for (size_t i = 0; i < num_mappings; ++i) {
                        const auto& mapping = atom_element_remappings[i];
                        ImGui::Text("%s -> %s (%s)", mapping.lbl, md_util_element_name(mapping.elem).ptr, md_util_element_symbol(mapping.elem).ptr);
                    }
                }
            }
    }

    void draw(ApplicationState& data) {
        if (!show_window) return;

        ImGui::SetNextWindowSize(ImVec2(560, 700), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("System", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::IsWindowHovered()) {
                md_bitfield_clear(&data.selection.highlight_mask);
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
                    draw_atom_types(data);
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
