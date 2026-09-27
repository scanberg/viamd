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

#include <imgui.h>
#include <imgui_widgets.h>
#include <implot_widgets.h>

#include <string>
#include <algorithm>

namespace dataset {

// ## Residue codes

// The class a residue is colored by in a sequence. Amino acids follow the Clustal X grouping,
// nucleotides are colored by their base.
enum ResidueClass : uint8_t {
    ResidueClass_Unknown = 0,
    ResidueClass_Hydrophobic,
    ResidueClass_Cysteine,
    ResidueClass_Positive,
    ResidueClass_Negative,
    ResidueClass_Polar,
    ResidueClass_Glycine,
    ResidueClass_Proline,
    ResidueClass_Aromatic,
    ResidueClass_Adenine,
    ResidueClass_Cytosine,
    ResidueClass_Guanine,
    ResidueClass_Thymine,
    ResidueClass_Uracil,
    ResidueClass_Count
};

static inline bool residue_class_is_nucleotide(ResidueClass cls) {
    return cls >= ResidueClass_Adenine && cls < ResidueClass_Count;
}

static const char* residue_class_label(ResidueClass cls) {
    switch (cls) {
    case ResidueClass_Hydrophobic: return "Hydrophobic";
    case ResidueClass_Cysteine:    return "Cysteine";
    case ResidueClass_Positive:    return "Positive";
    case ResidueClass_Negative:    return "Negative";
    case ResidueClass_Polar:       return "Polar";
    case ResidueClass_Glycine:     return "Glycine";
    case ResidueClass_Proline:     return "Proline";
    case ResidueClass_Aromatic:    return "Aromatic";
    case ResidueClass_Adenine:     return "Adenine";
    case ResidueClass_Cytosine:    return "Cytosine";
    case ResidueClass_Guanine:     return "Guanine";
    case ResidueClass_Thymine:     return "Thymine";
    case ResidueClass_Uracil:      return "Uracil";
    default:                       return "Unknown";
    }
}

// Pastel fills which dark text reads well on, in light and dark themes alike
static uint32_t residue_class_color(ResidueClass cls) {
    switch (cls) {
    case ResidueClass_Hydrophobic: return IM_COL32(160, 190, 245, 255);
    case ResidueClass_Cysteine:    return IM_COL32(245, 170, 170, 255);
    case ResidueClass_Positive:    return IM_COL32(245, 125, 115, 255);
    case ResidueClass_Negative:    return IM_COL32(215, 150, 215, 255);
    case ResidueClass_Polar:       return IM_COL32(140, 215, 140, 255);
    case ResidueClass_Glycine:     return IM_COL32(245, 185, 130, 255);
    case ResidueClass_Proline:     return IM_COL32(225, 225, 110, 255);
    case ResidueClass_Aromatic:    return IM_COL32(120, 205, 205, 255);
    case ResidueClass_Adenine:     return IM_COL32(252, 105, 122, 255);
    case ResidueClass_Cytosine:    return IM_COL32(248, 238,  92, 255);
    case ResidueClass_Guanine:     return IM_COL32(213, 179, 239, 255);
    case ResidueClass_Thymine:     return IM_COL32(159, 243, 160, 255);
    case ResidueClass_Uracil:      return IM_COL32(254, 157,  45, 255);
    default:                       return IM_COL32(190, 190, 190, 255);
    }
}

struct ResidueCode {
    const char*  name;
    char         code;
    ResidueClass cls;
};

#define RC(name, code, cls) {name, code, ResidueClass_##cls}

// Residue names as they appear in structure files and force fields, mapped to their single letter code.
// Besides the standard names this covers the protonation and disulfide variants of the common MD force
// fields (AMBER, CHARMM, GROMOS/OPLS as named in GROMACS), which is what simulated systems usually carry.
static const ResidueCode residue_codes[] = {
    RC("ALA",'A',Hydrophobic),
    RC("ARG",'R',Positive), RC("ARN",'R',Positive),
    RC("ASN",'N',Polar),
    RC("ASP",'D',Negative), RC("ASH",'D',Negative), RC("ASPH",'D',Negative), RC("ASPP",'D',Negative),
    RC("CYS",'C',Cysteine), RC("CYX",'C',Cysteine), RC("CYM",'C',Cysteine), RC("CYN",'C',Cysteine), RC("CYS2",'C',Cysteine), RC("CYSH",'C',Cysteine),
    RC("GLN",'Q',Polar),
    RC("GLU",'E',Negative), RC("GLH",'E',Negative), RC("GLUH",'E',Negative), RC("GLUP",'E',Negative),
    RC("GLY",'G',Glycine),
    RC("HIS",'H',Aromatic), RC("HID",'H',Aromatic), RC("HIE",'H',Aromatic), RC("HIP",'H',Aromatic),
    RC("HSD",'H',Aromatic), RC("HSE",'H',Aromatic), RC("HSP",'H',Aromatic),
    RC("HISD",'H',Aromatic), RC("HISE",'H',Aromatic), RC("HISH",'H',Aromatic), RC("HISA",'H',Aromatic), RC("HISB",'H',Aromatic), RC("HIS1",'H',Aromatic), RC("HIS2",'H',Aromatic),
    RC("ILE",'I',Hydrophobic),
    RC("LEU",'L',Hydrophobic),
    RC("LYS",'K',Positive), RC("LYN",'K',Positive), RC("LSN",'K',Positive), RC("LYSH",'K',Positive), RC("LYP",'K',Positive),
    RC("MET",'M',Hydrophobic), RC("MSE",'M',Hydrophobic),
    RC("PHE",'F',Hydrophobic),
    RC("PRO",'P',Proline),
    RC("SER",'S',Polar),
    RC("THR",'T',Polar),
    RC("TRP",'W',Hydrophobic),
    RC("TYR",'Y',Aromatic),
    RC("VAL",'V',Hydrophobic),
    RC("SEC",'U',Cysteine),
    RC("PYL",'O',Positive),

    // DNA
    RC("DA",'A',Adenine), RC("DC",'C',Cytosine), RC("DG",'G',Guanine), RC("DT",'T',Thymine), RC("DU",'U',Uracil),
    // RNA
    RC("A",'A',Adenine),  RC("C",'C',Cytosine),  RC("G",'G',Guanine),  RC("T",'T',Thymine),  RC("U",'U',Uracil),
    RC("RA",'A',Adenine), RC("RC",'C',Cytosine), RC("RG",'G',Guanine), RC("RU",'U',Uracil),
    // CHARMM
    RC("ADE",'A',Adenine), RC("CYT",'C',Cytosine), RC("GUA",'G',Guanine), RC("THY",'T',Thymine), RC("URA",'U',Uracil),
};

#undef RC

static int find_residue_code_exact(str_t name) {
    for (int i = 0; i < (int)ARRAY_SIZE(residue_codes); ++i) {
        if (str_eq_cstr(name, residue_codes[i].name)) return i;
    }
    return -1;
}

// Index into residue_codes for a component, -1 if it has no single letter code.
// The terminal variants are only accepted when the component flags agree on what the component is,
// so that a ligand which happens to be named like one (e.g. CN) is not taken for a terminal nucleotide.
static int find_residue_code(str_t name, md_flags_t comp_flags) {
    name = str_trim(name);
    int idx = find_residue_code_exact(name);
    if (idx != -1) return idx;

    // AMBER terminal amino acids: NALA, CALA, ...
    if ((comp_flags & MD_FLAG_AMINO_ACID) && name.len == 4 && (name.ptr[0] == 'N' || name.ptr[0] == 'C')) {
        idx = find_residue_code_exact(str_substr(name, 1));
        if (idx != -1 && !residue_class_is_nucleotide(residue_codes[idx].cls)) return idx;
    }

    // AMBER terminal nucleotides: DA5, DA3, DAN, RA5, A5, ...
    if ((comp_flags & (MD_FLAG_NUCLEOTIDE | MD_FLAG_NUCLEIC_ACID)) && name.len >= 2) {
        const char last = name.ptr[name.len - 1];
        if (last == '5' || last == '3' || last == 'N') {
            idx = find_residue_code_exact(str_substr(name, 0, name.len - 1));
            if (idx != -1 && residue_class_is_nucleotide(residue_codes[idx].cls)) return idx;
        }
    }

    return -1;
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

// What an entity is, as far as presenting it goes
enum EntityKind : uint8_t {
    EntityKind_Molecule = 0,    // Anything not covered below, e.g. lipids or solvents other than water
    EntityKind_Protein,
    EntityKind_NucleicAcid,
    EntityKind_Polymer,         // Connected chain of components which are neither amino acids nor nucleotides
    EntityKind_Water,
    EntityKind_Ion,
    EntityKind_Ligand,          // Marked as hetero by the file
};

static const char* entity_kind_label(EntityKind kind) {
    switch (kind) {
    case EntityKind_Protein:     return "Protein";
    case EntityKind_NucleicAcid: return "Nucleic acid";
    case EntityKind_Polymer:     return "Polymer";
    case EntityKind_Water:       return "Water";
    case EntityKind_Ion:         return "Ion";
    case EntityKind_Ligand:      return "Ligand";
    default:                     return "Molecule";
    }
}

// What the components of an entity's instances are called
static const char* entity_kind_comp_noun(EntityKind kind, size_t count, bool capitalized = false) {
    const bool one = count == 1;
    switch (kind) {
    case EntityKind_Protein:
    case EntityKind_NucleicAcid: return capitalized ? (one ? "Residue"  : "Residues")  : (one ? "residue"  : "residues");
    case EntityKind_Water:       return capitalized ? (one ? "Molecule" : "Molecules") : (one ? "molecule" : "molecules");
    case EntityKind_Ion:         return capitalized ? (one ? "Ion"      : "Ions")      : (one ? "ion"      : "ions");
    default:                     return capitalized ? (one ? "Component": "Components"): (one ? "component": "components");
    }
}

// Instances of these are chains, shown as a sequence
static inline bool entity_kind_is_sequence(EntityKind kind) {
    return kind == EntityKind_Protein || kind == EntityKind_NucleicAcid || kind == EntityKind_Polymer;
}

// The entity flags only carry the chain level flags (polypeptide, nucleic acid), so the flags of the
// components are combined into them before classifying: a derived entity of amino acids which lack
// the backbone bonds to be recognized as a polypeptide is still a protein.
static EntityKind classify_entity(md_flags_t flags) {
    if (flags & MD_FLAG_WATER)                                  return EntityKind_Water;
    if (flags & MD_FLAG_ION)                                    return EntityKind_Ion;
    if (flags & (MD_FLAG_POLYPEPTIDE  | MD_FLAG_AMINO_ACID))    return EntityKind_Protein;
    if (flags & (MD_FLAG_NUCLEIC_ACID | MD_FLAG_NUCLEOTIDE))    return EntityKind_NucleicAcid;
    if (flags & MD_FLAG_POLYMER)                                return EntityKind_Polymer;
    if (flags & MD_FLAG_HETERO)                                 return EntityKind_Ligand;
    return EntityKind_Molecule;
}

// Per entity summary, computed once when the system is initialized
struct EntityInfo {
    md_array(int) instances = 0;            // Indices of the instances of the entity
    size_t   num_comps = 0;                 // Summed over all instances
    size_t   num_atoms = 0;                 // Summed over all instances
    uint32_t min_comps = 0, max_comps = 0;  // Components per instance
    uint32_t min_atoms = 0, max_atoms = 0;  // Atoms per instance
    md_flags_t flags = MD_FLAG_NONE;        // Entity flags combined with the flags of all its components
    EntityKind kind  = EntityKind_Molecule;
    bool derived = false;                   // Inferred upon load rather than defined by the file
    char name[64] = "";                     // Description to show
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

    md_array(EntityInfo) entities = 0;  // One per entity of the system
    md_array(int16_t) comp_code = 0;    // Per component: index into residue_codes, -1 if it has no single letter code

    // Settings menu
    bool use_single_letter_codes = true;    // Show the residues of proteins and nucleic acids by their single letter code
    bool color_residues          = true;    // Color the residues of sequences by their class

    char atom_type_filter[64] = "";
    md_array(int) atom_type_order = 0;      // Display order of the atom types as sorted by the table, allocated in `arena`
    bool atom_type_resort     = false;      // The order is to be sorted again, e.g. after the atom types changed
    int  atom_type_selected   = -1;         // Atom type selected in the table, edited below it
    int  element_popup_z = -1;              // Element whose defaults are being edited

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
        atom_type_order = 0;
        atom_type_selected = -1;
        entities   = 0;
        comp_code  = 0;
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

        // Single letter codes of the components, looked up once rather than every frame
        if (comp_count > 0) {
            md_array_resize(comp_code, comp_count, arena);
            for (size_t i = 0; i < comp_count; ++i) {
                comp_code[i] = (int16_t)find_residue_code(md_component_name(&sys.component, i), md_component_flags(&sys.component, i));
            }
        }

        // Entities: which instances belong to each, and what those add up to
        const size_t ent_count = md_system_entity_count(&sys);
        if (ent_count > 0) {
            md_array_resize(entities, ent_count, arena);
            for (size_t e = 0; e < ent_count; ++e) {
                entities[e] = EntityInfo{};
                entities[e].flags = md_entity_flags(&sys.entity, e);
            }

            for (size_t i = 0; i < inst_count; ++i) {
                const int e = md_instance_entity_idx(&sys.instance, i);
                if (e < 0 || (size_t)e >= ent_count) continue;
                EntityInfo& info = entities[e];

                const md_urange_t comp_range = md_instance_component_range(&sys.instance, i);
                const uint32_t num_comps = comp_range.end - comp_range.beg;
                const uint32_t num_atoms = (uint32_t)md_system_instance_atom_count(&sys, i);
                if (md_array_size(info.instances) == 0) {
                    info.min_comps = info.max_comps = num_comps;
                    info.min_atoms = info.max_atoms = num_atoms;
                } else {
                    info.min_comps = MIN(info.min_comps, num_comps);
                    info.max_comps = MAX(info.max_comps, num_comps);
                    info.min_atoms = MIN(info.min_atoms, num_atoms);
                    info.max_atoms = MAX(info.max_atoms, num_atoms);
                }
                info.num_comps += num_comps;
                info.num_atoms += num_atoms;
                for (uint32_t c = comp_range.beg; c < comp_range.end; ++c) {
                    info.flags |= md_component_flags(&sys.component, c);
                }
                md_array_push(info.instances, (int)i, arena);
            }

            for (size_t e = 0; e < ent_count; ++e) {
                EntityInfo& info = entities[e];
                info.kind = classify_entity(info.flags);

                // Entities the file does not define are derived by mdlib (md_util_system_infer_entity_and_instance)
                str_t desc = str_trim(md_entity_description(&sys.entity, e));
                info.derived = (md_entity_flags(&sys.entity, e) & MD_FLAG_DERIVED) != 0;
                if (info.derived && (entity_kind_is_sequence(info.kind) || info.kind == EntityKind_Water)) {
                    // The generated description only restates the kind (polypeptide, water, ...), which is shown anyway
                    desc = STR_LIT("");
                }
                if (str_empty(desc)) {
                    desc = str_from_cstr(entity_kind_label(info.kind));
                }
                str_copy_to_char_buf(info.name, sizeof(info.name), desc);
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

    // ## Contents: what the system holds

    void draw_contents(ApplicationState& data) {
        if (!ImGui::CollapsingHeader("Contents", ImGuiTreeNodeFlags_DefaultOpen)) return;
        const md_system_t& sys = data.mold.sys;

        const ImGuiTableFlags flags = ImGuiTableFlags_SizingFixedFit;
        if (ImGui::BeginTable("##contents", 2, flags)) {
            auto count_row = [](const char* label, size_t count) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted(label);
                ImGui::TableSetColumnIndex(1);
                ImGui::Text("%zu", count);
            };
            count_row("Entities",   sys.entity.count);
            count_row("Instances",  sys.instance.count);
            count_row("Components", sys.component.count);
            count_row("Atoms",      sys.atom.count);
            count_row("Bonds",      sys.bond.count);
            ImGui::EndTable();
        }

        const md_unitcell_t& cell = data.mold.state.unitcell;
        if (cell.flags) {
            char unit_buf[32];
            const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), md_unit_angstrom());
            const bool ortho = cell.flags & MD_UNITCELL_ORTHO;
            const bool tricl = cell.flags & MD_UNITCELL_TRICLINIC;
            const bool px = cell.flags & MD_UNITCELL_PBC_X, py = cell.flags & MD_UNITCELL_PBC_Y, pz = cell.flags & MD_UNITCELL_PBC_Z;
            ImGui::Text("Box: %s, periodic in %s%s%s%s", ortho ? "orthorhombic" : tricl ? "triclinic" : "-",
                px ? "x" : "", py ? "y" : "", pz ? "z" : "", (px || py || pz) ? "" : "no direction");
            ImGui::Indent();
            if (ortho) {
                ImGui::Text("%.4g x %.4g x %.4g %s", cell.x * scl, cell.y * scl, cell.z * scl, unit_buf);
            } else if (tricl) {
                ImGui::Text("x %.4g, y %.4g, z %.4g %s", cell.x * scl, cell.y * scl, cell.z * scl, unit_buf);
                ImGui::Text("xy %.4g, xz %.4g, yz %.4g %s", cell.xy * scl, cell.xz * scl, cell.yz * scl, unit_buf);
            }
            ImGui::Unindent();
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

    static const char* fmt_count_range(char* buf, size_t cap, size_t lo, size_t hi) {
        char a[32], b[32];
        if (lo == hi) {
            snprintf(buf, cap, "%s", fmt_count(a, sizeof(a), lo));
        } else {
            snprintf(buf, cap, "%s - %s", fmt_count(a, sizeof(a), lo), fmt_count(b, sizeof(b), hi));
        }
        return buf;
    }

    static const char* fmt_mass(char* buf, size_t cap, double mass) {
        if (mass >= 1.0e6) {
            snprintf(buf, cap, "%.3g MDa", mass * 1.0e-6);
        } else if (mass >= 1.0e3) {
            snprintf(buf, cap, "%.2f kDa", mass * 1.0e-3);
        } else {
            snprintf(buf, cap, "%.2f Da", mass);
        }
        return buf;
    }

    static double atom_range_mass(const md_system_t& sys, md_urange_t atoms) {
        double mass = 0;
        for (uint32_t i = atoms.beg; i < atoms.end; ++i) {
            mass += md_atom_mass(&sys.atom, i);
        }
        return mass;
    }

    static int element_symbol_cmp(md_atomic_number_t a, md_atomic_number_t b) {
        char sa[8], sb[8];
        str_copy_to_char_buf(sa, sizeof(sa), md_atomic_number_symbol(a));
        str_copy_to_char_buf(sb, sizeof(sb), md_atomic_number_symbol(b));
        return strcmp(sa, sb);
    }

    // Chemical formula in Hill order: C and H first when there is carbon, the rest alphabetically.
    // Empty if any atom lacks an element, such as coarse grained beads.
    static void hill_formula(char* buf, size_t cap, const md_system_t& sys, md_urange_t atoms) {
        buf[0] = '\0';
        uint32_t count[MD_Z_Count] = {};
        for (uint32_t i = atoms.beg; i < atoms.end; ++i) {
            const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, i);
            if (z == 0 || z >= MD_Z_Count) return;
            count[z] += 1;
        }

        md_atomic_number_t order[MD_Z_Count];
        int num = 0;
        const bool organic = count[MD_Z_C] > 0;
        if (organic) {
            order[num++] = MD_Z_C;
            if (count[MD_Z_H]) order[num++] = MD_Z_H;
        }
        const int first_sorted = num;
        for (int z = 1; z < MD_Z_Count; ++z) {
            if (!count[z]) continue;
            if (organic && (z == MD_Z_C || z == MD_Z_H)) continue;
            order[num++] = (md_atomic_number_t)z;
        }
        for (int i = first_sorted + 1; i < num; ++i) {
            const md_atomic_number_t key = order[i];
            int j = i - 1;
            while (j >= first_sorted && element_symbol_cmp(order[j], key) > 0) {
                order[j + 1] = order[j];
                --j;
            }
            order[j + 1] = key;
        }

        size_t len = 0;
        for (int i = 0; i < num; ++i) {
            const md_atomic_number_t z = order[i];
            const str_t sym = md_atomic_number_symbol(z);
            const int n = count[z] > 1 ? snprintf(buf + len, cap - len, STR_FMT "%u", STR_ARG(sym), count[z])
                                       : snprintf(buf + len, cap - len, STR_FMT, STR_ARG(sym));
            if (n < 0 || (size_t)n >= cap - len) break; // Truncated
            len += (size_t)n;
        }
    }

    // The single letter sequence of an instance, X for components without a code, cut short with ...
    void sequence_string(char* buf, size_t cap, const md_system_t& sys, size_t inst_idx) const {
        const md_urange_t range = md_system_instance_comp_range(&sys, inst_idx);
        size_t len = 0;
        for (uint32_t c = range.beg; c < range.end; ++c) {
            if (len + 4 >= cap) {
                buf[len++] = '.'; buf[len++] = '.'; buf[len++] = '.';
                break;
            }
            const int code = c < md_array_size(comp_code) ? comp_code[c] : -1;
            buf[len++] = code >= 0 ? residue_codes[code].code : 'X';
        }
        buf[len] = '\0';
    }

    // Draws disabled text flush against right_x on the line of the previous item, when there is room for it
    static void right_aligned_disabled(float right_x, const char* text) {
        const float w = ImGui::CalcTextSize(text).x;
        ImGui::SameLine();
        const ImVec2 pos = ImGui::GetCursorScreenPos();
        if (right_x - w > pos.x) {
            ImGui::SetCursorScreenPos(ImVec2(right_x - w, pos.y));
            ImGui::TextDisabled("%s", text);
        } else {
            ImGui::NewLine();
        }
    }

    // Key value row for the tooltip tables
    static void info_row(const char* key) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::TextDisabled("%s", key);
        ImGui::TableSetColumnIndex(1);
    }

    // ## Selection targets: entities, instances and components

    enum Target {
        Target_Entity,
        Target_Instance,
        Target_Component,
    };

    void target_mask(md_bitfield_t* bf, const md_system_t& sys, Target target, size_t idx) const {
        switch (target) {
        case Target_Entity:
            if (idx < md_array_size(entities)) {
                const EntityInfo& info = entities[idx];
                for (size_t k = 0; k < md_array_size(info.instances); ++k) {
                    const md_urange_t range = md_system_instance_atom_range(&sys, info.instances[k]);
                    md_bitfield_set_range(bf, range.beg, range.end);
                }
            }
            break;
        case Target_Instance: {
            const md_urange_t range = md_system_instance_atom_range(&sys, idx);
            md_bitfield_set_range(bf, range.beg, range.end);
            break;
        }
        case Target_Component: {
            const md_urange_t range = md_system_component_atom_range(&sys, idx);
            md_bitfield_set_range(bf, range.beg, range.end);
            break;
        }
        }
    }

    // Hovering the target highlights it, Shift + click (de)selects it and a right click opens a menu to do the same.
    // Called right after the item which represents the target, with whether it is hovered.
    void target_interaction(ApplicationState& data, Target target, size_t idx, bool hovered) {
        if (hovered) {
            md_bitfield_clear(&data.selection.highlight_mask);
            target_mask(&data.selection.highlight_mask, data.mold.sys, target, idx);
            handle_item_click(data);
            // Shift + right click is taken by deselect
            if (!ImGui::IsKeyDown(ImGuiMod_Shift) && ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                ImGui::OpenPopup("##target_menu");
            }
        }
        if (ImGui::BeginPopup("##target_menu")) {
            md_bitfield_t* sel = &data.selection.selection_mask;
            if (ImGui::MenuItem("Select")) {
                md_bitfield_clear(sel);
                target_mask(sel, data.mold.sys, target, idx);
            }
            if (ImGui::MenuItem("Add to selection")) {
                target_mask(sel, data.mold.sys, target, idx);
            }
            if (ImGui::MenuItem("Remove from selection")) {
                md_temp_scope_t temp = md_temp_begin();
                md_bitfield_t mask = {};
                md_bitfield_init(&mask, md_temp_allocator(temp));
                target_mask(&mask, data.mold.sys, target, idx);
                md_bitfield_andnot_inplace(sel, &mask);
                md_temp_end(temp);
            }
            ImGui::EndPopup();
        }
    }

    // ## Tooltips

    void entity_tooltip(const ApplicationState& data, size_t e) const {
        const md_system_t& sys = data.mold.sys;
        const EntityInfo& info = entities[e];
        const size_t num_inst  = md_array_size(info.instances);
        const str_t  id        = md_entity_id(&sys.entity, e);
        char buf[128], num[32];

        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);

        ImGui::Text("Entity " STR_FMT ": %s", STR_ARG(id), info.name);
        if (!str_eq_cstr_ignore_case(str_from_cstr(info.name), entity_kind_label(info.kind))) {
            ImGui::SameLine();
            ImGui::TextDisabled("%s", entity_kind_label(info.kind));
        }

        ImGui::Separator();
        if (info.derived) {
            // How md_util_system_infer_entity_and_instance arrives at the instances of this kind of entity
            bool has_chain_ids = false;
            for (size_t k = 0; k < num_inst; ++k) {
                has_chain_ids |= !str_empty(md_instance_auth_id(&sys.instance, info.instances[k]));
            }
            const char* grouping = "";
            if (entity_kind_is_sequence(info.kind)) {
                grouping = has_chain_ids ? "Each instance is a run of consecutive residues sharing a chain ID of the file."
                                         : "Each instance is a run of consecutive residues joined by bonds.";
            } else if (info.kind == EntityKind_Water || info.kind == EntityKind_Ion) {
                grouping = has_chain_ids ? "Consecutive molecules of the same name and chain ID are grouped into one instance."
                                         : "Consecutive molecules of the same name are grouped into one instance.";
            } else {
                grouping = has_chain_ids ? "Each instance is a run of consecutive components sharing a chain ID of the file."
                                         : "Each instance is a molecule: consecutive components joined by bonds, or sharing a residue number.";
            }

            ImGui::TextUnformatted("Derived on load");
            ImGui::TextDisabled("The file does not define entities, so they were derived from the components. %s "
                                "Instances made up of the same components form one entity. Entity and instance IDs are generated.", grouping);
        } else {
            ImGui::TextUnformatted("Defined by the file");
        }
        ImGui::Spacing();

        if (ImGui::BeginTable("##entity_info", 2, ImGuiTableFlags_SizingFixedFit)) {
            // IDs of the first few instances, followed by the author (chain) IDs where the file provided those
            {
                const size_t max_ids = 8;
                size_t len = 0;
                bool has_auth = false;
                buf[0] = '\0';
                for (size_t k = 0; k < MIN(num_inst, max_ids); ++k) {
                    const str_t inst_id = md_instance_id(&sys.instance, info.instances[k]);
                    const int n = snprintf(buf + len, sizeof(buf) - len, "%s" STR_FMT, k ? ", " : "", STR_ARG(inst_id));
                    if (n < 0 || (size_t)n >= sizeof(buf) - len) break;
                    len += (size_t)n;
                    const str_t auth_id = md_instance_auth_id(&sys.instance, info.instances[k]);
                    has_auth |= !str_empty(auth_id) && !str_eq(auth_id, inst_id);
                }
                info_row("Instances");
                ImGui::Text("%s  (%s%s)", fmt_count(num, sizeof(num), num_inst), buf, num_inst > max_ids ? ", ..." : "");

                if (has_auth) {
                    len = 0;
                    buf[0] = '\0';
                    for (size_t k = 0; k < MIN(num_inst, max_ids); ++k) {
                        const str_t auth_id = md_instance_auth_id(&sys.instance, info.instances[k]);
                        const int n = snprintf(buf + len, sizeof(buf) - len, "%s" STR_FMT, k ? ", " : "", STR_ARG(auth_id));
                        if (n < 0 || (size_t)n >= sizeof(buf) - len) break;
                        len += (size_t)n;
                    }
                    info_row("Author chains");
                    ImGui::Text("%s%s", buf, num_inst > max_ids ? ", ..." : "");
                }
            }

            char comps[64], atoms[64];
            fmt_count_range(comps, sizeof(comps), info.min_comps, info.max_comps);
            fmt_count_range(atoms, sizeof(atoms), info.min_atoms, info.max_atoms);
            info_row(num_inst > 1 ? "Per instance" : "Size");
            ImGui::Text("%s %s, %s atoms", comps, entity_kind_comp_noun(info.kind, info.max_comps), atoms);

            if (num_inst > 1) {
                const size_t total = md_system_atom_count(&sys);
                info_row("Total");
                ImGui::Text("%s atoms (%.1f%% of the system)", fmt_count(num, sizeof(num), info.num_atoms), total ? 100.0 * info.num_atoms / (double)total : 0.0);
            }

            // Mass, formula and contents of the first instance, the others are alike by construction
            if (num_inst > 0) {
                const int first = info.instances[0];
                const md_urange_t atom_range = md_system_instance_atom_range(&sys, first);
                const char* per = num_inst > 1 ? " each" : "";

                info_row("Mass");
                ImGui::Text("%s%s", fmt_mass(buf, sizeof(buf), atom_range_mass(sys, atom_range)), per);

                hill_formula(buf, sizeof(buf), sys, atom_range);
                if (buf[0]) {
                    info_row("Formula");
                    ImGui::Text("%s%s", buf, per);
                }

                if (info.kind == EntityKind_Protein || info.kind == EntityKind_NucleicAcid) {
                    sequence_string(buf, 64, sys, first);
                    info_row("Sequence");
                    ImGui::TextUnformatted(buf);
                } else {
                    // Distinct component names
                    const md_urange_t comp_range = md_system_instance_comp_range(&sys, first);
                    const size_t max_names = 6;
                    str_t names[max_names];
                    size_t num_names = 0;
                    bool more = false;
                    for (uint32_t c = comp_range.beg; c < comp_range.end; ++c) {
                        const str_t name = md_component_name(&sys.component, c);
                        bool found = false;
                        for (size_t k = 0; k < num_names; ++k) found |= str_eq(names[k], name);
                        if (found) continue;
                        if (num_names == max_names) { more = true; break; }
                        names[num_names++] = name;
                    }
                    if (num_names > 0) {
                        size_t len = 0;
                        buf[0] = '\0';
                        for (size_t k = 0; k < num_names; ++k) {
                            const int n = snprintf(buf + len, sizeof(buf) - len, "%s" STR_FMT, k ? ", " : "", STR_ARG(names[k]));
                            if (n < 0 || (size_t)n >= sizeof(buf) - len) break;
                            len += (size_t)n;
                        }
                        info_row(num_names > 1 ? "Components" : "Component");
                        ImGui::Text("%s%s", buf, more ? ", ..." : "");
                    }
                }
            }

            if (info.flags & (MD_FLAG_ISOMER_L | MD_FLAG_ISOMER_D)) {
                info_row("Chirality");
                ImGui::TextUnformatted((info.flags & MD_FLAG_ISOMER_L) ? "L" : "D");
            }
            if (info.flags & MD_FLAG_COARSE_GRAINED) {
                info_row("Model");
                ImGui::TextUnformatted("Coarse grained");
            }
            ImGui::EndTable();
        }

        ImGui::Separator();
        ImGui::TextDisabled("Shift + click to select, Shift + right click to deselect, right click for options");
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }

    void instance_tooltip(const ApplicationState& data, size_t e, size_t inst_idx) const {
        const md_system_t& sys = data.mold.sys;
        const EntityInfo& info = entities[e];
        const str_t id      = md_instance_id(&sys.instance, inst_idx);
        const str_t auth_id = md_instance_auth_id(&sys.instance, inst_idx);
        const md_urange_t comp_range = md_system_instance_comp_range(&sys, inst_idx);
        const md_urange_t atom_range = md_system_instance_atom_range(&sys, inst_idx);
        const size_t num_comps = comp_range.end - comp_range.beg;
        char buf[128], num[32];

        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
        ImGui::Text("Instance " STR_FMT, STR_ARG(id));
        if (!str_empty(auth_id) && !str_eq(auth_id, id)) {
            ImGui::SameLine();
            ImGui::TextDisabled("(author chain " STR_FMT ")", STR_ARG(auth_id));
        }
        ImGui::TextDisabled("Entity " STR_FMT ": %s", STR_ARG(md_entity_id(&sys.entity, e)), info.name);
        ImGui::Separator();

        if (ImGui::BeginTable("##instance_info", 2, ImGuiTableFlags_SizingFixedFit)) {
            info_row(entity_kind_comp_noun(info.kind, num_comps, true));
            if (num_comps > 1) {
                ImGui::Text("%s  (%d - %d)", fmt_count(num, sizeof(num), num_comps),
                    md_component_seq_id(&sys.component, comp_range.beg), md_component_seq_id(&sys.component, comp_range.end - 1));
            } else {
                ImGui::Text("%s", fmt_count(num, sizeof(num), num_comps));
            }

            info_row("Atoms");
            ImGui::Text("%s", fmt_count(num, sizeof(num), atom_range.end - atom_range.beg));

            info_row("Mass");
            ImGui::TextUnformatted(fmt_mass(buf, sizeof(buf), atom_range_mass(sys, atom_range)));

            hill_formula(buf, sizeof(buf), sys, atom_range);
            if (buf[0]) {
                info_row("Formula");
                ImGui::TextUnformatted(buf);
            }

            if (info.kind == EntityKind_Protein || info.kind == EntityKind_NucleicAcid) {
                sequence_string(buf, 64, sys, inst_idx);
                info_row("Sequence");
                ImGui::TextUnformatted(buf);
            }
            ImGui::EndTable();
        }
        if (info.derived) {
            ImGui::TextDisabled("Grouped on load, the instance ID is generated");
        }
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }

    void component_tooltip(const ApplicationState& data, size_t e, size_t inst_idx, size_t comp_idx) const {
        const md_system_t& sys = data.mold.sys;
        const EntityInfo& info = entities[e];
        const str_t name = md_component_name(&sys.component, comp_idx);
        const md_flags_t flags = md_component_flags(&sys.component, comp_idx);
        const int code = comp_idx < md_array_size(comp_code) ? comp_code[comp_idx] : -1;
        const md_urange_t atom_range = md_system_component_atom_range(&sys, comp_idx);

        ImGui::BeginTooltip();
        ImGui::Text(STR_FMT " %d", STR_ARG(name), md_component_seq_id(&sys.component, comp_idx));
        if (code >= 0) {
            ImGui::SameLine();
            ImGui::TextDisabled("%c, %s", residue_codes[code].code, residue_class_label(residue_codes[code].cls));
        }

        const bool nucleic = info.kind == EntityKind_NucleicAcid;
        char num[32];
        ImGui::Text("%s atoms%s%s", fmt_count(num, sizeof(num), atom_range.end - atom_range.beg),
            (flags & MD_FLAG_TERMINAL_BEG) ? (nucleic ? ", 5' terminus" : ", N-terminus") : "",
            (flags & MD_FLAG_TERMINAL_END) ? (nucleic ? ", 3' terminus" : ", C-terminus") : "");
        ImGui::TextDisabled("Instance " STR_FMT ", entity " STR_FMT, STR_ARG(md_instance_id(&sys.instance, inst_idx)), STR_ARG(md_entity_id(&sys.entity, e)));
        ImGui::EndTooltip();
    }

    // ## Settings

    static void residue_color_legend() {
        ImGui::TextDisabled("Amino acids by Clustal X class, nucleotides by base");
        const float sz = ImGui::GetTextLineHeight();
        for (int cls = ResidueClass_Hydrophobic; cls < ResidueClass_Count; ++cls) {
            if (cls == ResidueClass_Adenine) ImGui::Separator();

            // The single letter codes of the class
            char codes[32] = "";
            size_t len = 0;
            for (size_t i = 0; i < ARRAY_SIZE(residue_codes) && len + 2 < sizeof(codes); ++i) {
                if (residue_codes[i].cls != cls) continue;
                if (memchr(codes, residue_codes[i].code, len)) continue;
                if (len) codes[len++] = ' ';
                codes[len++] = residue_codes[i].code;
                codes[len] = '\0';
            }

            const ImVec2 p = ImGui::GetCursorScreenPos();
            ImGui::GetWindowDrawList()->AddRectFilled(p, ImVec2(p.x + sz, p.y + sz), residue_class_color((ResidueClass)cls), 2.0f);
            ImGui::Dummy(ImVec2(sz, sz));
            ImGui::SameLine();
            if (residue_class_is_nucleotide((ResidueClass)cls)) {
                ImGui::TextUnformatted(residue_class_label((ResidueClass)cls));
            } else {
                ImGui::Text("%s (%s)", residue_class_label((ResidueClass)cls), codes);
            }
        }
    }

    void draw_menu_bar() {
        if (!ImGui::BeginMenuBar()) return;
        if (ImGui::BeginMenu("Settings")) {
            ImGui::SeparatorText("Sequences");
            ImGui::MenuItem("Single letter residue codes", nullptr, &use_single_letter_codes);
            ImGui::SetItemTooltip("Show the amino acids and nucleotides of proteins and nucleic acids by their single letter code.\n"
                                  "Residues without a code keep their name.");
            ImGui::MenuItem("Color residues by type", nullptr, &color_residues);
            if (ImGui::BeginItemTooltip()) {
                residue_color_legend();
                ImGui::EndTooltip();
            }
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    // ## Entities: entity > instance > component

    // The components of a chain as a wrapped sequence of chips, single letter codes in groups of ten
    void draw_sequence(ApplicationState& data, size_t e, size_t inst_idx) {
        const md_system_t& sys = data.mold.sys;
        const md_urange_t range = md_system_instance_comp_range(&sys, inst_idx);
        if (range.beg >= range.end) return;

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float font      = ImGui::GetFontSize();
        const float pad_x     = floorf(font * 0.2f);
        const float pad_y     = 1.0f;
        const float gap       = 1.0f;
        const float group_gap = floorf(font * 0.5f);
        const float h         = ImGui::GetTextLineHeight() + 2 * pad_y;
        const float letter_w  = ImGui::CalcTextSize("W").x + 2 * pad_x;

        const vec4_t sel = data.selection.color.selection.visible;
        const vec4_t hl  = data.selection.color.highlight.visible;
        const uint32_t sel_col    = ImGui::ColorConvertFloat4ToU32(ImVec4(sel.x, sel.y, sel.z, 1.0f));
        const uint32_t hl_col     = ImGui::ColorConvertFloat4ToU32(ImVec4(hl.x, hl.y, hl.z, 0.6f));
        const uint32_t plain_fill = ImGui::GetColorU32(ImGuiCol_FrameBg);
        const uint32_t text_col   = color_residues ? IM_COL32(20, 20, 20, 255) : ImGui::GetColorU32(ImGuiCol_Text);

        const ImVec2 origin = ImGui::GetCursorScreenPos();
        const float  max_x  = origin.x + ImGui::GetContentRegionAvail().x;
        ImVec2 p = origin;

        for (uint32_t c = range.beg; c < range.end; ++c) {
            const uint32_t n = c - range.beg;
            const int code = c < md_array_size(comp_code) ? comp_code[c] : -1;

            char  lbl[16];
            float w;
            if (use_single_letter_codes && code >= 0) {
                lbl[0] = residue_codes[code].code;
                lbl[1] = '\0';
                w = letter_w;
                if (n > 0 && n % 10 == 0) p.x += group_gap;
            } else {
                str_copy_to_char_buf(lbl, sizeof(lbl), md_component_name(&sys.component, c));
                w = ImGui::CalcTextSize(lbl).x + 2 * pad_x;
            }

            if (p.x + w > max_x && p.x > origin.x) {
                p.x  = origin.x;
                p.y += h + gap;
            }
            const ImVec2 a = p;
            const ImVec2 b = ImVec2(p.x + w, p.y + h);
            p.x += w + gap;

            if (!ImGui::IsRectVisible(a, b)) continue;

            ImGui::SetCursorScreenPos(a);
            ImGui::PushID((int)c);
            ImGui::InvisibleButton("##res", ImVec2(w, h));
            const bool hovered = ImGui::IsItemHovered();
            const bool tooltip = ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip);
            target_interaction(data, Target_Component, c, hovered);
            if (tooltip) component_tooltip(data, e, inst_idx, c);
            ImGui::PopID();

            // Reflects the selection and what is hovered in the viewport
            const md_urange_t atoms = md_system_component_atom_range(&sys, c);
            const bool has_atoms   = atoms.beg < atoms.end;
            const bool selected    = has_atoms && md_bitfield_test_bit(&data.selection.selection_mask, atoms.beg);
            const bool highlighted = hovered || (has_atoms && md_bitfield_test_bit(&data.selection.highlight_mask, atoms.beg));

            const uint32_t fill = color_residues ? residue_class_color(code >= 0 ? residue_codes[code].cls : ResidueClass_Unknown) : plain_fill;
            dl->AddRectFilled(a, b, fill, 2.0f);
            if (highlighted) dl->AddRectFilled(a, b, hl_col, 2.0f);
            if (selected)    dl->AddRect(a, b, sel_col, 2.0f, 0, 2.0f);
            dl->AddText(ImVec2(a.x + pad_x, a.y + pad_y), text_col, lbl);
        }

        // Claim the space the sequence occupies
        ImGui::SetCursorScreenPos(ImVec2(origin.x, p.y + h));
        ImGui::Dummy(ImVec2(0.0f, gap));
    }

    // Components of an instance which is not a chain (e.g. a group of water molecules), one row each
    void draw_component_list(ApplicationState& data, size_t e, size_t inst_idx) {
        const md_system_t& sys = data.mold.sys;
        const md_urange_t range = md_system_instance_comp_range(&sys, inst_idx);

        ImGuiListClipper clipper;
        clipper.Begin((int)(range.end - range.beg));
        while (clipper.Step()) {
            for (int k = clipper.DisplayStart; k < clipper.DisplayEnd; ++k) {
                const size_t c = range.beg + (size_t)k;
                ImGui::PushID((int)c);
                const float right_x = ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
                const str_t name = md_component_name(&sys.component, c);
                ImGui::TreeNodeEx("##comp", ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen,
                    STR_FMT " %d", STR_ARG(name), md_component_seq_id(&sys.component, c));
                const bool hovered = ImGui::IsItemHovered();
                const bool tooltip = ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip);
                target_interaction(data, Target_Component, c, hovered);
                if (tooltip) component_tooltip(data, e, inst_idx, c);

                char num[32], txt[64];
                const size_t num_atoms = md_component_atom_count(&sys.component, c);
                snprintf(txt, sizeof(txt), "%s %s", fmt_count(num, sizeof(num), num_atoms), num_atoms == 1 ? "atom" : "atoms");
                right_aligned_disabled(right_x, txt);
                ImGui::PopID();
            }
        }
        clipper.End();
    }

    void draw_instance(ApplicationState& data, size_t e, int inst_idx) {
        const md_system_t& sys = data.mold.sys;
        const EntityInfo& info = entities[e];

        ImGui::PushID(inst_idx);
        defer { ImGui::PopID(); };

        const str_t id      = md_instance_id(&sys.instance, inst_idx);
        const str_t auth_id = md_instance_auth_id(&sys.instance, inst_idx);
        const md_urange_t comp_range = md_system_instance_comp_range(&sys, inst_idx);
        const size_t num_comps = comp_range.end - comp_range.beg;
        const size_t num_atoms = md_system_instance_atom_count(&sys, inst_idx);
        const bool   leaf      = num_comps <= 1;

        char label[128];
        size_t len = (size_t)MAX(0, snprintf(label, sizeof(label), STR_FMT, STR_ARG(id)));
        if (!str_empty(auth_id) && !str_eq(auth_id, id) && len < sizeof(label)) {
            len += (size_t)MAX(0, snprintf(label + len, sizeof(label) - len, " (auth " STR_FMT ")", STR_ARG(auth_id)));
        }
        if (leaf && num_comps == 1 && len < sizeof(label)) {
            snprintf(label + len, sizeof(label) - len, "  " STR_FMT, STR_ARG(md_component_name(&sys.component, comp_range.beg)));
        }

        ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_SpanAvailWidth;
        if (leaf) {
            flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
        } else if (md_array_size(info.instances) == 1 && entity_kind_is_sequence(info.kind)) {
            flags |= ImGuiTreeNodeFlags_DefaultOpen;
        }

        const float right_x = ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
        const bool open    = ImGui::TreeNodeEx("##inst", flags, "%s", label);
        const bool hovered = ImGui::IsItemHovered();
        const bool tooltip = ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip);
        target_interaction(data, Target_Instance, inst_idx, hovered);
        if (tooltip) instance_tooltip(data, e, inst_idx);

        char n0[32], n1[32], txt[96];
        if (leaf) {
            snprintf(txt, sizeof(txt), "%s atoms", fmt_count(n0, sizeof(n0), num_atoms));
        } else {
            snprintf(txt, sizeof(txt), "%s %s, %s atoms", fmt_count(n0, sizeof(n0), num_comps), entity_kind_comp_noun(info.kind, num_comps), fmt_count(n1, sizeof(n1), num_atoms));
        }
        right_aligned_disabled(right_x, txt);

        if (open && !leaf) {
            if (entity_kind_is_sequence(info.kind)) {
                draw_sequence(data, e, inst_idx);
            } else {
                draw_component_list(data, e, inst_idx);
            }
            ImGui::TreePop();
        }
    }

    void draw_entities(ApplicationState& data) {
        const md_system_t& sys = data.mold.sys;
        const size_t num_entities = md_system_entity_count(&sys);
        if (num_entities == 0 || md_array_size(entities) != num_entities) return;

        char header[64];
        snprintf(header, sizeof(header), "Entities (%zu)###Entities", num_entities);
        if (!ImGui::CollapsingHeader(header, ImGuiTreeNodeFlags_DefaultOpen)) return;

        for (size_t e = 0; e < num_entities; ++e) {
            const EntityInfo& info = entities[e];
            const size_t num_inst = md_array_size(info.instances);
            if (num_inst == 0) continue;

            ImGui::PushID((int)e);
            defer { ImGui::PopID(); };

            const str_t id = md_entity_id(&sys.entity, e);
            ImGuiTreeNodeFlags flags = ImGuiTreeNodeFlags_SpanAvailWidth;
            if (num_entities == 1) flags |= ImGuiTreeNodeFlags_DefaultOpen;

            const float right_x = ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
            const bool open    = ImGui::TreeNodeEx("##entity", flags, STR_FMT "  %s", STR_ARG(id), info.name);
            const bool hovered = ImGui::IsItemHovered();
            const bool tooltip = ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip);
            target_interaction(data, Target_Entity, e, hovered);
            if (tooltip) entity_tooltip(data, e);

            // The kind next to the name, unless the name already is the kind
            {
                const char* kind = entity_kind_label(info.kind);
                const bool name_is_kind = str_eq_cstr_ignore_case(str_from_cstr(info.name), kind);
                char txt[64] = "";
                if (!name_is_kind) {
                    snprintf(txt, sizeof(txt), "%s%s", kind, info.derived ? ", derived" : "");
                } else if (info.derived) {
                    snprintf(txt, sizeof(txt), "derived");
                }
                if (txt[0]) {
                    ImGui::SameLine();
                    ImGui::TextDisabled("%s", txt);
                }
            }

            {
                char n0[32], n1[32], txt[128];
                size_t len = 0;
                if (num_inst > 1) {
                    len += (size_t)MAX(0, snprintf(txt + len, sizeof(txt) - len, "%s instances, ", fmt_count(n0, sizeof(n0), num_inst)));
                } else if (!entity_kind_is_sequence(info.kind) && info.num_comps > 1) {
                    len += (size_t)MAX(0, snprintf(txt + len, sizeof(txt) - len, "%s %s, ", fmt_count(n0, sizeof(n0), info.num_comps), entity_kind_comp_noun(info.kind, info.num_comps)));
                }
                if (len < sizeof(txt)) {
                    snprintf(txt + len, sizeof(txt) - len, "%s atoms", fmt_count(n1, sizeof(n1), info.num_atoms));
                }
                right_aligned_disabled(right_x, txt);
            }

            if (open) {
                if (info.max_comps <= 1) {
                    // Every instance is a single component (ions, lipids, small molecules), possibly thousands of them
                    ImGuiListClipper clipper;
                    clipper.Begin((int)num_inst);
                    while (clipper.Step()) {
                        for (int k = clipper.DisplayStart; k < clipper.DisplayEnd; ++k) {
                            draw_instance(data, e, info.instances[k]);
                        }
                    }
                    clipper.End();
                } else {
                    for (size_t k = 0; k < num_inst; ++k) {
                        draw_instance(data, e, info.instances[k]);
                    }
                }
                ImGui::TreePop();
            }
        }
        ImGui::Spacing();
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

        item.use_defaults = !(load.flags & MD_FLAG_COARSE_GRAINED);
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
        type.flags[i]  = (type.flags[i] & ~MD_FLAG_COARSE_GRAINED) | (load.flags & MD_FLAG_COARSE_GRAINED);
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
        const bool cg = type.flags[i] & MD_FLAG_COARSE_GRAINED;
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
        const bool cg_a = type.flags[a] & MD_FLAG_COARSE_GRAINED;
        const bool cg_b = type.flags[b] & MD_FLAG_COARSE_GRAINED;
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

        bool coarse_grained = type.flags[i] & MD_FLAG_COARSE_GRAINED;
        if (ImGui::Checkbox("Coarse grained", &coarse_grained)) {
            if (coarse_grained) {
                type.flags[i] |=  MD_FLAG_COARSE_GRAINED;
            } else {
                type.flags[i] &= ~MD_FLAG_COARSE_GRAINED;
            }
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

    void draw_atom_types(ApplicationState& data, bool& radius_changed, bool& color_changed, bool& mass_changed) {
        md_system_t& sys = data.mold.sys;
        md_atom_type_data_t& type = sys.atom.type;
        const size_t num_types = md_system_atom_type_count(&sys);
        if (num_types == 0 || md_array_size(atom_types) != num_types) return;

        // The sentinel "unknown" type is only listed when used
        const size_t first = atom_types[0].count == 0 ? 1 : 0;

        char header[64];
        snprintf(header, sizeof(header), "Atom Types (%zu)###AtomTypes", num_types - first);
        if (!ImGui::CollapsingHeader(header)) return;

        if (atom_type_selected >= (int)num_types) atom_type_selected = -1;

        size_t num_modified = 0;
        for (size_t i = first; i < num_types; ++i) {
            num_modified += compute_atom_type_delta(type, atom_types[i], i) ? 1 : 0;
        }

        const bool show_filter = num_types - first > 8;
        if (show_filter || num_modified) {
            char reset_lbl[48] = "";
            float reset_w = 0;
            if (num_modified) {
                snprintf(reset_lbl, sizeof(reset_lbl), "Reset %zu modified", num_modified);
                reset_w = ImGui::CalcTextSize(reset_lbl).x + ImGui::GetStyle().FramePadding.x * 2;
            }
            if (show_filter) {
                ImGui::SetNextItemWidth(num_modified ? -(reset_w + ImGui::GetStyle().ItemSpacing.x) : -FLT_MIN);
                ImGui::InputTextWithHint("##atom_type_filter", "Filter by name or element", atom_type_filter, sizeof(atom_type_filter));
                if (num_modified) ImGui::SameLine();
            }
            if (num_modified) {
                if (ImGui::Button(reset_lbl)) {
                    for (size_t i = first; i < num_types; ++i) {
                        if (compute_atom_type_delta(type, atom_types[i], i)) {
                            reset_atom_type(data, i, radius_changed, color_changed, mass_changed);
                        }
                    }
                }
                ImGui::SetItemTooltip("Restore the modified atom types (marked *) to what was loaded");
            }
        }
        const str_t filter = show_filter ? str_trim(str_from_cstr(atom_type_filter)) : STR_LIT("");

        auto passes_filter = [&](size_t i) {
            if (i < first) return false;
            if (str_empty(filter)) return true;
            const bool cg = type.flags[i] & MD_FLAG_COARSE_GRAINED;
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
            md_bitfield_clear(&data.selection.highlight_mask);
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
                const bool cg       = type.flags[i] & MD_FLAG_COARSE_GRAINED;
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
                        md_bitfield_clear(&data.selection.highlight_mask);
                        md_bitfield_set_indices_u32(&data.selection.highlight_mask, (uint32_t*)item.indices, md_array_size(item.indices));
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
            md_bitfield_clear(&data.selection.highlight_mask);
            for (size_t i = 0; i < num_types; ++i) {
                const DatasetItem& item = atom_types[i];
                if (type.z[i] == table_res.z) {
                    md_bitfield_set_indices_u32(&data.selection.highlight_mask, (uint32_t*)item.indices, md_array_size(item.indices));
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

    void draw(ApplicationState& data) {
        if (!show_window) return;

        ImGui::SetNextWindowSize(ImVec2(500, 600), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("System", &show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
            draw_menu_bar();
            draw_files(data);
            draw_contents(data);

            if (ImGui::IsWindowHovered()) {
                md_bitfield_clear(&data.selection.highlight_mask);
            }

            draw_series(data);
            draw_entities(data);

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
        ImGui::End();
    }
};

static Dataset instance;

}  // namespace dataset
