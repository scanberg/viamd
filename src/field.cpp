#include <field.h>

#include <viamd.h>
#include <imgui.h>

#include <core/md_common.h>

// The fields of the electronic structure: each a source of the representation (the density properties apart, a leaf
// each under es/density_property)
struct EsField {
    const char* path;
    const char* label;
    ElectronicStructureSource source;
    const char* keywords;
    const char* description;    // what it is, as the picker and the field's row show it on hover
};

static const EsField es_fields[] = {
    { "es/orbital", "Molecular Orbital", ElectronicStructureSource::MolecularOrbital, "mo homo lumo",
      "A molecular orbital psi_i of the calculation, of either spin: its value, the positive and negative lobes apart (their isosurfaces at +-0.05 au show its shape and its nodes), or its square, the density of an electron in it." },
    { "es/nto", "Natural Transition Orbital", ElectronicStructureSource::NaturalTransitionOrbital, "excited state particle hole",
      "A natural transition orbital of an excited state: the particle (where the excited electron goes) or the hole (where it leaves), of the pairs that together are the transition, each pair with its weight lambda. One pair of large weight: the transition is that one excitation." },
    { "es/transition_density", "Transition Density", ElectronicStructureSource::TransitionDensity, "excited state attachment detachment",
      "The densities of an excitation: the attachment and the detachment, where the excited state has more and less electron density than the ground state, or the transition density between the two states (whose dipole is the transition dipole)." },
    { "es/density", "Electron Density", ElectronicStructureSource::ElectronDensity, "rho spin ground state",
      "The ground state electron density rho of the QM region, in all or of a spin (alpha, beta, or the spin density, their difference). Its isosurface at 0.001-0.002 au is the usual molecular surface, what the energies and potentials below are read on." },
    { "es/elf", "Electron Localization Function", ElectronicStructureSource::ElectronLocalization, "elf basins",
      "The electron localization function (Becke and Edgecombe), 0 to 1: near 1 where the electrons are paired and held (cores, bonds, lone pairs), 1/2 in an electron gas of the same density. From the kinetic energy density of the orbitals against the gas's; its basins are the shells, the bonds and the lone pairs." },
    { "es/lol", "Localized Orbital Locator", ElectronicStructureSource::LocalizedOrbitalLocator, "lol",
      "The localized orbital locator (Schmider and Becke), 0 to 1: like ELF large where the electrons are localized, from the kinetic energy density alone, t / (1 + t) of t = tau_gas / tau. Sharper than ELF between the atoms." },
    { "es/esp", "Electrostatic Potential", ElectronicStructureSource::ElectrostaticPotential, "esp mep molecular",
      "The electrostatic potential of the QM region: its nuclei and its ground state electron density, exact over the basis. Negative (red) where the electrons dominate - lone pairs, pi systems, where an electrophile is drawn - and positive (blue) where the nuclei are less screened. On an isodensity surface (0.001-0.002 au) it is the molecular electrostatic potential. An embedding's environment is not part of it." },
    { "es/embedding_potential", "Embedding Potential", ElectronicStructureSource::EmbeddingPotential, "polarizable pe environment solvent",
      "The electrostatic potential of the environment: its charges, the dipoles and quadrupoles where it has them, and for a polarizable embedding the dipoles the calculation induced, where it kept them - the potential the QM region was embedded in. The QM region is left out, whatever charges a topology gives its atoms." },
    { "es/alie", "Average Local Ionization Energy", ElectronicStructureSource::AverageLocalIonizationEnergy, "alie",
      "The average local ionization energy of the QM region: the orbital energies of its occupied orbitals, weighted by how much each contributes to the electron density there, sum n |e| |phi|^2 / sum n |phi|^2. Lowest (red) where the least tightly held electrons are - the sites an electrophile or a radical attacks - and never below the HOMO's |e|. Meant for an isodensity surface (0.001-0.002 au). Kohn-Sham orbital energies give lower values than Hartree-Fock." },
    { "es/leae", "Local Electron Attachment Energy", ElectronicStructureSource::LocalElectronAttachmentEnergy, "leae",
      "The local electron attachment energy of the QM region: the energies of its virtual orbitals below zero, the ones that bind an added electron, weighted by their densities against the electron density. Never positive; lowest where a nucleophile attacks or a halogen bond forms (the sigma hole). Defined for hybrid DFT orbitals and read on the 0.004 au isodensity surface; it depends strongly on the method and basis, and is 0 everywhere with Hartree-Fock virtuals." },
    { "es/density_scale_length", "Density Scale Length", ElectronicStructureSource::DensityScaleLength, "dsl",
      "The density scale length rho / |grad rho| of the QM region's ground state density: the distance over which the density changes by a factor e. On an isodensity surface it is how far the surface moves along its normal per relative change of the density: long where the boundary of the electron cloud is soft, short where it is stiff, as when the density is compressed under pressure (GOSTSHYP). About 0.3-1 bohr on the 0.001 au surface." },
    { "es/gostshyp/wall_potential", "Model Potential", ElectronicStructureSource::GostshypWallPotential, "pressure cavity wall gaussian v1",
      "The GOSTSHYP model potential, Σ_j G_j(r) = Σ_j p_j G~_j(r): on each tessera j of the cavity a normalized Gaussian G~_j = (ω_j/π)^(3/2) exp(-ω_j |r - r_j|²), ω_j = π ln2 / a_j, of the amplitude p_j = p_inp a_j / f~_j the SCF settled, that holds the electron density in against the input pressure. Positive, largest at the cavity's surface, about a tessera's size thick. Its expectation value is the GOSTSHYP energy E^GOST = Σ_j p_j g~_j." },
    { "es/gostshyp/effective_potential", "Fock Contribution", ElectronicStructureSource::GostshypEffectivePotential, "pressure cavity effective potential v1 v2",
      "The potential whose matrix is the GOSTSHYP Fock matrix contribution V^GOST: what the electrons feel, Σ_j p_j G~_j(r) (1 + 2 ω_j λ_j n_j·(r - r_j)), n_j the outward normal, λ_j = g~_j / f~_j the decay length. The second term comes of the amplitudes depending on the density; to first order each Gaussian moved outward by its decay length. Higher than the model potential outside the cavity's surface, lower inside, and negative a little way in." },
    { "es/gostshyp/energy_density", "Energy Density", ElectronicStructureSource::GostshypEnergyDensity, "pressure cavity wall energy",
      "The GOSTSHYP energy per volume, ρ(r) Σ_j G_j(r): the electron density times the model potential. Where the GOSTSHYP energy is; its integral over space is E^GOST = Σ_j p_j g~_j, the cavity's GOSTSHYP Energy." },
};

static const md_str_t density_property_prefix = MD_STR_LIT("es/density_property/");

struct GroupLabel {
    const char* group;
    const char* label;
};

static const GroupLabel group_labels[] = {
    { "es",                  "Electronic Structure" },
    { "es/gostshyp",         "GOSTSHYP" },
    { "es/density_property", "Density Properties" },
};

md_str_t field_group_label(md_str_t group) {
    for (const GroupLabel& g : group_labels) {
        if (md_str_eq_cstr(group, g.group)) return md_str_from_cstr(g.label);
    }
    const size_t slash = md_str_rfind_char(group, '/');
    return slash == MD_STR_NPOS ? group : md_str_substr(group, slash + 1);
}

// A density property's path below its group in the attributes (vlx/density_property/<this>)
static md_str_t density_property_rel(md_str_t attr_path) {
    md_str_t rel = md_str_substr(attr_path, es_path::density_property.len);
    if (rel.len && rel.ptr[0] == '/') rel = md_str_substr(rel, 1);
    return rel;
}

static void entry_path(FieldEntry* e, const char* path) {
    snprintf(e->path, sizeof(e->path), "%s", path);
}

size_t field_catalog(md_array(FieldEntry)* out, const md_system_t& sys, md_allocator_i* alloc) {
    const ElectronicStructureSourceFlags mask = es_source_mask(sys);
    size_t count = 0;
    for (const EsField& f : es_fields) {
        if (!electronic_structure_source_supported(mask, f.source)) continue;
        FieldEntry e = {};
        entry_path(&e, f.path);
        e.label = md_str_from_cstr(f.label);
        e.keywords = f.keywords;
        e.description = md_str_from_cstr(f.description);
        e.source = f.source;
        // the density colours as the total density, whatever its spin shown
        e.field_kind = f.source == ElectronicStructureSource::ElectronDensity ? SurfaceFieldKind::ElectronDensity : electronic_structure_source_field_kind(f.source);
        md_array_push(*out, e, alloc);
        count += 1;
    }
    // A leaf each, after the densities they are of
    for (md_attribute_iter_t it = md_attributes_iter(&sys.attributes, es_path::density_property); md_attributes_next(&it);) {
        const md_attribute_t* attr = it.attr;
        const md_str_t rel = density_property_rel(attr->path);
        FieldEntry e = {};
        snprintf(e.path, sizeof(e.path), "%.*s%.*s", (int)density_property_prefix.len, density_property_prefix.ptr, (int)rel.len, rel.ptr);
        e.label = md_str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
        e.keywords = "density property";
        e.description = attr->description;
        e.source = ElectronicStructureSource::DensityProperty;
        e.density_property_key = attr->id;
        md_array_push(*out, e, alloc);
        count += 1;
    }
    return count;
}

size_t field_path_of_source(char* buf, size_t cap, ElectronicStructureSource source, md_attribute_id_t key, const md_system_t& sys) {
    if (source == ElectronicStructureSource::DensityProperty) {
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, key);
        if (!attr) return 0;
        const md_str_t rel = density_property_rel(attr->path);
        return (size_t)snprintf(buf, cap, "%.*s%.*s", (int)density_property_prefix.len, density_property_prefix.ptr, (int)rel.len, rel.ptr);
    }
    for (const EsField& f : es_fields) {
        if (f.source == source) return (size_t)snprintf(buf, cap, "%s", f.path);
    }
    return 0;
}

size_t field_path_of_kind(char* buf, size_t cap, SurfaceFieldKind kind) {
    for (const EsField& f : es_fields) {
        if (electronic_structure_source_field_kind(f.source) == kind) return (size_t)snprintf(buf, cap, "%s", f.path);
    }
    // the electron density as a surface field is the density's field
    if (kind == SurfaceFieldKind::ElectronDensity) return (size_t)snprintf(buf, cap, "es/density");
    return 0;
}

bool field_resolve(ElectronicStructureSource* source, md_attribute_id_t* key, md_str_t path) {
    if (md_str_begins_with(path, density_property_prefix)) {
        char attr[FIELD_PATH_MAX];
        const md_str_t rel = md_str_substr(path, density_property_prefix.len);
        const int n = snprintf(attr, sizeof(attr), "%.*s/%.*s", (int)es_path::density_property.len, es_path::density_property.ptr, (int)rel.len, rel.ptr);
        if (n <= 0 || n >= (int)sizeof(attr)) return false;
        *source = ElectronicStructureSource::DensityProperty;
        *key = md_attributes_id_from_path(md_str_t{attr, (size_t)n});
        return true;
    }
    for (const EsField& f : es_fields) {
        if (md_str_eq_cstr(path, f.path)) {
            *source = f.source;
            *key = MD_ATTRIBUTE_INVALID;
            return true;
        }
    }
    return false;
}

bool field_resolve_kind(SurfaceFieldKind* kind, md_str_t path) {
    if (md_str_eq_cstr(path, "es/density")) {
        *kind = SurfaceFieldKind::ElectronDensity;
        return true;
    }
    for (const EsField& f : es_fields) {
        if (md_str_eq_cstr(path, f.path)) {
            const SurfaceFieldKind k = electronic_structure_source_field_kind(f.source);
            if (k == SurfaceFieldKind::Count) return false;
            *kind = k;
            return true;
        }
    }
    return false;
}

static bool contains_ignore_case(md_str_t haystack, md_str_t needle) {
    if (needle.len == 0) return true;
    for (size_t i = 0; i + needle.len <= haystack.len; ++i) {
        if (md_str_eq_ignore_case(md_str_substr(haystack, i, needle.len), needle)) return true;
    }
    return false;
}

// Every whitespace separated word of the filter in the path, the label, the keywords or the group's label
static bool field_matches(md_str_t filter, const FieldEntry& e, md_str_t group_label) {
    const md_str_t path = md_str_from_cstr(e.path);
    const md_str_t kw = md_str_from_cstr(e.keywords);
    size_t i = 0;
    while (i < filter.len) {
        while (i < filter.len && (filter.ptr[i] == ' ' || filter.ptr[i] == '\t')) ++i;
        size_t j = i;
        while (j < filter.len && filter.ptr[j] != ' ' && filter.ptr[j] != '\t') ++j;
        if (j > i) {
            const md_str_t w = {filter.ptr + i, j - i};
            if (!contains_ignore_case(path, w) && !contains_ignore_case(e.label, w) && !contains_ignore_case(kw, w) &&
                !contains_ignore_case(group_label, w)) {
                return false;
            }
        }
        i = j;
    }
    return true;
}

// The segments of a path before its leaf
static md_str_t path_group(const char* path) {
    const md_str_t p = md_str_from_cstr(path);
    const size_t slash = md_str_rfind_char(p, '/');
    return slash == MD_STR_NPOS ? md_str_t{} : md_str_substr(p, 0, slash);
}

// "Electronic Structure / GOSTSHYP": the labels of a group and those it is in
static void group_title(char* buf, size_t cap, md_str_t group) {
    size_t len = 0;
    buf[0] = '\0';
    for (size_t end = 0; end <= group.len; ++end) {
        if (end < group.len && group.ptr[end] != '/') continue;
        const md_str_t label = field_group_label(md_str_substr(group, 0, end));
        const int n = snprintf(buf + len, cap - len, "%s%.*s", len ? " / " : "", (int)label.len, label.ptr);
        if (n < 0 || (size_t)n >= cap - len) break;
        len += (size_t)n;
    }
}

int field_picker(const char* id, const char* preview, const FieldEntry* entries, size_t count, md_str_t selected, bool (*pickable)(const FieldEntry&)) {
    int picked = -1;
    ImGui::PushID(id);
    ImGui::SetNextWindowSizeConstraints(ImVec2(0, 0), ImVec2(FLT_MAX, ImGui::GetTextLineHeightWithSpacing() * 24));
    if (ImGui::BeginCombo("##field", preview, ImGuiComboFlags_HeightLargest)) {
        // one popup is open at a time: one filter, cleared as it opens
        static char filter[128];
        if (ImGui::IsWindowAppearing()) {
            filter[0] = '\0';
            ImGui::SetKeyboardFocusHere();
        }
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputTextWithHint("##filter", "Search: words of the name or the path", filter, sizeof(filter));
        const md_str_t f = md_str_from_cstr(filter);

        md_str_t last_group = {};
        bool any = false;
        for (size_t i = 0; i < count; ++i) {
            const FieldEntry& e = entries[i];
            if (pickable && !pickable(e)) continue;
            const md_str_t group = path_group(e.path);
            char title[160];
            group_title(title, sizeof(title), group);
            if (!field_matches(f, e, md_str_from_cstr(title))) continue;
            if (!any || !md_str_eq(group, last_group)) {
                if (any) ImGui::Spacing();
                ImGui::TextDisabled("%s", title);
                last_group = group;
            }
            any = true;
            const bool is_selected = md_str_eq_cstr(selected, e.path);
            char label[256];
            snprintf(label, sizeof(label), "%.*s##%s", (int)e.label.len, e.label.ptr, e.path);
            ImGui::Indent();
            if (ImGui::Selectable(label, is_selected)) picked = (int)i;
            if (ImGui::IsItemHovered()) field_tooltip(e);
            if (is_selected && ImGui::IsWindowAppearing()) ImGui::SetScrollHereY();
            ImGui::Unindent();
        }
        if (!any) ImGui::TextDisabled("Nothing matches");
        ImGui::EndCombo();
    } else if (ImGui::IsItemHovered(ImGuiHoveredFlags_ForTooltip)) {
        // the one picked, closed: what it is
        for (size_t i = 0; i < count; ++i) {
            if (md_str_eq_cstr(selected, entries[i].path)) { field_tooltip(entries[i]); break; }
        }
    }
    ImGui::PopID();
    return picked;
}

const char* field_description_of_kind(SurfaceFieldKind kind) {
    for (const EsField& f : es_fields) {
        if (electronic_structure_source_field_kind(f.source) == kind) return f.description;
    }
    return "";
}

void field_tooltip(const FieldEntry& e) {
    if (!ImGui::BeginTooltip()) return;
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 40.0f);
    if (e.description.len) {
        ImGui::TextUnformatted(e.description.ptr, e.description.ptr + e.description.len);
        ImGui::Spacing();
    }
    ImGui::TextDisabled("%s", e.path);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
}
