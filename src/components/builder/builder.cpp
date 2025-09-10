#include <event.h>

#ifdef VIAMD_ENABLE_BUILDER

// Include RDKit headers first to avoid naming conflicts
#ifdef VIAMD_ENABLE_RDKIT
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#endif

#include <viamd.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <md_molecule.h>
#include <md_util.h>

#include <imgui_widgets.h>
#include <imgui.h>
#include <loader.h>

#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <map>
#include <cmath>

namespace builder {

struct MoleculeBuilder : viamd::EventHandler {
    bool show_window = false;
    bool rdkit_available = false;
    
    char smiles_input[256] = "CCO";  // Default to ethanol
    char error_message[512] = "";
    char info_message[512] = "";
    
    ApplicationState* app_state = nullptr;
    md_allocator_i* arena = nullptr;
    
    // Built molecule data
    struct BuiltMolecule {
        md_molecule_t mol = {};
        bool valid = false;
        int num_atoms = 0;
        int num_bonds = 0;
        std::string formula;
    } built_molecule;

    MoleculeBuilder() { 
        viamd::event_system_register_handler(*this);
        
#ifdef VIAMD_ENABLE_RDKIT
        rdkit_available = true;
        MD_LOG_INFO("RDKit support enabled for molecule builder");
#else
        rdkit_available = false;
        MD_LOG_INFO("RDKit not available - molecule builder will show informational message");
#endif
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& e = events[i];

            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                app_state = (ApplicationState*)e.payload;
                arena = md_arena_allocator_create(app_state->allocator.persistent, MEGABYTES(1));
                MD_LOG_INFO("Molecule Builder component initialized");
                break;
            }
            case viamd::EventType_ViamdShutdown:
                if (arena) {
                    cleanup_built_molecule();
                    md_arena_allocator_destroy(arena);
                }
                break;
            case viamd::EventType_ViamdFrameTick:
                draw_window();
                break;
            case viamd::EventType_ViamdWindowDrawMenu: {
                if (ImGui::BeginMenu("Builder")) {
                    ImGui::Checkbox("Molecule Builder", &show_window);
                    ImGui::EndMenu();
                }
                break;
            }
            default:
                break;
            }
        }
    }

    void cleanup_built_molecule() {
        if (built_molecule.valid) {
            md_molecule_free(&built_molecule.mol, arena);
            built_molecule = {};
        }
    }

#ifdef VIAMD_ENABLE_RDKIT
    bool build_molecule_from_smiles(const char* smiles) {
        if (!smiles || strlen(smiles) == 0) {
            strcpy(error_message, "Please enter a SMILES string");
            return false;
        }

        try {
            // Parse SMILES
            std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));
            if (!mol) {
                snprintf(error_message, sizeof(error_message), "Invalid SMILES: %s", smiles);
                return false;
            }

            // Add hydrogens
            RDKit::MolOps::addHs(*mol);
            
            // Generate 3D coordinates
            auto confId = RDKit::DGeomHelpers::EmbedMolecule(*mol);
            if (confId == -1) {
                strcpy(error_message, "Failed to generate 3D coordinates");
                return false;
            }

            // Optimize geometry with UFF
            try {
                RDKit::UFF::UFFOptimizeMolecule(*mol);
            } catch (...) {
                // UFF optimization failed, but we can still use the molecule
                MD_LOG_DEBUG("UFF optimization failed, using unoptimized geometry");
            }

            // Convert RDKit molecule to VIAMD format
            return convert_rdkit_to_viamd(*mol);

        } catch (const std::exception& e) {
            snprintf(error_message, sizeof(error_message), "RDKit error: %s", e.what());
            return false;
        } catch (...) {
            strcpy(error_message, "Unknown error in RDKit processing");
            return false;
        }
    }

    bool convert_rdkit_to_viamd(const RDKit::RWMol& rdkit_mol) {
        cleanup_built_molecule();

        // Initialize VIAMD molecule structure - zero initialize
        built_molecule.mol = {};

        const auto& conf = rdkit_mol.getConformer();
        unsigned int num_atoms = rdkit_mol.getNumAtoms();
        unsigned int num_bonds = rdkit_mol.getNumBonds();

        // Allocate arrays for atoms
        md_array_resize(built_molecule.mol.atom.x, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.y, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.z, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.element, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.type, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.radius, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.mass, num_atoms, arena);
        md_array_resize(built_molecule.mol.atom.flags, num_atoms, arena);

        // Convert atoms
        for (unsigned int i = 0; i < num_atoms; ++i) {
            const auto* atom = rdkit_mol.getAtomWithIdx(i);
            const auto& pos = conf.getAtomPos(i);

            // Position (convert from Angstrom to nanometers)
            built_molecule.mol.atom.x[i] = (float)(pos.x * 0.1);
            built_molecule.mol.atom.y[i] = (float)(pos.y * 0.1);
            built_molecule.mol.atom.z[i] = (float)(pos.z * 0.1);

            // Element
            built_molecule.mol.atom.element[i] = (uint8_t)atom->getAtomicNum();

            // Atom type (use element symbol)
            str_t element_name = md_util_element_symbol(atom->getAtomicNum());
            // Convert str_t to md_label_t by copying the string data
            size_t copy_len = MIN(element_name.len, sizeof(built_molecule.mol.atom.type[i].buf) - 1);
            memcpy(built_molecule.mol.atom.type[i].buf, element_name.ptr, copy_len);
            built_molecule.mol.atom.type[i].buf[copy_len] = '\0';

            // Properties
            built_molecule.mol.atom.radius[i] = md_util_element_vdw_radius(atom->getAtomicNum());
            built_molecule.mol.atom.mass[i] = (float)md_util_element_atomic_mass(atom->getAtomicNum());
            built_molecule.mol.atom.flags[i] = 0;
        }

        built_molecule.mol.atom.count = num_atoms;

        // Convert bonds
        if (num_bonds > 0) {
            md_array_resize(built_molecule.mol.bond.pairs, num_bonds, arena);
            md_array_resize(built_molecule.mol.bond.order, num_bonds, arena);

            for (unsigned int i = 0; i < num_bonds; ++i) {
                const auto* bond = rdkit_mol.getBondWithIdx(i);
                
                built_molecule.mol.bond.pairs[i].idx[0] = bond->getBeginAtomIdx();
                built_molecule.mol.bond.pairs[i].idx[1] = bond->getEndAtomIdx();
                
                // Convert bond order
                switch (bond->getBondType()) {
                case RDKit::Bond::SINGLE:
                    built_molecule.mol.bond.order[i] = 1;
                    break;
                case RDKit::Bond::DOUBLE:
                    built_molecule.mol.bond.order[i] = 2;
                    break;
                case RDKit::Bond::TRIPLE:
                    built_molecule.mol.bond.order[i] = 3;
                    break;
                case RDKit::Bond::AROMATIC:
                    built_molecule.mol.bond.order[i] = 1 | MD_BOND_FLAG_AROMATIC;
                    break;
                default:
                    built_molecule.mol.bond.order[i] = 1;
                    break;
                }
            }
            built_molecule.mol.bond.count = num_bonds;
        }

        // Set molecule info
        built_molecule.num_atoms = num_atoms;
        built_molecule.num_bonds = num_bonds;
        built_molecule.formula = RDKit::Descriptors::calcMolFormula(rdkit_mol);
        built_molecule.valid = true;

        snprintf(info_message, sizeof(info_message), 
                 "Built molecule: %s (%d atoms, %d bonds)", 
                 built_molecule.formula.c_str(), num_atoms, num_bonds);

        error_message[0] = '\0';  // Clear any previous errors
        return true;
    }
#endif

    void load_molecule_into_viamd() {
        if (!built_molecule.valid || !app_state) {
            strcpy(error_message, "No valid molecule to load");
            return;
        }

        // Free existing molecule
        md_molecule_free(&app_state->mold.mol, app_state->allocator.persistent);
        
        // Copy the built molecule to the app state
        app_state->mold.mol = built_molecule.mol;
        
        // Transfer ownership by clearing our copy without freeing
        built_molecule.mol = {};
        built_molecule.valid = false;

        // Update app state flags
        app_state->mold.dirty_buffers = MolBit_DirtyPosition | MolBit_DirtyRadius | MolBit_DirtyBonds;
        
        // Reset animation
        app_state->animation.frame = 0.0;
        
        // Clear trajectory - check if exists first
        if (app_state->mold.traj) {
            ::load::traj::close(app_state->mold.traj);
            app_state->mold.traj = nullptr;
        }

        snprintf(info_message, sizeof(info_message), 
                 "Molecule loaded successfully: %s", built_molecule.formula.c_str());
        
        MD_LOG_INFO("Molecule builder: loaded molecule into VIAMD");
    }

    void draw_example_buttons() {
        const char* examples[][2] = {
            {"Water", "O"},
            {"Methane", "C"},
            {"Ethanol", "CCO"},
            {"Benzene", "c1ccccc1"},
            {"Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
            {"Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"},
            {"Glucose", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"},
        };

        ImGui::Text("Quick Examples:");
        
        for (size_t i = 0; i < sizeof(examples) / sizeof(examples[0]); ++i) {
            if (i > 0 && i % 2 == 0) ImGui::NewLine();
            if (i % 2 == 1) ImGui::SameLine();
            
            if (ImGui::Button(examples[i][0])) {
                strcpy(smiles_input, examples[i][1]);
            }
        }
    }

    void draw_window() {
        if (!show_window) return;

        ImGui::SetNextWindowSize({400, 500}, ImGuiCond_FirstUseEver);
        
        if (ImGui::Begin("Molecule Builder", &show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (!rdkit_available) {
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.8f, 0.0f, 1.0f));
                ImGui::TextWrapped("RDKit is not available. To use the molecule builder, please install RDKit development libraries:");
                ImGui::PopStyleColor();
                
                ImGui::Separator();
                ImGui::Text("Installation instructions:");
                ImGui::BulletText("Ubuntu/Debian: sudo apt install librdkit-dev librdkit1");
                ImGui::BulletText("Conda: conda install -c conda-forge rdkit");
                ImGui::BulletText("Then rebuild VIAMD with -DVIAMD_ENABLE_RDKIT=ON");
                
                ImGui::Separator();
                ImGui::TextWrapped("The molecule builder allows you to:");
                ImGui::BulletText("Build molecules from SMILES strings");
                ImGui::BulletText("Generate 3D coordinates automatically");
                ImGui::BulletText("Load molecules directly into VIAMD");
                ImGui::BulletText("Use common molecule examples");
            } else {
                // SMILES input
                ImGui::Text("SMILES String:");
                ImGui::SetNextItemWidth(-1);
                ImGui::InputText("##smiles", smiles_input, sizeof(smiles_input));
                
                ImGui::Separator();
                
                // Example buttons
                draw_example_buttons();
                
                ImGui::Separator();
                
                // Build button
#ifdef VIAMD_ENABLE_RDKIT
                if (ImGui::Button("Build Molecule", ImVec2(-1, 0))) {
                    build_molecule_from_smiles(smiles_input);
                }
#else
                ImGui::BeginDisabled();
                ImGui::Button("Build Molecule (RDKit required)", ImVec2(-1, 0));
                ImGui::EndDisabled();
#endif
                
                // Load button (only enabled if we have a valid molecule)
                ImGui::BeginDisabled(!built_molecule.valid);
                if (ImGui::Button("Load into VIAMD", ImVec2(-1, 0))) {
                    load_molecule_into_viamd();
                }
                ImGui::EndDisabled();
                
                ImGui::Separator();
                
                // Status messages
                if (strlen(error_message) > 0) {
                    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.4f, 0.4f, 1.0f));
                    ImGui::TextWrapped("Error: %s", error_message);
                    ImGui::PopStyleColor();
                }
                
                if (strlen(info_message) > 0) {
                    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.4f, 1.0f, 0.4f, 1.0f));
                    ImGui::TextWrapped("Info: %s", info_message);
                    ImGui::PopStyleColor();
                }
                
                // Molecule info
                if (built_molecule.valid) {
                    ImGui::Separator();
                    ImGui::Text("Built Molecule:");
                    ImGui::BulletText("Formula: %s", built_molecule.formula.c_str());
                    ImGui::BulletText("Atoms: %d", built_molecule.num_atoms);
                    ImGui::BulletText("Bonds: %d", built_molecule.num_bonds);
                }
            }
        }
        ImGui::End();
    }
};

// Global instance
static MoleculeBuilder component;

} // namespace builder

#endif // VIAMD_ENABLE_BUILDER