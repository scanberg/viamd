#include <event.h>

#ifdef VIAMD_ENABLE_OPENMM

#include <viamd.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <md_molecule.h>

#include <imgui_widgets.h>
#include <imgui.h>

#include <OpenMM.h>

#include <memory>
#include <string>
#include <vector>
#include <cstring>

namespace openmm {

struct SimulationContext {
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::Integrator> integrator;
    
    // Force field parameters
    std::string force_field = "amber14-all.xml";
    std::string solvent = "amber14/tip3pfb.xml";
};

class OpenMMComponent : public viamd::EventHandler {
private:
    SimulationContext sim_context;
    md_allocator_i* allocator = nullptr;

public:
    OpenMMComponent() { 
        viamd::event_system_register_handler(*this); 
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                initialize(state);
                break;
            }
            case viamd::EventType_ViamdShutdown:
                shutdown();
                break;
            case viamd::EventType_ViamdFrameTick: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                update(state);
                draw_ui(state);
                break;
            }
            case viamd::EventType_ViamdTopologyInit: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                on_topology_init(state);
                break;
            }
            case viamd::EventType_ViamdTopologyFree:
                on_topology_free();
                break;
            case viamd::EventType_ViamdWindowDrawMenu: {
                draw_menu();
                break;
            }
            default:
                break;
            }
        }
    }

    void initialize(ApplicationState& state) {
        allocator = state.allocator.persistent;
        MD_LOG_INFO("OpenMM component initialized");
    }

    void shutdown() {
        // Note: cleanup_simulation needs an ApplicationState, so we can't clean up here
        // Cleanup will happen when topology is freed
        MD_LOG_INFO("OpenMM component shutdown");
    }

    void update(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (state.simulation.running && !state.simulation.paused && state.simulation.initialized) {
            // Run simulation steps
            run_simulation_step(state);
        }
#endif
    }

    void draw_ui(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (state.simulation.show_window) {
            draw_simulation_window(state);
        }
#endif
    }

    void draw_menu() {
#ifdef VIAMD_ENABLE_OPENMM
        if (ImGui::BeginMenu("Simulation")) {
            if (ImGui::MenuItem("OpenMM Simulation", nullptr, false)) {
                // This will be set via the ApplicationState
            }
            ImGui::EndMenu();
        }
#endif
    }

    void on_topology_init(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (state.mold.mol.atom.count > 0) {
            MD_LOG_INFO("Topology loaded, OpenMM system can be initialized");
            // Reset simulation state when new topology is loaded
            cleanup_simulation(state);
        }
#endif
    }

    void on_topology_free() {
        // Note: We can't call cleanup_simulation here without ApplicationState
        // The ApplicationState will handle clearing the simulation flags
        MD_LOG_INFO("Topology freed, OpenMM simulation stopped");
    }

private:
    void setup_system(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        try {
            if (state.mold.mol.atom.count == 0) {
                MD_LOG_ERROR("No atoms available for simulation setup");
                return;
            }

            // Create OpenMM system
            sim_context.system = std::make_unique<OpenMM::System>();
            
            // Add particles (atoms) to the system
            for (size_t i = 0; i < state.mold.mol.atom.count; ++i) {
                // Get mass from atomic number - simplified approach
                uint8_t atomic_number = get_atomic_number_from_atom_type(state.mold.mol.atom.type[i]);
                double mass = get_atomic_mass(atomic_number);
                sim_context.system->addParticle(mass);
            }
            
            // Set up basic force field (simplified)
            setup_force_field(state);
            
            // Create integrator
            sim_context.integrator = std::make_unique<OpenMM::LangevinIntegrator>(
                state.simulation.temperature, state.simulation.friction, state.simulation.timestep);
            
            // Create context
            sim_context.context = std::make_unique<OpenMM::Context>(
                *sim_context.system, *sim_context.integrator);
            
            // Set initial positions
            set_positions(state);
            
            // Perform energy minimization to prevent simulation explosion
            minimize_energy(state);
            
            state.simulation.initialized = true;
            state.simulation.current_frame = 0;
            state.simulation.simulation_time = 0.0;
            
            MD_LOG_INFO("OpenMM simulation system initialized with %zu atoms", 
                       state.mold.mol.atom.count);
                       
        } catch (const std::exception& e) {
            MD_LOG_ERROR("Failed to initialize OpenMM system: %s", e.what());
            cleanup_simulation(state);
        }
#endif
    }

    void setup_force_field(ApplicationState& state) {
        // Add harmonic bond forces (simplified)
        if (state.mold.mol.bond.count > 0) {
            auto* bondForce = new OpenMM::HarmonicBondForce();
            
            for (size_t i = 0; i < state.mold.mol.bond.count; ++i) {
                const auto& bond = state.mold.mol.bond.pairs[i];
                
                // Calculate actual bond length from current coordinates
                uint32_t atom1 = bond.idx[0];
                uint32_t atom2 = bond.idx[1];
                
                double dx = state.mold.mol.atom.x[atom1] - state.mold.mol.atom.x[atom2];
                double dy = state.mold.mol.atom.y[atom1] - state.mold.mol.atom.y[atom2];
                double dz = state.mold.mol.atom.z[atom1] - state.mold.mol.atom.z[atom2];
                
                double length_angstrom = sqrt(dx*dx + dy*dy + dz*dz);
                double length_nm = length_angstrom * 0.1; // Convert to nm
                
                // Use much softer force constant to prevent explosion
                // Reduced from 462750.4 to 30000 kJ/mol/nm^2 (about 10x softer)
                double k = 30000.0; // kJ/mol/nm^2 - softer bond strength
                
                // Sanity check: ensure reasonable bond length (0.5 - 3.0 Angstroms)
                if (length_angstrom < 0.5 || length_angstrom > 3.0) {
                    MD_LOG_ERROR("Unusual bond length detected: %.3f Å between atoms %u and %u", 
                               length_angstrom, atom1, atom2);
                    // Use default reasonable bond length if detected bond length is unrealistic
                    length_nm = 0.15; // 1.5 Angstroms in nm
                }
                
                bondForce->addBond(atom1, atom2, length_nm, k);
            }
            
            sim_context.system->addForce(bondForce);
        }
        
        // Add non-bonded forces (very simplified but safer)
        auto* nonbondedForce = new OpenMM::NonbondedForce();
        
        for (size_t i = 0; i < state.mold.mol.atom.count; ++i) {
            // Get element-specific parameters for better stability
            uint8_t atomic_number = get_atomic_number_from_atom_type(state.mold.mol.atom.type[i]);
            
            double charge = 0.0; // Neutral for simplicity - prevents electrostatic explosion
            double sigma, epsilon;
            
            // Element-specific van der Waals parameters (softer than before)
            switch (atomic_number) {
            case 1: // Hydrogen
                sigma = 0.24;  // nm - reduced from 0.3
                epsilon = 0.15; // kJ/mol - reduced from 0.5
                break;
            case 6: // Carbon  
                sigma = 0.34;
                epsilon = 0.35;
                break;
            case 7: // Nitrogen
                sigma = 0.32;
                epsilon = 0.30;
                break;
            case 8: // Oxygen
                sigma = 0.30;
                epsilon = 0.25;
                break;
            default:
                sigma = 0.32;  // Default values
                epsilon = 0.30;
                break;
            }
            
            nonbondedForce->addParticle(charge, sigma, epsilon);
        }
        
        nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
        nonbondedForce->setCutoffDistance(1.2); // nm - increased cutoff for smoother potential
        
        sim_context.system->addForce(nonbondedForce);
    }

    void set_positions(ApplicationState& state) {
        std::vector<OpenMM::Vec3> positions;
        positions.reserve(state.mold.mol.atom.count);
        
        for (size_t i = 0; i < state.mold.mol.atom.count; ++i) {
            // Convert from Angstroms to nanometers
            double x = state.mold.mol.atom.x[i] * 0.1;
            double y = state.mold.mol.atom.y[i] * 0.1;
            double z = state.mold.mol.atom.z[i] * 0.1;
            positions.emplace_back(x, y, z);
        }
        
        sim_context.context->setPositions(positions);
    }

    void minimize_energy(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (!sim_context.context) {
            return;
        }
        
        try {
            MD_LOG_INFO("Performing energy minimization to stabilize system...");
            
            // Perform local energy minimization to relax bad contacts
            // This is crucial to prevent simulation explosion
            OpenMM::LocalEnergyMinimizer::minimize(*sim_context.context, 1e-4, 1000);
            
            // Get minimized positions and update VIAMD coordinates
            OpenMM::State minimizedState = sim_context.context->getState(OpenMM::State::Positions | OpenMM::State::Energy);
            const std::vector<OpenMM::Vec3>& positions = minimizedState.getPositions();
            
            // Update VIAMD atom positions with minimized coordinates
            for (size_t i = 0; i < state.mold.mol.atom.count && i < positions.size(); ++i) {
                state.mold.mol.atom.x[i] = static_cast<float>(positions[i][0] * 10.0);
                state.mold.mol.atom.y[i] = static_cast<float>(positions[i][1] * 10.0);
                state.mold.mol.atom.z[i] = static_cast<float>(positions[i][2] * 10.0);
            }
            
            // Mark buffers as dirty for visualization update
            state.mold.dirty_buffers |= MolBit_DirtyPosition;
            
            double energy = minimizedState.getPotentialEnergy();
            MD_LOG_INFO("Energy minimization completed. Final energy: %.3f kJ/mol", energy);
            
        } catch (const std::exception& e) {
            MD_LOG_ERROR("Energy minimization failed: %s. Proceeding without minimization.", e.what());
        }
#endif
    }

    void run_simulation_step(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (!state.simulation.initialized || !sim_context.context) {
            return;
        }
        
        try {
            // Run simulation steps
            sim_context.integrator->step(state.simulation.steps_per_update);
            
            // Get updated positions and check for explosion
            OpenMM::State openmmState = sim_context.context->getState(OpenMM::State::Positions | OpenMM::State::Energy);
            const std::vector<OpenMM::Vec3>& positions = openmmState.getPositions();
            
            // Check for simulation explosion (coordinates too large)
            bool explosion_detected = false;
            double max_coord = 0.0;
            
            for (size_t i = 0; i < positions.size(); ++i) {
                double x = positions[i][0] * 10.0; // Convert to Angstroms
                double y = positions[i][1] * 10.0;
                double z = positions[i][2] * 10.0;
                
                double coord_magnitude = sqrt(x*x + y*y + z*z);
                max_coord = std::max(max_coord, coord_magnitude);
                
                // Check for explosion: coordinates > 100 Angstroms from origin
                if (coord_magnitude > 100.0) {
                    explosion_detected = true;
                    break;
                }
            }
            
            // Check for NaN or infinite values in energy
            double energy = openmmState.getPotentialEnergy();
            if (std::isnan(energy) || std::isinf(energy)) {
                explosion_detected = true;
            }
            
            if (explosion_detected) {
                MD_LOG_ERROR("Simulation explosion detected! Max coordinate: %.3f Å, Energy: %.3f kJ/mol", 
                            max_coord, energy);
                MD_LOG_ERROR("Stopping simulation to prevent further instability");
                state.simulation.running = false;
                state.simulation.paused = false;
                return;
            }
            
            // Update VIAMD atom positions (convert from nm to Angstroms)
            for (size_t i = 0; i < state.mold.mol.atom.count && i < positions.size(); ++i) {
                state.mold.mol.atom.x[i] = static_cast<float>(positions[i][0] * 10.0);
                state.mold.mol.atom.y[i] = static_cast<float>(positions[i][1] * 10.0);
                state.mold.mol.atom.z[i] = static_cast<float>(positions[i][2] * 10.0);
            }
            
            // Mark buffers as dirty for visualization update
            state.mold.dirty_buffers |= MolBit_DirtyPosition;
            
            // Update simulation state
            state.simulation.current_frame++;
            state.simulation.simulation_time += state.simulation.timestep * state.simulation.steps_per_update;
            
        } catch (const std::exception& e) {
            MD_LOG_ERROR("Simulation step failed: %s", e.what());
            state.simulation.running = false;
        }
#endif
    }

    void draw_simulation_window(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (!ImGui::Begin("OpenMM Simulation", &state.simulation.show_window)) {
            ImGui::End();
            return;
        }

        ImGui::Text("OpenMM Molecular Dynamics Simulation");
        ImGui::Separator();

        // System information
        if (state.mold.mol.atom.count > 0) {
            ImGui::Text("System: %zu atoms, %zu bonds", 
                       state.mold.mol.atom.count, state.mold.mol.bond.count);
        } else {
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "No molecular system loaded");
        }

        ImGui::Text("Status: %s", state.simulation.initialized ? "Initialized" : "Not initialized");
        
        if (state.simulation.initialized) {
            ImGui::Text("Frame: %d", state.simulation.current_frame);
            ImGui::Text("Time: %.3f ps", state.simulation.simulation_time);
        }

        ImGui::Separator();

        // Simulation parameters
        if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
            float temp = static_cast<float>(state.simulation.temperature);
            if (ImGui::SliderFloat("Temperature (K)", &temp, 250.0f, 400.0f)) {
                // Bounds checking for temperature
                temp = std::max(250.0f, std::min(400.0f, temp));
                state.simulation.temperature = temp;
            }
            
            float timestep = static_cast<float>(state.simulation.timestep);
            if (ImGui::SliderFloat("Timestep (ps)", &timestep, 0.0005f, 0.002f)) {
                // Much smaller timestep range to prevent instability
                // Reduced max from 0.005 to 0.002 ps
                timestep = std::max(0.0005f, std::min(0.002f, timestep));
                state.simulation.timestep = timestep;
            }
            
            float friction = static_cast<float>(state.simulation.friction);
            if (ImGui::SliderFloat("Friction (ps^-1)", &friction, 0.5f, 5.0f)) {
                // Constrain friction to reasonable range
                friction = std::max(0.5f, std::min(5.0f, friction));
                state.simulation.friction = friction;
            }
            
            ImGui::SliderInt("Steps per update", &state.simulation.steps_per_update, 1, 50);
            
            // Add helpful tooltips for safety
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Lower values = smoother animation, higher values = faster simulation");
            }
        }

        ImGui::Separator();

        // Control buttons
        if (!state.simulation.initialized) {
            if (ImGui::Button("Initialize System")) {
                if (state.mold.mol.atom.count > 0) {
                    setup_system(state);
                } else {
                    MD_LOG_ERROR("Load a molecular structure first");
                }
            }
        } else {
            if (!state.simulation.running) {
                if (ImGui::Button("Start Simulation")) {
                    state.simulation.running = true;
                    state.simulation.paused = false;
                }
            } else {
                if (!state.simulation.paused) {
                    if (ImGui::Button("Pause")) {
                        state.simulation.paused = true;
                    }
                } else {
                    if (ImGui::Button("Resume")) {
                        state.simulation.paused = false;
                    }
                }
                ImGui::SameLine();
                if (ImGui::Button("Stop")) {
                    state.simulation.running = false;
                    state.simulation.paused = false;
                }
            }

            ImGui::SameLine();
            if (ImGui::Button("Reset")) {
                reset_simulation(state);
            }
        }

        ImGui::End();
#endif
    }

    void reset_simulation(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        if (state.simulation.initialized) {
            state.simulation.running = false;
            state.simulation.paused = false;
            state.simulation.current_frame = 0;
            state.simulation.simulation_time = 0.0;
            
            // Reset positions to initial state
            set_positions(state);
            
            MD_LOG_INFO("Simulation reset");
        }
#endif
    }

    void cleanup_simulation(ApplicationState& state) {
#ifdef VIAMD_ENABLE_OPENMM
        state.simulation.running = false;
        state.simulation.paused = false;
        state.simulation.initialized = false;
        
        sim_context.context.reset();
        sim_context.integrator.reset();
        sim_context.system.reset();
#endif
    }

    uint8_t get_atomic_number_from_atom_type(const md_label_t& atom_type) {
        // Convert atom type label to atomic number
        const char* type_str = atom_type.buf;
        
        // Simple mapping based on common element symbols
        if (strncmp(type_str, "H", 1) == 0) return 1;
        if (strncmp(type_str, "C", 1) == 0) return 6;
        if (strncmp(type_str, "N", 1) == 0) return 7;
        if (strncmp(type_str, "O", 1) == 0) return 8;
        if (strncmp(type_str, "P", 1) == 0) return 15;
        if (strncmp(type_str, "S", 1) == 0) return 16;
        if (strncmp(type_str, "Cl", 2) == 0) return 17;
        if (strncmp(type_str, "Na", 2) == 0) return 11;
        if (strncmp(type_str, "Mg", 2) == 0) return 12;
        if (strncmp(type_str, "Ca", 2) == 0) return 20;
        if (strncmp(type_str, "Fe", 2) == 0) return 26;
        if (strncmp(type_str, "Zn", 2) == 0) return 30;
        
        // Default to carbon if unknown
        return 6;
    }

    double get_atomic_mass(uint8_t atomic_number) {
        // Simplified atomic masses in amu
        static const double masses[] = {
            0.0,    // 0 (dummy)
            1.008,  // H
            4.003,  // He
            6.941,  // Li
            9.012,  // Be
            10.811, // B
            12.011, // C
            14.007, // N
            15.999, // O
            18.998, // F
            20.180, // Ne
            22.990, // Na
            24.305, // Mg
            26.982, // Al
            28.086, // Si
            30.974, // P
            32.065, // S
            35.453, // Cl
            39.948, // Ar
        };
        
        if (atomic_number < sizeof(masses) / sizeof(masses[0])) {
            return masses[atomic_number];
        }
        return 12.011; // Default to carbon mass
    }
};

static OpenMMComponent g_openmm_component;

} // namespace openmm

#endif // VIAMD_ENABLE_OPENMM