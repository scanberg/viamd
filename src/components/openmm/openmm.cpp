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
#include <map>
#include <tuple>
#include <algorithm>
#include <cmath>

namespace openmm {

// AMBER force field parameter structures
struct AmberAtomType {
    std::string name;
    double mass;
    double sigma;    // LJ sigma parameter (nm)
    double epsilon;  // LJ epsilon parameter (kJ/mol)
    double charge;   // Partial charge (default, can be overridden)
};

struct AmberBondType {
    double k;        // Force constant (kJ/mol/nm^2)
    double r0;       // Equilibrium length (nm)
};

struct AmberAngleType {
    double k;        // Force constant (kJ/mol/rad^2)
    double theta0;   // Equilibrium angle (radians)
};

struct SimulationContext {
    std::unique_ptr<OpenMM::System> system;
    std::unique_ptr<OpenMM::Context> context;
    std::unique_ptr<OpenMM::Integrator> integrator;
    
    // AMBER force field identifier
    std::string force_field_name = "AMBER14";
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

    void gradually_scale_charges(ApplicationState& state, double scale_factor) {
#ifdef VIAMD_ENABLE_OPENMM
        if (!sim_context.system || !state.simulation.initialized) {
            return;
        }
        
        try {
            // Find the NonbondedForce in the system
            for (int i = 0; i < sim_context.system->getNumForces(); ++i) {
                OpenMM::NonbondedForce* nbForce = dynamic_cast<OpenMM::NonbondedForce*>(&sim_context.system->getForce(i));
                if (nbForce) {
                    // Scale all charges
                    for (int j = 0; j < nbForce->getNumParticles(); ++j) {
                        double charge, sigma, epsilon;
                        nbForce->getParticleParameters(j, charge, sigma, epsilon);
                        
                        // Apply the scaling
                        charge *= scale_factor;
                        nbForce->setParticleParameters(j, charge, sigma, epsilon);
                    }
                    
                    // Update the context with new parameters
                    nbForce->updateParametersInContext(*sim_context.context);
                    MD_LOG_INFO("Scaled all charges by factor %.2f for improved stability", scale_factor);
                    break;
                }
            }
        } catch (const std::exception& e) {
            MD_LOG_ERROR("Failed to scale charges: %s", e.what());
        }
#endif
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
            
            // Add particles (atoms) to the system with AMBER masses
            std::vector<std::string> amber_types(state.mold.mol.atom.count);
            for (size_t i = 0; i < state.mold.mol.atom.count; ++i) {
                uint8_t atomic_number = get_atomic_number_from_atom_type(state.mold.mol.atom.type[i]);
                
                // Create a temporary connectivity for atom type assignment
                std::vector<uint32_t> temp_connectivity;
                amber_types[i] = map_to_amber_type(state.mold.mol.atom.type[i], atomic_number, temp_connectivity, state);
                
                // Get AMBER mass for this atom type
                AmberAtomType atom_params = get_amber_atom_params(amber_types[i]);
                double mass = atom_params.mass;
                
                sim_context.system->addParticle(mass);
                
                MD_LOG_DEBUG("Atom %zu: type=%s, AMBER_type=%s, mass=%.3f amu", 
                           i, state.mold.mol.atom.type[i].buf, amber_types[i].c_str(), mass);
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
        MD_LOG_INFO("Setting up AMBER force field...");
        
        // Step 1: Build connectivity information and assign AMBER atom types
        std::vector<std::string> amber_types(state.mold.mol.atom.count);
        std::vector<std::vector<uint32_t>> connectivity(state.mold.mol.atom.count);
        
        // Build connectivity matrix
        for (size_t i = 0; i < state.mold.mol.bond.count; ++i) {
            const auto& bond = state.mold.mol.bond.pairs[i];
            uint32_t atom1 = bond.idx[0];
            uint32_t atom2 = bond.idx[1];
            connectivity[atom1].push_back(atom2);
            connectivity[atom2].push_back(atom1);
        }
        
        // Assign AMBER atom types
        for (size_t i = 0; i < state.mold.mol.atom.count; ++i) {
            uint8_t atomic_number = get_atomic_number_from_atom_type(state.mold.mol.atom.type[i]);
            amber_types[i] = map_to_amber_type(state.mold.mol.atom.type[i], atomic_number, connectivity[i], state);
        }
        
        // Step 2: Set up bond forces with AMBER parameters
        if (state.mold.mol.bond.count > 0) {
            auto* bondForce = new OpenMM::HarmonicBondForce();
            
            for (size_t i = 0; i < state.mold.mol.bond.count; ++i) {
                const auto& bond = state.mold.mol.bond.pairs[i];
                uint32_t atom1 = bond.idx[0];
                uint32_t atom2 = bond.idx[1];
                
                // Get AMBER bond parameters
                AmberBondType bond_params = get_amber_bond_params(amber_types[atom1], amber_types[atom2]);
                
                // Calculate current bond length as starting point
                double dx = state.mold.mol.atom.x[atom1] - state.mold.mol.atom.x[atom2];
                double dy = state.mold.mol.atom.y[atom1] - state.mold.mol.atom.y[atom2];
                double dz = state.mold.mol.atom.z[atom1] - state.mold.mol.atom.z[atom2];
                double current_length = sqrt(dx*dx + dy*dy + dz*dz) * 0.1; // Convert to nm
                
                // Use AMBER equilibrium length, but if current length is reasonable, use it
                double equilibrium_length = bond_params.r0;
                if (current_length > 0.05 && current_length < 0.30) {  // Reasonable bond length range
                    equilibrium_length = current_length;
                }
                
                bondForce->addBond(atom1, atom2, equilibrium_length, bond_params.k);
                
                MD_LOG_DEBUG("Bond %zu-%zu: %s-%s, k=%.1f, r0=%.4f nm", 
                           atom1, atom2, amber_types[atom1].c_str(), amber_types[atom2].c_str(),
                           bond_params.k, equilibrium_length);
            }
            
            sim_context.system->addForce(bondForce);
            MD_LOG_INFO("Added %zu bond forces", state.mold.mol.bond.count);
        }
        
        // Step 3: Set up angle forces (missing from original implementation)
        auto* angleForce = new OpenMM::HarmonicAngleForce();
        size_t angle_count = 0;
        
        for (size_t center = 0; center < state.mold.mol.atom.count; ++center) {
            const auto& bonded = connectivity[center];
            
            // For each pair of atoms bonded to the center atom, create an angle
            for (size_t i = 0; i < bonded.size(); ++i) {
                for (size_t j = i + 1; j < bonded.size(); ++j) {
                    uint32_t atom1 = bonded[i];
                    uint32_t atom2 = static_cast<uint32_t>(center);
                    uint32_t atom3 = bonded[j];
                    
                    // Get AMBER angle parameters
                    AmberAngleType angle_params = get_amber_angle_params(
                        amber_types[atom1], amber_types[atom2], amber_types[atom3]);
                    
                    angleForce->addAngle(atom1, atom2, atom3, angle_params.theta0, angle_params.k);
                    angle_count++;
                    
                    MD_LOG_DEBUG("Angle %u-%u-%u: %s-%s-%s, k=%.1f, theta0=%.3f rad", 
                               atom1, atom2, atom3, 
                               amber_types[atom1].c_str(), amber_types[atom2].c_str(), amber_types[atom3].c_str(),
                               angle_params.k, angle_params.theta0);
                }
            }
        }
        
        if (angle_count > 0) {
            sim_context.system->addForce(angleForce);
            MD_LOG_INFO("Added %zu angle forces", angle_count);
        } else {
            delete angleForce;
        }
        
        // Step 4: Set up non-bonded forces with AMBER parameters
        auto* nonbondedForce = new OpenMM::NonbondedForce();
        
        for (size_t i = 0; i < state.mold.mol.atom.count; ++i) {
            AmberAtomType atom_params = get_amber_atom_params(amber_types[i]);
            
            // Use AMBER parameters for charge, sigma, epsilon
            double charge = atom_params.charge;
            double sigma = atom_params.sigma;
            double epsilon = atom_params.epsilon;
            
            // For initial stability, reduce charges by factor of 0.2
            // This is more conservative than 0.5 and prevents electrostatic explosions
            charge *= 0.2;
            
            nonbondedForce->addParticle(charge, sigma, epsilon);
            
            MD_LOG_DEBUG("Atom %zu (%s): charge=%.3f, sigma=%.3f nm, epsilon=%.3f kJ/mol", 
                       i, amber_types[i].c_str(), charge, sigma, epsilon);
        }
        
        nonbondedForce->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
        nonbondedForce->setCutoffDistance(1.0); // nm - standard cutoff for non-periodic
        
        sim_context.system->addForce(nonbondedForce);
        MD_LOG_INFO("Added non-bonded forces with AMBER parameters (charges scaled by 0.2 for stability)");
        
        // Step 5: Add hydrogen bond constraints for stability
        // This allows larger timesteps and prevents high-frequency H vibrations
        size_t constraint_count = 0;
        
        for (size_t i = 0; i < state.mold.mol.bond.count; ++i) {
            const auto& bond = state.mold.mol.bond.pairs[i];
            uint32_t atom1 = bond.idx[0];
            uint32_t atom2 = bond.idx[1];
            
            // Check if this bond involves hydrogen
            bool atom1_is_H = (amber_types[atom1][0] == 'H');
            bool atom2_is_H = (amber_types[atom2][0] == 'H');
            
            if (atom1_is_H || atom2_is_H) {
                // Calculate bond length for constraint
                double dx = state.mold.mol.atom.x[atom1] - state.mold.mol.atom.x[atom2];
                double dy = state.mold.mol.atom.y[atom1] - state.mold.mol.atom.y[atom2];
                double dz = state.mold.mol.atom.z[atom1] - state.mold.mol.atom.z[atom2];
                double bond_length = sqrt(dx*dx + dy*dy + dz*dz) * 0.1; // Convert to nm
                
                // Use reasonable constraint length based on AMBER parameters
                AmberBondType bond_params = get_amber_bond_params(amber_types[atom1], amber_types[atom2]);
                double constraint_length = bond_params.r0;
                
                // Ensure reasonable constraint length
                if (bond_length > 0.05 && bond_length < 0.25) {
                    constraint_length = bond_length;
                }
                
                sim_context.system->addConstraint(atom1, atom2, constraint_length);
                constraint_count++;
                
                MD_LOG_DEBUG("H-bond constraint %u-%u: %s-%s, length=%.4f nm", 
                           atom1, atom2, amber_types[atom1].c_str(), amber_types[atom2].c_str(), constraint_length);
            }
        }
        
        if (constraint_count > 0) {
            MD_LOG_INFO("Added %zu hydrogen bond constraints for stability", constraint_count);
        }
        
        MD_LOG_INFO("AMBER force field setup complete: %zu atoms, %zu bonds, %zu angles, %zu constraints", 
                   state.mold.mol.atom.count, state.mold.mol.bond.count, angle_count, constraint_count);
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
            MD_LOG_INFO("Performing thorough energy minimization to stabilize AMBER 14 system...");
            
            // Step 1: Initial coarse minimization with higher tolerance
            // This quickly removes the worst contacts
            OpenMM::LocalEnergyMinimizer::minimize(*sim_context.context, 1e-2, 500);
            
            // Step 2: Fine minimization with stricter tolerance  
            // This ensures the system is well-relaxed before dynamics
            OpenMM::LocalEnergyMinimizer::minimize(*sim_context.context, 1e-6, 2000);
            
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
            
            // Initialize velocities from Maxwell-Boltzmann distribution
            // This is crucial for stable dynamics
            sim_context.context->setVelocitiesToTemperature(state.simulation.temperature);
            MD_LOG_INFO("Velocities initialized to %.1f K temperature", state.simulation.temperature);
            
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
            OpenMM::State openmmState = sim_context.context->getState(OpenMM::State::Positions | OpenMM::State::Energy | OpenMM::State::Forces);
            const std::vector<OpenMM::Vec3>& positions = openmmState.getPositions();
            const std::vector<OpenMM::Vec3>& forces = openmmState.getForces();
            
            // Enhanced stability monitoring for AMBER 14
            bool explosion_detected = false;
            double max_coord = 0.0;
            double max_force = 0.0;
            
            for (size_t i = 0; i < positions.size(); ++i) {
                double x = positions[i][0] * 10.0; // Convert to Angstroms
                double y = positions[i][1] * 10.0;
                double z = positions[i][2] * 10.0;
                
                double coord_magnitude = sqrt(x*x + y*y + z*z);
                max_coord = std::max(max_coord, coord_magnitude);
                
                // Check for explosion: coordinates > 50 Angstroms from origin (reduced threshold)
                if (coord_magnitude > 50.0) {
                    explosion_detected = true;
                    break;
                }
                
                // Check force magnitude for early explosion detection
                if (i < forces.size()) {
                    double fx = forces[i][0];
                    double fy = forces[i][1]; 
                    double fz = forces[i][2];
                    double force_magnitude = sqrt(fx*fx + fy*fy + fz*fz);
                    max_force = std::max(max_force, force_magnitude);
                    
                    // Check for excessive forces (> 10000 kJ/mol/nm)
                    if (force_magnitude > 10000.0) {
                        explosion_detected = true;
                        MD_LOG_ERROR("Excessive force detected on atom %zu: %.1f kJ/mol/nm", i, force_magnitude);
                        break;
                    }
                }
            }
            
            // Check for NaN or infinite values in energy
            double energy = openmmState.getPotentialEnergy();
            if (std::isnan(energy) || std::isinf(energy) || energy > 1e6) {
                explosion_detected = true;
            }
            
            if (explosion_detected) {
                MD_LOG_ERROR("AMBER 14 simulation explosion detected!");
                MD_LOG_ERROR("Max coordinate: %.3f Å, Max force: %.1f kJ/mol/nm, Energy: %.3f kJ/mol", 
                            max_coord, max_force, energy);
                MD_LOG_ERROR("Stopping simulation to prevent further instability");
                MD_LOG_INFO("Try: 1) Load a better minimized structure, 2) Reduce timestep, 3) Lower temperature");
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
        ImGui::Text("Force Field: %s", sim_context.force_field_name.c_str());
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
            
            // Show stability information
            ImGui::Text("Stability: Charges scaled to %.1f%%, H-bonds constrained", 
                       20.0);  // Initial 0.2 scaling = 20%
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
            if (ImGui::SliderFloat("Timestep (ps)", &timestep, 0.0001f, 0.001f)) {
                // Conservative timestep range for AMBER 14 stability
                // Reduced max from 0.002 to 0.001 ps for better stability
                timestep = std::max(0.0001f, std::min(0.001f, timestep));
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
                    // Safety checks before starting AMBER 14 simulation
                    if (state.simulation.temperature > 500.0) {
                        MD_LOG_INFO("Temperature %.1f K is very high for AMBER 14! Consider using lower temperature for stability.", 
                                   state.simulation.temperature);
                    }
                    if (state.simulation.timestep > 0.001) {
                        MD_LOG_INFO("Timestep %.4f ps is large for AMBER 14! Consider using smaller timestep for stability.", 
                                   state.simulation.timestep);
                    }
                    
                    state.simulation.running = true;
                    state.simulation.paused = false;
                    MD_LOG_INFO("Starting AMBER 14 simulation: T=%.1f K, dt=%.4f ps", 
                               state.simulation.temperature, state.simulation.timestep);
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
            
            // Advanced stability controls
            if (ImGui::CollapsingHeader("Stability Controls")) {
                ImGui::Text("AMBER 14 stability enhancement tools:");
                
                if (ImGui::Button("Scale Charges +20%")) {
                    gradually_scale_charges(state, 1.2);
                }
                ImGui::SameLine();
                if (ImGui::Button("Scale Charges -20%")) {
                    gradually_scale_charges(state, 0.8);
                }
                
                if (ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Gradually adjust electrostatic interactions for stability.\n"
                                    "Reduce charges if simulation is unstable.\n"
                                    "Increase charges gradually once system is stable.");
                }
                
                if (ImGui::Button("Re-minimize Energy")) {
                    if (state.simulation.initialized) {
                        bool was_running = state.simulation.running;
                        state.simulation.running = false;
                        minimize_energy(state);
                        state.simulation.running = was_running;
                    }
                }
                if (ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Perform additional energy minimization if system becomes unstable");
                }
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

    std::string map_to_amber_type(const md_label_t& atom_type, uint8_t atomic_number, const std::vector<uint32_t>& bonded_atoms, ApplicationState& state) {
        // Convert VIAMD atom type to AMBER atom type based on chemical environment
        const char* type_str = atom_type.buf;
        
        // First try direct mapping for common AMBER types
        std::string type_upper = type_str;
        std::transform(type_upper.begin(), type_upper.end(), type_upper.begin(), ::toupper);
        
        // Direct AMBER type mapping - check if this is already a known AMBER type
        AmberAtomType test_params = get_amber_atom_params(type_str);
        if (test_params.name == type_str) {  // Found exact match
            return type_str;
        }
        
        // Element-based mapping with chemical environment consideration
        switch (atomic_number) {
        case 1: // Hydrogen
            return "H";  // Generic hydrogen, could be refined based on bonding
        case 6: // Carbon
            if (bonded_atoms.size() == 4) return "CA";  // sp3 carbon (approximation)
            if (bonded_atoms.size() == 3) return "C";   // sp2 carbon (approximation)
            return "C*";  // Generic carbon
        case 7: // Nitrogen
            return "N";   // Generic nitrogen
        case 8: // Oxygen
            if (bonded_atoms.size() == 1) return "O";   // Carbonyl oxygen
            if (bonded_atoms.size() == 2) return "OH";  // Hydroxyl oxygen
            return "O*";  // Generic oxygen
        case 16: // Sulfur
            return "SH";  // Sulfur (thiol)
        default:
            // Create generic type based on element
            switch (atomic_number) {
            case 1: return "H*";
            case 6: return "C*";
            case 7: return "N*";
            case 8: return "O*";
            default: return "C*";  // Fallback to carbon
            }
        }
    }

    AmberAtomType get_amber_atom_params(const std::string& amber_type) {
        // Enhanced AMBER14 force field parameters for better stability
        static const std::map<std::string, AmberAtomType> amber_atom_types = {
            // Protein backbone
            {"C",   {"C", 12.011, 0.339967, 0.359824, 0.5973}},     // Carbonyl carbon
            {"CA",  {"CA", 12.011, 0.339967, 0.359824, 0.0337}},    // Alpha carbon  
            {"CB",  {"CB", 12.011, 0.339967, 0.359824, -0.0875}},   // Beta carbon
            {"CT",  {"CT", 12.011, 0.339967, 0.4577, 0.0000}},      // Aliphatic carbon
            {"N",   {"N", 14.007, 0.325000, 0.711280, -0.4157}},    // Backbone nitrogen
            {"N3",  {"N3", 14.007, 0.325000, 0.711280, -0.3000}},   // Amino nitrogen
            {"O",   {"O", 15.999, 0.296000, 0.878640, -0.5679}},    // Carbonyl oxygen
            {"O2",  {"O2", 15.999, 0.296000, 0.878640, -0.8000}},   // Carboxyl oxygen
            {"H",   {"H", 1.008,  0.247135, 0.065270, 0.2719}},     // Backbone hydrogen
            {"HA",  {"HA", 1.008,  0.264953, 0.065270, 0.0337}},    // Alpha hydrogen
            {"HB",  {"HB", 1.008,  0.264953, 0.065270, 0.0295}},    // Beta hydrogen
            {"HC",  {"HC", 1.008,  0.264953, 0.065270, 0.0000}},    // Aliphatic hydrogen
            
            // Common side chains  
            {"OH",  {"OH", 15.999, 0.306647, 0.880314, -0.6546}},   // Hydroxyl oxygen
            {"HO",  {"HO", 1.008,  0.100000, 0.046000, 0.4275}},    // Hydroxyl hydrogen (fixed parameters)
            {"SH",  {"SH", 32.065, 0.356359, 1.046000, -0.3119}},   // Sulfur
            {"HS",  {"HS", 1.008,  0.106908, 0.065270, 0.1933}},    // Sulfur hydrogen
            
            // Aromatic carbons for better coverage
            {"CW",  {"CW", 12.011, 0.339967, 0.359824, -0.0275}},   // Aromatic carbon
            {"CC",  {"CC", 12.011, 0.339967, 0.359824, 0.0000}},    // Aromatic carbon
            {"HP",  {"HP", 1.008,  0.247135, 0.065270, 0.1000}},    // Aromatic hydrogen
            
            // Generic fallbacks with improved parameters
            {"C*",  {"C*", 12.011, 0.339967, 0.359824, 0.0000}},    // Generic carbon
            {"N*",  {"N*", 14.007, 0.325000, 0.711280, 0.0000}},    // Generic nitrogen  
            {"O*",  {"O*", 15.999, 0.296000, 0.878640, 0.0000}},    // Generic oxygen
            {"H*",  {"H*", 1.008,  0.247135, 0.065270, 0.0000}},    // Generic hydrogen
            {"S*",  {"S*", 32.065, 0.356359, 1.046000, 0.0000}},    // Generic sulfur
        };
        
        auto it = amber_atom_types.find(amber_type);
        if (it != amber_atom_types.end()) {
            return it->second;
        }
        
        // Fallback to generic parameters based on first character
        if (amber_type[0] == 'H') return amber_atom_types.at("H*");
        if (amber_type[0] == 'C') return amber_atom_types.at("C*");
        if (amber_type[0] == 'N') return amber_atom_types.at("N*");
        if (amber_type[0] == 'O') return amber_atom_types.at("O*");
        if (amber_type[0] == 'S') return amber_atom_types.at("S*");
        
        return amber_atom_types.at("C*");  // Ultimate fallback
    }

    AmberBondType get_amber_bond_params(const std::string& type1, const std::string& type2) {
        // Enhanced AMBER bond parameters (kJ/mol/nm^2, nm) for better stability
        static const std::map<std::pair<std::string, std::string>, AmberBondType> amber_bond_types = {
            // Protein backbone bonds
            {{"C", "O"},   {502080.0, 0.1229}},   // C=O bond
            {{"C", "N"},   {351456.0, 0.1335}},   // C-N amide bond  
            {{"C", "CA"},  {317984.0, 0.1522}},   // C-C bond
            {{"CA", "N"},  {337648.0, 0.1449}},   // CA-N bond
            {{"CA", "HA"}, {307105.6, 0.1090}},   // C-H bond
            {{"CA", "CB"}, {317984.0, 0.1526}},   // CA-CB bond
            {{"N", "H"},   {363171.2, 0.1010}},   // N-H bond
            
            // Side chain bonds
            {{"OH", "HO"}, {462750.4, 0.0974}},   // O-H bond
            {{"SH", "HS"}, {274887.2, 0.1336}},   // S-H bond
            {{"CT", "HC"}, {307105.6, 0.1090}},   // Aliphatic C-H
            {{"CT", "CT"}, {317984.0, 0.1526}},   // Aliphatic C-C
            {{"N3", "H"},  {363171.2, 0.1010}},   // Amino N-H
            {{"O2", "C"},  {469870.4, 0.1250}},   // Carboxyl C-O
            
            // Aromatic bonds
            {{"CW", "HP"}, {307105.6, 0.1080}},   // Aromatic C-H
            {{"CC", "HP"}, {307105.6, 0.1080}},   // Aromatic C-H
            {{"CW", "CW"}, {418400.0, 0.1371}},   // Aromatic C-C
            {{"CC", "CC"}, {418400.0, 0.1371}},   // Aromatic C-C
            
            // Generic fallback bonds (more conservative force constants)
            {{"C*", "C*"}, {250000.0, 0.1540}},   // Generic C-C (reduced force)
            {{"C*", "H*"}, {284512.0, 0.1090}},   // Generic C-H (reduced force)
            {{"N*", "H*"}, {320000.0, 0.1010}},   // Generic N-H (reduced force)
            {{"O*", "H*"}, {400000.0, 0.0960}},   // Generic O-H (reduced force)
            {{"S*", "H*"}, {250000.0, 0.1336}},   // Generic S-H
            {{"S*", "C*"}, {200000.0, 0.1810}},   // Generic S-C
        };
        
        // Try direct lookup
        auto key1 = std::make_pair(type1, type2);
        auto key2 = std::make_pair(type2, type1);  // Try reverse order
        
        auto it = amber_bond_types.find(key1);
        if (it != amber_bond_types.end()) {
            return it->second;
        }
        
        it = amber_bond_types.find(key2);
        if (it != amber_bond_types.end()) {
            return it->second;
        }
        
        // Fallback to generic types
        std::string gen1 = std::string(1, type1[0]) + "*";
        std::string gen2 = std::string(1, type2[0]) + "*";
        
        key1 = std::make_pair(gen1, gen2);
        key2 = std::make_pair(gen2, gen1);
        
        it = amber_bond_types.find(key1);
        if (it != amber_bond_types.end()) {
            return it->second;
        }
        
        it = amber_bond_types.find(key2);
        if (it != amber_bond_types.end()) {
            return it->second;
        }
        
        // Ultimate fallback (more conservative parameters)
        return {250000.0, 0.1540};  // Conservative C-C bond parameters
    }

    AmberAngleType get_amber_angle_params(const std::string& type1, const std::string& type2, const std::string& type3) {
        // Enhanced AMBER angle parameters (kJ/mol/rad^2, radians) for better stability
        static const std::map<std::tuple<std::string, std::string, std::string>, AmberAngleType> amber_angle_types = {
            // Protein backbone angles
            {std::make_tuple("N", "CA", "C"),   {527.184, 1.9373}},   // 111.0 degrees
            {std::make_tuple("CA", "C", "O"),   {568.518, 2.0944}},   // 120.0 degrees
            {std::make_tuple("CA", "C", "N"),   {585.760, 2.0246}},   // 116.0 degrees
            {std::make_tuple("C", "N", "CA"),   {418.400, 2.0246}},   // 116.0 degrees
            {std::make_tuple("H", "N", "CA"),   {418.400, 2.0595}},   // 118.0 degrees
            {std::make_tuple("HA", "CA", "N"),  {418.400, 1.9373}},   // 111.0 degrees
            {std::make_tuple("HA", "CA", "C"),  {418.400, 1.9373}},   // 111.0 degrees
            {std::make_tuple("CB", "CA", "N"),  {527.184, 1.9373}},   // 111.0 degrees
            {std::make_tuple("CB", "CA", "C"),  {527.184, 1.9373}},   // 111.0 degrees
            
            // Side chain angles
            {std::make_tuple("CA", "CB", "HC"), {418.400, 1.9373}},   // 111.0 degrees
            {std::make_tuple("HC", "CT", "HC"), {276.144, 1.8762}},   // 107.8 degrees (tetrahedral)
            {std::make_tuple("CT", "CT", "HC"), {418.400, 1.9373}},   // 111.0 degrees
            {std::make_tuple("N3", "CT", "HC"), {418.400, 1.9373}},   // 111.0 degrees
            {std::make_tuple("H", "N3", "CT"),  {418.400, 1.9373}},   // 111.0 degrees
            {std::make_tuple("HO", "OH", "CT"), {460.240, 1.8849}},   // 107.8 degrees
            {std::make_tuple("HS", "SH", "CT"), {368.192, 1.6232}},   // 93.0 degrees
            
            // More conservative generic fallback angles (reduced force constants)
            {std::make_tuple("C*", "C*", "C*"), {400.000, 1.9373}},   // Generic C-C-C (reduced)
            {std::make_tuple("C*", "C*", "H*"), {350.000, 1.9373}},   // Generic C-C-H (reduced)
            {std::make_tuple("H*", "C*", "H*"), {250.000, 1.8762}},   // Generic H-C-H (reduced)
            {std::make_tuple("N*", "C*", "C*"), {400.000, 1.9373}},   // Generic N-C-C
            {std::make_tuple("O*", "C*", "C*"), {400.000, 1.9373}},   // Generic O-C-C
            {std::make_tuple("S*", "C*", "C*"), {350.000, 1.9373}},   // Generic S-C-C
        };
        
        auto key = std::make_tuple(type1, type2, type3);
        auto it = amber_angle_types.find(key);
        if (it != amber_angle_types.end()) {
            return it->second;
        }
        
        // Try reverse order
        key = std::make_tuple(type3, type2, type1);
        it = amber_angle_types.find(key);
        if (it != amber_angle_types.end()) {
            return it->second;
        }
        
        // Fallback to generic types
        std::string gen1 = std::string(1, type1[0]) + "*";
        std::string gen2 = std::string(1, type2[0]) + "*";
        std::string gen3 = std::string(1, type3[0]) + "*";
        
        key = std::make_tuple(gen1, gen2, gen3);
        it = amber_angle_types.find(key);
        if (it != amber_angle_types.end()) {
            return it->second;
        }
        
        key = std::make_tuple(gen3, gen2, gen1);
        it = amber_angle_types.find(key);
        if (it != amber_angle_types.end()) {
            return it->second;
        }
        
        // Ultimate fallback (more conservative parameters)
        return {400.000, 1.9373};  // Conservative tetrahedral angle
    }
};

static OpenMMComponent g_openmm_component;

} // namespace openmm

#endif // VIAMD_ENABLE_OPENMM