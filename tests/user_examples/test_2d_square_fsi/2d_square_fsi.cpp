#include "case.h"
#include "sphinxsys.h"

using namespace SPH;

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** handle command line arguments. */
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    /** Set the starting time to zero. */
    GlobalStaticVariables::physical_time_ = 0.0;
    /** The water block, body, material and particles container. */
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    /** The wall boundary, body and particles container. */
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    /** The baffle, body and particles container. */
    SolidBody shell_baffle(sph_system, makeShared<DefaultShape>("ShellBaffle"));
    shell_baffle.defineParticlesAndMaterial<ShellParticles, LinearElasticSolid>(rho0_s, Youngs_modulus, poisson);
    shell_baffle.generateParticles<ShellBaffleParticleGenerator>();

    /** Particle and body creation of baffle observer.*/
    ObserverBody fluid_observer(sph_system, "FluidObserver");
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_locations);
    /** topology */
    InnerRelation water_inner(water_block);
    InnerRelation baffle_inner(shell_baffle);

    ContactRelation water_wall_contact(water_block, {&wall_boundary});

    ComplexRelation water_shell_complex(water_inner, {&shell_baffle});
    ContactRelation shell_water_contact(shell_baffle, {&water_block});

    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});
    /** Emitter setup. */
    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    /** Emitter buffer inflow condition. */
    BodyAlignedBoxByCell emitter_buffer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_buffer_translation)), emitter_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<FreeStreamVelocity>> emitter_buffer_inflow_condition(emitter_buffer);

    BodyAlignedBoxByCell disposer(water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, 0);
    /** Free-stream BC. */
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex> free_stream_surface_indicator(water_wall_contact, water_shell_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density(water_wall_contact, water_shell_complex);
    /** Algorithms for Fluid dynamics. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d(0)));
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    SimpleDynamics<fluid_dynamics::FreeStreamVelocityCorrection<FreeStreamVelocity>> velocity_boundary_condition_constraint(water_block);
    Dynamics1Level<fluid_dynamics::FluidShellandWallIntegration1stHalfRiemann> pressure_relaxation(water_wall_contact, water_shell_complex);
    pressure_relaxation.post_processes_.push_back(&velocity_boundary_condition_constraint);
    Dynamics1Level<fluid_dynamics::FluidShellandWallIntegration2ndHalfRiemann> density_relaxation(water_wall_contact, water_shell_complex);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithShellandWall> viscous_acceleration(water_wall_contact, water_shell_complex);
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_wall_contact, water_shell_complex);
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_inner);
    /** Algorithms for solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell_baffle);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(baffle_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(baffle_inner, 5, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(baffle_inner);
    /** Algorithms of FSI. */
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_shell(shell_water_contact);
    InteractionDynamics<solid_dynamics::FluidForceOnShellUpdate> fluid_force_on_shell_update(shell_water_contact, viscous_force_on_shell);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_baffle);
    /** constraint and damping */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    BoundaryGeometry shell_boundary_geometry(shell_baffle, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> baffle_constrain(shell_boundary_geometry);
    /** IO */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("PreviousSurfaceIndicator");
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
    shell_baffle.addBodyStateForRecording<Vecd>("PseudoNormal");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    RestartIO restart_io(io_environment, sph_system.real_bodies_);
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", io_environment, fluid_observer_contact);
    /** Prepare the simulation with cell linked list, configuration and case specified initial condition if necessary */
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    shell_corrected_configuration.exec();
    write_real_body_states.writeToFile();
    /** Time parameters. */
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;

    Real end_time = 400.0;
    Real output_interval = end_time / 200.0;
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    /** Main loop starts here. */
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            free_stream_surface_indicator.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_on_shell.exec();

            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                /** Fluid dynamics and force on solid. */
                pressure_relaxation.exec(dt);
                fluid_force_on_shell_update.exec();
                density_relaxation.exec(dt);

                /** Solid dynamics time stepping. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = shell_time_step_size.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    shell_stress_relaxation_first.exec(dt_s);
                    baffle_constrain.exec();
                    shell_stress_relaxation_second.exec(dt_s);
                    dt_s_sum += dt_s;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                emitter_buffer_inflow_condition.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            water_block.updateCellLinkedListWithParticleSort(100);
            shell_baffle.updateCellLinkedList();

            water_shell_complex.updateConfiguration();
            water_wall_contact.updateConfiguration();
            shell_water_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        write_real_body_states.writeToFile();
        fluid_observer_contact.updateConfiguration();
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}