/**
 * @file dambreak_sycl.cpp
 * @brief 2D dambreak example using SYCL.
 * @author Xiangyu Hu
 */
#include "sphinxsys_sycl.h"
using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;                    /**< Water tank length. */
Real DH = 5.366;                    /**< Water tank height. */
Real LL = 2.0;                      /**< Water column length. */
Real LH = 1.0;                      /**< Water column height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Reference density of fluid. */
Real gravity_g = 1.0;                    /**< Gravity. */
Real U_ref = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_ref;                 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;   // translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    TransformShape<GeometricShapeBox> initial_water_block(Transform(water_block_translation), water_block_halfsize, "WaterBody");
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
    fluid_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    UpdateCellLinkedList<execution::ParallelDevicePolicy, CellLinkedList> water_cell_linked_list(water_block);
    UpdateCellLinkedList<execution::ParallelDevicePolicy, CellLinkedList> wall_cell_linked_list(wall_boundary);
    DiscreteVariable<UnsignedInt> *dv_particle_index_water =
        water_block.getBaseParticles().getVariableByName<UnsignedInt>("ParticleIndex");
    DiscreteVariable<UnsignedInt> *dv_particle_index_wall =
        wall_boundary.getBaseParticles().getVariableByName<UnsignedInt>("ParticleIndex");

    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation fluid_observer_contact(fluid_observer, {&water_block});

    SequencedCombination<UpdateRelation<
        execution::ParallelDevicePolicy, BodyRelationUpdate<Inner<>, Contact<>>>>
        water_block_update_complex_relation(water_block_inner, water_wall_contact);
    UpdateRelation<execution::SequencedPolicy, BodyRelationUpdate<Contact<>>>
        fluid_observer_contact_relation(fluid_observer_contact);
    DiscreteVariable<UnsignedInt> *dv_particle_offset_water_inner = water_block_inner.getParticleOffset();
    StdVec<DiscreteVariable<UnsignedInt> *> dv_particle_offset_water_contact = water_wall_contact.getContactParticleOffset();
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_wall_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxillary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    StateDynamics<execution::ParallelDevicePolicy, GravityForceCK<Gravity>> constant_gravity(water_block, gravity);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_boundary_normal_direction(wall_boundary);
    StateDynamics<execution::ParallelDevicePolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<execution::ParallelDevicePolicy, fluid_dynamics::AdvectionStepClose> water_advection_step_close(water_block);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_wall_contact);

    InteractionDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AcousticStep1stHalfWithWallRiemannCK>
        fluid_acoustic_step_1st_half(water_block_inner, water_wall_contact);
    InteractionDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AcousticStep2ndHalfWithWallRiemannCK>
        fluid_acoustic_step_2nd_half(water_block_inner, water_wall_contact);

    InteractionDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::DensityRegularizationComplexFreeSurface>
        fluid_density_by_summation_ck(water_block_inner, water_wall_contact);
    DiscreteVariable<Real> *dv_density = water_block.getBaseParticles().getVariableByName<Real>("Density");

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> fluid_advection_time_step(water_block, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);

    ReduceDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AdvectionTimeStepCK> fluid_advection_time_step_ck(water_block, U_ref);
    ReduceDynamicsCK<execution::ParallelDevicePolicy, fluid_dynamics::AcousticTimeStepCK> fluid_acoustic_time_step_ck(water_block);
    DiscreteVariable<Vecd> *dv_force_prior = water_block.getBaseParticles().getVariableByName<Vecd>("ForcePrior");
    DiscreteVariable<Vecd> *dv_force = water_block.getBaseParticles().getVariableByName<Vecd>("Force");
    DiscreteVariable<Vecd> *dv_velocity = water_block.getBaseParticles().getVariableByName<Vecd>("Velocity");
    DiscreteVariable<Real> *dv_rho = water_block.getBaseParticles().getVariableByName<Real>("Density");
    DiscreteVariable<Real> *dv_mass = water_block.getBaseParticles().getVariableByName<Real>("Mass");
    DiscreteVariable<Real> *dv_p = water_block.getBaseParticles().getVariableByName<Real>("Pressure");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalMechanicalEnergy>> write_water_mechanical_energy(water_block, gravity);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_recorded_water_pressure("Pressure", fluid_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<execution::ParallelDevicePolicy, Real>>
        fluid_observer_pressure("Pressure", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();

    constant_gravity.exec();
    dv_force_prior->synchronizeWithDevice();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    dv_particle_index_water->synchronizeWithDevice();
    dv_particle_index_wall->synchronizeWithDevice();

    water_block_update_complex_relation.exec();
    dv_particle_offset_water_inner->synchronizeWithDevice();
    for (size_t k = 0; k != dv_particle_offset_water_contact.size(); ++k)
    {
        dv_particle_offset_water_contact[k]->synchronizeWithDevice();
    }

    fluid_density_by_summation_ck.exec();
    dv_density->synchronizeWithDevice();

    fluid_observer_contact_relation.exec();
    fluid_observer_pressure.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    SingularVariable<Real> *sv_physical_time = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    if (sph_system.RestartStep() != 0)
    {
        sv_physical_time->setValue(restart_io.readRestartFiles(sph_system.RestartStep()));
        water_block.updateCellLinkedList();
        water_wall_complex.updateConfiguration();
        fluid_observer_contact.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 20.0;
    Real output_interval = 0.1;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    write_recorded_water_pressure.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (sv_physical_time->getValue() < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();

            Real advection_dt = fluid_advection_time_step.exec();
            if (number_of_iterations % restart_output_interval == 0)
            {
                dv_force_prior->synchronizeToDevice();
                dv_force->synchronizeToDevice();
                dv_velocity->synchronizeToDevice();
                Real advection_dt_ck = fluid_advection_time_step_ck.exec();
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                          << "	advection_dt = " << advection_dt << "   advection_dt_ck = " << advection_dt_ck << "\n";
                if (ABS(advection_dt - advection_dt_ck) > 1.0e-6)
                {
                    std::cout << "Error: the advection time step is not consistent with the CK time step." << std::endl;
                    exit(1);
                }
            }

            fluid_density_by_summation.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                acoustic_dt = fluid_acoustic_time_step.exec();

                if (number_of_iterations % restart_output_interval == 0)
                {
                    dv_rho->synchronizeToDevice();
                    dv_mass->synchronizeToDevice();
                    dv_p->synchronizeToDevice();
                    dv_force_prior->synchronizeToDevice();
                    dv_force->synchronizeToDevice();
                    dv_velocity->synchronizeToDevice();
                    Real acoustic_dt_ck = fluid_acoustic_time_step_ck.exec();
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations
                              << "	acoustic_dt = " << acoustic_dt << "   acoustic_dt_ck = " << acoustic_dt_ck << "\n";
                    if (ABS(acoustic_dt - acoustic_dt_ck) > 1.0e-6)
                    {
                        std::cout << "Error: the acoustic time step is not consistent with the CK time step." << std::endl;
                        exit(1);
                    }
                }
                fluid_pressure_relaxation.exec(acoustic_dt);
                fluid_density_relaxation.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                sv_physical_time->incrementValue(acoustic_dt);
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body observables and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << sv_physical_time->getValue()
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_water_mechanical_energy.writeToFile(number_of_iterations);
                    write_recorded_water_pressure.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedList();
            water_wall_complex.updateConfiguration();
            fluid_observer_contact.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        body_states_recording.writeToFile();
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_water_mechanical_energy.generateDataBase(1.0e-3);
        write_recorded_water_pressure.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_water_mechanical_energy.testResult();
        write_recorded_water_pressure.testResult();
    }

    return 0;
};
