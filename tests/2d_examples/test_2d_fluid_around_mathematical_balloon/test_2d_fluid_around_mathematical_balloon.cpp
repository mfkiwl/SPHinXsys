/**
 * @file 	test_2d_fluid_around_balloon_shell.cpp
 * @brief 	Test on fluid-shell interaction when 2 shell particles are close to each other
 * @details This is a case to test fluid-shell interaction.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    Kernel *kernel_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Real> W_ijV_j_ttl_contact;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, std::vector<BaseContactRelation *> &contact_relations)
        : particles_(&inner_relation.base_particles_),
          kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
          inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t n = 0; n < contact_relations.size(); n++)
        {
            auto &contact_relation = *contact_relations[n];
            for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
            {
                contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
                contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
            }
        }
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl, "TotalKernel");
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl_contact, "TotalKernelContact");
        inner_relation.base_particles_.registerVariable(dW_ijV_je_ij_ttl, "TotalKernelGrad");
        inner_relation.base_particles_.registerVariable(number_of_inner_neighbor, "InnerNeighborNumber");
        inner_relation.base_particles_.registerVariable(number_of_contact_neighbor, "ContactNeighborNumber");
    }

    inline void exec()
    {
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                int N_inner_number = 0;
                int N_contact_number = 0;
                Real W_ijV_j_ttl_i = particles_->Vol_[index_i] * kernel_->W(0, ZeroVecd);
                Vecd dW_ijV_je_ij_ttl_i = Vecd::Zero();
                const Neighborhood &inner_neighborhood = (*inner_configuration_)[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    W_ijV_j_ttl_i += inner_neighborhood.W_ij_[n] * particles_->Vol_[index_j];
                    dW_ijV_je_ij_ttl_i += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                    N_inner_number++;
                }

                double W_ijV_j_ttl_contact_i = 0;
                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    const SPH::Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        W_ijV_j_ttl_contact_i += wall_neighborhood.W_ij_[n] * contact_particles_[k]->Vol_[index_j];
                        dW_ijV_je_ij_ttl_i += wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
                        N_contact_number++;
                    }
                }
                W_ijV_j_ttl[index_i] = W_ijV_j_ttl_i + W_ijV_j_ttl_contact_i;
                W_ijV_j_ttl_contact[index_i] = W_ijV_j_ttl_contact_i;
                dW_ijV_je_ij_ttl[index_i] = dW_ijV_je_ij_ttl_i;
                number_of_inner_neighbor[index_i] = N_inner_number;
                number_of_contact_neighbor[index_i] = N_contact_number;
            });
    }
};

// Make sure the initial normal direction points from shell to fluid
class CorrectNormalDirection
{
  private:
    BaseParticles *particles_;
    Vecd direction_point_;
    bool if_reverse_;

  public:
    CorrectNormalDirection(SPHBody &shell, Vecd direction_point, bool if_reverse = false)
        : particles_(&shell.getBaseParticles()), direction_point_(direction_point), if_reverse_(if_reverse){};
    inline void exec()
    {
        auto &n = *particles_->getVariableByName<Vecd>("NormalDirection");
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                Vecd displacement = particles_->pos_[index_i] - direction_point_;
                if (if_reverse_)
                    displacement *= -1;
                if (n[index_i].dot(displacement) > 0)
                    n[index_i] *= -1;
            });
    }
};
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real scale = 0.001;
const Real DH = 12.0 * scale; /**< Reference and the height of main channel. */
const Real DL = 3.0 * DH;     /**< Reference length. */
const Real DL_balloon = 16.0 * scale;
const Real radius_balloon = 3.5 * scale; // radius of the mid surface

const Real resolution_ref = 0.6 * scale;
const Real resolution_shell = 0.3 * scale;  // outer surface radius
const Real thickness_balloon = 1.2 * scale; // thickness of the balloon

const Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
const Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */

const Real radius_balloon_outer = radius_balloon + 0.5 * resolution_shell; // radius of the outer surface
const Real radius_balloon_inner = radius_balloon - 0.5 * resolution_shell; // radius of the outer surface
const Real level_set_refinement_ratio = resolution_shell / (0.1 * resolution_shell);

const BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -0.5 * DH - BW), Vec2d(DL + BW, 0.5 * DH + BW));
const Vec2d balloon_center(0.5 * DL, 0);
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1056.0;               /**< Reference density of fluid. */
const Real mu_f = 3.5e-3;                 /**< Dynamics viscosity. */
const Real Re = 1000.0;                   /**< Reynolds number. */
const Real U_f = Re * mu_f / rho0_f / DH; /**< Characteristic velocity. */
const Real U_max = 1.5 * U_f * DH / (DH - 2 * radius_balloon_outer);
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
const Real c_f = 10.0 * U_max;
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const Real rho0_s = 1250.0; /**< Reference density.*/
const Real hardness = 50;   // Durometer hardnes: 50A
const Real youngs_modulus =
    std::pow(10, 0.0235 * hardness - 0.6403) * 1e3; // actual: 1e6, ref: https://doi.org/10.5254/1.3547752 eq. 12A
const Real poisson_ratio = 0.495;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * thickness_balloon;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_shape;
    water_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));
    water_shape.push_back(Vecd(-DL_sponge, 0.5 * DH));
    water_shape.push_back(Vecd(DL, 0.5 * DH));
    water_shape.push_back(Vecd(DL, -0.5 * DH));
    water_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));

    return water_shape;
}
/** create a wall outer shape */
std::vector<Vecd> createWallOuterShape()
{
    // geometry
    std::vector<Vecd> outer_shape;
    outer_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
    outer_shape.push_back(Vecd(-DL_sponge - BW, 0.5 * DH + BW));
    outer_shape.push_back(Vecd(DL + BW, 0.5 * DH + BW));
    outer_shape.push_back(Vecd(DL + BW, -0.5 * DH - BW));
    outer_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

    return outer_shape;
}
/** create a wall inner shape */
std::vector<Vecd> createWallInnerShape()
{
    // geometry
    std::vector<Vecd> inner_shape;
    inner_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
    inner_shape.push_back(Vecd(-DL_sponge - BW, 0.5 * DH));
    inner_shape.push_back(Vecd(DL + BW, 0.5 * DH));
    inner_shape.push_back(Vecd(DL + BW, -0.5 * DH));
    inner_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));

    return inner_shape;
}
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
const Vec2d circle_center_1 = balloon_center - Vec2d(0.5 * DL_balloon, 0);
const Vec2d circle_center_2 = balloon_center + Vec2d(0.5 * DL_balloon, 0);
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addABox(Transform(balloon_center), Vec2d(0.5 * DL_balloon, radius_balloon_outer), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_1, radius_balloon_outer, 100, ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_2, radius_balloon_outer, 100, ShapeBooleanOps::sub);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWallOuterShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWallInnerShape(), ShapeBooleanOps::sub);
    }
};
class Shell : public MultiPolygonShape
{
  public:
    explicit Shell(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(balloon_center), Vec2d(0.5 * DL_balloon, radius_balloon_outer), ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center_1, radius_balloon_outer, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center_2, radius_balloon_outer, 100, ShapeBooleanOps::add);
        multi_polygon_.addABox(Transform(balloon_center), Vec2d(0.5 * DL_balloon, radius_balloon_inner), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_1, radius_balloon_inner, 100, ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_2, radius_balloon_inner, 100, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(0.5),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        // Real run_time = GlobalStaticVariables::physical_time_;
        // Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        target_velocity[0] = 1.5 * u_ref_ * SMAX(0.0, 1.0 - 4 * position[1] * position[1] / DH / DH);
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define the boundary geometry
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry(){};

  private:
    void tagManually(size_t index_i)
    {
        if (std::abs(base_particles_.pos_[index_i][1]) < resolution_shell)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class ForceBoundaryGeometry : public BodyPartByParticle
{
  public:
    ForceBoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&ForceBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (fabs(base_particles_.pos_[index_i][0] - balloon_center.x()) < DL_balloon)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
const Real t_ref = 0.1;
const Real force_prior = 5.0;
class BalloonForce : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    explicit BalloonForce(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part),
          force_prior_(particles_->force_prior_),
          pos0_(particles_->pos0_),
          n_(particles_->n_){};

  protected:
    StdLargeVec<Vecd> &force_prior_;
    StdLargeVec<Vecd> &pos0_;
    StdLargeVec<Vecd> &n_;
    void update(size_t index_i, Real dt = 0.0)
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real x_ratio = (pos0_[index_i].x() - (balloon_center.x() - DL_balloon / 2.0)) / DL_balloon;
        Real force_avg = 0.5 * force_prior * (1 - cos(Pi * (2.0 * run_time / t_ref + x_ratio)));
        force_prior_[index_i] = -force_avg * n_[index_i];
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    std::cout << "U_max = " << U_max << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody shell(sph_system, makeShared<Shell>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    // shell.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        shell.generateParticles<ParticleGeneratorReload>(io_environment, shell.getName());
    }
    else
    {
        shell.defineBodyLevelSetShape(level_set_refinement_ratio);
        shell.generateParticles<ThickSurfaceParticleGeneratorLattice>(resolution_shell);
    }

    if (!sph_system.RunParticleRelaxation() && !sph_system.ReloadParticles())
    {
        std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
        return 0;
    }
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation shell_inner(shell);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for wall boundary.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> shell_random_particles(shell);
        relax_dynamics::ShellRelaxationStep relaxation_step_shell_inner(shell_inner);
        relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(shell_inner, resolution_shell, cos(Pi / 2.0));
        shell.addBodyStateForRecording<int>("UpdatedIndicator");
        CorrectNormalDirection correct_normal_direction(shell, balloon_center, true);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        MeshRecordingToPlt write_mesh_cell_linked_list(io_environment, shell.getCellLinkedList());
        ReloadParticleIO write_particle_reload_files(io_environment, {&shell});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        shell_random_particles.exec(0.25);

        relaxation_step_shell_inner.MidSurfaceBounding().exec();
        write_relaxed_particles.writeToFile(0);
        shell.updateCellLinkedList();
        write_mesh_cell_linked_list.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            for (int k = 0; k < 2; ++k)
                relaxation_step_shell_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_relaxed_particles.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        correct_normal_direction.exec();
        shell_normal_prediction.exec();
        correct_normal_direction.exec();

        shell.setNewlyUpdated();
        write_relaxed_particles.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    auto update_shell_thickness = [&]()
    {
        auto &thickness = *shell.getBaseParticles().getVariableByName<Real>("Thickness");
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                thickness[index_i] = thickness_balloon;
            });
    };
    update_shell_thickness();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    InnerRelation shell_inner(shell);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelationToShell water_shell_contact(water_block, {&shell}, true);
    ContactRelationFromShell shell_water_contact(shell, {&water_block}, true);
    ComplexRelation water_block_complex(water_inner, {&water_wall_contact, &water_shell_contact});
    // inner relation to compute curvature
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, water_block);
    // shell self contact
    ShellSelfContactRelation shell_self_contact(shell);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeStream>, Contact<>, Contact<>>> update_fluid_density_by_summation(water_inner, water_wall_contact, water_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> fluid_pressure_relaxation(water_inner, water_wall_contact, water_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, NoRiemannSolver>> fluid_density_relaxation(water_inner, water_wall_contact, water_shell_contact);
    InteractionDynamics<ComplexInteraction<fluid_dynamics::ViscousAcceleration<Inner<>, Contact<Wall>, Contact<Wall>>>> viscous_acceleration(water_inner, water_wall_contact, water_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>>>> inlet_outlet_surface_particle_indicator(water_inner, water_wall_contact, water_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SingleResolution>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_velocity_correction(water_inner, water_wall_contact, water_shell_contact);
    /** Algorithm for in-/outlet. */
    Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH + 2 * resolution_ref);
    Vec2d emitter_translation = Vec2d(-DL_sponge + 0.5 * BW, 0.0);
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);

    Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH + 2 * resolution_ref);
    Vec2d inlet_buffer_translation = Vec2d(-0.5 * DL_sponge - 2 * resolution_ref, 0.0);
    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(emitter_buffer);

    Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.6 * DH);
    Vec2d disposer_translation = Vec2d(DL, 0.0);
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, xAxis);
    /** Algorithm for solid dynamics. */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<TimeStepInitialization> shell_initialize_timestep(shell);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell);
    /** Algorithms for shell self contact. */
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> shell_curvature(shell_inner);
    InteractionDynamics<solid_dynamics::ShellSelfContactDensitySummation> shell_self_contact_density(shell_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> shell_self_contact_forces(shell_self_contact);
    auto update_shell_volume = [&]()
    {
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                shell.getBaseParticles().Vol_[index_i] = shell.getBaseParticles().mass_[index_i] / shell.getBaseParticles().rho_[index_i];
            });
    };
    /** FSI */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell);
    /** constraint and damping */
    BoundaryGeometry shell_boundary_geometry(shell, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> shell_constrain(shell_boundary_geometry);
    ForceBoundaryGeometry force_bc_geometry(shell, "ForceBcGeometry");
    SimpleDynamics<BalloonForce> balloon_force(force_bc_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<int>("Indicator");
    shell.addBodyStateForRecording<Real>("Thickness");
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("SelfContactDensity");
    shell.addBodyStateForRecording<Vecd>("PriorForce");
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    shell_corrected_configuration.exec();
    shell_curvature.compute_initial_curvature();
    // shell_average_curvature.exec();
    water_block_complex.updateConfiguration();
    shell_water_contact.updateConfiguration();

    //   Check dWijVjeij
    std::vector<BaseContactRelation *> contact_relations = {&water_shell_contact, &water_wall_contact};
    CheckKernelCompleteness check_kernel_completeness(water_inner, contact_relations);
    check_kernel_completeness.exec();
    water_block.addBodyStateForRecording<Real>("TotalKernel");
    water_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 1.0;
    Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    Real dt_s = 0.0;                         /**< Default acoustic time step sizes for solid. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    const double Dt_ref = fluid_advection_time_step.exec();
    const double dt_ref = fluid_acoustic_time_step.exec();
    const double dt_s_ref = shell_time_step_size.exec();

    // auto run_self_contact = [&]()
    // {
    //     while (GlobalStaticVariables::physical_time_ < end_time)
    //     {
    //         Real integration_time = 0.0;
    //         /** Integrate time (loop) until the next output time. */
    //         while (integration_time < output_interval)
    //         {
    //             shell_initialize_timestep.exec();

    //             balloon_force.exec();

    //             shell_self_contact_density.exec();
    //             shell_self_contact_forces.exec();

    //             Real dt_s = 0.5 * shell_time_step_size.exec();
    //             if (dt_s < dt_s_ref / 100)
    //             {
    //                 std::cout << "dt_s = " << dt_s << ", dt_s_ref = " << dt_s_ref << std::endl;
    //                 std::cout << "Shell time step decreased too much!" << std::endl;
    //                 throw std::runtime_error("Shell time step decreased too much!");
    //             }

    //             shell_stress_relaxation_first.exec(dt_s);
    //             shell_constrain.exec();
    //             shell_position_damping.exec(dt_s);
    //             shell_rotation_damping.exec(dt_s);
    //             shell_constrain.exec();
    //             shell_stress_relaxation_second.exec(dt_s);

    //             shell_update_normal.exec();
    //             update_shell_volume();
    //             // shell_curvature.exec();
    //             shell.updateCellLinkedList();
    //             shell_self_contact.updateConfiguration();

    //             number_of_iterations++;
    //             integration_time += dt_s;
    //             GlobalStaticVariables::physical_time_ += dt_s;

    //             if (number_of_iterations % screen_output_interval == 0)
    //             {
    //                 std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
    //                           << GlobalStaticVariables::physical_time_
    //                           << "  dt_s = " << dt_s << std::endl;
    //             }
    //         }

    //         TickCount t2 = TickCount::now();
    //         write_body_states.writeToFile();
    //         TickCount t3 = TickCount::now();
    //         interval += t3 - t2;
    //     }
    //     TickCount t4 = TickCount::now();

    //     TimeInterval tt;
    //     tt = t4 - t1 - interval;
    //     std::cout << "Total wall time for computation: " << tt.seconds()
    //               << " seconds." << std::endl;
    // };

    // run_self_contact();
    // exit(0);

    auto run_simulation = [&]()
    {
        std::cout << "Simulation starts here" << std::endl;
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                fluid_step_initialization.exec();
                Real Dt = fluid_advection_time_step.exec();
                if (Dt < Dt_ref / 20)
                {
                    std::cout << "Dt = " << Dt << ", Dt_ref = " << Dt_ref << std::endl;
                    std::cout << "Advective time step decreased too much!" << std::endl;
                    throw std::runtime_error("Advective time step decreased too much!");
                }
                inlet_outlet_surface_particle_indicator.exec();
                update_fluid_density_by_summation.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec();

                /** Dynamics including pressure relaxation. */
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt_temp = fluid_acoustic_time_step.exec();
                    if (dt_temp < dt_ref / 20)
                    {
                        std::cout << "dt = " << dt_temp << ", dt_ref = " << dt_ref << std::endl;
                        std::cout << "Acoustic time step decreased too much!" << std::endl;
                        throw std::runtime_error("Acoustic time step decreased too much!");
                    }
                    dt = SMIN(dt_temp, Dt - relaxation_time);
                    fluid_pressure_relaxation.exec(dt);
                    emitter_buffer_inflow_condition.exec();

                    // check_kernel_completeness.exec();
                    // water_block.setNewlyUpdated();
                    // shell.setNewlyUpdated();
                    // write_body_states.writeToFile();

                    fluid_density_relaxation.exec(dt);

                    /** Solid dynamics time stepping. */
                    Real dt_s_sum = 0.0;
                    average_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        shell_initialize_timestep.exec();

                        balloon_force.exec();

                        shell_self_contact_density.exec();
                        shell_self_contact_forces.exec();

                        Real dt_s_temp = shell_time_step_size.exec();
                        if (dt_s_temp < dt_s_ref / 100)
                        {
                            std::cout << "dt_s = " << dt_s_temp << ", dt_s_ref = " << dt_s_ref << std::endl;
                            std::cout << "Shell time step decreased too much!" << std::endl;
                            throw std::runtime_error("Shell time step decreased too much!");
                        }
                        dt_s = std::min(dt_s_temp, dt - dt_s_sum);

                        shell_stress_relaxation_first.exec(dt_s);
                        shell_constrain.exec();
                        shell_position_damping.exec(dt_s);
                        shell_rotation_damping.exec(dt_s);
                        shell_constrain.exec();
                        shell_stress_relaxation_second.exec(dt_s);

                        update_shell_volume();
                        shell_update_normal.exec();
                        shell_curvature.exec();
                        shell.updateCellLinkedList();
                        shell_self_contact.updateConfiguration();

                        dt_s_sum += dt_s;
                    }
                    average_velocity_and_acceleration.update_averages_.exec(dt);

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                }
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << "	Dt = " << Dt << "	dt = " << dt;
                    std::cout << "  dt_s = " << dt_s;
                    std::cout << "\n";
                }
                number_of_iterations++;

                /** inflow injection*/
                emitter_inflow_injection.exec();
                disposer_outflow_deletion.exec();

                /** Update cell linked list and configuration. */
                water_block.updateCellLinkedList();

                update_shell_volume();
                shell_update_normal.exec();
                shell.updateCellLinkedList();
                shell_curvature_inner.updateConfiguration();
                shell_average_curvature.exec();
                shell_water_contact.updateConfiguration();

                water_block_complex.updateConfiguration();
            }
            // exit(0);
            TickCount t2 = TickCount::now();
            check_kernel_completeness.exec();
            write_body_states.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds()
                  << " seconds." << std::endl;
    };

    try
    {
        run_simulation();
    }
    catch (const std::exception &e)
    {
        std::cout << "Error catched..." << std::endl;
        water_block.setNewlyUpdated();
        shell.setNewlyUpdated();
        write_body_states.writeToFile(1e8);
    }
    return 0;
}