/**
 * @file 	test_3d_incompressible_channel_flow.h
 * @brief 	This is a test to show inviscid incompressible channel flow using .msh files from ICEM and FLUENT.
 * @author 	Yash Mandaokar, Zhentong Wang and Xiangyu Hu
 */

#ifndef TEST_3D_INCOMPRESSIBLE_CHANNEL_FLOW_H
#define TEST_3D_INCOMPRESSIBLE_CHANNEL_FLOW_H
#include "unstructured_mesh_3d.h"             
#include "common_weakly_compressible_FVM_classes_3d.h"

using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1;                                        /**< Computation domain length. */
Real DH = 0.6494805454;                           /**< Computation domain height. */
Real DW = 0.038968832;                          /**< Computation domain width. */
Real particle_spacing_ref = 1.0 / 500.0;        /**< Initial reference particle spacing. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-0.3 / 0.769846, 0.0, 0.0), Vec3d(0.469846 / 0.769846, DH, DW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                        /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real P_ref = 101325;                                     /*operating pressure fluent [Pa]*/
Real gamma = 1.4;
//Real c_f = 10.0 * U_f;
Real c_f = sqrt(P_ref * gamma / (117.6655 * rho0_f)) * U_f;                           /**< Speed of sound. */
Real mu_f = 0.0;                                          /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_fullpath = "./input/Channel_ICEM.msh";

//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------

class WaveBody : public ComplexShape
{
public:
    explicit WaveBody(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wave(0.5 * DH, 0.5 * DL, 0.5 * DW);
        Transform translation_wave(halfsize_wave);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wave), halfsize_wave);
    }
};
///----------------------------------------------------------------------
//	Initialization
//----------------------------------------------------------------------
class InvCFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
    public:
    explicit InvCFInitialCondition(SPHBody& sph_body)
        : FluidInitialCondition(sph_body), rho_(particles_->rho_),
        p_(*particles_->getVariableByName<Real>("Pressure")), mom_(*particles_->getVariableByName<Vecd>("Momentum")),
        vel_(particles_->vel_) {};

protected:
    StdLargeVec<Real>& rho_, & p_;
    StdLargeVec<Vecd>& mom_;
    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 50 / 117.6655;
        vel_[index_i][0] = 1.0;
        vel_[index_i][1] = 0.0;
        vel_[index_i][2] = 0.0;
        
    }

protected:
    StdLargeVec<Vecd>& vel_;
};
///----------------------------------------------------------------------
//	InvCFBoundaryConditionSetup
//----------------------------------------------------------------------
class InvCFBoundaryConditionSetup : public BoundaryConditionSetupInFVM_3d
{
    public:
        InvCFBoundaryConditionSetup(BaseInnerRelationInFVM_3d & inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
            vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index)
            : BoundaryConditionSetupInFVM_3d(inner_relation, each_boundary_type_with_all_ghosts_index,
                each_boundary_type_with_all_ghosts_eij_, each_boundary_type_contact_real_index),
            fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())) {};
        virtual ~InvCFBoundaryConditionSetup() {};

        void applyReflectiveWallBoundary(size_t ghost_index, size_t index_i, Vecd e_ij) override
        {
            vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij)) + (-e_ij.dot(vel_[index_i]) * (e_ij));
            p_[ghost_index] = p_[index_i];
            rho_[ghost_index] = rho_[index_i];
        }
        void applyVelocityInletFlow(size_t ghost_index, size_t index_i) override
        {
            Vecd far_field_velocity(1.0, 0.0, 0.0);
            vel_[ghost_index] = far_field_velocity;
            p_[ghost_index] = p_[index_i];
            rho_[ghost_index] = rho_[index_i];
        }
        void applyPressureOutletBC(size_t ghost_index, size_t index_i) override
        {
            vel_[ghost_index] = vel_[index_i];
            p_[ghost_index] = 100.0 / 117.6655;
            rho_[ghost_index] = rho_[index_i];
        }
        void applysymmetry(size_t ghost_index, size_t index_i, Vecd e_ij) 
        {
            vel_[ghost_index] = (vel_[index_i] - 2 * e_ij.dot(vel_[index_i]) * e_ij);
            rho_[ghost_index] = rho_[index_i];
            p_[ghost_index] = p_[index_i];
        }

    protected:
        Fluid& fluid_;
};
#endif // TEST_3D_INCOMPRESSIBLE_CHANNEL_FLOW_H