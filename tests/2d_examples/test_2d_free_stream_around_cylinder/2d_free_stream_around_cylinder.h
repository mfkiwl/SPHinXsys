/**
 * @file 	freestream_flow_around_cylinder.h
 * @brief 	This is the case file for the test of free-stream boundary condition.
 * @details  We consider a flow pass the cylinder with freestream boundary condition in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 30.0;                               /**< Channel length. */
Real DH = 16.0;                               /**< Channel height. */
Real particle_spacing_ref = 0.2;              /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose emitter. */
Real BW = 4.0 * particle_spacing_ref;         /**< Sponge region to impose injection. */
Vec2d insert_circle_center(10.0, 0.5 * DH);   /**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;              /**< Radius of the cylinder. */
// Observation locations
Vec2d point_coordinate_1(3.0, 5.0);
Vec2d point_coordinate_2(4.0, 5.0);
Vec2d point_coordinate_3(5.0, 5.0);
StdVec<Vecd> observation_locations = {point_coordinate_1, point_coordinate_2, point_coordinate_3};
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                            /**< Density. */
Real U_f = 1.0;                                               /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                                        /**< Speed of sound. */
Real Re = 100.0;                                              /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * insert_circle_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
//	water block shape
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(-DL_sponge, 0.0)};
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;

//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    WaterBlock() : MultiPolygonShape()
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
    }
};
class CylinderShape : public MultiPolygonShape
{
  public:
    CylinderShape() : MultiPolygonShape()
    {
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
struct FreeStreamVelocity
{
    Real u_ref_, t_ref_;

    template <class BoundaryConditionType>
    FreeStreamVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();
        target_velocity[0] = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        return target_velocity;
    }
};
