/**
 * @file 	two_phase_dambreak_static_confinement.h
 * @brief 	Numerical parameters and body definition for 2D two-phase dambreak flow.
 * @author 	Yongchuan Yu and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 6;						  /**< Reference length. */
Real DH = 0.2;						  /**< Reference and the height of main channel. */
Real resolution_ref = 0.01;			  /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;		  /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
//-------------------------------------------------------
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + BW));
/** Observation locations*/
Real x_observe = 0.90 * DL;
Real x_observe_start = 0.90 * DL;
Real observe_spacing_x = 0.02 * DL;
int num_observer_points_x = 6;
int num_observer_points = 20;
Real observe_spacing = DH / num_observer_points;
StdVec<Vecd> observation_locations;
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;	   /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
//Real Re = 25000.0;					/**< Reynolds number. */
Real Re = 100.0;
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the emitter block. */
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d inlet_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize;

/** the water block . */
std::vector<Vecd> water_block_shape
{
	Vecd(-DL_sponge, 0.0),
	Vecd(-DL_sponge, DH),
	Vecd(DL, DH),
	Vecd(DL, 0.0),
	Vecd(-DL_sponge, 0.0)
};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape
{
	Vecd(-DL_sponge - 2.0 * BW, -BW), //1
	Vecd(-DL_sponge - 2.0 * BW, DH + BW), //2
	Vecd(DL + 2.0 * BW , DH + BW), //3
	Vecd(DL + 2.0 * BW , -BW), //4
	Vecd(-DL_sponge - 2.0 * BW, -BW), //1
};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape
{
	Vecd(-DL_sponge - 3.0 * BW, 0.0), //1
	Vecd(-DL_sponge - 3.0 * BW, DH), //2
	Vecd(DL + 3.0 * BW  , DH), //3
	Vecd(DL + 3.0 * BW , 0.0), //4
	Vecd(-DL_sponge - 3.0 * BW, 0.0), //1
};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};

class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
	Real u_ref_, t_ref_;
	AlignedBoxShape& aligned_box_;
	Vecd halfsize_;

	template <class BoundaryConditionType>
	InflowVelocity(BoundaryConditionType& boundary_condition)
		: u_ref_(U_f), t_ref_(2.0),
		aligned_box_(boundary_condition.getAlignedBox()),
		halfsize_(aligned_box_.HalfSize()) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = velocity;
		Real run_time = GlobalStaticVariables::physical_time_;
		Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		if (aligned_box_.checkInBounds(0, position))
		{
			target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
		}
		return target_velocity;
	}
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
	Real du_ave_dt_, u_ref_, t_ref_;

public:
	explicit TimeDependentAcceleration(Vecd gravity_vector)
		: Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_f), du_ave_dt_(0) {}

	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		Real run_time_ = GlobalStaticVariables::physical_time_;
		du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
		//if (position[1] <0.75*DH && position[1] > 0.25*DH)
		//{
		return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
		//}
		//else
		//{
		//	return global_acceleration_;
		//}
	}
};
//----------------------------------------------------------------------
//	Define turbulent inflow boundary condition
//----------------------------------------------------------------------
//class TurbulentEmitterBufferInflowCondition : public  fluid_dynamics::TurbulentInflowCondition
//{
//public:
//	TurbulentEmitterBufferInflowCondition(FluidBody& fluid_body, BodyAlignedBoxByCell& aligned_box_part)
//		: TurbulentInflowCondition(fluid_body, aligned_box_part) {}
//
//	Real getTurbulentInflowK(Vecd& position, Vecd& velocity, Real& turbu_k) override
//	{
//		Real u = velocity[0];
//		Real temp_in_turbu_k = 1.5 * pow((TurbulentIntensity * u), 2);
//		//Real temp_in_turbu_k = 1.5*pow((TurbulentIntensity * U_f), 2);
//		Real turbu_k_original = turbu_k;
//		//std::cout << "temp_in_turbu_k=" << temp_in_turbu_k << endl;
//		//std::cout << "turbu_k_original=" << turbu_k_original << endl;
//		if (position[0] < 0.0)
//		{
//			turbu_k_original = temp_in_turbu_k;
//		}
//		return turbu_k_original;
//	}
//	Real getTurbulentInflowE(Vecd& position, Real& turbu_k, Real& turbu_E) override
//	{
//		//Real temp_in_turbu_E = C_mu * pow(turbu_k, 1.5) / (0.1*getTurbulentLength());
//		Real temp_in_turbu_E = pow(turbu_k, 1.5) / (getTurbulentLength());
//		Real turbu_E_original = turbu_E;
//		//std::cout << "temp_in_turbu_k=" << temp_in_turbu_E << endl;
//		//std::cout << "turbu_k_original=" << turbu_E_original << endl;
//		if (position[0] < 0.0)
//		{
//			turbu_E_original = temp_in_turbu_E;
//		}
//		return turbu_E_original;
//	}
//	Real getTurbulentLength() override
//	{
//		return 0.07 * DH / pow(C_mu, 0.75); //Accoding to FLUNT Guide
//		//return 0.0045*DH;
//	}
//};
