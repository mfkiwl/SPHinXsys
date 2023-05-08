/**
 * @file 	 stfb.h
 * @brief 	 This is the case file for 2D still floating body.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 10.0; 				/**< TOTAL SIMULATION TIME*/
Real DL = 3.0;									/**< Tank length. */
Real DH = 2.5;									/**< Tank height. */
Real WH = 2.0;									/**< Water block height. */
Real L = 1.0;									/**< Base of the floating body. */
Real particle_spacing_ref = L /20;	
Real BW = particle_spacing_ref * 4.0;			/**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-DL -BW, -DH -BW), Vec2d(DL + BW, DH + BW));
Vec2d offset = Vec2d::Zero();
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;					   	/**< Reference density of fluid. */
Real gravity_g = 9.81;				       	/**< Value of gravity. */
Real U_f = 2.0 * sqrt(0.79 * gravity_g); 	/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                     	/**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Structure Properties G and Inertia
//----------------------------------------------------------------------
/* Weight of the solid structure*/
Real StructureMass = 700;
/**< Area of the solid structure*/
Real FlStA =L*L;
/**< Density of the solid structure*/
Real rho_s=StructureMass/FlStA; 
/* Equilibrium position of the solid structure*/
Real H = -(rho_s/rho0_f*L-L/2);		/**< Strart placemnt of Flt Body*/

Real bcmx=0;
Real bcmy=H+0;
Vec2d G(bcmx,bcmy);
Real Ix = L*L*L*L/3;
Real Iy = L*L*L*L/3;
Real Iz = StructureMass/12*(L*L+L*L);
/**
 *  
 * Structure observer position
 * 
 * */
Vec2d obs = G;
//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------

MultiPolygon createFltStr()
{	
		/** Geometry definition. */
		std::vector<Vecd> sructure;
		sructure.push_back(Vec2d(-L/2,H-L/2));
		sructure.push_back(Vec2d(-L/2,H+L/2));
		sructure.push_back(Vec2d(L/2,H+L/2));
		sructure.push_back(Vec2d(L/2,H-L/2));
		sructure.push_back(Vec2d(-L/2,H-L/2));

		MultiPolygon multi_polygon_;

		multi_polygon_.addAPolygon(sructure, ShapeBooleanOps::add);

		return multi_polygon_;
}

class FloatingStructure : public MultiPolygonShape
{
public:
	explicit FloatingStructure(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
			/** Geometry definition. */
		std::vector<Vecd> sructure;
		sructure.push_back(Vec2d(-L/2,H-L/2));
		sructure.push_back(Vec2d(-L/2,H+L/2));
		sructure.push_back(Vec2d(L/2,H+L/2));
		sructure.push_back(Vec2d(L/2,H-L/2));
		sructure.push_back(Vec2d(-L/2,H-L/2));

		multi_polygon_.addAPolygon(sructure, ShapeBooleanOps::add);

	}
};

class StructureSystemForSimbody : public SolidBodyPartForSimbody
{
public:
	StructureSystemForSimbody(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
		: SolidBodyPartForSimbody(sph_body, shape_ptr)
	{
		//Vec2d mass_center(G[0], G[1]);
		//initial_mass_center_ = SimTK::Vec3(mass_center[0], mass_center[1], 0.0);
		body_part_mass_properties_ =
			mass_properties_ptr_keeper_
				.createPtr<SimTK::MassProperties>(StructureMass, SimTK::Vec3(0.0), SimTK::UnitInertia(Ix,Iy,Iz));
	}
};
//----------------------------------------------------------------------
//	Water block
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_block_shape;
		water_block_shape.push_back(Vec2d(-DL/2,-DH/2));
		water_block_shape.push_back(Vec2d(-DL/2,0));
		water_block_shape.push_back(Vec2d(DL/2,0));
		water_block_shape.push_back(Vec2d(DL/2,-DH/2));
		water_block_shape.push_back(Vec2d(-DL/2,-DH/2));

		/*Structure substract*/

		std::vector<Vecd> sructure;
		sructure.push_back(Vec2d(-L/2,H-L/2));
		sructure.push_back(Vec2d(-L/2,H+L/2));
		sructure.push_back(Vec2d(L/2,H+L/2));
		sructure.push_back(Vec2d(L/2,H-L/2));
		sructure.push_back(Vec2d(-L/2,H-L/2));

		multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(sructure, ShapeBooleanOps::sub);

	}
};
//----------------------------------------------------------------------
//	Wall cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vec2d(-DL/2, -DH/2)+Vec2d(-BW,-BW));
		outer_wall_shape.push_back(Vec2d(-DL/2, DH/2)+Vec2d(-BW,0));
		outer_wall_shape.push_back(Vec2d(DL/2, DH/2)+Vec2d(+BW,0));
		outer_wall_shape.push_back(Vec2d(DL/2, -DH/2)+Vec2d(+BW,-BW));
		outer_wall_shape.push_back(Vec2d(-DL/2, -DH/2)+Vec2d(-BW,-BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vec2d(-DL/2, -DH/2));
		inner_wall_shape.push_back(Vec2d(-DL/2, DH/2));
		inner_wall_shape.push_back(Vec2d(DL/2, DH/2));
		inner_wall_shape.push_back(Vec2d(DL/2, -DH/2));
		inner_wall_shape.push_back(Vec2d(-DL/2, -DH/2));


		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);

	}
};
//----------------------------------------------------------------------
//	create mesuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
MultiPolygon createFreeSurfaceGauge()
{	
		/** Geometry definition. */
		std::vector<Vecd> point;
		point.push_back(Vecd(DL/3-h,0));
		point.push_back(Vecd(DL/3-h,DH));
		point.push_back(Vecd(DL/3+h,DH));
		point.push_back(Vecd(DL/3+h,DH));
		point.push_back(Vecd(DL/3-h,0));

		MultiPolygon multi_polygon_;

		multi_polygon_.addAPolygon(point, ShapeBooleanOps::add);

		return multi_polygon_;
}