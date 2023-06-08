#include "external_force.h"

namespace SPH {
//=================================================================================================//
	ExternalForce::ExternalForce() {}
//=================================================================================================//
	Gravity::Gravity(Vecd global_acceleration, Vecd reference_position)
		: ExternalForce(), global_acceleration_(global_acceleration),
		zero_potential_reference_(reference_position) {}
	//=================================================================================================//
	Vecd Gravity::InducedAcceleration(Vecd& position)
	{
        return InducedAcceleration_Device(position);
	}
	//=================================================================================================//
	Real Gravity::getPotential(Vecd& position)
	{
		return InducedAcceleration(position).dot(zero_potential_reference_ - position);
	}
	//=================================================================================================//
}
//=================================================================================================//