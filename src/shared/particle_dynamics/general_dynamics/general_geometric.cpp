#include "general_geometric.h"
#include "base_particles.hpp"

namespace SPH
{
//=============================================================================================//
NormalDirectionFromBodyShape::NormalDirectionFromBodyShape(SPHBody &sph_body, Shape &shape)
    : LocalDynamics(sph_body), shape_(shape),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      n_(particles_->registerStateVariable<Vecd>("NormalDirection")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      phi_(particles_->registerStateVariable<Real>("SignedDistance")),
      phi0_(particles_->registerStateVariable<Real>("InitialSignedDistance")) {}
//=============================================================================================//
void NormalDirectionFromBodyShape::update(size_t index_i, Real dt)
{
    Vecd normal_direction = shape_.findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
    Real signed_distance = shape_.findSignedDistance(pos_[index_i]);
    phi_[index_i] = signed_distance;
    phi0_[index_i] = signed_distance;
}
//=============================================================================================//
NormalDirectionFromSubShapeAndOp::
    NormalDirectionFromSubShapeAndOp(SPHBody &sph_body, SubShapeAndOp *shape_and_op)
    : LocalDynamics(sph_body), shape_and_op_(shape_and_op), shape_(shape_and_op_->first),
      switch_sign_(shape_and_op_->second == ShapeBooleanOps::add ? 1.0 : -1.0),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      n_(particles_->registerStateVariable<Vecd>("NormalDirection")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      phi_(particles_->registerStateVariable<Real>("SignedDistance")),
      phi0_(particles_->registerStateVariable<Real>("InitialSignedDistance")) {}
//=============================================================================================//
void NormalDirectionFromSubShapeAndOp::update(size_t index_i, Real dt)
{
    Vecd normal_direction = switch_sign_ * shape_->findNormalDirection(pos_[index_i]);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
    Real signed_distance = switch_sign_ * shape_->findSignedDistance(pos_[index_i]);
    phi_[index_i] = signed_distance;
    phi0_[index_i] = signed_distance;
}
//=============================================================================================//
NormalDirectionFromParticles::
    NormalDirectionFromParticles(BaseInnerRelation &inner_relation, Shape &shape)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      shape_(shape), pos_(particles_->getVariableDataByName<Vecd>("Position")),
      n_(particles_->registerStateVariable<Vecd>("NormalDirection")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      phi_(particles_->registerStateVariable<Real>("SignedDistance")),
      phi0_(particles_->registerStateVariable<Real>("InitialSignedDistance")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
void NormalDirectionFromParticles::interaction(size_t index_i, Real dt)
{
    Vecd normal_direction = ZeroData<Vecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        normal_direction -= inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
    }
    normal_direction = normal_direction / (normal_direction.norm() + TinyReal);
    n_[index_i] = normal_direction;
    n0_[index_i] = normal_direction;
    Real signed_distance = shape_.findSignedDistance(pos_[index_i]);
    phi_[index_i] = signed_distance;
    phi0_[index_i] = signed_distance;
}
//=================================================================================================//
} // namespace SPH
