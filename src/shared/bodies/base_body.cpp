#include "base_body.h"

#include "base_body_relation.h"
#include "base_particles.hpp"
#include "sph_system.h"

namespace SPH
{
//=================================================================================================//
SPHBody::SPHBody(SPHSystem &sph_system, const std::string &name)
    : sph_system_(sph_system), body_name_(name), newly_updated_(true),
      base_particles_(nullptr), is_bound_set_(false), total_body_parts_(0),
      sph_adaptation_(sph_adaptation_ptr_keeper_.createPtr<SPHAdaptation>(sph_system.ReferenceResolution())),
      base_material_(base_material_ptr_keeper_.createPtr<BaseMaterial>())
{
    sph_system_.sph_bodies_.push_back(this);
}
//=================================================================================================//
BoundingBox SPHBody::getSPHSystemBounds()
{
    return sph_system_.system_domain_bounds_;
}
//=================================================================================================//
int SPHBody::getNewBodyPartID()
{
    total_body_parts_++;
    return total_body_parts_;
};
//=================================================================================================//
SPHSystem &SPHBody::getSPHSystem()
{
    return sph_system_;
}
//=================================================================================================//
BaseParticles &SPHBody::getBaseParticles()
{
    if (base_particles_ == nullptr)
    {
        std::cout << "\n Error: BaseParticle of body: " << body_name_ << " not constructed yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *base_particles_;
};
//=================================================================================================//
BaseMaterial &SPHBody::getBaseMaterial()
{
    if (base_material_ == nullptr)
    {
        std::cout << "\n Error: BaseMaterial not generated yet! \n";
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    return *base_material_;
};
//=================================================================================================//
void SPHBody::setSPHBodyBounds(const BoundingBox &bound)
{
    bound_ = bound;
    is_bound_set_ = true;
}
//=================================================================================================//
void SPHBody::registerComputingKernel(execution::Implementation<Base> *implementation)
{
    all_simple_reduce_computing_kernels_.push_back(implementation);
}
//=================================================================================================//
void SPHBody::defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
{
    sph_adaptation_->resetAdaptationRatios(h_spacing_ratio, new_system_refinement_ratio);
}
//=================================================================================================//
void SPHBody::writeParticlesToXmlForRestart(std::string &filefullpath)
{
    base_particles_->writeParticlesToXmlForRestart(filefullpath);
}
//=================================================================================================//
void SPHBody::readParticlesFromXmlForRestart(std::string &filefullpath)
{
    base_particles_->readParticleFromXmlForRestart(filefullpath);
}
//=================================================================================================//
void SPHBody::writeToXmlForReloadParticle(std::string &filefullpath)
{
    base_particles_->writeToXmlForReloadParticle(filefullpath);
}
//=================================================================================================//
BaseCellLinkedList &RealBody::getCellLinkedList()
{
    if (!cell_linked_list_created_)
    {
        cell_linked_list_ptr_ =
            sph_adaptation_->createCellLinkedList(getSPHSystemBounds(), *base_particles_);
        cell_linked_list_created_ = true;
    }
    return *cell_linked_list_ptr_.get();
}
//=================================================================================================//
void RealBody::updateCellLinkedList()
{
    getCellLinkedList().UpdateCellLists(*base_particles_);
}
//=================================================================================================//
} // namespace SPH
