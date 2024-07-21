#ifndef PARTICLE_SORTING_HPP
#define PARTICLE_SORTING_HPP

#include "particle_sorting.h"

namespace SPH
{
//=================================================================================================//
template <typename SequenceMethodType, class ExecutionPolicy>
ParticleSorting<SequenceMethodType, ExecutionPolicy>::
    ComputingSequences::ComputingSequences(RealBody &real_body)
    : LocalDynamics(real_body), DataDelegateSimple(real_body),
      sequence_method_(particles_, &real_body.getCellLinkedList()),
      pos_(particles_->getVariableByName<Vecd>("Position")),
      sequence_(particles_->registerSharedVariable<SizeT>("Sequence")) {}
//=================================================================================================//
template <typename SequenceMethodType, class ExecutionPolicy>
SizeT ParticleSorting<SequenceMethodType, ExecutionPolicy>::
    ComputingSequences::update(size_t index_i, Real dt)
{
    sequence_[index_i] = sequence_method_(pos_[index_i]);
}
//=================================================================================================//
template <typename SequenceMethodType, class ExecutionPolicy>
ParticleSorting<SequenceMethodType, ExecutionPolicy>::ParticleSorting(RealBody &real_body)
    : LocalDynamics(real_body), DataDelegateSimple(real_body),
      sequence_(particles_->registerSharedVariable<SizeT>("Sequence")),
      swap_sortable_particle_data_(*particles_), compare_(),
      quick_sort_particle_range_(sequence_, 0, compare_, swap_sortable_particle_data_),
      quick_sort_particle_body_() {}
//=================================================================================================//
void ParticleSorting::sortingParticleData(size_t *begin, size_t size)
{
    quick_sort_particle_range_.begin_ = sequence_;
    quick_sort_particle_range_.size_ = particles_->TotalRealParticles();
    parallel_for(quick_sort_particle_range_, quick_sort_particle_body_, ap);
}
//=================================================================================================//
template <typename SequenceMethodType, class ExecutionPolicy>
ParticleSorting<SequenceMethodType, ExecutionPolicy>::
    UpdateSortedId::UpdateSortedId(RealBody &real_body)
    : LocalDynamics(real_body), DataDelegateSimple(real_body),
      original_id_(particles_->getVariableDataByName<SizeT>("OriginalID")),
      sorted_id_(particles_->getVariableDataByNameFrom<SizeT>("SortedID", "OriginalID")) {}
//=================================================================================================//
template <typename SequenceMethodType, class ExecutionPolicy>
void ParticleSorting<SequenceMethodType, ExecutionPolicy>::update(size_t index_i, Real dt)
{
    sorted_id_[original_id_[index_i]] = index_i;
}
//=================================================================================================//
} // namespace SPH
#endif // PARTICLE_SORTING_HPP