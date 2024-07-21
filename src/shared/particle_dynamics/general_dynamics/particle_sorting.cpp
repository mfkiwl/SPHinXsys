#include "particle_sorting.h"

#include "base_body.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
SwapSortableParticleData::SwapSortableParticleData(BaseParticles &base_particles)
    : sequence_(base_particles.ParticleSequences()),
      original_id_(base_particles.ParticleOriginalIds()),
      sortable_data_(base_particles.SortableParticleData()),
      swap_particle_data_value_(sortable_data_) {}
//=================================================================================================//
void SwapSortableParticleData::operator()(size_t *a, size_t *b)
{
    std::swap(*a, *b);

    size_t index_a = a - sequence_;
    size_t index_b = b - sequence_;
    std::swap(original_id_[index_a], original_id_[index_b]);
    swap_particle_data_value_(index_a, index_b);
}
//=================================================================================================//
SingleResolutionSequence::
    SingleResolutionSequence(BaseParticles *base_particles, BaseCellLinkedList *cell_linked_list)
    : cell_linked_list_(cell_linked_list) {}
//=================================================================================================//
SizeT SingleResolutionSequence::operator()(const Vecd &position)
{
    return cell_linked_list->transferMeshIndexToMortonOrder(
        cell_linked_list->CellIndexFromPosition(position));
}
//=================================================================================================//
MultiResolutionSequence::
    MultiResolutionSequence(BaseParticles *base_particles, BaseCellLinkedList *cell_linked_list)
    : multi_level_cell_linked_list_(DynamicCast<MultilevelCellLinkedList>(this, cell_linked_list)),
      mesh_levels_(multi_level_cell_linked_list_->getMeshLevels()),
      kernel_(multi_level_cell_linked_list_->getKernel()),
      h_ratio_(base_particles->getVariableDataByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
SizeT MultiResolutionSequence::operator()(const Vecd &position)
{
    size_t level = multi_level_cell_linked_list_->getMeshLevel(kernel_->CutOffRadius(h_ratio_[i]));
    return mesh_levels_[level]->transferMeshIndexToMortonOrder(
        mesh_levels_[level]->CellIndexFromPosition(position));
}
//=================================================================================================//
} // namespace SPH
