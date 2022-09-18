/**
 * @file 	cell_linked_list.hpp
 * @brief 	Here, template functions in CellLinkedList are defined.
 * @author	Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "base_particles.h"
#include "cell_linked_list.h"

namespace SPH
{
	//=================================================================================================//
	template <class DynamicsRange, typename GetSearchDepth, typename GetNeighborRelation>
	void CellLinkedList::searchNeighborsByParticles(DynamicsRange &dynamics_range,
													ParticleConfiguration &particle_configuration,
													GetSearchDepth &get_search_depth,
													GetNeighborRelation &get_neighbor_relation)
	{
		parallel_for(
			blocked_range<size_t>(0, dynamics_range.SizeOfLoopRange()),
			[&](const blocked_range<size_t> &r)
			{
				StdLargeVec<Vecd> &pos = dynamics_range.getBaseParticles().pos_;
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					size_t index_i = dynamics_range.getParticleIndex(num);
					int search_depth = get_search_depth(index_i);
					Vecu target_cell_index = CellIndexFromPosition(pos[index_i]);
					int i = (int)target_cell_index[0];
					int j = (int)target_cell_index[1];

					Neighborhood &neighborhood = particle_configuration[index_i];
					for (int l = SMAX(i - search_depth, 0); l <= SMIN(i + search_depth, int(number_of_cells_[0]) - 1); ++l)
						for (int m = SMAX(j - search_depth, 0); m <= SMIN(j + search_depth, int(number_of_cells_[1]) - 1); ++m)
						{
							ListDataVector &target_particles = cell_data_lists_[l][m];
							for (const ListData &list_data : target_particles)
							{
								// displacement pointing from neighboring particle to origin particle
								Vecd displacement = pos[index_i] - list_data.second;
								get_neighbor_relation(neighborhood, displacement, index_i, list_data.first);
							}
						}
				}
			},
			ap);
	}
	//=================================================================================================//
}
