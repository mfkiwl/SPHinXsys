/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	particle_generator_network.h
 * @brief 	This is a class of particle generator, which generates particles
 * 			with in network or tree form. 
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_GENERATOR_NETWORK_H
#define PARTICLE_GENERATOR_NETWORK_H


#include "sph_data_conainers.h"
#include "base_particle_generator.h"
#include "base_geometry.h"

namespace SPH 
{
	class Tree;
	class BaseLevelSet;
	class BaseMeshCellLinkedList;
	class ComplexShape;

	/**
	 * @class ParticleGeneratorNetwork
	 * @brief Generate a tree-shape network for the conduction system of a heart with particles.
	 */
	class ParticleGeneratorNetwork : public ParticleGenerator
	{
	public:
		ParticleGeneratorNetwork(Vecd starting_pnt, Vecd second_pnt, int iterator, Real grad_factor);
		virtual ~ParticleGeneratorNetwork() {};
		/**
		 *@brief Parameters initialization. 
		 *@param[in] sph_body*(SPHBody) SPHBody to whom it generate particles.
		 */
		virtual void initialize(SPHBody* sph_body) override;
		/**
		 *@brief Created base particles based on edges in branch.
		 *@param[in] base_particles(BaseParticles) Pointer to baseparticle link to a SPHBody.
		 */
		virtual void createBaseParticles(BaseParticles* base_particles) override;
	protected:
		Vecd starting_pnt_;			/**< Starting point for net work. */
		Vecd second_pnt_;				/**< Second point, approximate the growing direction. */
		size_t n_it_; 					/**< Number of iterations (generations of branch. */
		bool fascicles_;				/**< Create fascicles? */
		size_t	segments_in_branch_;	/**< approximated number of segments in a branch. */
		Real segment_length_;			/**< segment length of the branch. */
		Real angle_ = 0.3;				/**< angle with respect to the direction of the previous edge and the new edge. */
		Real repulsivity_ = 0.175; 		/**< repulsivity parameter. */
		Real grad_factor_;				/**< Factor for computing gradient from nearest node. */
		std::vector<Real> fascicle_angles_ = {-1.25, 0.75}; 	/**< angles with respect to the initial edge of the fascicles.*/
		Real fascicle_ratio_ = 15.0; 						/**< ratio of length  of the fascicles. Include one per fascicle to include.*/
		ComplexShape* body_shape_;
		RealBody* real_body_;
		/**
		 *@brief Get the gradient from nearest points, for imposing repulsive force. 
		 *@param[in] pt(Vecd) Inquiry point.
		 *@param[in] delta(Real) parameter for gradient calculation.
		 *@param[in] mesh_cell_linked_list(BaseMeshCellLinkedList) CLL for nearest point searching. 
		 */
		Vecd getGradientFromNearestPoints(Vecd pt, Real delta, BaseMeshCellLinkedList* mesh_cell_linked_list);
		/**
		 *@brief Create a new branch if it is valid.
		 *@param[in] sph_body(SPHBody) The SPHBody to whom the tree belongs.
		 *@param[in] parent_id(size_t) Id of parent branch.
		 *@param[in] angle(Real) The angle for growing new points.
		 *@param[in] repulsivity(Real) The repulsivity for creating new points.
		 *@param[in] number_segments(size_t) Number of segments in this branch.
		 *@param[in] tree(Tree) The tree to whom the branch belongs.
		 */
		bool createABranchIfValid(SPHBody* sph_body, size_t parent_id, Real angle,
			Real repulsivity, size_t number_segments, Tree* tree);
		/**
		 *@brief Functions that creates a new node in the mesh surface and it to the queue is it lies in the surface.
		 *@param[in] init_node vector that contains the coordinates of the last node added in the branch.
		 * 			 vector that contains the coordinates of the last node added in the branch.
		 *@param[in] dir a vector that contains the direction from the init_node to the node to project.
		 *@param[out] end point of the created segment.
		 */
		Vecd createATentativeNewBranchPoint(Vecd init_point, Vecd dir);
		/**
		 *@brief Check if the new point has collision with the existing points.
		 *@param[in] new_point(Vecd) The enquiry point.
		 *@param[in] nearest_neighbor(ListData) The nearest point of the existing points.
		 *@param[in] parent_id(size_t)  Id of parent branch
		 *@param[in] tree(Tree) The tree to whom the branch belongs.
		 */
		bool isCollision(Vecd& new_point, ListData& nearest_neighbor, size_t parent_id, Tree* tree);
		/**
		 *@brief Check if the new point is valid according to extra constraint.
		 *@param[in] new_point(Vecd) The enquiry point.
		 */
		virtual bool extraCheck(Vecd& new_point){return false;};
	};
}
#endif //PARTICLE_GENERATOR_NETWORK_H