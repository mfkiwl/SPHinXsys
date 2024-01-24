#pragma once
#include "continuum_dynamics_inner.h"
namespace SPH
{
	namespace continuum_dynamics
	{
		template <class FluidDynamicsType>
		BaseIntegration1stHalf<FluidDynamicsType>::BaseIntegration1stHalf(BaseInnerRelation& inner_relation)
			: FluidDynamicsType(inner_relation),
			acc_shear_(*this->particles_->template getVariableByName<Vecd>("AccelerationByShear")) {}
		//=================================================================================================//
		template <class FluidDynamicsType>
		void BaseIntegration1stHalf<FluidDynamicsType>::update(size_t index_i, Real dt)
		{
			this->vel_[index_i] += ((this->force_prior_[index_i] + this->force_[index_i]) / this->mass_[index_i] + this->acc_shear_[index_i]) * dt;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseStressRelaxation1stHalf<RiemannSolverType>::BaseStressRelaxation1stHalf(BaseInnerRelation& inner_relation)
			: BaseRelaxationPlastic(inner_relation), riemann_solver_(plastic_continuum_, plastic_continuum_),
			velocity_gradient_(particles_->velocity_gradient_) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseStressRelaxation1stHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			p_[index_i] = -stress_tensor_3D_[index_i].trace() / 3;
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		Vecd BaseStressRelaxation1stHalf<RiemannSolverType>::computeNonConservativeForce(size_t index_i)
		{
			Vecd force = force_prior_[index_i] * rho_[index_i];
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];

				force += mass_[index_i] * (p_[index_i] - p_[index_j]) * dW_ijV_j * e_ij;
			}
			return force / rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseStressRelaxation1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Vecd force = Vecd::Zero();
			Real rho_dissipation(0);
			Real rho_i = rho_[index_i];
			Matd stress_tensor_i = reduceTensor(stress_tensor_3D_[index_i]);
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];

			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				Matd stress_tensor_j = reduceTensor(stress_tensor_3D_[index_j]);
				force += mass_[index_i] * rho_[index_j] * ((stress_tensor_i + stress_tensor_j) / (rho_i * rho_[index_j])) * nablaW_ijV_j;
				rho_dissipation += riemann_solver_.DissipativeUJump(p_[index_i] - p_[index_j]) * dW_ijV_j;
			}
			force_[index_i] += force;
			drho_dt_[index_i] = rho_dissipation * rho_[index_i];
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseStressRelaxation1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		BaseStressRelaxation2ndHalf<RiemannSolverType>::BaseStressRelaxation2ndHalf(BaseInnerRelation& inner_relation)
			: BaseRelaxationPlastic(inner_relation), riemann_solver_(plastic_continuum_, plastic_continuum_, 20.0 * (Real)Dimensions),
			velocity_gradient_(particles_->velocity_gradient_), acc_deviatoric_plastic_strain_(particles_->acc_deviatoric_plastic_strain_),
			vertical_stress_(particles_->vertical_stress_), Vol_(particles_->Vol_), mass_(particles_->mass_),
			E_(plastic_continuum_.getYoungsModulus()), nu_(plastic_continuum_.getPoissonRatio()) {}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseStressRelaxation2ndHalf<RiemannSolverType>::initialization(size_t index_i, Real dt)
		{
			pos_[index_i] += vel_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseStressRelaxation2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
		{
			Real density_change_rate(0);
			Vecd p_dissipation = Vecd::Zero();
			Matd velocity_gradient = Matd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				const Vecd& e_ij = inner_neighborhood.e_ij_[n];
				Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
				Real u_jump = (vel_[index_i] - vel_[index_j]).dot(e_ij);
				density_change_rate += u_jump * dW_ijV_j;
				p_dissipation += mass_[index_i] * riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
				velocity_gradient -= (vel_[index_i] - vel_[index_j]) * dW_ijV_j * e_ij.transpose();
			}
			drho_dt_[index_i] += density_change_rate * rho_[index_i];
			force_[index_i] = p_dissipation / rho_[index_i];
			velocity_gradient_[index_i] = velocity_gradient;
		}
		//=================================================================================================//
		template <class RiemannSolverType>
		void BaseStressRelaxation2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
		{
			rho_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_[index_i];
			//update stress
			Mat3d velocity_gradient = increaseTensor(velocity_gradient_[index_i]);
			Mat3d stress_tensor_rate_3D_ = plastic_continuum_.ConstitutiveRelation(velocity_gradient, stress_tensor_3D_[index_i]);
			stress_rate_3D_[index_i] += stress_tensor_rate_3D_;
			stress_tensor_3D_[index_i] += stress_rate_3D_[index_i] * dt;
			//For plasticity
			Mat3d stress_tensor_ = plastic_continuum_.ReturnMapping(stress_tensor_3D_[index_i]);
			stress_tensor_3D_[index_i] = stress_tensor_;
			vertical_stress_[index_i] = stress_tensor_3D_[index_i](1, 1);
			strain_rate_3D_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
			strain_tensor_3D_[index_i] += strain_rate_3D_[index_i] * dt;
			//calculate elastic strain
			Mat3d deviatoric_stress = stress_tensor_3D_[index_i] - (1.0 / 3.0) * stress_tensor_3D_[index_i].trace() * Mat3d::Identity();
			Real hydrostatic_pressure = (1.0 / 3.0) * stress_tensor_3D_[index_i].trace();
			elastic_strain_tensor_3D_[index_i] = deviatoric_stress / (2.0 * plastic_continuum_.getShearModulus(E_, nu_)) + hydrostatic_pressure * Mat3d::Identity() / (9.0 * plastic_continuum_.getBulkModulus(E_, nu_));
			Mat3d plastic_strain_tensor_3D = strain_tensor_3D_[index_i] - elastic_strain_tensor_3D_[index_i];
			acc_deviatoric_plastic_strain_[index_i] = particles_->getDeviatoricPlasticStrain(plastic_strain_tensor_3D);
		}
	}
} // namespace SPH