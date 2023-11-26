#include "continuum_dynamics_inner.hpp"
namespace SPH
{
//=================================================================================================//
    namespace continuum_dynamics
    {
        ContinuumInitialCondition::ContinuumInitialCondition(SPHBody &sph_body)
            : LocalDynamics(sph_body), PlasticContinuumDataSimple(sph_body),
              pos_(particles_->pos_), vel_(particles_->vel_), stress_tensor_3D_(particles_->stress_tensor_3D_) {}
        //=================================================================================================//
        ContinuumAcousticTimeStepSize::ContinuumAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
            : fluid_dynamics::AcousticTimeStepSize(sph_body, acousticCFL) {}
        //=================================================================================================//
        Real ContinuumAcousticTimeStepSize::reduce(size_t index_i, Real dt)
        {
            return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
        }
        //=================================================================================================//
        Real ContinuumAcousticTimeStepSize::outputResult(Real reduced_value)
        {
            return acousticCFL_ * smoothing_length_min_ / (fluid_.ReferenceSoundSpeed() + TinyReal);
        }
        //=================================================================================================//
        BaseRelaxation::BaseRelaxation(BaseInnerRelation &inner_relation)
            : LocalDynamics(inner_relation.getSPHBody()), ContinuumDataInner(inner_relation),
              continuum_(particles_->continuum_), rho_(particles_->rho_),
              p_(*particles_->getVariableByName<Real>("Pressure")), drho_dt_(*particles_->getVariableByName<Real>("DensityChangeRate")),
              pos_(particles_->pos_), vel_(particles_->vel_),
              acc_(particles_->acc_), acc_prior_(particles_->acc_prior_) {}
        //=================================================================================================//
        //===============================ArtificialStressRelaxation======================================//
        //=================================================================================================//
        BaseArtificialStressRelaxation ::
            BaseArtificialStressRelaxation(BaseInnerRelation &inner_relation, Real epsilon)
            : BaseRelaxation(inner_relation), smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
              reference_spacing_(sph_body_.sph_adaptation_->ReferenceSpacing()), epsilon_(epsilon) {}
        Matd BaseArtificialStressRelaxation::repulsiveForce(Matd stress_tensor_i, Real rho_i)
        {
            Real sigma_xx = stress_tensor_i(0, 0);
            Real sigma_xy = stress_tensor_i(0, 1);
            Real sigma_yy = stress_tensor_i(1, 1);
            Real tiny_real(0);
            sigma_xx - sigma_yy > 0 ? tiny_real = TinyReal : tiny_real = -TinyReal;
            Real tan_sita_2 = 2 * sigma_xy / (sigma_xx - sigma_yy + tiny_real);
            Real sita_2 = atan(tan_sita_2);
            Real sita = sita_2 / 2.0;
            Real sigma_xx_dot = cos(sita) * cos(sita) * sigma_xx + 2 * cos(sita) * sin(sita) * sigma_xy + sin(sita) * sin(sita) * sigma_yy;
            Real sigma_yy_dot = sin(sita) * sin(sita) * sigma_xx + 2 * cos(sita) * sin(sita) * sigma_xy + cos(sita) * cos(sita) * sigma_yy;
            Real R_xx_dot = 0;
            Real R_yy_dot = 0;
            if (sigma_xx_dot > 0)
            {
                R_xx_dot = -epsilon_ * sigma_xx_dot / (rho_i * rho_i);
            }
            if (sigma_yy_dot > 0)
            {
                R_yy_dot = -epsilon_ * sigma_yy_dot / (rho_i * rho_i);
            }
            Matd R = Matd::Zero();
            R(0, 0) = R_xx_dot * cos(sita) * cos(sita) + R_yy_dot * sin(sita) * sin(sita);
            R(1, 1) = R_xx_dot * sin(sita) * sin(sita) + R_yy_dot * cos(sita) * cos(sita);
            R(0, 1) = (R_xx_dot - R_yy_dot) * cos(sita) * sin(sita);
            R(1, 0) = R(0, 1);
            return R;
        }
        ArtificialNormalShearStressRelaxation ::
            ArtificialNormalShearStressRelaxation(BaseInnerRelation &inner_relation, Real exponent)
            : BaseArtificialStressRelaxation(inner_relation), shear_stress_(particles_->shear_stress_), acc_shear_(particles_->acc_shear_), exponent_(exponent) {}
        void ArtificialNormalShearStressRelaxation::interaction(size_t index_i, Real dt)
        {
            Vecd acceleration = Vecd::Zero();
            Real rho_i = rho_[index_i];
            Matd stress_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            Real W_ini = sph_body_.sph_adaptation_->getKernel()->W_2D(reference_spacing_ / smoothing_length_);
            Matd R_i = repulsiveForce(stress_i, rho_i);
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

                Real W_ij = sph_body_.sph_adaptation_->getKernel()->W_2D(r_ij / smoothing_length_);
                Real f_ij = W_ij / W_ini;
                Matd stress_j = shear_stress_[index_j] - p_[index_j] * Matd::Identity();
                Matd R_j = repulsiveForce(stress_j, rho_[index_j]);
                Matd repulsive_force = pow(f_ij, exponent_) * (R_i + R_j);

                acceleration += rho_[index_j] * repulsive_force * nablaW_ijV_j;
            }
            acc_shear_[index_i] = acceleration;
        }

        //=================================================================================================//
        //===============================ShearStressRelaxation1stHalf======================================//
        //=================================================================================================/
        ShearStressRelaxation1stHalf ::
            ShearStressRelaxation1stHalf(BaseInnerRelation &inner_relation)
            : BaseRelaxation(inner_relation),
              shear_stress_(particles_->shear_stress_), shear_stress_rate_(particles_->shear_stress_rate_),
              acc_shear_(particles_->acc_shear_) {}
        //=================================================================================================//
        void ShearStressRelaxation1stHalf::initialization(size_t index_i, Real dt)
        {
            shear_stress_[index_i] += shear_stress_rate_[index_i] * dt * 0.5;
        }
        //=================================================================================================//
        void ShearStressRelaxation1stHalf::interaction(size_t index_i, Real dt)
        {
            Real rho_i = rho_[index_i];
            Matd shear_stress_i = shear_stress_[index_i];
            Vecd acceleration = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                acceleration += rho_[index_j] * (shear_stress_i / (rho_i * rho_i) + shear_stress_[index_j] / (rho_[index_j] * rho_[index_j])) * nablaW_ijV_j;
            }
            acc_shear_[index_i] += acceleration; // for with artificial stress
            //acc_shear_[index_i] = acceleration; // for without artificial stress
        }

        void ShearStressRelaxation1stHalf::update(size_t index_i, Real dt)
        {
            vel_[index_i] += acc_shear_[index_i] * dt;
        }
        //=================================================================================================//
        //===============================ShearStressRelaxation2ndHalf======================================//
        //=================================================================================================//
        ShearStressRelaxation2ndHalf ::
            ShearStressRelaxation2ndHalf(BaseInnerRelation &inner_relation)
            : BaseRelaxation(inner_relation),
              shear_stress_(particles_->shear_stress_), shear_stress_rate_(particles_->shear_stress_rate_),
              velocity_gradient_(particles_->velocity_gradient_), strain_tensor_(particles_->strain_tensor_),
              strain_tensor_rate_(particles_->strain_tensor_rate_), von_mises_stress_(particles_->von_mises_stress_) {}
        //=================================================================================================//
        void ShearStressRelaxation2ndHalf::interaction(size_t index_i, Real dt)
        {
            Matd velocity_gradient = Matd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

                Matd velocity_gradient_ij = -(vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
                velocity_gradient += velocity_gradient_ij;
            }
            velocity_gradient_[index_i] = velocity_gradient;
        }

        void ShearStressRelaxation2ndHalf::update(size_t index_i, Real dt)
        {
            shear_stress_rate_[index_i] = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
            shear_stress_[index_i] += shear_stress_rate_[index_i] * dt * 0.5;
            Matd stress_tensor_i = shear_stress_[index_i] + p_[index_i] * Matd::Identity();
            von_mises_stress_[index_i] = getVonMisesStressFromMatrix(stress_tensor_i);
        }
        //=================================================================================================//
        //===================================Non-hourglass formulation=====================================//
        //=================================================================================================//
        ShearAccelerationRelaxation::ShearAccelerationRelaxation(BaseInnerRelation& inner_relation)
            : BaseRelaxation(inner_relation),
            G_(continuum_.getShearModulus(continuum_.getYoungsModulus(), continuum_.getPoissonRatio())),
            smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()), shear_stress_(particles_->shear_stress_),
            B_(*this->particles_->template registerSharedVariable<Matd>("KernelCorrectionMatrix", Matd::Identity())), acc_shear_(particles_->acc_shear_) {}
        void ShearAccelerationRelaxation::interaction(size_t index_i, Real dt)
        {
            Real rho_i = rho_[index_i];
            Vecd acceleration = Vecd::Zero();
            Vecd vel_i = vel_[index_i];
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

                Vecd v_ij = (vel_i - vel_[index_j]);
                acceleration += 2 * v_ij * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
            }
            acc_shear_[index_i] += G_ * acceleration * dt / rho_i;
        }
        //=================================================================================================//
        void AngularConservativeShearAccelerationRelaxation::interaction(size_t index_i, Real dt)
        {
            Real rho_i = rho_[index_i];
            Vecd acceleration = Vecd::Zero();
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Vecd& e_ij = inner_neighborhood.e_ij_[n];
                Real eta_ij = 2 * (0.7*(Real)Dimensions+2.1) * (vel_[index_i] - vel_[index_j]).dot(e_ij) / (r_ij + TinyReal);
                acceleration += eta_ij * dW_ijV_j * e_ij;
            }
            acc_shear_[index_i] += G_ * acceleration * dt / rho_i;
        }
        //=================================================================================================//
        ShearStressRelaxation ::
            ShearStressRelaxation(BaseInnerRelation& inner_relation)
            : BaseRelaxation(inner_relation),
            shear_stress_(particles_->shear_stress_), shear_stress_rate_(particles_->shear_stress_rate_),
            velocity_gradient_(particles_->velocity_gradient_), strain_tensor_(particles_->strain_tensor_),
            strain_tensor_rate_(particles_->strain_tensor_rate_), von_mises_stress_(particles_->von_mises_stress_),
            von_mises_strain_(particles_->von_mises_strain_), Vol_(particles_->Vol_), 
            B_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")) {}
        void ShearStressRelaxation::initialization(size_t index_i, Real dt)
        {
            strain_tensor_[index_i] += strain_tensor_rate_[index_i] * 0.5 * dt;
            shear_stress_[index_i] += shear_stress_rate_[index_i] * 0.5 * dt;
        }
        void ShearStressRelaxation::interaction(size_t index_i, Real dt)
        {
            Matd velocity_gradient = Matd::Zero();
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j_ = inner_neighborhood.dW_ijV_j_[n];
                Vecd& e_ij = inner_neighborhood.e_ij_[n];

                //Matd velocity_gradient_ij = - (vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
                Vecd v_ij = vel_[index_i] - vel_[index_j];
                Matd velocity_gradient_ij = -v_ij * (B_[index_i] * e_ij * dW_ijV_j_).transpose();
                velocity_gradient += velocity_gradient_ij;
            }
            velocity_gradient_[index_i] = velocity_gradient;
            //calculate strain
            Matd strain_rate = 0.5 * (velocity_gradient + velocity_gradient.transpose());
            strain_tensor_rate_[index_i] = strain_rate;
            strain_tensor_[index_i] += strain_tensor_rate_[index_i] * 0.5 * dt;
            Matd strain_i = strain_tensor_[index_i];
            von_mises_strain_[index_i] = getVonMisesStressFromMatrix(strain_i);
        }
        void ShearStressRelaxation::update(size_t index_i, Real dt)
        {
            shear_stress_rate_[index_i] = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
            shear_stress_[index_i] += shear_stress_rate_[index_i] * dt * 0.5;

            //VonMises Stress
            Matd stress_tensor_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
            von_mises_stress_[index_i] = getVonMisesStressFromMatrix(stress_tensor_i);
        }
        //=================================================================================================//
        FixedInAxisDirection::FixedInAxisDirection(BodyPartByParticle& body_part, Vecd constrained_axises)
            : BaseMotionConstraint<BodyPartByParticle>(body_part), constrain_matrix_(Matd::Identity())
        {
            for (int k = 0; k != Dimensions; ++k)
                constrain_matrix_(k, k) = constrained_axises[k];
        };
        //=================================================================================================//
        void FixedInAxisDirection::update(size_t index_i, Real dt)
        {
            vel_[index_i] = constrain_matrix_ * vel_[index_i];
        };
        //=================================================================================================//
        ConstrainSolidBodyMassCenter::
            ConstrainSolidBodyMassCenter(SPHBody& sph_body, Vecd constrain_direction)
            : LocalDynamics(sph_body), ContinuumDataSimple(sph_body),
            correction_matrix_(Matd::Identity()), vel_(particles_->vel_),
            compute_total_momentum_(sph_body, "Velocity")
        {
            for (int i = 0; i != Dimensions; ++i)
                correction_matrix_(i, i) = constrain_direction[i];
            ReduceDynamics<QuantitySummation<Real>> compute_total_mass_(sph_body, "MassiveMeasure");
            total_mass_ = compute_total_mass_.exec();
        }
        //=================================================================================================//
        void ConstrainSolidBodyMassCenter::setupDynamics(Real dt)
        {
            velocity_correction_ =
                correction_matrix_ * compute_total_momentum_.exec(dt) / total_mass_;
        }
        //=================================================================================================//
        void ConstrainSolidBodyMassCenter::update(size_t index_i, Real dt)
        {
            vel_[index_i] -= velocity_correction_;
        }
        //=================================================================================================//
        //================================================Plastic==========================================//
        //=================================================================================================//
        BaseRelaxationPlastic::BaseRelaxationPlastic(BaseInnerRelation& inner_relation)
            : LocalDynamics(inner_relation.getSPHBody()), PlasticContinuumDataInner(inner_relation),
            plastic_continuum_(particles_->plastic_continuum_), rho_(particles_->rho_),
            p_(*particles_->getVariableByName<Real>("Pressure")), drho_dt_(*particles_->registerSharedVariable<Real>("DensityChangeRate")), pos_(particles_->pos_), vel_(particles_->vel_), acc_(particles_->acc_), acc_prior_(particles_->acc_prior_),
            stress_tensor_3D_(particles_->stress_tensor_3D_), strain_tensor_3D_(particles_->strain_tensor_3D_),
            stress_rate_3D_(particles_->stress_rate_3D_), strain_rate_3D_(particles_->strain_rate_3D_),
            elastic_strain_tensor_3D_(particles_->elastic_strain_tensor_3D_), elastic_strain_rate_3D_(particles_->elastic_strain_rate_3D_) {}
        Matd BaseRelaxationPlastic::reduceTensor(Mat3d tensor_3d)
        {
            Matd tensor_2d;
            for (int i = 0; i < (Real)Dimensions; i++)
            {
                for (int j = 0; j < (Real)Dimensions; j++)
                {
                    tensor_2d(i, j) = tensor_3d(i, j);
                }
            }
            return tensor_2d;
        }
        Mat3d BaseRelaxationPlastic::increaseTensor(Matd tensor_2d)
        {
            Mat3d tensor_3d = Mat3d::Zero();
            for (int i = 0; i < (Real)Dimensions; i++)
            {
                for (int j = 0; j < (Real)Dimensions; j++)
                {
                    tensor_3d(i, j) = tensor_2d(i, j);
                }
            }
            return tensor_3d;
        }
        //====================================================================================//
        //===============================StressDiffusion======================================//
        //====================================================================================//
        StressDiffusion::StressDiffusion(BaseInnerRelation& inner_relation)
            : BaseRelaxationPlastic(inner_relation), fai_(DynamicCast<PlasticContinuum>(this, plastic_continuum_).getFrictionAngle()), smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
            sound_speed_(plastic_continuum_.ReferenceSoundSpeed()) {}
        void StressDiffusion::interaction(size_t index_i, Real dt)
        {
            Real gravity = abs(acc_prior_[index_i](1, 0));
            Real density = plastic_continuum_.getDensity();
            Mat3d diffusion_stress_rate_ = Mat3d::Zero();
            Mat3d diffusion_stress_ = Mat3d::Zero();

            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Real y_ij = pos_[index_i](1, 0) - pos_[index_j](1, 0);
                diffusion_stress_ = stress_tensor_3D_[index_i] - stress_tensor_3D_[index_j];
                diffusion_stress_(0, 0) = diffusion_stress_(0, 0) - (1 - sin(fai_)) * density * gravity * y_ij;
                diffusion_stress_(1, 1) = diffusion_stress_(1, 1) - density * gravity * y_ij;
                diffusion_stress_(2, 2) = diffusion_stress_(2, 2) - (1 - sin(fai_)) * density * gravity * y_ij;
                diffusion_stress_rate_ += 2 * zeta_ * smoothing_length_ * sound_speed_ * diffusion_stress_ * r_ij * dW_ijV_j / (r_ij * r_ij + 0.01 * smoothing_length_);
            }
            stress_rate_3D_[index_i] = diffusion_stress_rate_;
        }
        //=================================================================================================//
        //==============================ShearStressRelaxationHourglassControl==============================//
        //=================================================================================================//
        ShearStressRelaxationHourglassControl ::
            ShearStressRelaxationHourglassControl(BaseInnerRelation& inner_relation, int hourglass_control)
            : BaseRelaxation(inner_relation),
            shear_stress_(particles_->shear_stress_), shear_stress_rate_(particles_->shear_stress_rate_),
            velocity_gradient_(particles_->velocity_gradient_), acc_shear_(particles_->acc_shear_),
            von_mises_stress_(particles_->von_mises_stress_),
            //B_(*this->particles_->template registerSharedVariable<Matd>("CorrectionMatrix", Matd::Identity())),
            B_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")),
            hourglass_control_(hourglass_control),
            pos0_(particles_->pos0_),
            strain_tensor_rate_(particles_->strain_tensor_rate_), strain_tensor_(particles_->strain_tensor_)
        {
            particles_->registerVariable(acc_hourglass_, "AccelerationHourglass");
            particles_->registerSortableVariable<Vecd>("AccelerationHourglass");
            particles_->addVariableToWrite<Vecd>("AccelerationHourglass");
            //
            particles_->registerVariable(scale_coef_, "ScaleCoefficient");
            particles_->registerSortableVariable<Matd>("ScaleCoefficient");
            particles_->addVariableToWrite<Matd>("ScaleCoefficient");
            //
            particles_->registerVariable(shear_strain_, "ShearStrain");
            particles_->registerSortableVariable<Matd>("ShearStrain");
        }

        void ShearStressRelaxationHourglassControl::interaction(size_t index_i, Real dt)
        {
            Matd velocity_gradient = Matd::Zero();
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            Vecd vel_i = vel_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Vecd& e_ij = inner_neighborhood.e_ij_[n];
                Matd velocity_gradient_ij;
                if (hourglass_control_)
                {
                    velocity_gradient_ij = -(vel_i - vel_[index_j]) * (B_[index_i] * e_ij * dW_ijV_j).transpose();
                    //std::cout << B_[index_i] << std::endl;
                }
                else
                {
                    velocity_gradient_ij = -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
                }
                velocity_gradient += velocity_gradient_ij;
            }
            velocity_gradient_[index_i] = velocity_gradient;
            shear_stress_rate_[index_i] = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);

            shear_stress_[index_i] += shear_stress_rate_[index_i] * dt;
            Matd stress_tensor_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
            von_mises_stress_[index_i] = getVonMisesStressFromMatrix(stress_tensor_i);
            //
            strain_tensor_rate_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
            strain_tensor_[index_i] += strain_tensor_rate_[index_i] * dt;
            shear_strain_[index_i] = strain_tensor_[index_i] - Matd::Identity() * strain_tensor_[index_i].trace() / (Real)Dimensions;
        }
        void ShearStressRelaxationHourglassControl::update(size_t index_i, Real dt)
        {
            Real G_ = continuum_.getShearModulus(continuum_.getYoungsModulus(), continuum_.getPoissonRatio());
            const Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd);

            Real rho_i = rho_[index_i];
            Matd shear_stress_i = shear_stress_[index_i];
            Vecd acceleration = Vecd::Zero();
            Vecd acceleration_hourglass = Vecd::Zero();
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Vecd& e_ij = inner_neighborhood.e_ij_[n];
                acceleration += ((shear_stress_i + shear_stress_[index_j]) / rho_i) * dW_ijV_j * e_ij;
                if (hourglass_control_)
                {
                    Vecd v_ij = vel_[index_i] - vel_[index_j];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd v_ij_correction = v_ij - 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_j]) * r_ij * e_ij;

                    Real coef = 3.0 * (Real)Dimensions;   //oscillating_beam_2D when PH = 0.02
                    //acceleration_hourglass += coef * G_ * v_ij_correction * gamma_limiter * dW_ijV_j * dt / (rho_i * r_ij);
                    acceleration_hourglass += coef * G_ * v_ij_correction * dW_ijV_j * dt / (rho_i * r_ij);

                    // the inverse method
                    //Matd shear_strain__i_inverse = shear_strain_[index_i].completeOrthogonalDecomposition().pseudoInverse();
                    //Matd shear_strain__j_inverse = shear_strain_[index_j].completeOrthogonalDecomposition().pseudoInverse();
                    //scale_coef_[index_i] = (shear_stress_[index_i] * shear_strain__i_inverse + shear_stress_[index_j] * shear_strain__j_inverse) / (4 * G_);

                    //acceleration_hourglass += coef * scale_coef_[index_i] * G_ * v_ij_correction * dW_ijV_j * dt / (rho_i * r_ij);
                }
            }
            if (hourglass_control_)
            {
                acc_hourglass_[index_i] += acceleration_hourglass;
            }
            acc_shear_[index_i] = acceleration + acc_hourglass_[index_i];
        }
        //=============================================================================================//
        //====================================J2Plasticity=============================================//
        //=============================================================================================//
        BaseRelaxationJ2Plasticity::BaseRelaxationJ2Plasticity(BaseInnerRelation& inner_relation)
            : LocalDynamics(inner_relation.getSPHBody()), J2PlasticicityDataInner(inner_relation),
            J2_plasticity_(particles_->J2_plasticity_), rho_(particles_->rho_),
            p_(*particles_->getVariableByName<Real>("Pressure")), drho_dt_(*particles_->registerSharedVariable<Real>("DensityChangeRate")), pos_(particles_->pos_), vel_(particles_->vel_), acc_(particles_->acc_), acc_prior_(particles_->acc_prior_) {}
        Matd BaseRelaxationJ2Plasticity::reduceTensor(Mat3d tensor_3d)
        {
            Matd tensor_2d;
            for (int i = 0; i < (Real)Dimensions; i++)
            {
                for (int j = 0; j < (Real)Dimensions; j++)
                {
                    tensor_2d(i, j) = tensor_3d(i, j);
                }
            }
            return tensor_2d;
        }
        Mat3d BaseRelaxationJ2Plasticity::increaseTensor(Matd tensor_2d)
        {
            Mat3d tensor_3d = Mat3d::Zero();
            for (int i = 0; i < (Real)Dimensions; i++)
            {
                for (int j = 0; j < (Real)Dimensions; j++)
                {
                    tensor_3d(i, j) = tensor_2d(i, j);
                }
            }
            return tensor_3d;
        }
        //=============================================================================================//
        ShearStressRelaxationHourglassControlJ2Plasticity ::
            ShearStressRelaxationHourglassControlJ2Plasticity(BaseInnerRelation& inner_relation, int hourglass_control)
            : BaseRelaxationJ2Plasticity(inner_relation),
            shear_stress_3D_(particles_->shear_stress_3D_), shear_stress_rate_3D_(particles_->shear_stress_rate_3D_),
            shear_strain_3D_(particles_->shear_strain_3D_), shear_strain_rate_3D_(particles_->shear_strain_rate_3D_),
            plastic_indicator_(particles_->plastic_indicator_),
            velocity_gradient_(particles_->velocity_gradient_), acc_shear_(particles_->acc_shear_),
            von_mises_stress_(particles_->von_mises_stress_),
            B_(*particles_->getVariableByName<Matd>("KernelCorrectionMatrix")),
            strain_tensor_3D_(particles_->strain_tensor_3D_), strain_rate_3D_(particles_->strain_rate_3D_),
            hourglass_control_(hourglass_control), E_(J2_plasticity_.getYoungsModulus()), nu_(J2_plasticity_.getPoissonRatio())
        {
            particles_->registerVariable(acc_hourglass_, "AccelerationHourglass");
            particles_->registerSortableVariable<Vecd>("AccelerationHourglass");

            particles_->registerVariable(scale_coef_, "ScaleCoefficient");
            particles_->registerSortableVariable<Matd>("ScaleCoefficient");
            particles_->addVariableToWrite<Matd>("ScaleCoefficient");

            particles_->registerVariable(plastic_strain_tensor_3D_, "PlasticStrainTensor3D");
            particles_->registerSortableVariable<Mat3d>("PlasticStrainTensor3D");

            particles_->registerVariable(hardening_parameter_, "HardeningParameter");
            particles_->registerSortableVariable<Real>("HardeningParameter");
            particles_->addVariableToWrite<Real>("HardeningParameter");
        }
        void ShearStressRelaxationHourglassControlJ2Plasticity::initialization(size_t index_i, Real dt)
        {
            Matd velocity_gradient = Matd::Zero();
            const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            Vecd vel_i = vel_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                const Vecd& e_ij = inner_neighborhood.e_ij_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                Matd velocity_gradient_ij = -(vel_i - vel_[index_j]) * (B_[index_i] * e_ij * dW_ijV_j).transpose();
                velocity_gradient += velocity_gradient_ij;
            }
            velocity_gradient_[index_i] = velocity_gradient;
        }
        void ShearStressRelaxationHourglassControlJ2Plasticity::interaction(size_t index_i, Real dt)
        {
            Real G_ = J2_plasticity_.getShearModulus(J2_plasticity_.getYoungsModulus(), J2_plasticity_.getPoissonRatio());
            Mat3d velocity_gradient = increaseTensor(velocity_gradient_[index_i]);
            shear_stress_rate_3D_[index_i] = J2_plasticity_.ConstitutiveRelationShearStress(velocity_gradient, shear_stress_3D_[index_i]);
            shear_stress_3D_[index_i] += shear_stress_rate_3D_[index_i] * dt;

            //For plasticity
            plastic_indicator_[index_i] = J2_plasticity_.PlasticIndicator(shear_stress_3D_[index_i]);
            //shear_stress_3D_[index_i] = J2_plasticity_.ReturnMappingShearStress(shear_stress_3D_[index_i]);

            Mat3d stress_tensor_i = shear_stress_3D_[index_i] - p_[index_i] * Mat3d::Identity();
            von_mises_stress_[index_i] = getVonMisesStressFromMatrix(reduceTensor(stress_tensor_i));
            // strain tensor
            strain_rate_3D_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
            strain_tensor_3D_[index_i] += strain_rate_3D_[index_i] * dt;
            // shear strain
            Mat3d strain_tensor_i = strain_tensor_3D_[index_i];
            shear_strain_3D_[index_i] = strain_tensor_i - Mat3d::Identity() * strain_tensor_i.trace() / 3.0;
            // calculate scale matrix
            if (plastic_indicator_[index_i])
            {
                scale_coef_[index_i] = reduceTensor(shear_stress_3D_[index_i]) * reduceTensor(shear_strain_3D_[index_i]).inverse() / (2.0 * G_);
            }
            else
            {
                scale_coef_[index_i] = Matd::Identity();
            }
            //calculate elastic strain
            Mat3d deviatoric_stress = shear_stress_3D_[index_i];
            Real hydrostatic_pressure = - p_[index_i];
            Mat3d elastic_strain_tensor_3D = deviatoric_stress / (2 * J2_plasticity_.getShearModulus(E_, nu_)) + hydrostatic_pressure * Mat3d::Identity() / (9 * J2_plasticity_.getBulkModulus(E_, nu_));
            plastic_strain_tensor_3D_[index_i] = strain_tensor_3D_[index_i] - elastic_strain_tensor_3D;
            Real J2_plastic_strain = 0.5 * (plastic_strain_tensor_3D_[index_i].cwiseProduct(plastic_strain_tensor_3D_[index_i].transpose())).sum();
            hardening_parameter_[index_i] = sqrt(2.0 / 3.0) * sqrt(2.0 * J2_plastic_strain);

        }
        void ShearStressRelaxationHourglassControlJ2Plasticity::update(size_t index_i, Real dt)
        {
            Real G_ = J2_plasticity_.getShearModulus(J2_plasticity_.getYoungsModulus(), J2_plasticity_.getPoissonRatio());
            Vecd acceleration_hourglass = Vecd::Zero();

            Vecd acceleration = Vecd::Zero();
            Real rho_dissipation(0);
            Real rho_i = rho_[index_i];
            Matd shear_stress_i = reduceTensor(shear_stress_3D_[index_i]);
            const Neighborhood& inner_neighborhood = inner_configuration_[index_i];

            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                const Vecd& e_ij = inner_neighborhood.e_ij_[n];
                Matd shear_stress_j = reduceTensor(shear_stress_3D_[index_j]);
                acceleration += ((shear_stress_i + shear_stress_j) / rho_i) * nablaW_ijV_j;
                //hourglass control
                Vecd v_ij = vel_[index_i] - vel_[index_j];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd v_ij_correction = v_ij - 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_j]) * r_ij * e_ij;
                Real coef = 3.0;
                acceleration_hourglass += coef * 0.5 * (scale_coef_[index_i] + scale_coef_[index_j]) * G_ * v_ij_correction * dW_ijV_j * dt / (rho_i * r_ij);
            }
            acc_hourglass_[index_i] += acceleration_hourglass;
            acc_shear_[index_i] = acceleration + acc_hourglass_[index_i];
        }
    } // namespace continuum_dynamics
} // namespace SPH
