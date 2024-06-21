#include "smoothing_kernel.h"

namespace SPH
{
//=================================================================================================//
TabulatedFunction::TabulatedFunction(Real h, Real dq, KernelDataArray data)
    : inv_h_(1.0 / h), dq_(dq), data_(data)
{
    delta_q_[0] = (-1.0 * dq) * (-2.0 * dq) * (-3.0 * dq);
    delta_q_[1] = dq * (-1.0 * dq) * (-2.0 * dq);
    delta_q_[2] = (2.0 * dq) * dq * (-1.0 * dq);
    delta_q_[3] = (3.0 * dq) * (2.0 * dq) * dq;
}
//=================================================================================================//
Real TabulatedFunction::operator()(Real distance) const
{
    Real q = distance * inv_h_;
    int location = (int)floor(q / dq_);
    int i = location + 1;
    Real fraction_1 = q - Real(location) * dq_; // fraction_1 correspond to i
    Real fraction_0 = fraction_1 + dq_;         // fraction_0 correspond to i-1
    Real fraction_2 = fraction_1 - dq_;         // fraction_2 correspond to i+1
    Real fraction_3 = fraction_1 - 2 * dq_;     ////fraction_3 correspond to i+2

    return (fraction_1 * fraction_2 * fraction_3) / delta_q_[0] * data_[i - 1] +
           (fraction_0 * fraction_2 * fraction_3) / delta_q_[1] * data_[i] +
           (fraction_0 * fraction_1 * fraction_3) / delta_q_[2] * data_[i + 1] +
           (fraction_0 * fraction_1 * fraction_2) / delta_q_[3] * data_[i + 2];
}
//=================================================================================================//
BaseKernel::BaseKernel(const std::string &name, Real h, Real kernel_size, Real truncation)
    : name_(name), h_(h), kernel_size_(kernel_size),
      cutoff_radius_(truncation * h * kernel_size) {}
//=================================================================================================//
WendlandC2::WendlandC2(Real h) : BaseKernel("WendlandC2", h, 2.0)
{
    factor_w1d_ = 3.0 / 4.0 / h_;
    factor_w2d_ = 7.0 / (4.0 * Pi) / h_ / h_;
    factor_w3d_ = 21.0 / (16.0 * Pi) / h_ / h_ / h_;
}
//=================================================================================================//
Real WendlandC2::W(const Real q) const
{
    return pow(1.0 - 0.5 * q, 4) * (1.0 + 2.0 * q);
}
//=================================================================================================//
Real WendlandC2::dW(const Real q) const
{
    return 0.625 * pow(q - 2.0, 3) * q;
}
//=================================================================================================//
Real SmoothingKernel::computeLatticeNumberDensity(Vec2d zero, Real particle_spacing)
{
    Real sigma(0);
    int search_depth = int(cutoff_radius_ / particle_spacing) + 1;
    for (int j = -search_depth; j <= search_depth; ++j)
        for (int i = -search_depth; i <= search_depth; ++i)
        {
            Vec2d particle_location(i * particle_spacing, j * particle_spacing);
            Real distance = particle_location.norm();
            if (distance < cutoff_radius_)
                sigma += w2d_(distance);
        }
    return sigma;
}
//=================================================================================================//
Real SmoothingKernel::computeLatticeNumberDensity(Vec3d zero, Real particle_spacing)
{
    Real sigma(0);
    int search_depth = int(cutoff_radius_ / particle_spacing) + 1;
    for (int k = -search_depth; k <= search_depth; ++k)
        for (int j = -search_depth; j <= search_depth; ++j)
            for (int i = -search_depth; i <= search_depth; ++i)
            {
                Vec3d particle_location(i * particle_spacing,
                                        j * particle_spacing, k * particle_spacing);
                Real distance = particle_location.norm();
                if (distance < cutoff_radius_)
                    sigma += w3d_(distance);
            }
    return sigma;
}
//=================================================================================================//
} // namespace SPH
