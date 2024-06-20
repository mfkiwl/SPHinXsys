#include "smoothing_kernel.h"

namespace SPH
{
//=================================================================================================//
TabulatedFunction::TabulatedFunction(Real dq, std::array<Real, 4> delta_q, KernelDataSize data)
    : dq_(dq), delta_q_(delta_q), data_(data) {}
//=================================================================================================//
Real TabulatedFunction::operator()(Real q) const
{
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
Real TabulatedFunction::operator()(Real h_ratio, Real q) const
{
    return operator()(q * h_ratio);
}
//=================================================================================================//
bool WithinCutOff::operator()(Vecd &displacement) const
{
    return displacement.squaredNorm() < rc_ref_sqr_ ? true : false;
}
//=================================================================================================//
bool WithinCutOff::operator()(Real h_ratio, Vecd &displacement) const
{
    return (h_ratio * displacement).squaredNorm() < rc_ref_sqr_ ? true : false;
}
//=================================================================================================//
BaseKernel::BaseKernel(const std::string &name, Real h, Real kernel_size, Real truncation)
    : name_(name), h_(h), kernel_size_(kernel_size),
      cutoff_radius_(truncation * h * kernel_size) {}
//=================================================================================================//
WendlandC2::KernelWendlandC2(Real h)
    : BaseKernel("WendlandC2", h, 2.0)
{
    factor_w1d_ = 3.0 / 4.0 / h_;
    factor_w2d_ = 7.0 / (4.0 * Pi) / h_ / h_;
    factor_w3d_ = 21.0 / (16.0 * Pi) / h_ / h_ / h_;
    setDerivativeParameters();
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
} // namespace SPH
