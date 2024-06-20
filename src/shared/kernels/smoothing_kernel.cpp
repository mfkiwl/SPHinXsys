#include "smoothing_kernel.h"

namespace SPH
{
//=================================================================================================//
TabulatedFunction::TabulatedFunction(Real dq, std::array<Real, 4> delta_q, std::array<Real, 32> data)
    : dq_(dq), delta_q_(delta_q), data_(data) {}
//=================================================================================================//
Real TabulatedRadialFunction::operator()(Real q)
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
BaseKernel::BaseKernel(const std::string &name, Real h, Real kernel_size, Real truncation)
    : name_(name), h_(h), kernel_size_(kernel_size),
      cutoff_radius_(truncation * h * kernel_size) {}
//=================================================================================================//
void BaseKernel::setDerivativeParameters()
{
    factor_dW_1D_ = factor_w1d_ / h_;
    factor_dW_2D_ = factor_w2d_ / h_;
    factor_dW_3D_ = factor_w3d_ / h_;
    factor_d2W_1D_ = factor_dW_1D_ / h_;
    factor_d2W_2D_ = factor_dW_2D_ / h_;
    factor_d2W_3D_ = factor_dW_3D_ / h_;
}
//=================================================================================================//
KernelWendlandC2::KernelWendlandC2(Real h)
    : BaseKernel("Wendland2CKernel", h, 2.0)
{
    factor_w1d_ = 3.0 / 4.0 / h_;
    factor_w2d_ = 7.0 / (4.0 * Pi) / h_ / h_;
    factor_w3d_ = 21.0 / (16.0 * Pi) / h_ / h_ / h_;
    setDerivativeParameters();
}
//=================================================================================================//
Real KernelWendlandC2::W(const Real q) const
{
    return pow(1.0 - 0.5 * q, 4) * (1.0 + 2.0 * q);
}
//=================================================================================================//
Real KernelWendlandC2::dW(const Real q) const
{
    return 0.625 * pow(q - 2.0, 3) * q;
}
//=================================================================================================//
} // namespace SPH
