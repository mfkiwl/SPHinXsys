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
bool WithinCutOff::operator()(Vecd &displacement) const
{
    return displacement.squaredNorm() < cutoff_radius_sqr_ ? true : false;
}
//=================================================================================================//
bool WithinCutOff::operator()(Real h_ratio, Vecd &displacement) const
{
    return (h_ratio * displacement).squaredNorm() < cutoff_radius_sqr_ ? true : false;
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
Real WendlandC2::W(Real q) const
{
    return pow(1.0 - 0.5 * q, 4) * (1.0 + 2.0 * q);
}
//=================================================================================================//
Real WendlandC2::dW(Real q) const
{
    return 0.625 * pow(q - 2.0, 3) * q;
}
//=================================================================================================//
LaguerreGauss::LaguerreGauss(Real h) : BaseKernel("LaguerreGauss", h, 2.0)
{
    factor_w1d_ = 8.0 / (5.0 * pow(Pi, 0.5)) / h_;
    factor_w2d_ = 3.0 / Pi / h_ / h_;
    factor_w3d_ = 8.0 / pow(Pi, 1.5) / h_ / h_ / h_;
}
//=================================================================================================//
Real LaguerreGauss::W(Real q) const
{
    return (1.0 - pow(q, 2) + pow(q, 4) / 6.0) * exp(-pow(q, 2));
}
//=================================================================================================//
Real LaguerreGauss::dW(Real q) const
{
    return (-pow(q, 5) / 3.0 + 8.0 * pow(q, 3) / 3.0 - 4.0 * q) * exp(-pow(q, 2));
}
//=================================================================================================//
} // namespace SPH
