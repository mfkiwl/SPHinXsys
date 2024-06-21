#ifndef SMOOTHING_KERNEL_HPP
#define SMOOTHING_KERNEL_HPP

#include "smoothing_kernel.h"

namespace SPH
{
//=================================================================================================//
template <typename T>
Real TabulatedFunction::operator()(Real distance, const T &displacement) const
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
template <typename T>
Real TabulatedFunction::operator()(Real h_ratio, Real distance, const T &displacement) const
{
    return operator()(h_ratio *distance, displacement);
}
//=================================================================================================//
template <typename KernelType>
SmoothingKernel::SmoothingKernel(const KernelType &kernel)
    : BaseKernel(kernel)
{
    Real dq = cutoff_radius_ / Real(KernelResolution);
    Real support = kernel_size_ * h_;

    KernelDataArray w1d_data, w2d_data, w3d_data;
    KernelDataArray dw1d_data, dw2d_data, dw3d_data;
    for (int i = 0; i < KernelDataSize; i++)
    {
        Real distance = ABS(Real(i - 1)) * dq; // zero distance value at i=1

        Real value = distance < support ? kernel.W(distance) : 0.0;
        w1d_data[i] = factor_w1d_ * value;
        w2d_data[i] = factor_w2d_ * value;
        w3d_data[i] = factor_w3d_ * value;

        Real derivative = distance < support ? kernel.dW(distance) : 0.0;
        dw1d_data[i] = factor_w1d_ * derivative / h_;
        dw2d_data[i] = factor_w2d_ * derivative / h_;
        dw3d_data[i] = factor_w3d_ * derivative / h_;
    }

    w1d_ = TabulatedFunction(h_, dq, w1d_data);
    w2d_ = TabulatedFunction(h_, dq, w2d_data);
    w3d_ = TabulatedFunction(h_, dq, w3d_data);
    dw1d_ = TabulatedFunction(h_, dq, dw1d_data);
    dw2d_ = TabulatedFunction(h_, dq, dw2d_data);
    dw3d_ = TabulatedFunction(h_, dq, dw3d_data);
}
//=================================================================================================//
} // namespace SPH
#endif // SMOOTHING_KERNEL_HPP
