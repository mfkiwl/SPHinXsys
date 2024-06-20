#ifndef SMOOTHING_KERNEL_HPP
#define SMOOTHING_KERNEL_HPP

#include "smoothing_kernel.h"

namespace SPH
{
//=================================================================================================//
template <typename KernelType>
SmoothingKernel::SmoothingKernel(const KernelType &kernel)
    : BaseKernel(kernel), within_cutoff_(cutoff_radius_)
{
    Real dq = cutoff_radius_ / Real(KernelResolution);
    Real support = kernel_size_ * h_;

    std::array<Real, 4> delta_q; // interpolation coefficients
    delta_q[0] = (-1.0 * dq) * (-2.0 * dq) * (-3.0 * dq);
    delta_q[1] = dq * (-1.0 * dq) * (-2.0 * dq);
    delta_q[2] = (2.0 * dq) * dq * (-1.0 * dq);
    delta_q[3] = (3.0 * dq) * (2.0 * dq) * dq;

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

    w1d_ = TabulatedFunction(dq, delta_q, w1d_data);
    w2d_ = TabulatedFunction(dq, delta_q, w2d_data);
    w3d_ = TabulatedFunction(dq, delta_q, w3d_data);
    dw1d_ = TabulatedFunction(dq, delta_q, dw1d_data);
    dw2d_ = TabulatedFunction(dq, delta_q, dw2d_data);
    dw3d_ = TabulatedFunction(dq, delta_q, dw3d_data);
}
//=================================================================================================//
} // namespace SPH
#endif // SMOOTHING_KERNEL_HPP
