/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file smoothing_kernel.h
 * @brief This is the classes of kernel functions.
 * The kernal function define the relevance between two neighboring particles.
 * Basically, the further the two particles, the less relevance they have.
 * @author	Xiangyu Hu
 */

#ifndef SMOOTHING_KERNELS_H
#define SMOOTHING_KERNELS_H

#include "base_data_package.h"
#include "sph_data_containers.h"

namespace SPH
{
class BaseKernel
{
  protected:
    std::string name_;
    Real h_;             /**< reference smoothing length and its inverse **/
    Real kernel_size_;   /**<kernel_size_ *  h_ gives the zero kernel value */
    Real cutoff_radius_; /** reference cut off radius, beyond this kernel value is neglected. **/
    /** Normalization factors for the kernel function  **/
    Real factor_w1d_, factor_w2d_, factor_w3d_;

  public:
    /** Truncation gives the case that not the entire support of kernel is used */
    BaseKernel(const std::string &name, Real h, Real kernel_size, Real truncation = 1.0);
};

class WendlandC2 : public BaseKernel
{
  public:
    explicit WendlandC2(Real h);
    virtual ~WendlandC2(){};

    Real W(const Real q) const;
    Real dW(const Real q) const;
};

const int KernelResolution = 28;
static constexpr int KernelDataSize = KernelResolution + 4;
using KernelDataArray = std::array<Real, KernelDataSize>;

/**
 * @class TabulatedFunction
 * @brief Four-point Lagrangian interpolation is used to obtain kernel values.
 */
class TabulatedFunction;
{
  public:
    TabulatedFunction(Real dq, std::array<Real, 4> delta_q, KernelDataArray data);
    Real operator()(Real q) const;
    Real operator()(Real h_ratio, Real q) const
    {
        return operator()(q * h_ratio);
    };

  protected:
    Real dq_;
    std::array<Real, 4> delta_q_;
    KernelDataArray data_;
};

class WithinCutOff
{
  public:
    WithinCutOff(Real rc_ref) : rc_ref_sqr_(rc_ref * rc_ref){};
    bool operator()(Vecd &displacement) const
    {
        return displacement.squaredNorm() < rc_ref_sqr_ ? true : false;
    };
    bool operator()(Real h_ratio, Vecd &displacement) const
    {
        return (h_ratio * displacement).squaredNorm() < rc_ref_sqr_ ? true : false;
    };

  protected:
    Real rc_ref_sqr_;
};

class SmoothingKernel : public BaseKernel
{
    WithinCutOff within_cutoff_;           /**< functor to check if particles are within cut off radius */
    TabulatedFunction w1d_, w2d_, w3d_;    /**< kernel value for 1, 2 and 3d **/
    TabulatedFunction dw1d_, dw2d_, dw3d_; /**< kernel derivative for 1, 2 and 3d **/

  public:
    template <typename KernelType>
    explicit SmoothingKernel(const KernelType &kernel)
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

    WithinCutOff getWithinCutOff() const { return within_cutoff_; };

    Real KernelAtOrigin(TypeIdentity<Vec2d> empty_Vec2d) const { return factor_w2d_; };
    TabulatedFunction KernelFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return w2d_; };
    TabulatedFunction KernelDerivativeFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return dw2d_; };

    Real KernelAtOrigin(TypeIdentity<Vec3d> empty_Vec3d) const { return factor_w3d_; };
    TabulatedFunction KernelFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return w3d_; };
    TabulatedFunction KernelDerivativeFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return dw3d_; };

    Real SurfaceKernelAtOrigin(TypeIdentity<Vec2d> empty_Vec2d) const { return factor_w1d_; };
    TabulatedFunction SurfaceKerneFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return w1d_; };
    TabulatedFunction SurfaceKernelDerivativeFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return dw1d_; };

    Real SurfaceKernelAtOrigin(TypeIdentity<Vec3d> empty_Vec3d) const { return factor_w2d_; };
    TabulatedFunction SurfaceKerneFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return w2d_; };
    TabulatedFunction SurfaceKernelDerivativeFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return dw2d_; };

    Real LinearKernelAtOrigin() const { return factor_w1d_; };                  // only for 3D application
    TabulatedFunction LinearKernelFunction() const { return w1d_; };            // only for 3D application
    TabulatedFunction LinearKernelDerivativeFunction() const { return dw1d_; }; // only for 3D application
};
} // namespace SPH
#endif // SMOOTHING_KERNELS_H
