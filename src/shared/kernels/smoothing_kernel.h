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
    /** Auxiliary factors for the derivative of kernel function  **/
    Real factor_dW_1D_, factor_dW_2D_, factor_dW_3D_;
    /** Auxiliary factors for the second order derivative of kernel function  **/
    Real factor_d2W_1D_, factor_d2W_2D_, factor_d2W_3D_;

    void setDerivativeParameters();

  public:
    BaseKernel(const std::string &name, Real h, Real kernel_size, Real truncation = 1.0);
};

class KernelWendlandC2 : public BaseKernel
{
  public:
    explicit KernelWendlandC2(Real h);
    virtual ~KernelWendlandC2(){};

    Real W(const Real q) const;
    Real dW(const Real q) const;
};

/**
 * @class TabulatedRadialFunction
 * @brief Four-point Lagrangian interpolation is
 * used to obtain kernel values.
 */
class TabulatedRadialFunction;
{
  public:
    static const int resolution_ = 28;
    static constexpr int data_size_ = resolution_ + 4;

    TabulatedFunction(Real dq, std::array<Real, 4> delta_q, std::array<Real, 32> data);
    Real operator()(Real q) const;
    Real operator()(Real h_ratio, Real q) const
    {
        return operator()(q * h_ratio);
    };

  protected:
    Real dq_;
    std::array<Real, 4> delta_q_;
    std::array<Real, data_size_> data_;
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
    Real truncation_; /**< fraction of full support to obtain cut off radius */

    WithinCutOff within_cutoff_;                 /**< functor to check if particles are within cut off radius */
    TabulatedRadialFunction w1d_, w2d_, w3d_;    /**< kernel value for 1, 2 and 3d **/
    TabulatedRadialFunction dw1d_, dw2d_, dw3d_; /**< kernel derivative for 1, 2 and 3d **/

  public:
    template <typename KernelType>
    explicit SmoothingKernel(const KernelType &kernel)
        : BaseKernel(kernel), within_cutoff_(cutoff_radius_)
    {
        using resolution = TabulatedRadialFunction::resolution_;
        using data_size = TabulatedRadialFunction::data_size_;

        Real dq = KernelSize() / Real(resolution);

        std::array<Real, 4> delta_q; // interpolation coefficients
        delta_q[0] = (-1.0 * dq) * (-2.0 * dq) * (-3.0 * dq);
        delta_q[1] = dq * (-1.0 * dq) * (-2.0 * dq);
        delta_q[2] = (2.0 * dq) * dq * (-1.0 * dq);
        delta_q[3] = (3.0 * dq) * (2.0 * dq) * dq;

        std::array<Real, data_size> w1d_data, w2d_data, w3d_data;
        std::array<Real, data_size> dw1d_data, dw2d_data, dw3d_data;
        for (int i = 0; i < data_size; i++)
        {
            w1d_data[i] = factor_w1d_ * kernel.W(Real(i - 1) * dq_);
            w2d_data[i] = factor_w2d_ * kernel.W(Real(i - 1) * dq_);
            w3d_data[i] = factor_w3d_ * kernel.W(Real(i - 1) * dq_);

            dw1d_data[i] = factor_w1d_ * kernel.dW(Real(i - 1) * dq_) / h_;
            dw2d_data[i] = factor_w2d_ * kernel.dW(Real(i - 1) * dq_) / h_;
            dw3d_data[i] = factor_w3d_ * kernel.dW(Real(i - 1) * dq_) / h_;
        }

        w1d_ = TabulatedRadialFunction(dq, delta_q, w1d_data);
        w2d_ = TabulatedRadialFunction(dq, delta_q, w2d_data);
        w3d_ = TabulatedRadialFunction(dq, delta_q, w3d_data);
        dw1d_ = TabulatedRadialFunction(dq, delta_q, dw1d_data);
        dw2d_ = TabulatedRadialFunction(dq, delta_q, dw2d_data);
        dw3d_ = TabulatedRadialFunction(dq, delta_q, dw3d_data);
    }

    WithinCutOff getWithinCutOff() const { return within_cutoff_; };

    Real KernelAtOrigin(TypeIdentity<Vec2d> empty_Vec2d) const { return factor_w2d_; };
    TabulatedRadialFunction KernelFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return w2d_; };
    TabulatedRadialFunction KernelDerivativeFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return dw2d_; };

    Real KernelAtOrigin(TypeIdentity<Vec3d> empty_Vec3d) const { return factor_w3d_; };
    TabulatedRadialFunction KernelFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return w3d_; };
    TabulatedRadialFunction KernelDerivativeFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return dw3d_; };

    Real SurfaceKernelAtOrigin(TypeIdentity<Vec2d> empty_Vec2d) const { return factor_w1d_; };
    TabulatedRadialFunction SurfaceKerneFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return w1d_; };
    TabulatedRadialFunction SurfaceKernelDerivativeFunction(TypeIdentity<Vec2d> empty_Vec2d) const { return dw1d_; };

    Real SurfaceKernelAtOrigin(TypeIdentity<Vec3d> empty_Vec3d) const { return factor_w2d_; };
    TabulatedRadialFunction SurfaceKerneFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return w2d_; };
    TabulatedRadialFunction SurfaceKernelDerivativeFunction(TypeIdentity<Vec3d> empty_Vec3d) const { return dw2d_; };

    Real LinearKernelAtOrigin() const { return factor_w1d_; };
    TabulatedRadialFunction LinearKernelFunction() const { return w1d_; };            // only for 3D application
    TabulatedRadialFunction LinearKernelDerivativeFunction() const { return dw1d_; }; // only for 3D application
};
} // namespace SPH
#endif // SMOOTHING_KERNELS_H
