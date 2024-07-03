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
 * @brief This is a generalized implementation of SPH smoothing kernel.
 * The kernal function defines the relevance between two neighboring particles.
 * Basically, the further the two particles, the less relevance they have.
 * Note that the smoothing kernel here is defined on an isotropic uniform
 * space without deformation. For anisotropic space with adaptive resolutions,
 * the complication will be introduced in other classes.
 * Also note, in order to generalize for different kernel types,
 * the kernel profiles is mapped onto an array so that the evaluation
 * process is dependent on an interpolation only.
 * @author	Xiangyu Hu
 */

#ifndef SMOOTHING_KERNEL_H
#define SMOOTHING_KERNEL_H

#include "base_data_package.h"
#include "sph_data_containers.h"

namespace SPH
{
class BaseKernel
{
  protected:
    std::string name_;
    Real h_;             /**< reference smoothing length and its inverse **/
    Real kernel_size_;   /**< kernel_size_ *  h_ gives the zero kernel value */
    Real cutoff_radius_; /**< reference cut off radius, beyond this kernel value is neglected. **/
    /** Normalization factors for the kernel function  **/
    Real factor_w1d_, factor_w2d_, factor_w3d_;

  public:
    /** Truncation gives the case that not the entire support of kernel is used */
    BaseKernel(const std::string &name, Real h, Real kernel_size, Real truncation = 1.0);
    std::string Name() const { return name_; };
    Real SmoothingLength() const { return h_; };
    Real KernelSize() const { return kernel_size_; };
    Real CutOffRadius() const { return cutoff_radius_; };
};

class WendlandC2 : public BaseKernel
{
  public:
    explicit WendlandC2(Real h);
    virtual ~WendlandC2(){};

    Real W(Real q) const;
    Real dW(Real q) const;
};

class LaguerreGauss : public BaseKernel
{
  public:
    explicit LaguerreGauss(Real h);
    virtual ~LaguerreGauss(){};

    Real W(Real q) const;
    Real dW(Real q) const;
};

const int KernelResolution = 20;
constexpr int KernelDataSize = KernelResolution + 4;
using KernelDataArray = std::array<Real, KernelDataSize>;

/**
 * @class RadialFunction
 * @brief Four-point Lagrangian interpolation is used to obtain kernel values.
 */
class RadialFunction
{
  public:
    RadialFunction(Real h, Real dq, KernelDataArray data);
    template <typename T>
    Real operator()(Real distance, const T &displacement) const;
    template <typename ScalingType, typename T>
    Real operator()(const ScalingType &scaling, Real distance, const T &displacement) const;

  protected:
    Real inv_h_, dq_;
    std::array<Real, 4> delta_q_; // interpolation coefficients
    KernelDataArray data_;
};

class WithinCutOff
{
  public:
    WithinCutOff(Real cutoff_radius)
        : cutoff_radius_sqr_(cutoff_radius * cutoff_radius){};
    bool operator()(Vecd &displacement) const; // for constant smoothing length
    template <typename ScalingType>
    bool operator()(const ScalingType &scaling, Vecd &displacement) const; // for variable smoothing length

  protected:
    Real cutoff_radius_sqr_;
};

class SmoothingKernel : public BaseKernel
{
    RadialFunction w1d_, w2d_, w3d_;    /**< kernel value for 1, 2 and 3d **/
    RadialFunction dw1d_, dw2d_, dw3d_; /**< kernel derivative for 1, 2 and 3d **/

  public:
    template <typename KernelType>
    explicit SmoothingKernel(const KernelType &kernel);

    RadialFunction KernelFunction(Vec2d zero) const { return w2d_; };
    RadialFunction KernelDerivativeFunction(Vec2d zero) const { return dw2d_; };

    RadialFunction KernelFunction(Vec3d zero) const { return w3d_; };
    RadialFunction KernelDerivativeFunction(Vec3d zero) const { return dw3d_; };

    RadialFunction SurfaceKerneFunction(Vec2d zero) const { return w1d_; };
    RadialFunction SurfaceKernelDerivativeFunction(Vec2d zero) const { return dw1d_; };

    RadialFunction SurfaceKerneFunction(Vec3d zero) const { return w2d_; };
    RadialFunction SurfaceKernelDerivativeFunction(Vec3d zero) const { return dw2d_; };

    RadialFunction LinearKernelFunction() const { return w1d_; };            // only for 3D application
    RadialFunction LinearKernelDerivativeFunction() const { return dw1d_; }; // only for 3D application
};
} // namespace SPH
#endif // SMOOTHING_KERNEL_H
