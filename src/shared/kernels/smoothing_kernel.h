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
/**
 * @class TabulatedRadialFunction
 * @brief Four-point Lagrangian interpolation is
 * used to obtain kernel values.
 */
class TabulatedRadialFunction;
{
  public:
    template <typename OriginalRadialFunction>
    TabulatedFunction();
    Real operator()(Real q) const;
    Real operator()(Real h_ratio, Real q) const;

  protected:
    Real dq_, delta_q_0_, delta_q_1_, delta_q_2_, delta_q_3_;
    std::array<Real, 32> discreted_data_;
};

class withinCutOff
{
  public:
    withinCutOff(Real rc_ref) : rc_ref_sqr_(rc_ref * rc_ref){};
    bool operator()(Vecd &displacement) const
    {
        return displacement.squaredNorm() < rc_ref_sqr_ ? true : false;
    };
    bool operator()(Real h_ratio, Vecd &displacement) const;

  protected:
    Real rc_ref_sqr_;
};

struct SmoothingKernel
{
    std::string name_;
    Real h_;           /**< reference smoothing length **/
    Real kernel_size_; /**<kernel_size_ *  h_ gives the zero kernel value */
    Real truncation_;  /**< to obtain cut off radius */
    Real rc_ref_;      /** reference cut off radius, beyond this kernel value is neglected. **/

    withinCutoff within_cutoff_;                 /**< functor to check if particles are within cut off radius */
    TabulatedRadialFunction w1d_, w2d_, w3d_;    /**< kernel value for 1, 2 and 3d **/
    TabulatedRadialFunction dw1d_, dw2d_, dw3d_; /**< kernel derivative for 1, 2 and 3d **/
    TabulatedRadialFunction dw1d_, dw2d_, dw3d_; /**< kernel 2nd-order derivative for 1, 2 and 3d **/

    template <typename KernelType>
    explicit SmoothingKernel(KernelType kernel);

    //----------------------------------------------------------------------
    //		Below are for variable smoothing length.
    //		Note that we input the ratio between the reference smoothing length
    //		to the variable smoothing length.
    //----------------------------------------------------------------------
  protected:
    /** Functor for variable smoothing length. */
    typedef std::function<Real(const Real &)> FactorFunctor;
    FactorFunctor h_factor_W_1D_, h_factor_W_2D_, h_factor_W_3D_;
    FactorFunctor h_factor_dW_1D_, h_factor_dW_2D_, h_factor_dW_3D_;
    FactorFunctor h_factor_d2W_1D_, h_factor_d2W_2D_, h_factor_d2W_3D_;

    Real factorW1D(const Real &h_ratio) const { return h_ratio; };
    Real factorW2D(const Real &h_ratio) const { return h_ratio * h_ratio; };
    Real factorW3D(const Real &h_ratio) const { return h_ratio * h_ratio * h_ratio; };
    Real factordW1D(const Real &h_ratio) const { return factorW1D(h_ratio) * h_ratio; };
    Real factordW2D(const Real &h_ratio) const { return factorW2D(h_ratio) * h_ratio; };
    Real factordW3D(const Real &h_ratio) const { return factorW3D(h_ratio) * h_ratio; };
    Real factord2W1D(const Real &h_ratio) const { return factordW1D(h_ratio) * h_ratio; };
    Real factord2W2D(const Real &h_ratio) const { return factordW2D(h_ratio) * h_ratio; };
    Real factord2W3D(const Real &h_ratio) const { return factordW3D(h_ratio) * h_ratio; };

  public:
    Real CutOffRadius(Real h_ratio) const { return rc_ref_ / h_ratio; };
    Real CutOffRadiusSqr(Real h_ratio) const { return rc_ref_sqr_ / (h_ratio * h_ratio); };

    Real W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const;
    Real W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const;
    Real W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const;

    /** Calculates the kernel value at the origin **/
    Real W0(const Real &h_ratio, const Real &point_i) const;
    Real W0(const Real &h_ratio, const Vec2d &point_i) const;
    Real W0(const Real &h_ratio, const Vec3d &point_i) const;

    /** Calculates the kernel derivation for the given distance of two particles **/
    Real dW(const Real &h_ratio, const Real &r_ij, const Real &displacement) const;
    Real dW(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const;
    Real dW(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const;

    /** Calculates the kernel second order derivation for the given distance of two particles **/
    Real d2W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const;
    Real d2W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const;
    Real d2W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const;
    //----------------------------------------------------------------------
    //		Below are for reduced kernels.
    //----------------------------------------------------------------------
  public:
    void reduceOnce();  /** reduce for thin structures or films */
    void reduceTwice(); /** reduce for linear structures or filaments */
};
} // namespace SPH
#endif // SMOOTHING_KERNELS_H
