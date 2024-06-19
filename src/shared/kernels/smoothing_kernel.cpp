#include "smoothing_kernel.h"

namespace SPH
{
//=================================================================================================//
Real TabulatedRadialFunction::operator()(Real q)
{
    int location = (int)floor(q / dq_);
    int i = location + 1;
    Real fraction_1 = q - Real(location) * dq_; // fraction_1 correspond to i
    Real fraction_0 = fraction_1 + dq_;         // fraction_0 correspond to i-1
    Real fraction_2 = fraction_1 - dq_;         // fraction_2 correspond to i+1
    Real fraction_3 = fraction_1 - 2 * dq_;     ////fraction_3 correspond to i+2

    return (fraction_1 * fraction_2 * fraction_3) / delta_q_0_ * discreted_data_[i - 1] +
           (fraction_0 * fraction_2 * fraction_3) / delta_q_1_ * discreted_data_[i] +
           (fraction_0 * fraction_1 * fraction_3) / delta_q_2_ * discreted_data_[i + 1] +
           (fraction_0 * fraction_1 * fraction_2) / delta_q_3_ * discreted_data_[i + 2];
}
//=================================================================================================//
} // namespace SPH
