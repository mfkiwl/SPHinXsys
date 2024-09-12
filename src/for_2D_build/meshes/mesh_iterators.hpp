/**
 * @file 	mesh_iterators.hpp
 * @brief 	Here, Functions belong to mesh iterators.
 * @author	hi Zhang and Xiangyu Hu
 */

#pragma once

#include "mesh_iterators.h"

namespace SPH
{
//=================================================================================================//
template <int lower0, int upper0,
          int lower1, int upper1, typename FunctionOnEach>
inline void mesh_for_each2d(const FunctionOnEach &function)
{
    for (int l = lower0; l != upper0; ++l)
        for (int m = lower1; m != upper1; ++m)
        {
            function(l, m);
        }
}
//=================================================================================================//
template <int lower0, int upper0,
          int lower1, int upper1, typename CheckOnEach>
inline Array2i mesh_find_if2d(const CheckOnEach &function)
{
    for (int l = lower0; l != upper0; ++l)
        for (int m = lower1; m != upper1; ++m)
        {
            if (function(l, m))
                return Array2i(l, m);
        }
    return Array2i(upper0, upper1);
}
//=================================================================================================//
template <typename FunctionOnEach>
void mesh_for_each(const Array2i &lower, const Array2i &upper, const FunctionOnEach &function)
{
    for (int l = lower[0]; l != upper[0]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
        {
            function(Array2i(l, m));
        }
}
//=================================================================================================//
template <typename FunctionOnEach>
Array2i mesh_find_if(const Array2i &lower, const Array2i &upper, const FunctionOnEach &function)
{
    for (int l = lower[0]; l != upper[0]; ++l)
        for (int m = lower[1]; m != upper[1]; ++m)
        {
            if (function(l, m))
                return Array2i(l, m);
        }
    return upper;
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
{
    for (int i = (mesh_range.first)[0]; i != (mesh_range.second)[0]; ++i)
        for (int j = (mesh_range.first)[1]; j != (mesh_range.second)[1]; ++j)
        {
            local_function(Array2i(i, j));
        }
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_parallel_for(const MeshRange &mesh_range, const LocalFunction &local_function, Args &&...args)
{
    parallel_for(
        IndexRange2d((mesh_range.first)[0], (mesh_range.second)[0],
                     (mesh_range.first)[1], (mesh_range.second)[1]),
        [&](const IndexRange2d &r)
        {
            for (size_t i = r.rows().begin(); i != r.rows().end(); ++i)
                for (size_t j = r.cols().begin(); j != r.cols().end(); ++j)
                {
                    local_function(Array2i(i, j));
                }
        },
        ap);
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_stride_forward_for(const MeshRange &mesh_range, const Arrayi &stride, const LocalFunction &local_function, Args &&...args)
{
    for (int m = 0; m < stride[0]; m++)
        for (int n = 0; n < stride[1]; n++)
            for (auto i = (mesh_range.first)[0] + m; i < (mesh_range.second)[0]; i += stride[0])
                for (auto j = (mesh_range.first)[1] + n; j < (mesh_range.second)[1]; j += stride[1])
                {
                    local_function(Array2i(i, j));
                }
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_stride_forward_parallel_for(const MeshRange &mesh_range, const Arrayi &stride, const LocalFunction &local_function, Args &&...args)
{
    for (int m = 0; m < stride[0]; m++)
        for (int n = 0; n < stride[1]; n++)
        {
            parallel_for((mesh_range.first)[0] + m, (mesh_range.second)[0], stride[0], [&](auto i)
                         { parallel_for((mesh_range.first)[1] + n, (mesh_range.second)[1], stride[1], [&](auto j)
                                        { local_function(Array2i(i, j)); }, ap); },
                         ap);
        }
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_stride_backward_for(const MeshRange &mesh_range, const Arrayi &stride, const LocalFunction &local_function, Args &&...args)
{
    for (int m = stride[0] - 1; m >= 0; m--)
        for (int n = stride[1] - 1; n >= 0; n--)
            for (auto i = (mesh_range.first)[0] + m; i < (mesh_range.second)[0]; i += stride[0])
                for (auto j = (mesh_range.first)[1] + n; j < (mesh_range.second)[1]; j += stride[1])
                {
                    local_function(Array2i(i, j));
                }
}
//=================================================================================================//
template <typename LocalFunction, typename... Args>
void mesh_stride_backward_parallel_for(const MeshRange &mesh_range, const Arrayi &stride, const LocalFunction &local_function, Args &&...args)
{
    for (int m = stride[0] - 1; m >= 0; m--)
        for (int n = stride[1] - 1; n >= 0; n--)
        {
            parallel_for((mesh_range.first)[0] + m, (mesh_range.second)[0], stride[0], [&](auto i)
                         { parallel_for((mesh_range.first)[1] + n, (mesh_range.second)[1], stride[1], [&](auto j)
                                        { local_function(Array2i(i, j)); }, ap); },
                         ap);
        }
}
//=================================================================================================//
} // namespace SPH
