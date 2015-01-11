/*
Copyright 2014, Jernej Kovacic

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/**
 * @file
 * @author Jernej Kovacic
 *
 * An internal header file, it should not be included directly.
 * @headername{DiffGeneric.h}
 *
 * Declaration of functions within the namespace Diff with that
 * perform numerical differentiation of continuous functions.
 */


#ifndef _MATH_DIFF_GENERIC_HPP_
#define _MATH_DIFF_GENERIC_HPP_

#include <cstddef>

#include "../settings/calc_settings.h"
#include "util/IFunctionGeneric.hpp"
#include "exception/CalculusException.hpp"
#include "exception/FunctionException.hpp"


namespace math
{

/**
 * @brief Supported algorithms for numerical differentiation
 */
struct EDiffMethod
{
    enum method
    {
        FORWARD,             /// Forward difference method
        BACKWARD,            /// Backward difference method
        CENTRAL,             /// Central difference method
        FIVE_POINT,          /// 5 point method
    };
};  // struct EDiffMethod


/**
 * A namespace with several implemented algorithms for
 * numerical differentiation.
 * 
 * @note Numerical differentiation is typically not very accurate.
 *       Make sure that the step 'h' is not too small to avoid
 *       large floating point rounding errors.
 */
namespace Diff
{

    template <typename F>
    F diff(
               const IFunctionGeneric<F>& f,
               const F& x,
               const F& h = static_cast<F>(DIFF_STEP_NUM) / static_cast<F>(DIFF_STEP_DEN),
               EDiffMethod::method method = DIFF_DEFAULT_METHOD
             ) throw (CalculusException);


    template <typename F>
    F diff2(
               const IFunctionGeneric<F>& f,
               const F& x,
               const F& h = static_cast<F>(DIFF_STEP_NUM) / static_cast<F>(DIFF_STEP_DEN)
             ) throw (CalculusException);



}  // namespace Diff

}  // namespace math

// DEFINITION
#include "calculus/DiffGeneric.cpp"

#endif  // _MATH_DIFF_GENERIC_HPP_
