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
 * @headername{IntegGeneric.h}
 *
 * Declaration of the namespace Integ with functions that
 * perform numerical integration of continuous functions.
 */

#ifndef _MATH_INTEG_GENERIC_HPP_
#define _MATH_INTEG_GENERIC_HPP_

#include <cstddef>

#include "../settings/calc_settings.h"
#include "exception/CalculusException.hpp"
#include "util/IFunctionGeneric.hpp"


namespace math
{

/**
 * @brief Supported algorithms for numerical calculation
 *        of definite integrals
 */
struct EIntegAlg
{
    enum alg
    {
        RECTANGLE,              /// Rectangle method
        TRAPEZOIDAL,            /// Trapezoidal rule
        SIMPSON,                /// Simpson's rule
        SIMPSON_3_8,            /// Simpson's 3/8 rule
        BOOLE,                  /// Boole's rule
    };
};  // struct EIntegAlg


/**
 * A namespace with several implemented algorithms for
 * numerical integration.
 */
namespace Integ
{

    template <typename F>
    F integ(
               const IFunctionGeneric<F>& f,
               const F& a,
               const F& b,
               size_t n = INTEG_DEFAULT_STEPS,
               EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
             ) throw(CalculusException);


    template <typename F>
    F integH(
               const IFunctionGeneric<F>& f,
               const F& a,
               const F& b,
               const F& h = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
               EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
             ) throw(CalculusException);


    template <typename F>
    F integImpNegInf(
               const IFunctionGeneric<F>& f,
               const F& b,
               size_t nimp = INTEG_DEFAULT_STEPS,
               size_t nprop = INTEG_DEFAULT_STEPS,
               const F& bp = -static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
               EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
             ) throw(CalculusException);


    template <typename F>
    F integImpPosInf(
                   const IFunctionGeneric<F>& f,
                   const F& a,
                   size_t nimp = INTEG_DEFAULT_STEPS,
                   size_t prop = INTEG_DEFAULT_STEPS,
                   const F& bp = static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
                 ) throw(CalculusException);


    template <typename F>
    F integImp(
                   const IFunctionGeneric<F>& f,
                   size_t nimp = INTEG_DEFAULT_STEPS,
                   size_t nprop = INTEG_DEFAULT_STEPS,
                   const F& bpneg = -static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   const F& bppos = static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
                 ) throw(CalculusException);


    template <typename F>
    F integImpNegInfH(
                   const IFunctionGeneric<F>& f,
                   const F& b,
                   const F& himp = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
                   const F& hprop = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
                   const F& bp = -static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
                 ) throw(CalculusException);


    template <typename F>
    F integImpPosInfH(
                   const IFunctionGeneric<F>& f,
                   const F& a,
                   const F& himp = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
                   const F& hprop = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
                   const F& bp = static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
                 ) throw(CalculusException);


    template <typename F>
    F integImpH(
                   const IFunctionGeneric<F>& f,
                   const F& himp = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
                   const F& hprop = static_cast<F>(INTEG_STEP_NUM) / static_cast<F>(INTEG_STEP_DEN),
                   const F& bpneg = -static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   const F& bppos = static_cast<F>(INTEG_IMP_INT_BREAKPOINT_NUM) / static_cast<F>(INTEG_IMP_INT_BREAKPOINT_DEN),
                   EIntegAlg::alg algorithm = INTEG_DEFAULT_METHOD
                 ) throw(CalculusException);

}  // namespace Integ


}  // namespace math

// DEFINITION
#include "calculus/IntegGeneric.cpp"

#endif  // _MATH_INTEG_GENERIC_HPP_
