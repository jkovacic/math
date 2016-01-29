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
 * @headername{RootFindGeneric.h}
 *
 * Declaration of functions within the namespace RootFind with several
 * root finding algorithms.
 */

#ifndef _MATH_ROOT_FIND_GENERIC_HPP_
#define _MATH_ROOT_FIND_GENERIC_HPP_

#include <cstddef>

#include "../settings/rootfind_settings.h"
#include "exception/RootFindException.hpp"
#include "util/IFunctionGeneric.hpp"


namespace math
{

/**
 * @brief A namespace with functions that implement several root finding
 *        algorithms.
 */
namespace RootFind
{

    template <typename F>
    F bisection(
           const IFunctionGeneric<F>& f,
           const F& from,
           const F& to,
           const F& epsx = static_cast<F>(ROOTFIND_EPSX_NUM) / static_cast<F>(ROOTFIND_EPSX_DEN),
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN)
         ) throw (RootFindException);


    template <typename F>
    F regulaFalsi(
           const IFunctionGeneric<F>& f,
           const F& from,
           const F& to,
           const F& epsx = static_cast<F>(ROOTFIND_EPSX_NUM) / static_cast<F>(ROOTFIND_EPSX_DEN),
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN)
         ) throw (RootFindException);


    template <typename F>
    F secant(
           const IFunctionGeneric<F>& f,
           const F& r0,
           const F& r1,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <typename F>
    F newton(
           const IFunctionGeneric<F>& f,
           const IFunctionGeneric<F>& diff,
           const F& x0,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <typename F>
    F quasiNewton(
           const IFunctionGeneric<F>& f,
           const F& x0,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const F& h = static_cast<F>(ROOTFIND_DIFF_STEP_NUM) / static_cast<F>(ROOTFIND_DIFF_STEP_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <typename F>
    F halley(
           const IFunctionGeneric<F>& f,
           const IFunctionGeneric<F>& diff,
           const IFunctionGeneric<F>& diff2,
           const F& x0,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <typename F>
    F quasiHalley(
           const IFunctionGeneric<F>& f,
           const F& x0,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const F& h = static_cast<F>(ROOTFIND_DIFF_STEP_NUM) / static_cast<F>(ROOTFIND_DIFF_STEP_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <typename F>
    F halleyMod(
           const IFunctionGeneric<F>& f,
           const IFunctionGeneric<F>& diff,
           const IFunctionGeneric<F>& diff2,
           const F& x0,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <typename F>
    F quasiHalleyMod(
           const IFunctionGeneric<F>& f,
           const F& x0,
           const F& epsy = static_cast<F>(ROOTFIND_EPSY_NUM) / static_cast<F>(ROOTFIND_EPSY_DEN),
           const F& h = static_cast<F>(ROOTFIND_DIFF_STEP_NUM) / static_cast<F>(ROOTFIND_DIFF_STEP_DEN),
           const std::size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);

}  // namepsace RootFind

}  // namespace math


// DEFINITION
#include "root_find/RootFindGeneric.cpp"

#endif  // _MATH_ROOT_FIND_GENERIC_HPP_
