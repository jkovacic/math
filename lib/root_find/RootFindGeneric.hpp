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

    template <class T>
    T bisection(
           const IFunctionGeneric<T>& f,
           const T& from,
           const T& to,
           const T& epsx = static_cast<T>(ROOTFIND_EPSX_NUM) / static_cast<T>(ROOTFIND_EPSX_DEN),
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN)
         ) throw (RootFindException);


    template <class T>
    T regulaFalsi(
           const IFunctionGeneric<T>& f,
           const T& from,
           const T& to,
           const T& epsx = static_cast<T>(ROOTFIND_EPSX_NUM) / static_cast<T>(ROOTFIND_EPSX_DEN),
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN)
         ) throw (RootFindException);


    template <class T>
    T secant(
           const IFunctionGeneric<T>& f,
           const T& r0,
           const T& r1,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <class T>
    T newton(
           const IFunctionGeneric<T>& f,
           const IFunctionGeneric<T>& diff,
           const T& x0,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <class T>
    T quasiNewton(
           const IFunctionGeneric<T>& f,
           const T& x0,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           const T& h = static_cast<T>(ROOTFIND_DIFF_STEP_NUM) / static_cast<T>(ROOTFIND_DIFF_STEP_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <class T>
    T halley(
           const IFunctionGeneric<T>& f,
           const IFunctionGeneric<T>& diff,
           const IFunctionGeneric<T>& diff2,
           const T& x0,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <class T>
    T quasiHalley(
           const IFunctionGeneric<T>& f,
           const T& x0,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           const T& h = static_cast<T>(ROOTFIND_DIFF_STEP_NUM) / static_cast<T>(ROOTFIND_DIFF_STEP_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <class T>
    T halleyMod(
           const IFunctionGeneric<T>& f,
           const IFunctionGeneric<T>& diff,
           const IFunctionGeneric<T>& diff2,
           const T& x0,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);


    template <class T>
    T quasiHalleyMod(
           const IFunctionGeneric<T>& f,
           const T& x0,
           const T& epsy = static_cast<T>(ROOTFIND_EPSY_NUM) / static_cast<T>(ROOTFIND_EPSY_DEN),
           const T& h = static_cast<T>(ROOTFIND_DIFF_STEP_NUM) / static_cast<T>(ROOTFIND_DIFF_STEP_DEN),
           size_t Nmax = ROOTFIND_MAX_ITER
         ) throw (RootFindException);

}  // namepsace RootFind

}  // namespace math


// DEFINITION
#include "root_find/RootFindGeneric.cpp"

#endif  // _MATH_ROOT_FIND_GENERIC_HPP_
