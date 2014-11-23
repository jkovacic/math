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
 * Declaration of the class RootFindGeneric with several
 * root finding algorithms.
 */

#ifndef _MATH_ROOT_FIND_GENERIC_HPP_
#define _MATH_ROOT_FIND_GENERIC_HPP_

#include <cstddef>

#include "exception/RootFindException.hpp"
#include "util/IFunctionGeneric.hpp"
#include "util/NumericUtil.hpp"


namespace math
{

/**
 * @brief A class with implementation of several root finding
 *        algorithms.
 *
 * The whole class is static, so no instantiation is necessary.
 */
template <class T>
class RootFindGeneric
{

public:
    static T bisection(
           const IFunctionGeneric<T>& f,
           const T& from,
           const T& to,
           const T& EPSX = NumericUtil<T>::getEPS(),
           const T& EPSY = NumericUtil<T>::getEPS()
         ) throw(RootFindException);


    static T regulaFalsi(
           const IFunctionGeneric<T>& f,
           const T& from,
           const T& to,
           const T& EPSX = NumericUtil<T>::getEPS(),
           const T& EPSY = NumericUtil<T>::getEPS()
         ) throw(RootFindException);


    static T secant(
           const IFunctionGeneric<T>& f,
           const T& r0,
           const T& r1,
           const T& EPSY = NumericUtil<T>::getEPS(),
           size_t Nmax = 10000
         ) throw(RootFindException);


    static T newton(
           const IFunctionGeneric<T>& f,
           const IFunctionGeneric<T>& diff,
           const T& x0,
           const T& EPSY = NumericUtil<T>::getEPS(),
           size_t Nmax = 10000
         ) throw(RootFindException);


    static T quasiNewton(
           const IFunctionGeneric<T>& f,
           const T& x0,
           const T& EPSY = NumericUtil<T>::getEPS(),
           const T& h = static_cast<T>(1) / static_cast<T>(1000),
           size_t Nmax = 10000
         ) throw(RootFindException);


private:

    static T __newtonCommon(
           const IFunctionGeneric<T>& f,
           const IFunctionGeneric<T>& diff,
           const T& x0,
           const T& EPSY,
           const T& h,
           size_t Nmax,
           bool diffFunc
         ) throw(RootFindException);
};  // class RootFindGeneric

// Functions of type float, double and long double
// make most sense, therefore the following types are predefined
typedef RootFindGeneric<float>        FRootFind;
typedef RootFindGeneric<double>       RootFind;
typedef RootFindGeneric<long double>  LDRootFind;

}  // namespace math


// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "root_find/RootFindGeneric.cpp"

#endif  // _MATH_ROOT_FIND_GENERIC_HPP_
