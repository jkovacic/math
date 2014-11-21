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
 * Declaration of the class IntegGeneric and related classes that
 * perform numerical integration of continuous functions.
 */

#ifndef _MATH_INTEG_GENERIC_HPP_
#define _MATH_INTEG_GENERIC_HPP_

#include <cstddef>

#include "exception/CalculusException.hpp"
#include "exception/FunctionException.hpp"
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
 * A class with several implemented algorithms for
 * numerical integration.
 *
 * The class is static, so no instantiation is necessary.
 */
template <class T>
class IntegGeneric
{

public:

    static T integ(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n,
               EIntegAlg::alg algorithm=EIntegAlg::SIMPSON
             ) throw(CalculusException);

    static T integH(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               const T& h,
               EIntegAlg::alg algorithm=EIntegAlg::SIMPSON
             ) throw(CalculusException);


private:

    static T __rectangle(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(FunctionException);

    static T __trapezoidal(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(FunctionException);

    static T __simpson(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(FunctionException);

    static T __simpson38(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(FunctionException);

    static T __boole(
               const IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(FunctionException);

};  // class IntegGeneric

// Float, double and long double integrators make most sense,
// therefore the following types are predefined
typedef IntegGeneric<float>          FInteg;
typedef IntegGeneric<double>         Integ;
typedef IntegGeneric<long double>    LDInteg;

}

// DEFINITION
#include "calculus/IntegGeneric.cpp"

#endif  // _MATH_INTEG_GENERIC_HPP_
