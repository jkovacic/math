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
 * Declaration of the class DiffGeneric and related classes that
 * perform numerical differentiation of continuous functions.
 */


#ifndef _MATH_DIFF_GENERIC_HPP_
#define _MATH_DIFF_GENERIC_HPP_

#include <cstddef>

#include "exception/CalculusException.hpp"
#include "exception/FunctionException.hpp"
#include "util/IFunctionGeneric.hpp"


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
 * A class with several implemented algorithms for
 * numerical differentiation.
 *
 * The class is static, so no instantiation is necessary.
 */
template <class T>
class DiffGeneric
{

public:

    static T diff(
               const IFunctionGeneric<T>& f,
               const T& x,
               const T& h,
               EDiffMethod::method method=EDiffMethod::CENTRAL
             ) throw(CalculusException);


private:

    static T __forwardDiff(
               const IFunctionGeneric<T>& f,
               const T& x,
               const T& h
             ) throw(FunctionException);

    static T __backwardDiff(
               const IFunctionGeneric<T>& f,
               const T& x,
               const T& h
             ) throw(FunctionException);

    static T __centralDiff(
               const IFunctionGeneric<T>& f,
               const T& x,
               const T& h
             ) throw(FunctionException);

    static T __5pointDiff(
               const IFunctionGeneric<T>& f,
               const T& x,
               const T& h
             ) throw(FunctionException);

};  // class DiffGeneric

// Float, double and long double differentiators make most sense,
// therefore the following types are predefined
typedef DiffGeneric<float>          FDiff;
typedef DiffGeneric<double>         Diff;
typedef DiffGeneric<long double>    LDDiff;

}

// DEFINITION
#include "calculus/DiffGeneric.cpp"

#endif  // _MATH_DIFF_GENERIC_HPP_
