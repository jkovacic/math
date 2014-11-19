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

#include "exception/IntegException.hpp"


namespace math
{

/**
 * A base class with the interface of the function to
 * be numerically integrated.
 *
 * You must derive a class from this one and implement the
 * pure virtual function
 *
 *   T func(const T& x) const throw(IntegException)
 *
 *  that returns function's value for the given input argument 'x'.
 *  The function is expected to twrow IntegException::UNDEFINED
 *  if the value is not defined for the given input argument 'x'.
 *
 *  The instance of the derived class is then passed to integration
 *  algorithms that call the function 'func' where applicable.
 *
 *  @note 'func' should not be stateful, its return value should only
 *  depend on 'x'.
 *
 *  It is possible to parameterize the class by introducing additional
 *  properties that can be set via setter methods.
 */
template <class T>
class IIntegFunctionGeneric
{

public:

	/**
	 * An interface for the function to be numerically integrated.
	 * This is a pure virtual function and must be implemented
	 * in the derived class.
	 *
	 * The function should not be stateful, i.e. its output should
	 * only depend on 'x'.
	 *
	 * @param x - input argument
	 *
	 * @return func(x)
	 *
	 * @throw IntegExcpetion::UNDEFINED if the function is not defined at given 'x'
	 */
    virtual T func(const T& x) const throw(IntegException) = 0;

    virtual ~IIntegFunctionGeneric();
};  // class IIntegFunctionGeneric


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
               const IIntegFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n,
               EIntegAlg::alg algorithm=EIntegAlg::SIMPSON
             ) throw(IntegException);

    static T integH(
               const IIntegFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               const T& h,
               EIntegAlg::alg algorithm=EIntegAlg::SIMPSON
             ) throw(IntegException);


private:

    static T __rectangle(
               const IIntegFunctionGeneric<T>& f,
               const T&a,
               const T& b,
               size_t n
             ) throw(IntegException);

    static T __trapezoidal(
               const IIntegFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(IntegException);

    static T __simpson(
               const IIntegFunctionGeneric<T>& f,
               const T&a,
               const T& b,
               size_t n
             ) throw(IntegException);

    static T __simpson38(
               const IIntegFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n
             ) throw(IntegException);

};  // class IntegGeneric

// Float, double and long double integrators make most sense,
// therefore the following types are predefined
typedef IntegGeneric<float>          FInteg;
typedef IntegGeneric<double>         Integ;
typedef IntegGeneric<long double>    LDInteg;

typedef IIntegFunctionGeneric<float>         FIIntegFunction;
typedef IIntegFunctionGeneric<double>        IIntegFunction;
typedef IIntegFunctionGeneric<long double>   LDIIntegFunction;

}

// DEFINITION
#include "integ/IntegGeneric.cpp"

#endif  // _MATH_INTEG_GENERIC_HPP_
