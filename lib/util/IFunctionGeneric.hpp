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
 * Declaration of the base class IFunctionGeneric
 */

#ifndef _MATH_IFUNCTIONGENERIC_HPP_
#define _MATH_IFUNCTIONGENERIC_HPP_


#include "exception/FunctionException.hpp"


namespace math
{

/**
 * A base class with the interface of the function to
 * be numerically processed by calculus algorithms.
 *
 * You must derive a class from this one and implement the
 * pure virtual function
 *
 *   T func(const T& x) const throw(FunctionException)
 *
 *  that returns function's value for the given input argument 'x'.
 *  The function is expected to throw FuncException::UNDEFINED
 *  if the value is not defined for the given input argument 'x'.
 *
 *  The instance of the derived class is then passed to numerical
 *  algorithms that call the function 'func' where applicable.
 *
 *  @note 'func' should not be stateful, its return value should only
 *  depend on 'x'.
 *
 *  It is possible to parameterize the class by introducing additional
 *  properties that can be set via setter methods.
 */
template <class T>
class IFunctionGeneric
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
     * @throw FunctionExcpetion::UNDEFINED if the function is not defined at given 'x'
     */
    virtual T func(const T& x) const throw(FunctionException) = 0;

    virtual ~IFunctionGeneric();
};  // class IFunctionGeneric

// Float, double and long double functions make most sense,
// therefore the following types are predefined

typedef IFunctionGeneric<float>         FIFunction;
typedef IFunctionGeneric<double>        IFunction;
typedef IFunctionGeneric<long double>   LDIFunction;

}

// DEFINITION
#include "util/IFunctionGeneric.cpp"

#endif  // _MATH_IFUNCTIONGENERIC_HPP_
