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
 *
 * Declaration of auxiliary classes and functions that efficiently
 * evaluate continued fractions.
 */

#ifndef _MATH_CTDFRAC_HPP_
#define _MATH_CTDFRAC_HPP_


#include <cstddef>

#include "exception/FunctionException.hpp"
#include "exception/SpecFunException.hpp"


namespace math
{

/*
 * @brief A namespace with auxiliary functions for
 *        evaluation of continued fractions.
 */
namespace  CtdFrac
{


/*
 * A base class with the interface of the functions that return
 * the i^th coefficient 'a_i' and 'b_i' of a continued fraction.
 *
 * You must derive a class from this one and implement the
 * pure virtual functions
 *
 *   T fa(const T& x, size_t i) const throw(FunctionException)
 *
 * and
 *
 *   T fa(const T& x, size_t i) const throw(FunctionException)
 *
 *  that return values of 'a_i' and 'b_i', respectively, optionally
 *  depending on 'x'. The functions are expected to throw
 *  FuncException::UNDEFINED if the value is not defined for the given
 *  input argument 'x'.
 *
 *  The instance of the derived class is then passed to the function
 *  math::CtdFrac::ctdFrac() that actually evaluates the continued
 *  fraction and calls fa() and fb() as applicable.
 *
 *  @note 'fa' and 'fb' should not be stateful, their return values
 *  should only depend on 'x' and 'i'.
 *
 *  It is possible to parameterize the class by introducing additional
 *  properties that can be set via setter methods.
 */
template <class T>
class ICtdFracFuncGeneric
{

public:

    /*
     * An interface for the function that returns the i^th coefficient
     * 'a_i' (may also depend on 'x'). This is a pure virtual function
     * and must be implemented in the derived class.
     *
     * The function should not be stateful, i.e. its output should
     * only depend on 'x' and 'i'.
     *
     * @param x - input argument
     * @param i - number of the coefficient 'a'
     *
     * @return a_i(x)
     *
     * @throw FunctionExcpetion::UNDEFINED if the function is not defined at given 'x'
     */
    virtual T fa(const T& x, size_t i) const throw (FunctionException) = 0;

    /*
     * An interface for the function that returns the i^th coefficient
     * 'b_i' (may also depend on 'x'). This is a pure virtual function
     * and must be implemented in the derived class.
     *
     * The function should not be stateful, i.e. its output should
     * only depend on 'x' and 'i'.
     *
     * @param x - input argument
     * @param i - number of the coefficient 'b'
     *
     * @return b_i(x)
     *
     * @throw FunctionExcpetion::UNDEFINED if the function is not defined at given 'x'
     */
    virtual T fb(const T& x, size_t i) const throw (FunctionException) = 0;

    virtual ~ICtdFracFuncGeneric();

};  // class ICtdFracFuncGeneric




template <class T>
T ctdFrac(
           const ICtdFracFuncGeneric<T>& ctdf,
           const T& x,
           const T& tol = static_cast<T>(1)/static_cast<T>(1000000)
         ) throw(SpecFunException);

}  // namespace CtdFrac

}  // namespace math


// DEFINITION
#include "specfun/CtdFracGeneric.cpp"

#endif  // _MATH_CTDFRAC_HPP_
