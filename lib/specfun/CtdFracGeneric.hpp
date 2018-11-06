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

#ifndef _MATH_CTDFRACGENERIC_HPP_
#define _MATH_CTDFRACGENERIC_HPP_


#include <cstddef>

#include "../settings/specfun_settings.h"

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
 *   T fa(const std::size_t i) const throw(FunctionException)
 *
 * and
 *
 *   T fb(const std::size_t i) const throw(FunctionException)
 *
 * that return values of 'a_i' and 'b_i', respectively. 
 * The functions are expected to throw FunctionException::UNDEFINED if
 * the value is not defined at any combination of 'i' and/or 'x'.
 *
 * The instance of the derived class is then passed to the function
 * math::CtdFrac::ctdFrac() that actually evaluates the continued
 * fraction and calls fa() and fb() as applicable.
 *
 * @note 'fa' and 'fb' should not be stateful, their return values
 * should only depend on i'.
 *
 * It is possible to parameterize the class by introducing additional
 * properties that can be set via setter methods.
 */
template <class T>
class ICtdFracFuncGeneric
{

public:

    /*
     * An interface for the function that returns the i^th coefficient
     * 'a_i'. This is a pure virtual function and must be implemented
     * in the derived class.
     *
     * The function should not be stateful, i.e. its output should
     * only depend on 'i'.
     *
     * @param i - number of the coefficient 'a'
     *
     * @return a_i(i)
     *
     * @throw FunctionExcpetion::UNDEFINED if the function is not defined
     */
    virtual T fa(const std::size_t i) const = 0;

    /*
     * An interface for the function that returns the i^th coefficient
     * 'b_i'. This is a pure virtual function and must be implemented
     * in the derived class.
     *
     * The function should not be stateful, i.e. its output should
     * only depend on 'i'.
     *
     * @param i - number of the coefficient 'b'
     *
     * @return b_i(i)
     *
     * @throw FunctionExcpetion::UNDEFINED if the function is not defined
     */
    virtual T fb(const std::size_t i) const = 0;

    virtual ~ICtdFracFuncGeneric();

};  // class ICtdFracFuncGeneric




template <class T>
T ctdFrac(
           const ICtdFracFuncGeneric<T>& ctdf,
           const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
         );

}  // namespace CtdFrac

}  // namespace math


// DEFINITION
#include "specfun/CtdFracGeneric.cpp"

#endif  // _MATH_CTDFRACGENERIC_HPP_
