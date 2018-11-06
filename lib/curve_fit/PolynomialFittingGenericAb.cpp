/*
Copyright 2013, Jernej Kovacic

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
 * Implementation of the class PolynomialFittingGenericAb. This is an abstract 
 * class with implemented some functionality, common to all classes
 * that perform polynomial interpolation, regression, etc.
 */

// no #include "PolynomialFittingGenericAb.hpp" !!!
#include "polynomial/PolynomialGeneric.hpp"


/**
 * @return a copy of the polynomial that best fits entered points.
 * 
 * @throw CurveFittingException if the polynomial has not been generated yet 
 */
template <typename F>
math::PolynomialGeneric<F> math::PolynomialFittingGenericAb<F>::getPolynomial() const
{
    if ( false==this->m_curveGenerated )
    {
        throw math::CurveFittingException(math::CurveFittingException::CURVE_NOT_GENERATED);
    }

    return this->m_poly;
}


/**
 * Evaluates the best fitting polynomial at the given input value 'x'.
 * Sometimes extrapolation (when 'x' is not within bounds defined by the entered 
 * points) is not desired, so by default an exception is thrown on such occurrences.
 * 
 * @param x - input argument to the polynomial
 * @param strict - whether an exception is thrown if 'x' is out of bounds (default: true)
 * 
 * @return value of the polynomial for the given 'x'
 * 
 * @throw CurveFittingException if the polynomial has not been generated yet or 'x' is out of definition range bounds
 */
template <typename F>
F math::PolynomialFittingGenericAb<F>::valueAt(const F& x, const bool strict) const
{
    // the curve must be already generated
    if ( false==this->m_curveGenerated )
    {
        throw math::CurveFittingException(math::CurveFittingException::CURVE_NOT_GENERATED);
    }

    // check if 'x' is within the definition range
    if ( true==strict && (x<this->m_points.front().m_x || x>this->m_points.back().m_x) )
    {
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_BOUNDS);
    }

    // evaluate the polynomial using the PolynomialGeneric's member function
    return this->m_poly.value(x);
}
