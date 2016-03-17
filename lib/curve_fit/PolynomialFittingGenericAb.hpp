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
 * An internal header file, it should not be included directly.
 *
 * Declaration of the class PolynomialFittingGenericAb.h. This is an 
 * abstract class with some functionality common to all derived classes
 * that perform polynomial interpolation, regression, etc.
 */

#ifndef _MATH_POLYNOMIALFITTINGGENERICAB_HPP_
#define _MATH_POLYNOMIALFITTINGGENERICAB_HPP_

#include "CurveFittingGenericAb.hpp"
#include "polynomial/PolynomialGeneric.hpp"

namespace math
{

// Some compilers may require forward declaration
// of the templated base class:
template <typename F> class CurveFittingGenericAb;


/**
 * @brief Common functionality for curve fitting/interpolation algorithms that
 *        find the best fitting polynomial.
 * 
 * The class is abstract and cannot be instantiated. More specialized algorithms
 * are implemented by its derived classes.
 */
template <typename F>
class PolynomialFittingGenericAb : public CurveFittingGenericAb<F>
{

protected:
    // Curve fitting polynomial
    PolynomialGeneric<F> m_poly;

public:
    // get a copy of the best fitting polynomial
    PolynomialGeneric<F> getPolynomial() const throw (CurveFittingException);
    // inherited as an abstract function from the base class, value of the polynomial at the given abscissa
    F valueAt(const F& x, const bool strict=true) const throw (CurveFittingException);

    //generateCurve() remains a pure virtual method and must be implemented by derived classes
};

}  // namespace math

// DEFINITION
#include "curve_fit/PolynomialFittingGenericAb.cpp"

#endif  // _MATH_POLYNOMIALFITTINGGENERICAB_HPP_
