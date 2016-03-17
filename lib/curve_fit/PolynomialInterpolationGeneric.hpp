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
 * @headername{PolynomialInterpolation.h}
 *
 * Declaration of the class PolynomialInterpolationGeneric.h that calculates 
 * an interpolation polynomial that goes exactly through entered points.
 */

#ifndef _MATH_POLYNOMIALINTERPOLATIONGENERIC_HPP_
#define _MATH_POLYNOMIALINTERPOLATIONGENERIC_HPP_

#include <cstddef>

#include "PolynomialFittingGenericAb.hpp"
#include "exception/CurveFittingException.hpp"


namespace math
{

// Some compilers may require forward declaration
// of the templated base class:
template <typename F> class PolynomialFittingGenericAb;


/**
 * @brief A class that finds a polynomial that exactly fits an input series of points.
 * 
 * The algorithm takes N+1 points and calculates a N-degree polynomial exactly
 * fitting all input points.
 */
template <typename F>
class PolynomialInterpolationGeneric : public PolynomialFittingGenericAb<F>
{
public:
    // Constructor
    PolynomialInterpolationGeneric();

    // inherited as an abstract function from the base class
    void generateCurve(const std::size_t degree = 0) throw (CurveFittingException);
};

// Interpolation classes with elements of types float, double and long double make
// most sense so the following three types are predefined:
typedef PolynomialInterpolationGeneric<float>       FPolynomialInterpolation;
typedef PolynomialInterpolationGeneric<double>      PolynomialInterpolation;
typedef PolynomialInterpolationGeneric<long double> LDPolynomialInterpolation;

}  // namespace math


// DEFINITION
#include "curve_fit/PolynomialInterpolationGeneric.cpp"

#endif  // _MATH_POLYNOMIALINTERPOLATIONGENERIC_HPP_
