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
 * @headername{PolynomialRegression.h}
 *
 * Declaration of the class PolynomialFittingGenericAb.h that calculates 
 * a regression polynomial using the method of least squares.
 */

#ifndef _MATH_POLYNOMIALREGRESSIONGENERIC_HPP
#define _MATH_POLYNOMIALREGRESSIONGENERIC_HPP_

#include <cstddef>

#include "PolynomialFittingGenericAb.hpp"
#include "exception/CurveFittingException.hpp"
#include "polynomial/PolynomialGeneric.hpp"


namespace math
{

// Some compilers may require forward declaration
// of the templated base class:
template <typename F> class PolynomialFittingGenericAb;


/**
 * @brief A class implementing a polynomial regression algorithm, satisfying
 *        the least sum of squares criteria.
 * 
 * The algorithm accepts an arbitrary number of input points and finds a
 * polynomial of the desired degree that best fits to the points and whose
 * sum of square deviations is minimized.
 */
template <typename F>
class PolynomialRegressionGeneric : public PolynomialFittingGenericAb<F>
{

public:
    // Constructor
    PolynomialRegressionGeneric();
    
    // inherited as an abstract function from the base class
    void generateCurve(const std::size_t degree = 1) throw (CurveFittingException);

};

// Regression classes with elements of types float, double and long double make
// most sense so the following three types are predefined:
typedef PolynomialRegressionGeneric<float>       FPolynomialRegression;
typedef PolynomialRegressionGeneric<double>      PolynomialRegression;
typedef PolynomialRegressionGeneric<long double> LDPolynomialRegression;

}  // namespace math


// DEFINITION
#include "curve_fit/PolynomialRegressionGeneric.cpp"

#endif  // _MATH_POLYNOMIALREGRESSIONGENERIC_HPP_
