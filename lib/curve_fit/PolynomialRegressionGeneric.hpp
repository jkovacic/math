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
 * @headername{PolynomialRegressionGeneric.h}
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

/**
 * @brief A class implementing a polynomial regression algorithm, satisfying
 *        the least sum of squares criteria.
 * 
 * The algorithm accepts an arbitrary number of input points and finds a
 * polynomial of the desired degree that best fits to the points and whose
 * sum of square deviations is minimized.
 */
template<class T>
class PolynomialRegressionGeneric : public PolynomialFittingGenericAb<T>
{

public:
    // Constructor
    PolynomialRegressionGeneric();
    
    // inherited as an abstract function from the base class
    void generateCurve(size_t degree = 1) throw (CurveFittingException);

};

// Regression classes with elements of types float and double make most sense
// so the following two types are predefined:
typedef PolynomialRegressionGeneric<float>   FPolynomialRegression;
typedef PolynomialRegressionGeneric<double>  PolynomialRegression;

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math


// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.

#include "curve_fit/PolynomialRegressionGeneric.cpp"

#endif  // _MATH_POLYNOMIALREGRESSIONGENERIC_HPP_
