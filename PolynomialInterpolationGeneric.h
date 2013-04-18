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
 * @file PolynomialInterpolationGeneric.h
 *
 * Declaration of the class PolynomialInterpolationGeneric.h that calculates 
 * an interpolation polynomial that goes exactly through entered points.
 *
 * @author Jernej Kovacic
 */

#ifndef _MATH_POLYNOMIALINTERPOLATIONGENERIC_H_
#define _MATH_POLYNOMIALINTERPOLATIONGENERIC_H_

#include "PolynomialFittingGenericAb.h"


namespace math
{

/**
 * @brief A class that finds a polynomial that exactly fits an input series of points.
 * 
 * The algorithm takes N+1 points and calculates a N-degree polynomial exactly
 * fitting all input points.
 */
template<class T>
class PolynomialInterpolationGeneric : public PolynomialFittingGenericAb<T>
{
public:
    // Constructor
    PolynomialInterpolationGeneric();

    // inherited as an abstract function from the base class
    void generateCurve(size_t degree = 0) throw (CurveFittingException);
};

// Interpolation classes with elements of types float and double make most sense
// so the following two types are predefined:
typedef PolynomialInterpolationGeneric<float>   FPolynomialInterpolation;
typedef PolynomialInterpolationGeneric<double>  PolynomialInterpolation;

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math


// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.

#include "PolynomialInterpolationGeneric.cpp"

#endif // _MATH_POLYNOMIALINTERPOLATIONGENERIC_H_
