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
 * @file PolynomialFittingGenericAb.h
 *
 * Declaration of the class PolynomialFittingGenericAb.h. This is an 
 * abstract class with some functionality common to all derived classes
 * that perform polynomial interpolation, regression, etc.
 *
 * @author Jernej Kovacic
 */

#ifndef _MATH_POLYNOMIALFITTINGGENERICAB_H_
#define _MATH_POLYNOMIALFITTINGGENERICAB_H_

#include "CurveFittingGenericAb.h"
#include "PolynomialGeneric.h"

namespace math
{

/**
 * @brief Common functionality for curve fitting/interpolation algorithms that
 *        find the best fitting polynomial.
 * 
 * The class is abstract and cannot be instantiated. More specialized algorithms
 * are implemented by its derived classes.
 */
template<class T>
class PolynomialFittingGenericAb : public CurveFittingGenericAb<T>
{

protected:
    // Curve fitting polynomial
    PolynomialGeneric<T> poly;

public:
    // get a copy of the best fitting polynomial
    PolynomialGeneric<T> getPolynomial() const throw (CurveFittingException);
    // inherited as an abstract function from the base class, value of the polynomial at the given abscissa
    T valueAt(const T& x, bool strict=true) const throw (CurveFittingException);

    //generateCurve() remains a pure virtual method and must be implemented by derived classes
};

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "PolynomialFittingGenericAb.cpp"

#endif // _MATH_POLYNOMIALFITTINGGENERICAB_H_
