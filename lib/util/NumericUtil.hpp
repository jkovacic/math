/*
Copyright 2011, Jernej Kovacic

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
 * @headername{NumericUtil.h}
 *
 * Declaration of the class NumericUtil, a collection of some useful
 * numerical utilities.
 */

#ifndef _MATH_NUMERICUTIL_HPP_
#define _MATH_NUMERICUTIL_HPP_

#include "rational/Rational.hpp"

#include <cstddef>
#include <complex>

namespace math
{

/**
 * @brief A class with some common functionality, e.g. determination whether
 *        a numeric value is close enough to zero, efficient exponentiation, etc.
 * 
 * Note that all functions are static so instantiation of this class is not necessary.
 */
template<class T>
class NumericUtil
{
private:
    // Equivalent to eps in MATLAB (and its clones)
    // Its value depends on type T
    static T EPS;
public:

    // Does the value equal zero? Note that in numerical mathematics
    // (where mostly float or double values are in use), "equals" has a
    // different meaning than in discrete mathematics (int etc.)
    static bool isZero(const T& value);
    static bool isZero(const T& value, const T& eps);

    // Get or set the internal value of 'eps'
    static T getEPS();
    static void setEPS(const T& eps);
    
    // num's sign:
    static short int sign(const T& num);
    
    // Efficient calculation of base^n:
    static T power(const T& base, size_t n);
};

// Declaration of specialized methods inside the name space declaration
// is essential if implemented elsewhere:
template<> bool NumericUtil<float>::isZero(const float& value, const float& eps);
template<> bool NumericUtil<double>::isZero(const double& value, const double& eps);
template<> bool NumericUtil<long double>::isZero(const long double& value, const long double& eps);
template<> bool NumericUtil<Rational>::isZero(const Rational& value, const Rational& eps);
template<> bool NumericUtil<std::complex<float> >::isZero(const std::complex<float>& value, const std::complex<float>& eps);
template<> bool NumericUtil<std::complex<double> >::isZero(const std::complex<double>& value, const std::complex<double>& eps);
template<> bool NumericUtil<std::complex<long double> >::isZero(const std::complex<long double>& value, const std::complex<long double>& eps);

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
} // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.

#include "util/NumericUtil.cpp"

#endif	// _MATH_NUMERICUTIL_HPP_
