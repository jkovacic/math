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
@file NumericUtil.cpp

Implementation  of the class NumericUtil, a collection of some useful
numerical utilities. This is a templated class and must not be compiled.
Instead it must be included after the class declaration in the .h file

@author Jernej Kovacic
*/

// Delibarately there is no #include "NumericUtil.h"
#include "Rational.h"


// Note that the optimal EPS depends on application requirements
/*
 * Definition of EPS for type float
 */
template<>
float math::NumericUtil<float>::EPS = 1e-9f;

/*
 * Double is a more accurate type...
 */
template<>
double math::NumericUtil<double>::EPS = 1e-16;

/*
 * In int and other types, EPS doesn't make sense, so set it to 0
 */
template<class T>
T math::NumericUtil<T>::EPS = static_cast<T>(0);


/**
 * A constant value with the T's representation of zero (0)
 */
template<class T>
const T math::NumericUtil<T>::ZERO ( static_cast<T>(0) );

/**
 * A constant value with the T's representation of one (1)
 */
template<class T>
const T math::NumericUtil<T>::ONE ( static_cast<T>(1) );


/*
 * The implementation for integers et al. where the == operator
 * does make sense and no comparison to EPS is necessary.
 */
template<class T>
bool math::NumericUtil<T>::isZero(const T& value)
{
    bool retVal = ( static_cast<T>(0)==value ? true : false );

    return retVal;
}


/*
 * Float and double require specialized implementations of isZero().
 * In case of these two types, the equality operator (==) is useless.
 * In numerical mathematics, two numbers are considered "equal", when
 * absolute value of their difference does not exceed a reasonably set EPS.
 * Both specializations are very similar and only differ in types of an input value.
 * For easier maintainability, the specialization will be implemented
 * only once using a parameterized #define
 */

#define _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(FD) \
template<> \
bool math::NumericUtil<FD>::isZero(const FD& value) \
{ \
    bool retVal = false; \
    retVal = ( value>-EPS && value<EPS ? true : false ); \
    return retVal; \
}
// end od #define

// derive specialization for float:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(float)

// ... and for double:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(double)

// #definition of _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO not needed anymore, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO

/*
 * Implementation for Rational
 */
template<>
bool math::NumericUtil<math::Rational>::isZero(const math::Rational& value)
{
    // Rational already contains its own isZero()...
    return value.isZero();
}
