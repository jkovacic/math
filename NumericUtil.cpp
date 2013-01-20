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

#include <complex>


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
 * and long double is even more accuarte
 */
 template<>
long double math::NumericUtil<long double>::EPS = 1e-22L;

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


/**
 * Does the given value equal (or is close enough to) zero?
 * Implementation depends on the type T.
 * For floating point types (float, double, long double), it checks
 * whether its absolute value is less than a hardcoded constant 'eps'.
 *
 * @param value
 *
 @ @return true or false
 */
template<class T>
bool math::NumericUtil<T>::isZero(const T& value)
{
    /*
     * The implementation for integers et al. where the == operator
     * does make sense and no comparison to EPS is necessary.
     */
    bool retVal = ( ZERO==value ? true : false );

    return retVal;
}


/*
 * Float, double and long double require specialized implementations of isZero().
 * In case of these three types, the equality operator (==) is useless.
 * In numerical mathematics, two numbers are considered "equal", when
 * absolute value of their difference does not exceed a reasonably set EPS.
 * All specializations are very similar and only differ in types of an input value.
 * For easier maintainability, the specialization will be implemented
 * only once using a parameterized #define
 */

#define _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(FDL) \
template<> \
bool math::NumericUtil<FDL>::isZero(const FDL& value) \
{ \
    bool retVal = false; \
    retVal = ( value>-EPS && value<EPS ? true : false ); \
    return retVal; \
}
// end of #define

// derive specialization for float:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(float)

// ... for double:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(double)

// ... and long double:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(long double)

// #definition of _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO not needed anymore, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO


/*
 * Specialization for complex.
 * As complex is a templated class, again it must be implemented for each supported subtemplated
 * type. To facilitate this, a parameterized macro is introduced.
 * Note: norm() calculates a sum of both parts' squares. It is a bit more efficient to compare
 * it with EPS^2 than calculating its square root (the actual definition of complex abs. value).
 */
#define _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX(FDL) \
template<> \
bool math::NumericUtil<std::complex<FDL> >::isZero(const std::complex<FDL>& value) \
{ \
    bool retVal = false; \
    const FDL eps = math::NumericUtil<FDL>::getEPS(); \
    retVal = ( std::norm(value)<=eps*eps ? true : false ); \
    return retVal; \
}
// end of #define

// derive specialization for float:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX(float)

// ... for double:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX(double)

// ... and long double:
_MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX(long double)

// #definition of _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX not needed anymore, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX
/*
 * Implementation for Rational
 */
template<>
bool math::NumericUtil<math::Rational>::isZero(const math::Rational& value)
{
    // Rational already contains its own isZero()...
    return value.isZero();
}

/**
 * @return internally hardcoded value of 'eps'
 */
template<class T>
T math::NumericUtil<T>::getEPS()
{
    return EPS;
}
