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
 * Implementation  of functions of the namespace NumericUtil,
 * a collection of some useful numerical utilities.
 */

// no #include "NumericUtil.hpp" !!!
#include <cstddef>
#include <complex>
#include <limits>



/**
 * @return value of 'eps' for the desired type
 */
template <class T>
T math::NumericUtil::getEPS()
{
    // for int and most other types, eps does not
	// make any sense, in this case just return 0.

    return static_cast<T>(0);
}


/*
 * Specialization of getEPS() for float, double and long double.
 * In this case just return the value obtained from
 * std::numeric_limits<type>::epsilon().
 *
 * All specializations are very similar and only differ in types of the
 * returned value. For easier maintainability, the specialization will be
 * implemented only once using a parameterized #define
 */

namespace math
{
namespace NumericUtil
{

#define _MATH_NUMERICUTIL_SPECIALIZED_GETEPS(FDL) \
template <> \
FDL getEPS<FDL>() \
{ \
    return std::numeric_limits<FDL>::epsilon(); \
}
// end of #define

_MATH_NUMERICUTIL_SPECIALIZED_GETEPS(float)
_MATH_NUMERICUTIL_SPECIALIZED_GETEPS(double)
_MATH_NUMERICUTIL_SPECIALIZED_GETEPS(long double)

// #definition of _MATH_NUMERICUTIL_SPECIALIZED_GETEPS not needed anymore, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_GETEPS

}  // namespace NumericUtil
}  // namespace math


/**
 * Does the given value equal (or is close enough to) zero?
 * Implementation depends on the type T.
 * For floating point types (float, double, long double), it checks
 * whether its absolute value is less than the system epsilon
 * for the type T.
 *
 * @param value - value to be evaluated
 *
 * @return true or false
 */
template <class T>
bool math::NumericUtil::isZero(const T& value)
{
    return math::NumericUtil::isZero(value, math::NumericUtil::getEPS<T>() );
}


/**
 * Does the given value equal (or is close enough to) zero?
 * Implementation depends on the type T.
 * For floating point types (float, double, long double), it checks
 * whether its absolute value is less than a small value 'eps'.
 *
 * @param value - value to be evaluated
 * @param eps - a "threshold" to compare 'value' to, where applicable
 *
 * @return true or false
 */
template <class T>
bool math::NumericUtil::isZero(const T& value, const T& eps)
{
    /*
     * The implementation for integers et al. where the == operator
     * does make sense and no comparison to EPS is necessary.
     */

    return ( static_cast<T>(0)==value ? true : false );

    (void) eps;
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

namespace math
{

namespace NumericUtil
{

#define _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO(FDL) \
template <> \
bool isZero(const FDL& value, const FDL& eps) \
{ \
    return ( value>-eps && value<eps ? true : false ); \
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
 * As complex is a templated class, again it must be implemented for each supported
 * subtemplated type. To facilitate this, a parameterized macro is introduced.
 * Note: norm() calculates a sum of both parts' squares. It is compared to
 * the norm of 'eps'.
 */
#define _MATH_NUMERICUTIL_SPECIALIZED_IS_ZERO_COMPLEX(FDL) \
template <> \
bool isZero(const std::complex<FDL>& value, const std::complex<FDL>& eps) \
{ \
    return ( std::norm(value)<=std::norm(eps) ? true : false ); \
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

}  // namepsace NumericUtil
}  // namespace math


/**
 * Sign of the number
 * 
 * @param num - number
 * 
 * @return -1 if 'num' is negative, 0 if (close to) 0, 1 if positive 
 */
template <class T>
short int math::NumericUtil::sign(const T& num)
{
    // handle 0 first:
    if ( true==math::NumericUtil::isZero(num) )
    {
        return 0;
    }

    return ( num < static_cast<T>(0) ? -1 : 1 );
}


/**
 * Multiplication unit for the type T, typically 1 casted to T.
 * 
 * @param t - an arbitrary instance of T, ignored by the generic implementation,
 *            at specializations it might be useful to determine unit's dimensions etc.
 *
 * @return an instance of T acting as a multiplication unit
 */
template <class T>
T math::NumericUtil::unit(const T& t)
{
    (void) t;
    return T(static_cast<T>(1));
}
