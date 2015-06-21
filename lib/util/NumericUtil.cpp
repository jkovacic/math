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
#include <cmath>
#include <algorithm>

#include "../settings/numutil_settings.h"


namespace math { namespace NumericUtil { namespace __private
{

/*
 * Checks validity of the desired value of the application's default 'eps'.
 * 
 * @note The function is only applicable for floating point types.
 * 
 * @param val - desired default 'eps'
 * 
 * @return max( 'val', system specific machine epsilon )
 */
template <class T>
inline T __epsDirect(const T& val)
{
    return std::max<T>( val, std::numeric_limits<T>::epsilon() );
}


/*
 * Returns a valid value of the application's default 'eps' if a multiplier
 * by the system specific machine epsilon is given.
 * 
 * @note The function is only applicable for floating point types.
 * 
 * @param k - multiplier
 * 
 * @return max( k*epsilon(), epsilon() )
 */
template <class T>
inline T __epsMult(const T& k)
{
    return ( std::max<T>( k, static_cast<T>(1) ) * 
             std::numeric_limits<T>::epsilon() );
}


/*
 * A class (struct) that stores the default value of 'eps'
 * as a static member, hence the class is never instantiated.
 * To speedup the access, the member is accessed directly
 * as a public member, without getters and setters.
 */
template <class T>
struct Eps
{
public:
    static T eps;
};

// General initializer of the static member
template <class T>
T Eps<T>::eps = static_cast<T>(0);

// A macro for specialization of the initializer for floating point types:
#define _MATH_NUMERICUTIL_SPECIALIZED_INITEPS(FDL, DIR, VAL)                \
    template <>                                                             \
    FDL Eps<FDL>::eps = ( true == DIR ?                                     \
        math::NumericUtil::__private::__epsDirect<FDL>(VAL) :               \
        math::NumericUtil::__private::__epsMult<FDL>(VAL) );
// end of macro definition

// Apply the macro for each floating point type:
_MATH_NUMERICUTIL_SPECIALIZED_INITEPS(float, NUMUTIL_FLOAT_DIR, NUMUTIL_FLOAT_VAL);
_MATH_NUMERICUTIL_SPECIALIZED_INITEPS(double, NUMUTIL_DOUBLE_DIR, NUMUTIL_DOUBLE_VAL);
_MATH_NUMERICUTIL_SPECIALIZED_INITEPS(long double, NUMUTIL_LONGDOUBLE_DIR, NUMUTIL_LONGDOUBLE_VAL);

// the macro is no longer needed, #undef it
#undef _MATH_NUMERICUTIL_SPECIALIZED_INITEPS

}}}  // namespace math::NumericUtil::__private


/**
 * @return value of default 'eps' for the desired type
 */
template <class T>
T math::NumericUtil::getEPS()
{
    /*
     * for integral types, 'eps' does not
     * make any sense, in this case just return 0.
     */

    return static_cast<T>(0);
}


/*
 * Specialization of getEPS() for floating point types.
 * In this case just return the value, stored in Eps<T>::eps.
 *
 * All specializations are very similar and only differ in types of the
 * returned value. For easier maintainability, the specialization will be
 * implemented only once using a parameterized #define
 */

namespace math {  namespace NumericUtil
{

#define _MATH_NUMERICUTIL_SPECIALIZED_GETEPS(FDL)        \
template <>                                              \
FDL getEPS<FDL>()                                        \
{                                                        \
    return math::NumericUtil::__private::Eps<FDL>::eps;  \
}
// end of #define

_MATH_NUMERICUTIL_SPECIALIZED_GETEPS(float);
_MATH_NUMERICUTIL_SPECIALIZED_GETEPS(double);
_MATH_NUMERICUTIL_SPECIALIZED_GETEPS(long double);

// the macro is no longer needed, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_GETEPS

}}  // namespace math::NumericUtil


/**
 * Sets the application default value of 'eps', valid within
 * the entire unit of compilation.
 * 
 * @note For integral types, the function does not make any sense,
 *       hence it does not do anything in this case.
 * 
 * @note The actual default epsilon will be set to the system specific machine
 *       epsilon if 'eps' is less than this value.
 * 
 * @param eps - desired default value of 'epsilon' (default: 0)
 */
template <class T>
void math::NumericUtil::setEPS(const T& eps = static_cast<T>(0))
{
    // the "general" version (for integral types) of this
    // function does not do anything

    (void) eps;
}


/**
 * Sets the application default value of 'epsilon' as the specified multiple of
 * the system specific value of the machine epsilon.
 * 
 * The default 'epsilon' will be valid within the entire unit of compilation.
 * 
 * @note For integral types, the function does not make any sense,
 *       hence it does not do anything in this case.
 * 
 * @note The actual default epsilon will be set to the system specific machine
 *       epsilon if 'keps' is less 1.
 * 
 * @param keps - a factor to multiply the system specific machine epsilon
 */
template <class T>
void math::NumericUtil::setMultEPS(const T& keps)
{
    // the "general" version (for integral types) of this
    // function does not do anything

    (void) keps;
}
    

namespace math {  namespace NumericUtil
{

// Specializations of setEPS() for floating point types:

#define _MATH_NUMERICUTIL_SPECIALIZED_SETEPS(FDL)                  \
template <>                                                        \
void setEPS(const FDL& eps)                                        \
{                                                                  \
    math::NumericUtil::__private::Eps<FDL>::eps =                  \
       math::NumericUtil::__private::__epsDirect<FDL>(eps);        \
}
// end of macro definition

_MATH_NUMERICUTIL_SPECIALIZED_SETEPS(float);
_MATH_NUMERICUTIL_SPECIALIZED_SETEPS(double);
_MATH_NUMERICUTIL_SPECIALIZED_SETEPS(long double);

// the macro is no longer needed
#undef _MATH_NUMERICUTIL_SPECIALIZED_SETEPS

// Specializations of setMultEPS() for floating point types:

#define _MATH_NUMERICUTIL_SPECIALIZED_SETMULTEPS(FDL)         \
template <>                                                   \
void setMultEPS(const FDL& keps)                              \
{                                                             \
    math::NumericUtil::__private::Eps<FDL>::eps =             \
        math::NumericUtil::__private::__epsMult<FDL>(keps);   \
}
// end of macro definition

_MATH_NUMERICUTIL_SPECIALIZED_SETMULTEPS(float);
_MATH_NUMERICUTIL_SPECIALIZED_SETMULTEPS(double);
_MATH_NUMERICUTIL_SPECIALIZED_SETMULTEPS(long double);

// the macro is no longer needed
#undef _MATH_NUMERICUTIL_SPECIALIZED_SETMULTEPS

}}  // namespace math::NumericUtil


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
inline bool math::NumericUtil::isZero(const T& value, const T& eps)
{
    /*
     * The implementation for integral types where the == operator
     * does make sense and no comparison to EPS is necessary.
     */

    return ( static_cast<T>(0) == value );

    (void) eps;
}


/*
 * Floating point types require specialized implementations of isZero().
 * In case of these three types, the equality operator (==) is useless.
 * In numerical mathematics, two numbers are considered "equal", when
 * absolute value of their difference does not exceed a reasonably set EPS.
 * All specializations are very similar and only differ in types of an input value.
 * For easier maintainability, the specialization will be implemented
 * only once using a parameterized #define
 */

namespace math {  namespace NumericUtil
{

// specialization of isZero for floating point types:

#define _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_2ARG(FDL)       \
template <>                                                  \
inline bool isZero(const FDL& value, const FDL& eps)         \
{                                                            \
    return ( value>-eps && value<eps );                      \
}
// end of #define

_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_2ARG(float);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_2ARG(double);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_2ARG(long double);

// the macro is no longer needed, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_2ARG


/*
 * Specialization for complex.
 * As complex is a templated class, again it must be implemented for each supported
 * subtemplated type. To facilitate this, a parameterized macro is introduced.
 * Note: norm() calculates a sum of both parts' squares. It is compared to
 * the norm of 'eps'.
 */
#define _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_2ARG(FDL)                   \
template <>                                                                      \
inline bool isZero(const std::complex<FDL>& value, const std::complex<FDL>& eps) \
{                                                                                \
    return ( std::norm(value) <= std::norm(eps) );                               \
}
// end of #define

_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_2ARG(float);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_2ARG(double);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_2ARG(long double);

// the macro is no longer needed, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_2ARG

}}  // namepsace math::NumericUtil


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
inline bool math::NumericUtil::isZero(const T& value)
{
    return ( static_cast<T>(0) == value );
}

// Specialization for floating point types
namespace math {  namespace NumericUtil
{
    
#define _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_1ARG(FDL)     \
template <>                                                \
inline bool isZero(const FDL& value)                       \
{                                                          \
    return math::NumericUtil::isZero(value,                \
        math::NumericUtil::__private::Eps<FDL>::eps);      \
}
// end of macro definition

_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_1ARG(float);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_1ARG(double);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_1ARG(long double);

// the macro is no longer needed, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_1ARG

// and specialization for complex numbers
#define _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_1ARG(FDL)                  \
template <>                                                                     \
inline bool isZero(const std::complex<FDL>& value)                              \
{                                                                               \
    return math::NumericUtil::isZero(value,                                     \
           std::complex<FDL>(math::NumericUtil::__private::Eps<FDL>::eps) );    \
}
// end of macro definition

_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_1ARG(float);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_1ARG(double);
_MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_1ARG(long double);

// the macro is no longer needed, #undef it:
#undef _MATH_NUMERICUTIL_SPECIALIZED_ISZERO_COMPLEX_1ARG

}}  // namespace math::NumericUtil


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
 * Rounds a real number (of type F) to the nearest integer (of type I).
 * 
 * @note An attempt to convert a negative 'n' to an unsigned type will
 *       return 0.
 * 
 * @note If 'n' is out of I's range, an overflow will occur.
 * 
 * @param n - a real value to be rounded to the nearest integer
 * 
 * @return 'n' rounded to the nearest integer
 */
template <typename F, typename I>
inline I math::NumericUtil::intRound(const F& n)
{
    // if attempting to cast a negative 'n' to an unsigned type,
    // return 0 immediately:
    if ( n < static_cast<F>(0) && 
         false == std::numeric_limits<I>::is_signed )
    {
        return static_cast<I>(0);
    }

    /*
     * 0.5 is added to a positive 'n' or subtracted from a negative 'n'.
     * Casting to I will then cut off the fractional part.
     */
    return ( n >= static_cast<F>(0) ? 
             static_cast<I>(n + static_cast<F>(1)/static_cast<F>(2) ) : 
             static_cast<I>(n - static_cast<F>(1)/static_cast<F>(2) ) );
}


/**
 * Rounds the real value (of type F) 'n' downwards, returning the
 * largest integer value (of type I) that is not greater than 'n'.
 * 
 * @note An attempt to convert a negative 'n' to an unsigned type will
 *       return 0.
 * 
 * @note If 'n' is out of I's range, an overflow will occur.
 * 
 * @param n - a real value to be rounded to an integer
 * 
 * @return floor(n) casted to I
 */
template <typename F, typename I>
inline I math::NumericUtil::intFloor(const F& n)
{
    return math::NumericUtil::intRound<F, I>(std::floor(n));
}


/**
 * Rounds the real value (of type F) 'n' upwards, returning the
 * largest integer value (of type I) that is not less than 'n'.
 * 
 * @note An attempt to convert a negative 'n' to an unsigned type will
 *       return 0.
 * 
 * @note If 'n' is out of I's range, an overflow will occur.
 * 
 * @param n - a real value to be rounded to an integer
 * 
 * @return ceil(n) casted to I
 */
template <typename F, typename I>
inline I math::NumericUtil::intCeil(const F& n)
{
    return math::NumericUtil::intRound<F, I>(std::ceil(n));
}
