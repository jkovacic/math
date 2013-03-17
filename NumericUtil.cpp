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

// Deliberately there is no #include "NumericUtil.h"
#include "Rational.h"
#include "SqMatrixGeneric.h"
#include "PolynomialGeneric.h"

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
 * and long double is even more accurate
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
 * @return value of 'eps' for the desired type
 */
template<class T>
T math::NumericUtil<T>::getEPS()
{
    return EPS;
}

/**
 * Sets the new value of 'eps' for the desired type if the default one does
 * not meet application's requirements.
 * 
 * @note There is no check of input so make sure a sensible value
 *       (typically a very small positive number) is entered.
 * 
 * @param eps - new value of EPS
 */
template<class T>
void math::NumericUtil<T>::setEPS(const T& eps)
{
    EPS = eps;
}

/*
 * In C++, it is not possible to specialize a function of a templated class
 * if the specialization is also a generic type (e.g. T -> math::SqMatrixGeneric<T>).
 * However it is possible if a single function is templated.
 * For more details, see the discussion at:
 * http://www.cplusplus.com/forum/general/68298/
 * 
 * At the moment, getUnit() is only used by power(). Until it is changed,
 * getUnit() will be defined in this file as a templated standalone function.
 * Additionlly, implementations of the function are "hidden" in the namespace
 * math::getunit which is not supposed to be known to other members of math::
 */
namespace math
{
    namespace getunit
    {
        /*
         * @param t - an arbitrary instance of T, ignored by the generic implementation,
         *            at specializations it might be useful to determine unit's dimensions etc.
         * 
         * @return an instance of T acting as a multiplication unit 
         */
        template<class T>
        T getUnit(const T& t)
        {
            return T(math::NumericUtil<T>::ONE);
        }

        /*
         * Specialization for generic class math::SqMatrixGeneric<T>
         * 
         * @param t - an arbitrary instance of a square matrix to determine 
         *            dimensions of the returned unit matrix
         *
         * @return a unit nxn square matrix where n is a dimension of 't'
         */
        template<class T>
        math::SqMatrixGeneric<T> getUnit(const math::SqMatrixGeneric<T>& t)
        {
            const unsigned int N = t.nrRows();
            math::SqMatrixGeneric<T> retVal(N);
            retVal.setUnit();
            return retVal;
        }

        /*
         * Specialization for generic class math::PolynomialGeneric<T>
         * 
         * @param t - ignored
         * 
         * @return a unit polynomial p(x) = (T)1 
         */
        template<class T>
        math::PolynomialGeneric<T> getUnit(const math::PolynomialGeneric<T>& t)
        {
            math::PolynomialGeneric<T> retVal(1);
            retVal.set(0, math::NumericUtil<T>::ONE);
            return retVal;
        }
    } // namespace units
} // namespace math

/**
 * Efficient calculation of positive integer power. 
 * Complexity of the algorithm is O(log2 n).
 * 
 * @note T must have implemented operator*=
 * 
 * @param base - base of the exponentiation
 * @param n - exponent (a positive integer number)
 * 
 * @return base^n 
 */
template<class T>
T math::NumericUtil<T>::power(const T& base, unsigned int n)
{
    /*
     * "Exponentiation by squaring" algorithm will be applied.
     * 
     * Exponentiation can be expanded into:
     *        n    n%2       n/2
     *       a  = a   * (a^2)
     * 
     * where / and % are integer division and remainder operators.
     * 
     * The expression above can be further recursively derived to:
     *  n    n%2    (n/2)%2     (n/2)/2
     * a  = a   (a^2)       (a^4)       =
     * 
     *       n%2    (n/2)%2    (n/4)%2    ((n/2)/2)/2
     *    = a   (a^2)      (a^4)      (a^8)           = ....
     * 
     * It is simple to develop an iterative algorithm where factor
     * is iteratively increased by squaring itself and coefficients 
     * ai (i=0..n) are iteratively calculated by the remainder of
     * division by 2 
     */
    
    T retVal = math::getunit::getUnit(base);
    T factor = base;
    
    // Obtain coefficients ai from the exponent's binary form.
    // Note: "i>>=1" is a bitwise equivalent bitwise equivalent of "i/=2"
    for ( unsigned int i=n; i>0; i>>=1 )
    {
        // Check the coefficient ai (no need to multiply retVal by 1 if ai is 0)
        // Note: "i&1" is a bitwise equivalent of "i%2"
        if ( 0!=(i & static_cast<unsigned int>(1) ) )
        {
            retVal *= factor;
        }
        
        factor *= factor;
    }

    return retVal;
}
