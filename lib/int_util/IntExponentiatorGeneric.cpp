/*
Copyright 2014, Jernej Kovacic

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
 * Implementation of functionality in the namespace IntExponentiator
 * with functions for efficient calculation of powers of
 * positive integer exponents.
 */


// No #include "IntExponentiatorGeneric.hpp" !!!
#include "matrix/SqMatrixGeneric.hpp"
#include "polynomial/PolynomialGeneric.hpp"


/*
 * In C++, it is not possible to specialize a function of a templated class
 * if the specialization is also a generic type (e.g. T -> math::SqMatrixGeneric<T>).
 * However it is possible if a single function is templated.
 * For more details, see the discussion at:
 * http://www.cplusplus.com/forum/general/68298/
 *
 * At the moment, getUnit() is only used by power(). Until it is changed,
 * getUnit() will be defined in this file as a templated standalone function.
 * Additionally, implementations of the function are "hidden" in the namespace
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
            (void) t;
            return T(static_cast<T>(1));
        }

        /*
         * Specialization for generic class math::SqMatrixGeneric<T>
         *
         * @param t - an arbitrary instance of a square matrix to determine
         *            dimensions of the returned unit matrix
         *
         * @return a unit n x n square matrix where n is a dimension of 't'
         */
        template<class T>
        math::SqMatrixGeneric<T> getUnit(const math::SqMatrixGeneric<T>& t)
        {
            const size_t N = t.nrRows();
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
            (void) t;
            return math::PolynomialGeneric<T>(static_cast<T>(1));
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
T math::IntExponentiator::power(const T& base, size_t n)
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

    // Note: it is safe to use bitwise operators for arithmetic operations
    // on unsigned int values as they do not depend on endianess:
    // http://stackoverflow.com/questions/7184789/does-bit-shift-depends-on-endianness

    T retVal = math::getunit::getUnit(base);
    T factor = base;

    // Obtain coefficients ai from the exponent's binary form.
    // Note: "i>>=1" is a bitwise equivalent bitwise equivalent of "i/=2"
    for ( size_t i=n; i>0; i>>=1 )
    {
        // Check the coefficient ai (no need to multiply retVal by 1 if ai is 0)
        // Note: "i&1" is a bitwise equivalent of "i%2"
        if ( 0!=(i & static_cast<size_t>(1) ) )
        {
            retVal *= factor;
        }

        factor *= factor;
    }

    return retVal;
}
