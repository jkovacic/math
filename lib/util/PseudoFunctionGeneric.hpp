/*
Copyright 2015, Jernej Kovacic

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
 *
 * Declaration and implementation of various pseudo functions within namespace 
 * PseudoFunction that are monotonically increasing (or decreasing) function of
 * their official "counterparts", but more efficient to compute.
 * 
 * @note This header should only be used internally!
 */

#ifndef _MATH_PSEUDOFUNCTIONGENERIC_HPP_
#define _MATH_PSEUDOFUNCTIONGENERIC_HPP_

#include <complex>
#include <cmath>


namespace math
{

namespace PseudoFunction
{


/*
 * A "pseudo absolute value" of 'x'.
 * 
 * The function is used in situations when values must be ordered by
 * their absolute values and sometimes it is more convenient/more efficient to
 * obtain values that are a monotonically increasing function of
 * actual absolute values.
 * 
 * The general implementation returns the actual absolute value which
 * can be obtained efficiently for most scalar types.
 * 
 * @param x - value whose absolute value is returned
 * 
 * @return absolute value of 'x' 
 */
template <class T>
inline T pabs(const T& x)
{
    return std::abs(x);
}



/*
 * Partial "specialization" of 'pabs' for complex numbers.
 * 
 * This function returns a square of the actual absolute value and is
 * as such more efficient because no additional calculation of square root
 * (not a fast operation) is necessary.
 * 
 * @param x - a complex value
 * 
 * @return square of the absolute value of 'x', returned as a complex value with imag. part equal to 0
 */
template <class T>
inline std::complex<T> pabs(const std::complex<T>& x)
{
    return std::complex<T>( std::norm(x), static_cast<T>(0) );
}



/*
 * A convenience function that compares two (absolute) values.
 * 
 * @param a - first value
 * @param b - second value
 * 
 * @return true if a>b, false otherwise
 */
template <class T>
inline bool absgt(const T& a, const T& b)
{
    return (a > b);
}



/*
 * Partial "specialization" of 'absgt' for complex numbers.
 * The function compares real parts of 'a' and 'b'.
 * 
 * @param a - first complex value
 * @param b - second complex value
 * 
 * @return true if re(a)>re(b), false otherwise
 */
template <class T>
inline bool absgt(const std::complex<T>& a, const std::complex<T>& b)
{
    return ( std::real(a) > std::real(b) );
}



/*
 * Converts pseudo absolute value into the actual one.
 * 
 * @param x - pseudo absolute value
 *
 * @return actual absolute value
 */
template <class T>
inline T pabs2abs(const T& x)
{
    // for non-complex types, the 'x' already is the actual abs. value
    return x;
}



/*
 * Partial "specialization" of 'pabs2abs' for complex types.
 * 
 * @param x - pseudo absolute value
 *
 * @return actual absolute value
 */
template <class T>
inline std::complex<T> pabs2abs(const std::complex<T>& x)
{
    /*
     * For complex types, the real part of 'x' is actually a norm of the
     * absolute value. Hence its square root must be calculated.
     */

    return std::complex<T>(std::sqrt(std::real(x)) , static_cast<T>(0));
}


}  // namespace PseudoFunction

}  // namespace math


#endif  // _MATH_PSEUDOFUNCTIONGENERIC_HPP_
