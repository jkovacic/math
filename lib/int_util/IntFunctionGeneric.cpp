/*
Copyright 2016, Jernej Kovacic

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
 * Implementation of functionality in the namespace IntFunction with some
 * elementary functions, "specialized" for integer numbers.
 */


// no #include "IntFunctionGeneric.hpp" !!!
#include "int_util/IntUtilGeneric.hpp"
#include "exception/IntFunctionException.hpp"

#include <limits>



/**
 * Finds the highest integer that does not exceed sqrt(n).
 * If 'n' is a perfect int square, its exact square root will be returned.
 *
 * @param n - integer input argument whose "square root" will be calculated
 *
 * @return floor(sqrt(n))
 *
 * @throw IntFunctionException if 'n' is negative
 */
template <typename I>
I math::IntFunction::intSqrt(const I& n)
             throw (math::IntFunctionException)
{
    // sanity check
    if ( true == math::IntUtil::isNegative<I>(n) )
    {
        throw math::IntFunctionException(math::IntFunctionException::NEGATIVE_ARG);
    }

    /*
     * Integer square root can be efficiently calculated using the
     * C. Woo's algorithm, described and implemented at:
     * http://medialab.freaknet.org/martin/src/sqrt/
     */
    I res = static_cast<I>(0);
    I bit = static_cast<I>(1) << (static_cast<I>(8) * sizeof(I) - static_cast<I>(2));
    I sq = n;

    /*
     * Even if 'I' represents an unsigned type, 'bit' will never be
     * negative because its MSB will not be set.
     */

    while ( bit > sq )
    {
        bit >>= static_cast<I>(2);
    }

    while ( bit != static_cast<I>(0) )
    {
        if (sq >= res + bit)
        {
            sq -= (res + bit);
            res = (res>>static_cast<I>(1)) + bit;
        }
        else
        {
            res >>= static_cast<I>(1);
        }

        bit >>= static_cast<I>(2);
    }

    return res;
}


/**
 * Finds the highest integer that does not exceed log2(n).
 * If 'n' is a power of 2, its log2 will be returned.
 *
 * @param n - integer input argument whose "log base 2" will be calculated
 *
 * @return floor(log2(n))
 *
 * @throw IntFunctionException if 'n' is not strictly positive or if the type 'I' is not integral
 */
template <typename I>
I math::IntFunction::intLog2(const I& n) throw(math::IntFunctionException)
{
    // Sanity check:
    // - this function uses bitwise operators that might not be properly supported
    //   with non integral types.
    if ( false == std::numeric_limits<I>::is_integer )
    {
        throw math::IntFunctionException(math::IntFunctionException::INVALID_TYPE);
    }

    // - log of any base is only defined on strictly positive numbers
    if ( n <= static_cast<I>(0) )
    {
        throw math::IntFunctionException(math::IntFunctionException::NONPOSITIVE_ARG);
    }


    // Size of type 'I' in bits
    I l2 = sizeof(I) * static_cast<I>(8);

    /*
     * Extremely unlikely situation to occur but better to handle it
     * to prevent unusual results.
     */
    if ( static_cast<I>(0) == l2 )
    {
        throw math::IntFunctionException(math::IntFunctionException::INTERNAL_ERROR);
    }

    // Log2 will ALWAYS be less than the bit size of the integral type
    --l2;

    /*
     * Short description of the algorithm:
     * Set the mask's MSB to 1 (1<<l2),
     * check if the mask's bit is set in 'n',
     * if not, compare the next significant bit of 'n'
     * (and decrement 'l2) and so forth until a n's bit is set.
     *
     * Note: since 'n' must be strictly positive (and consequently
     * non-zero), a set bit of 'n' will always be found and the
     * algorithm will never be trapped into an infinite loop.
     *
     */
    for ( I mask = static_cast<I>(1) << l2;
         (n & mask) == static_cast<I>(0);
          --l2, mask >>= static_cast<I>(1) );

    return l2;
}


/*
 * Specialization of 'intLog2' for the type 'bool'.
 * Implementation of this type is vendor specific and
 * theoretically it is possible that the generic algorithm
 * above ends up in an infinite loop. Hence this specialization
 * for 'bool' will immediately throw an exception.
 */
namespace math {  namespace IntFunction
{

    template <>
    bool intLog2(const bool& n) throw(math::IntFunctionException)
    {
        (void) n;
        throw math::IntFunctionException(math::IntFunctionException::INVALID_TYPE);
    }
}}
