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
 * @file
 * @author Jernej Kovacic
 * 
 * Implementation of functionality in the namespace IntFactorization 
 * with functions for factorization of integers and some 
 * prime number utilities.
 */

// no #include "IntFactorizationGeneric.hpp" !!!
#include "int_util/IntUtilGeneric.hpp"
#include "exception/IntFactorizationException.hpp"

#include <limits>
#include <new>

// a namespace with "private" functions
namespace math {  namespace IntFactorization {  namespace __private
{

/*
 * Checks sign of 'n' and throws an exception if it is negative. 
 * 
 * @param n - integer value to check
 * 
 * @throw IntFactorizationException if 'n' is negative
 */
template <typename I>
void __checkSign(const I& n) throw (math::IntFactorizationException)
{
    // The actual implementation depends on signedness of 'I' and is
    // implemented in the class __checkSignImpl
    if ( true == math::IntUtil::isNegative<I>(n) )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::NEGATIVE_ARG);
    }
}

}}}  // namespace math::IntFactorization::__private


/**
 * @param n - integer to be checked
 * 
 * @return whether 'n' is a prime number
 */
template <typename I>
bool math::IntFactorization::isPrime(const I& n)
{
    /*
     * 0 is not a natural number and as such cannot be a prime.
     * By convention, 1 is also not a prime.
     */
    if ( n<=static_cast<I>(0) || static_cast<I>(1)==n )
    {
        return false;
    }
    
    /*
     * 2 and 3 are both primes.
     * As they are a sort of "specific" in comparison to other primes
     * (see below), they will be handled separately.
     */
    if ( static_cast<I>(2)==n || static_cast<I>(3)==n )
    {
        return true;
    }
    
    /*
     * Even numbers except 2 (already handled above) cannot be primes
     */
    if ( static_cast<I>(0) == n%static_cast<I>(2) )
    {
        return false;
    }
    
    /*
     * All prime numbers (p) greater than 3 satisfy at least one of the following two
     * criteria:
     * p = 6*k - 1   or  p = 6*k + 1  (where k is a natural number)
     * 
     * Proof:
     * In the following sequence: (p-1, p, p+1),
     * p cannot be even, so both of its "neighbours", p-1 and p+1, are even (divisible by 2).
     * Additionally only one member of the triple is divisible by 3. Again it cannot be p
     * (otherwise it wouldn't be a prime), so either p-1 or p+1 is divisible by 3 and 2
     * and thus also divisible by 6.
     * 
     * This criteria is applied as the first filter to determine whether 'n' is 
     * a candidate for a prime number.
     */
    const unsigned short int mod6 = static_cast<unsigned short int>(n % 6);
    if ( 1!=mod6 && 5!=mod6 )
    {
        return false;
    }
    
    /*
     * Finally apply the well known Eratosthenes algorithm.
     * Try dividing 'n' by all odd numbers (except 1) smaller than n. Odd numbers 
     * cannot be divisible by even divisors so it is sensible to skip them.
     * It also doesn't make sense to check divisibility by divisors 
     * that exceed sqrt(n).
     */

    // Integer sqrt(n) can be obtained quite efficiently.
    // Additionally it prevents any possibility of integer overflow
    // (theoretically possible for large 'n' if "i*i<=n" is used
    // as the for loop's condition)

    const I sqn = math::IntFactorization::intSqrt<I>(n);

    for ( I i=static_cast<I>(3); i<=sqn; i+=static_cast<I>(2) )
    {
        if ( static_cast<I>(0) == n%i )
        {
            // not a prime
            return false;
        }
    }
    
    return true;
}

/**
 * Greatest common divisor of 'first' and 'second'.
 *
 * @param first
 * @param second
 *
 * @return the greatest common divisor of 'first' and 'second'
 * 
 * @throw IntFactorizationException if any of input arguments is non-positive
 */
template <typename I>
I math::IntFactorization::greatestCommonDivisor(
                    const I& first, 
                    const I& second ) 
            throw(math::IntFactorizationException)
{
    // sanity check
    math::IntFactorization::__private::__checkSign<I>(first);
    math::IntFactorization::__private::__checkSign<I>(second);    

    /*
     * The well known Euclidean algorithm is utilized to find the greatest common divisor.
     * It is known to be efficient, more details about it are available, for instance,
     * at http://en.wikipedia.org/wiki/Euclidean_algorithm
     */

    // If any of both arguments is 0, the algorithm would "end up" in an infinite loop
    // or division by zero can occur. In such a case, throw an exception.
    if ( static_cast<I>(0) == first || static_cast<I>(0) == second )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::INVALID_INPUT);
    }

    // Now it is guaranteed to converge towards the GCD (or 1)
    I a = first;
    I b = second;
    I t;

    while ( static_cast<I>(0) != b )
    {
        t = b;
        b = a % b;
        a = t;
    }  // while

    return a;    
}

/**
 * Least common multiple of 'first' and 'second'.
 * It only makes sense when both arguments are not zero.
 *
 * @param first
 * @param second
 *
 * @return the least common multiple of 'first' and 'second'
 * 
 * @throw IntFactorizationException if any of input arguments is non-positive or if the LCM exceeds the I's range
 */
template <typename I>
I math::IntFactorization::leastCommonMultiple(
                    const I& first, 
                    const I& second) 
            throw(math::IntFactorizationException)
{
    /*
     * 'first' and 'second' can be expressed as:
     * 
     *     first = GCD * k1    and
     *    second = GCD * k2,
     * 
     * where GCD is their greatest common divisor and k1 and k2 are co-primes.
     * 
     * The largest common multiple (LCM) is defined as:
     *   
     *    LCM = GCD * k1 * k2
     * 
     * Since GCD cannot be less than 1, both sides of the equation can be
     * multiplied by GCD, resulting in:
     *
     *   GCD * LCM = (GCD * k1) * (GCD * k2) = first * second
     * 
     * Since an efficient algorithm for finding the GCD is known 
     * (the Euclidean algorithm), the LCM can be efficiently calculated as:
     * 
     *   LCM = first * second / GCD 
     */

    // sanity check will be performed by greatestCommonDivisor

    // At this point the GCD will be no less than 1, definitely not 0
    const I GCD = math::IntFactorization::greatestCommonDivisor<I>(first, second);

    // ...so it's safe to divide

    // Note that both 'first' and 'second' are divisible by GCD
    I temp = first / GCD;

    // make sure that the LCM will not fall out of I's range
    if ( std::numeric_limits<I>::max()/temp < second )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_RANGE);
    }

    return temp * second;
}

/**
 * Find the first prime number, greater than 'n'.
 * 
 * @param n
 * 
 * @return the first prime number p greater than n 
 * 
 * @throw IntFactorizationException if the next prime is out of I's range or 'n' is negative
 */
template <typename I>
I math::IntFactorization::nextPrime(const I& n) 
            throw(math::IntFactorizationException)
{
    // sanity check
    math::IntFactorization::__private::__checkSign<I>(n);

    I retVal = n;

    /*
     * integers <3 are a sort of "specific" in comparison to primes equal or
     * greater than 5 (see below), so they are handled separately.
     */
    if ( n < static_cast<I>(2) )
    {
        return static_cast<I>(2);
    }
    
    if ( static_cast<I>(2) == n )
    {
        return static_cast<I>(3);
    }
    
    /*
     * All primes equal or greater than 5 satisfy at least one of the conditions:
     *   p = 6*k-1  and  p = 6*k+1    where k is a positive integer number
     * (see comments of isPrime() for the proof). for that reason,
     * only integers satisfying these two criteria will be considered as
     * candidates for the next prime.
     */
    do
    {
        const unsigned short int mod6 = static_cast<unsigned short int>(retVal%6);
        unsigned short int summand = 0;

        switch (mod6)
        {
            case 0:
                // 6*k ==> 6*k+1
                summand = 1;
                break;

            case 1:
            case 2:
            case 3:
            case 4:
                // 6*k+n (1<=n<=4) ==> 6*k+5 = 6*(k+1)-1
                summand = 5 - mod6;
                break;

            case 5:
                // 6*k+5 = 6*(k+1)-1 ==> 6*k+7 = 6*(k+1)+1
                summand = 2;
                break;

            default:
                // should never occur, just in case:
                summand = 1;
        };

        // check of range:
        if ( std::numeric_limits<I>::max()-retVal < summand )
        {
            // the next candidate would be out of range
            throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_RANGE);
        }

        retVal += static_cast<I>(summand);
    }
    while ( false==math::IntFactorization::isPrime<I>(retVal) );
        
    return retVal;
}

/**
 * Finds the highest integer that does not exceed sqrt(n).
 * If 'n' is a perfect int square, its exact square root will be returned.
 * 
 * @param n - integer input argument whose "square root" will be calculated
 * 
 * @return floor(sqrt(n))
 * 
 * @throw IntFactorizationException if 'n' is negative
 */
template <typename I>
I math::IntFactorization::intSqrt(const I& n)
             throw (math::IntFactorizationException)
{
    // sanity check
    math::IntFactorization::__private::__checkSign<I>(n);

    /*
     * Integer square root can be efficiently calculated using the
     * C. Woo's algorithm, described and implemented at:
     * http://medialab.freaknet.org/martin/src/sqrt/
     */
    I res = static_cast<I>(0);
    I bit = static_cast<I>(1) << (static_cast<I>(8) * sizeof(I) - static_cast<I>(2));
    I sq = n;

    // Even if 'I' represents an unsigned type, 'bit' will never be
    // negative because its MSB will not be set.

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
 * Prime factorization.
 * 
 * @param n - an integer value to factorize
 * @param fac - a reference to a map to fill with prime factors and their powers 
 * 
 * @throw IntFactorizationException if allocation of memory fails or 'n' is negative
 */
template <typename I>
void math::IntFactorization::factor(
                    const I& n, 
                    std::map<I, I>& fac ) 
            throw(math::IntFactorizationException)
{
    // sanity check
    math::IntFactorization::__private::__checkSign(n);

    try
    {
        I comp = n;
        I pf = static_cast<I>(2);

        fac.clear();

        // "Specialized" handling of 0 and 1
        if ( static_cast<I>(0)==n || static_cast<I>(1)==n )
        {
            fac.insert(std::pair<I, I>(n, static_cast<I>(1)) );
            return;
        }
    
        /*
         * For all other integers, find consecutive primes and try dividing the
         * "remainder" by them. If it is divisible, create a new key of the map.
         * divide the "remainder" as long as possible.
         * For each successful division by a prime, update the value for the 
         * prime key and update the "remainder". 
         */
        for ( pf = static_cast<I>(2); 
              pf <= comp; 
              pf = math::IntFactorization::nextPrime<I>(pf) )
        {
            // Check divisibility by the current prime:
            if ( static_cast<I>(0) != comp%pf )
            {
                continue;
            }

            // If it is divisible, create a new map key first: 
            fac.insert( std::pair<I, I>(pf, static_cast<I>(0)) );
            
            // Divide by pf as many times as possible:
            while ( pf<=comp && static_cast<I>(0)==comp%pf )
            {
                // For each successful division, update the key's (pf) value:
                ++fac[pf];
                // and divide comp by pf:
                comp /= pf;
            }

        }

    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_MEMORY);
    }

}

/**
 * List of all integer numbers that divide 'n', including 1 and 'n' itself.
 * Divisors will be sorted in ascending order.
 * 
 * @param n - an integer whose divisors will be determined
 * @param div - a reference to a set to fill with integer divisors.
 * 
 * @throw IntFactorizationException if allocation of memory fails or 'n' is negative
 */
template <typename I>
void math::IntFactorization::divisors(
                    const I& n,
                    std::set<I>& div ) 
            throw(math::IntFactorizationException)
{
    // sanity check
    math::IntFactorization::__private::__checkSign(n);

    try
    {
        div.clear();

        // "Specialized" handling of 0 and 1
        if ( static_cast<I>(0)==n || static_cast<I>(1)==n )
        {
            div.insert(n);
            return;
        }
    
        /*
         * 'i' iterates from 1 to sqrt(n). For each 'i', try dividing
         * 'n' by it. If 'n' is divisible by 'i', insert 'i' and 'n/i'
         * into the set.
         */

        const I sqn = math::IntFactorization::intSqrt<I>(n);

        for ( I i=static_cast<I>(1); i<=sqn; ++i )
        {
            // Check divisibility by 'i':
            if ( static_cast<I>(0) == n%i )
            {
                // If successful, insert 'i' into the set:
                div.insert(i);
                // If 'i' is a perfect sqrt of 'n', do not insert it once again:
                const I ni = n/i;
                if ( i != ni )
                {
                    // If it is not a perfect square root, also insert 'n/i':
                    div.insert(ni);
                }
            }
        } // for i
        
        // Note that the set automatically sorts elements in ascending order
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_MEMORY);
    }

}
