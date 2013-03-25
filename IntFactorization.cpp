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
 * @file IntFactorization.cpp
 * 
 * Implementation of the class IntFactorization with static
 * functions for factorization of integers and some 
 * prime number utilities.
 * 
 * As all functions are static, no instantiation of this 
 * class is necessary.
 * 
 * @author Jernej Kovacic
 */

#include "IntFactorization.h"
#include "IntFactorizationException.h"

#include <climits>
#include <new>

    
/**
 * @param n - integer to be checked
 * 
 * @return whether 'n' is a prime number
 */
bool math::IntFactorization::isPrime(unsigned long int n)
{
    /*
     * 0 is not a natural number and as such cannot be a prime.
     * By convention, 1 is also not a prime.
     */
    if ( 0==n || 1==n )
    {
        return false;
    }
    
    /*
     * 2 and 3 are both primes.
     * As they are a sort of "specific" in comparison to other primes
     * (see below), they will be handled separately.
     */
    if ( 2==n || 3==n )
    {
        return true;
    }
    
    /*
     * Even numbers except 2 (already handled above) cannot be primes
     */
    if ( 0==n%2 )
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
    const unsigned short int mod6 = n % 6;
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
    for ( unsigned int i=3; i*i<=n; i+=2 )
    {
        if ( 0==n%i )
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
 * @throw IntFactorizationException if any of input arguments is zero
 */
unsigned long int math::IntFactorization::gcd(unsigned long int first, unsigned long int second) throw(math::IntFactorizationException)
{
    /*
     * The well known Euclidean algorithm is utilized to find the greatest common divisor.
     * It is known to be efficient, more details about it are available, for instance,
     * at http://en.wikipedia.org/wiki/Euclidean_algorithm
     */

    // If any of both arguments is 0, the algorithm would "end up" in an infinite loop
    // or division by zero can occur. In such a case, throw an exception.
    if ( 0 == first || 0 == second )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::INVALID_INPUT);
    }

    // Now it is guaranteed to converge towards the GCD (or 1)
    unsigned long int a = first;
    unsigned long int b = second;
    unsigned long int t;

    while ( 0 != b )
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
 * @throw IntFactorizationException if any of input arguments is zero or if the LCM exceeds the integer range
 */
unsigned long int math::IntFactorization::lcm(unsigned long int first, unsigned long int second) throw(math::IntFactorizationException)
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

    // If any of the two input values equals zero, neither GCD
    // nor LCM doesn't make any sense, throw an exception in this case.
    if ( 0 == first || 0 == second )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::INVALID_INPUT);
    }

    // At this point the GCD will be no less than 1, definitely not 0
    const unsigned long int GCD = gcd(first, second);

    // so it's safe to divide
    unsigned long long int retVal = first*second/GCD;
    // and check the range
    if ( retVal>ULONG_MAX )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_RANGE);
    }
        
    return static_cast<unsigned long int>(retVal);
}

/**
 * Find the first prime number, greater than 'n'.
 * 
 * @param n
 * 
 * @return the first prime number p greater than n 
 * 
 * @throw IntFactorizationException if the next prime is out of integer range
 */
unsigned long int math::IntFactorization::nextPrime(unsigned long int n) throw(math::IntFactorizationException)
{
    unsigned long int retVal = n;
    
    /*
     * integers <3 are a sort of "specific" in comparison to primes equal or
     * greater than 5 (see below), so they are handled separately.
     */
    if ( n < 2 )
    {
        return 2L;
    }
    
    if ( 2 == n )
    {
        return 3L;
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
        const unsigned short int mod6 = retVal%6;
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
        if ( UINT_MAX-retVal < summand )
        {
            // the next candidate would be out of range
            throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_RANGE);
        }
        
        retVal += summand;
    }
    while ( false==isPrime(retVal) );
        
    return retVal;
}

/**
 * Finds the highest integer that does not exceed sqrt(n).
 * If 'n' is a perfect int square, its exact square root will be returned.
 * 
 * @param n - integer input argument whose "square root" will be calculated
 * 
 * @return floor(sqrt(n)) 
 */
unsigned long int math::IntFactorization::intSqrt(unsigned long int n)
{
    /*
     * Integer square root can be efficiently calculated using the
     * C. Woo's algorithm, described and implemented at:
     * http://medialab.freaknet.org/martin/src/sqrt/
     */
    unsigned long int res = 0;
    unsigned long int bit = 1L << (8*sizeof(unsigned long int)-2);
    unsigned long sq = n;
    
    while ( bit > sq )
    {
        bit >>= 2;
    }
    
    while ( bit != 0 ) 
    {
        if (sq >= res + bit) 
        {
            sq -= res + bit;
            res = (res>>1) + bit;
        }
        else
        {
            res >>= 1;
        }
        
        bit >>= 2;
    }
    
    return res;
}

/**
 * Prime factorization.
 * 
 * @param n - an integer value to factorize
 * 
 * @return a map of all prime factors with their powers
 * 
 * @throw IntFactorizationException if allocation of memory fails
 */
std::map<unsigned long int, unsigned int> math::IntFactorization::factor(unsigned long int n) throw(math::IntFactorizationException)
{
    try
    {
        std::map<unsigned long int, unsigned int> retVal;
        retVal.clear();
        unsigned long int comp = n;
        unsigned long int pf = 2;

        // "Specialized" handling of 0 and 1
        if ( 0==n || 1==n )
        {
            retVal.insert(std::pair<unsigned long int, unsigned int>(n, 1) );
            return retVal;
        }
    
        /*
         * For all other integers, find consecutive primes and try dividing the
         * "remainder" by them. If it is divisible, create a new key of the map.
         * divide the "remainder" as long as possible.
         * For each successful division by a prime, update the value for the 
         * prime key and update the "remainder". 
         */
        for ( pf=2; pf<=comp; pf=nextPrime(pf) )
        {
            // Check divisibility by the current prime:
            if ( 0!=comp%pf )
            {
                continue;
            }

            // If it is divisible, create a new map key first: 
            retVal.insert( std::pair<unsigned long int, unsigned int>(pf, 0) );
            
            // Divide by pf as many times as possible:
            while ( pf<=comp && 0==comp%pf )
            {
                // For each successful division, update the key's (pf) value:
                retVal[pf]++;
                // and divide comp by pf:
                comp /= pf;
            }

        }
        
        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_MEMORY);
    }

}

/**
 * List of all integer numbers that divide 'n', including 1 and 'n' itself.
 *  
 * @param n - an integer whose divisors will be determined
 * 
 * @return a set of all integer divisors in ascending order
 * 
 * @throw IntFactorizationException if allocation of memory fails
 */
std::set<unsigned long int> math::IntFactorization::divisors(unsigned long int n) throw(math::IntFactorizationException)
{
    try
    {
        std::set<unsigned long int> retVal;
        retVal.clear();

        // "Specialized" handling of 0 and 1
        if ( 0==n || 1==n )
        {
            retVal.insert(n);
            return retVal;
        }
    
        /*
         * 'i' iterates from 1 to sqrt(n). for each 'i', try dividing
         * 'n' by it. If 'n' is divisible by 'i', insert 'i' and 'n/i'
         * into the set.
         */
        for ( unsigned long int i=1; i*i<=n; i++ )
        {
            // Check divisibility by 'i':
            if ( 0==n%i )
            {
                // If successful, insert 'i' into the set:
                retVal.insert(i);
                // If 'i' is a perfect sqrt of 'n', do not insert it once again:
                const unsigned long int ni = n/i;
                if ( i!=ni )
                {
                    // If it is not a perfect square root, also insert 'n/i':
                    retVal.insert(ni);
                }
            }
        } // for i
        
        // Note that the set automatically sorts elements in ascending order
        return retVal;
    }
    catch ( std::bad_alloc& ba )
    {
        throw math::IntFactorizationException(math::IntFactorizationException::OUT_OF_MEMORY);
    }

}
