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
 * @file IntCombinatorics.cpp
 * 
 * Implementation of the class IntCombinatorics with static
 * functions for calculation of factorials, binomial coefficients, etc.
 *
 * As all functions are static, no instantiation of this 
 * class is necessary.
 *  
 * @author Jernej Kovacic
 */

#include "IntCombinatorics.h"
#include "CombinatoricsException.h"

#include <climits>

#ifdef OPENMP
#    include <omp.h>
#endif


/**                                  -
 * Calculates the falling factorial power:
 * 
 *    K                           
 *    -                     
 *  N!  = N * (N-1) * ... * (N-K+1)
 * 
 * @param N - see above
 * @param K - see above
 * 
 * @return falling factorial
 * 
 * @throw Combinatorics exception if input arguments are invalid or the result would be out of integer range
 */
unsigned long long int math::IntCombinatorics::fallingFactorial(unsigned long long int N, unsigned long long int K) throw(math::CombinatoricsException)
{
    /*
     *  k
     *  -
     * n   =  n * (n-1) * (n-2) * ... * (n-k+1)
     */
    
    if ( K>N )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }
    
    if ( 0LL==K )
    {
        // TODO is this a definition or an exception should be thrown?
        return 1LL;
    }
    
    unsigned long long int retVal = 1LL;
    
    /*
     * If reaching this point, K is always greater than 1 (K==0 is handled above)
     * so N-K+1 can never exceed N (and thus ULLONG_MAX).
     * 
     * Since K cannot be greater than N, the loop will terminate before
     * i reaches 0 and an integer overflow is not possible.
     */
    for ( unsigned long long int i=N; i>=N-K+1; --i )
    {
        if ( ULLONG_MAX/retVal < i )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        retVal *= i;
    }
    
    return retVal;
}

/**
 * Calculates the rising factorial power:
 *    _                           
 *    K                     
 *  N!  = N * (N+1) * ... * (N+K-1)
 * 
 * @param N - see above
 * @param K - see above
 * 
 * @return rising factorial
 * 
 * @throw Combinatorics exception if input arguments are invalid or the result would be out of integer range
 */
unsigned long long int math::IntCombinatorics::risingFactorial(unsigned long long int N, unsigned long long int K) throw(math::CombinatoricsException)
{   
    if ( 0LL == N )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }
    
    if ( 0LL==K )
    {
        // TODO is this a definition or an exception should be thrown?
        return 1LL;
    }
    
    /*
     * Make sure, N-1+k does not cause an integer overflow.
     * 
     * N-1+K must be less or equal to ULLONG_MAX:
     *   N-1+K <= ULLONG_MAX  ==>  K <= ULLONG_MAX-N+1 
     */
    if ( K>(ULLONG_MAX-N+1) )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }
    
    // The algorithm is already implemented by fallingFactorial
    return fallingFactorial(N-1+K, K);
}

/**
 * Calculates the factorial of N:
 * 
 *    N! = 1 * 2 * 3 * ... * (N-1) * N
 * 
 * If 'from' is different than 1, only partial factorial is calculated:
 *   
 *    from * (from+1) * ... * (N-1) * N
 * 
 * where 0 < 'from' <= N
 * 
 * @param N - see above
 * @param from - start of factorization (default: 1)
 * 
 * @return factorial 
 * 
 * @throw Combinatorics exception if input arguments are invalid or the result would be out of integer range
 */
unsigned long long int math::IntCombinatorics::factorial(unsigned long long int N, unsigned long long int from) throw(math::CombinatoricsException)
{
    /*
     * Factorial of a positive integer 'n' is defined as:
     * 
     *    n! = 1 * 2 * 3 * ... * (n-2) * (n-1) * n
     * 
     * or
     *            n
     *          +---+
     *          |   |
     *     n! = |   |  i
     *          |   |
     *           i=1
     * 
     */

    // By convention, 0! = 1
    // 'from' is ignored in this case
    if ( 0LL == N )
    {
        return 1LL;
    }

    if ( 0LL==from || from>N )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }
    
    /*
     * The algorithm is actually already implemented by fallingFactorial.
     * 
     * 'from' is always greater than 0 and less than or equal to N.
     * This means, N-from+1 cannot be out of integer range.
     */
    return fallingFactorial(N, N-from+1);
}

/**
 * Calculates a general multifactorial:
 * 
 *           /
 *          |  1         <==   0 < N < K
 *   (K)    |       
 * N!    =  {
 *          |       (K)
 *          | (N-K)!     <==   N >= K
 *          \
 * 
 *      (K)
 *    N!    = N * (N-K) * (N-2*K) * ...
 * 
 * as long as (N-i*K) >= K
 *  
 * @param N - see above
 * @param K - see above
 * 
 * @return multi factorial
 * 
 * @throw Combinatorics exception if input arguments are invalid or the result would be out of integer range
 */
unsigned long long int math::IntCombinatorics::multiFactorial(unsigned long long int N, unsigned int K) throw(math::CombinatoricsException)
{    
    if ( N<K )
    {
        // TODO check definition of multifactorial
        return 1LL;
    }
    
    unsigned long long int retVal = 1LL;
    
    /*
     * Iterative implementation of the definition above.
     * 
     * 'i' is always greater than or equal to 0, so an integer overflow
     * cannot occur.
     */
    for ( unsigned long long int i=N; i>=K; i-=K )
    {
        if ( ULLONG_MAX/retVal < i )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }
        
        retVal *= i;
    }
    
    return retVal;
}
    
/**
 * calculates a double factorial:
 * 
 *    N!! = N * (N-2) * (N-4) * ... * (N-i*2)
 * 
 * as long as (N-i*2) >= 2
 * 
 * @param N - see above
 * 
 * @return double factorial
 * 
 * @throw Combinatorics exception if input arguments are invalid or the result would be out of integer range
 */
unsigned long long int math::IntCombinatorics::doubleFactorial(unsigned long long int N) throw(math::CombinatoricsException)
{
    // double factorial is a multi factorial with K=2
    return multiFactorial(N, 2);
}
    
/**
 * Calculates a binomial coefficient:
 * 
 *                    /     \
 *                    |  N  |         N!
 *      binom(N, K) = |     | = --------------
 *                    |  K  |     K! * (N-K)!
 *                    \     /
 * 
 * @param N - see above
 * @param K - see above
 * 
 * @return binomial coefficient
 * 
 * @throw Combinatorics exception if input arguments are invalid or the result would be out of integer range
 */ 
unsigned long long int math::IntCombinatorics::binom(unsigned long long int N, unsigned long long int K) throw(math::CombinatoricsException)
{    
    if ( 0LL == N )
    {
        return 0LL;
    }
    
    if ( 0LL==K || N==K )
    {
        return 1LL;
    }
    
    if ( K>N )
    {
        return 0LL;
    }
    
    /*
     * One possible algorithm to compute the binomial coefficient could be
     * the Pascal's triangle, however this method is memory inefficient.
     * 
     * Calculation of the binomial coefficient directly from its definition
     * (division of factorials) would be inefficient and would risk an integer 
     * overflow. If (n-k)! is eliminated from n!, the binomial coefficient
     * can be derived as a series of products:
     * 
     * 
     *       n!          (n-(k-1)) * (n-(k-2)) * ... * (n-2) * (n-1) * n
     *  ------------- = -------------------------------------------------- =
     *   k! * (n-k)!         1     *     2     * ... * (k-2) * (k-1) * k
     *
     * 
     *          n-(k-1)     n-(k-2)           n-2     n-1     n
     *      =  --------- * --------- * ... * ----- * ----- * ---  =
     *             1           2              k-2     k-1     k
     * 
     *           k
     *         +---+
     *         |   |  n-(k-i)
     *      =  |   | --------- 
     *         |   |     i
     *          i=1
     * 
     * 
     * To avoid errors because of truncation errors (due to integer division), it
     * is important that i runs from 1 to k and not vice versa. If 'r' consecutive
     * integers are multiplied, one of them is divisible by 'r', so the whole
     * product is also divisible by 'r' and no integer division truncation occurs. 
     */

    unsigned long long int retVal = 1LL;
    
    /*
     * Since binomial(n,k) == binomial(n,(n-k)), it is sensible to
     * to choose min(k, n-k) and thus reduce the number of multiplications.
     */
    const unsigned long long int k = ( K<=N/2 ? K : N-K );
    
    // Since k is min (K, N-K), it can never be equal to ULLONG_MAX,
    // so the integer overflow of 'i' is not possible.
    for ( unsigned long long int i=1; i<=k; ++i )
    {
        const unsigned long long int factor = N-k+i; 
        if ( ULLONG_MAX/retVal < factor )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }
        
        retVal *= factor;
        retVal /= i;
    }
    
    return retVal;
}
