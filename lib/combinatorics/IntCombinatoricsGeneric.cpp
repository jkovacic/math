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
 * Implementation of functionality in the namespace IntCombinatorics with
 * functions for calculation of factorials, binomial coefficients, etc.
 */

// no #include "IntCombinatorics.hpp" !!!
#include "int_util/IntUtilGeneric.hpp"
#include "exception/CombinatoricsException.hpp"

#include <limits>
#include <algorithm>


// A namespace with "private" functions
namespace math {  namespace IntCombinatorics {  namespace __private
{

/*
 * Checks sign of 'n' and throws an exception if it is negative. 
 * 
 * @param n - integer value to check
 * 
 * @throw CombinatoricsException if 'n' is negative
 */
template <typename I>
void __checkSign(const I& n) throw (math::CombinatoricsException)
{
    if ( true == math::IntUtil::isNegative<I>(n) )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }
}

}}}  // namespace math::IntCombinatorics::__private



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
template <typename I>
I math::IntCombinatorics::fallingFactorial(
                            const I& N, 
                            const I& K ) 
                        throw(math::CombinatoricsException)
{
    // sanity check
    math::IntCombinatorics::__private::__checkSign<I>(N);
    math::IntCombinatorics::__private::__checkSign<I>(K);

    /*
     *  k
     *  -
     * n   =  n * (n-1) * (n-2) * ... * (n-k+1)
     */

    if ( K > N )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }

    if ( static_cast<I>(0) == K )
    {
        // TODO is this a definition or an exception should be thrown?
        return static_cast<I>(1);
    }

    I fallFact = static_cast<I>(1);

    /*
     * If reaching this point, K is always greater than 1 (K==0 is handled above)
     * so N-K+1 can never exceed N (and thus max. I).
     * 
     * Since K cannot be greater than N, the loop will terminate before
     * i reaches 0 and an integer overflow is not possible.
     */
    const I IMAX = std::numeric_limits<I>::max();

    for ( I i=N; i>=N-K+static_cast<I>(1); --i )
    {
        if ( IMAX/fallFact < i )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }

        fallFact *= i;
    }

    return fallFact;
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
template <typename I>
I math::IntCombinatorics::risingFactorial( 
                            const I& N, 
                            const I& K ) 
                        throw(math::CombinatoricsException)
{
    // sanity check
    math::IntCombinatorics::__private::__checkSign<I>(N);
    math::IntCombinatorics::__private::__checkSign<I>(K);

    if ( static_cast<I>(0) == N )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }

    if ( static_cast<I>(0) == K )
    {
        // TODO is this a definition or an exception should be thrown?
        return static_cast<I>(1);
    }

    /*
     * Make sure, N-1+k does not cause an integer overflow.
     * 
     * N-1+K must be less or equal to I_MAX:
     *   N-1+K <= I_MAX  ==>  K <= I_MAX-N+1 
     */
    if ( K > (std::numeric_limits<I>::max() - N + static_cast<I>(1)) )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }

    // The algorithm is already implemented by fallingFactorial
    return math::IntCombinatorics::fallingFactorial<I>(N-static_cast<I>(1)+K, K);
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
template <typename I>
I math::IntCombinatorics::factorial(
                            const I& N, 
                            const I& from ) 
                        throw(math::CombinatoricsException)
{
    // sanity check
    math::IntCombinatorics::__private::__checkSign<I>(N);
    math::IntCombinatorics::__private::__checkSign<I>(from);

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
    if ( static_cast<I>(0) == N )
    {
        return static_cast<I>(1);
    }

    if ( static_cast<I>(0)==from || from>N )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
    }

    /*
     * The algorithm is actually already implemented by fallingFactorial.
     * 
     * 'from' is always greater than 0 and less than or equal to N.
     * This means, N-from+1 cannot be out of integer range.
     */
    return math::IntCombinatorics::fallingFactorial<I>(N, N-from+static_cast<I>(1));
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
template <typename I>
I math::IntCombinatorics::multiFactorial(
                            const I& N, 
                            const I& K ) 
                        throw(math::CombinatoricsException)
{    
    // sanity check
    math::IntCombinatorics::__private::__checkSign<I>(N);
    math::IntCombinatorics::__private::__checkSign<I>(K);

    if ( N < K )
    {
        // TODO check definition of multifactorial
        return static_cast<I>(1);
    }

    I multiFact = static_cast<I>(1);

    /*
     * Iterative implementation of the definition above.
     * 
     * 'i' is always greater than or equal to 0, so an integer overflow
     * cannot occur.
     */
    const I IMAX = std::numeric_limits<I>::max();

    for ( I i=N; i>=K; i-=K )
    {
        if ( IMAX/multiFact < i )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }

        multiFact *= i;
    }

    return multiFact;
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
template <typename I>
I math::IntCombinatorics::doubleFactorial(
                            const I& N) 
                        throw(math::CombinatoricsException)
{
    // sanity check is performed by multiFactorial())

    // double factorial is a multi factorial with K=2
    return math::IntCombinatorics::multiFactorial<I>(N, static_cast<I>(2));
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
template <typename I>
I math::IntCombinatorics::binom(
                            const I& N, 
                            const I& K ) 
                        throw(math::CombinatoricsException)
{
    // sanity check
    math::IntCombinatorics::__private::__checkSign<I>(N);
    math::IntCombinatorics::__private::__checkSign<I>(K);

    if ( static_cast<I>(0) == N )
    {
        return static_cast<I>(0);
    }

    if ( static_cast<I>(0)==K || N==K )
    {
        return static_cast<I>(1);
    }

    if ( K > N )
    {
        return static_cast<I>(0);
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

    I bn = static_cast<I>(1);

    /*
     * Since binomial(n,k) == binomial(n,(n-k)), it is sensible to
     * to choose min(k, n-k) and thus reduce the number of multiplications.
     */
    const I k = std::min(K, N-K); 

    // Since k is min (K, N-K), it can never be equal to I_MAX,
    // so the integer overflow of 'i' is not possible.
    const I IMAX = std::numeric_limits<I>::max();

    for ( I i=static_cast<I>(1); i<=k; ++i )
    {
        const I factor = N - k + i; 
        if ( IMAX/bn < factor )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }

        bn *= factor;
        bn /= i;
    }

    return bn;
}
