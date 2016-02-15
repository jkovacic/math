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
 * Implementation of functions within the namespace BinomDist
 * that perform various binomial distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability mass function, etc.
 */


/*
 * More details about the binomial distribution:
 * https://en.wikipedia.org/wiki/Binomial_distribution
 */


// no #include "BinomialDistGeneric.h" !!!
#include <cmath>
#include <algorithm>
#include <limits>

# include "../settings/probdist_settings.h"

#include "util/NumericUtil.hpp"
#include "int_util/IntUtilGeneric.hpp"
#include "combinatorics/IntCombinatoricsGeneric.hpp"
#include "int_util/IntExponentiatorGeneric.hpp"
#include "specfun/SpecFunGeneric.hpp"

#include "exception/StatisticsException.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/CombinatoricsException.hpp"



namespace math {  namespace BinomDist {  namespace __private {


/*
 * Checks validity of given parameters.
 *
 * 'n' and 'k' must be non-negative and 'k' must not
 * be greater than 'n', 'p' must be between 0 and 1.
 *
 * @param k - number of successes to check
 * @param n - number of trials to check
 * @param p - trial's success probability to check
 *
 * @throw StatisticsException if any parameter is not valid
 */
template <typename F, typename I>
void __checkParams(const I& k, const I& n, const F& p) throw(math::StatisticsException)
{

    /*
     * n >= 0
     * 0 <= k <= n
     */
    if ( true == math::IntUtil::isNegative<I>(k) ||
         true == math::IntUtil::isNegative<I>(n) ||
         k>n )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    /*
     * 0 <= p <= 1
     */
    if ( p < static_cast<F>(0) ||
         p > static_cast<F>(1) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }
}

}}}  // namespace math::BinomDist::__private



/**
 * Mean or expectation of the given binomial distribution
 *
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 *
 * @return binomial distribution's mean (or expected value)
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::BinomDist::mean(
          const I& n,
          const F& p
        ) throw(math::StatisticsException)
{
    math::BinomDist::__private::__checkParams<F, I>(
        static_cast<F>(0), n, p);

    /*
     * Binomial distribution's mean value is defined as:
     *
     *   mu = n * p
     */

    return static_cast<F>(n) * p;
}


/**
 * Variance of the given binomial distribution
 *
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 *
 * @return binomial distribution's variance
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::BinomDist::var(
          const I& n,
          const F& p
        ) throw(math::StatisticsException)
{
    math::BinomDist::__private::__checkParams<F, I>(
        static_cast<F>(0), n, p);

    /*
     * Binomial distribution's variance is defined as:
     *
     *   var = n * p * (1-p)
     */

    return static_cast<F>(n) * p * (static_cast<F>(1)-p);
}


/**
 * Standard deviation of the given binomial distribution
 *
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 *
 * @return binomial distribution's standard deviation
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::BinomDist::stdev(
          const I& n,
          const F& p
        ) throw(math::StatisticsException)
{
    /*
     * Standard deviation is defined as a square root of the variance:
     *
     *   s = sqrt(var)
     */

    return std::sqrt( math::BinomDist::var<F, I>(n, p) );
}


/**
 * Can the given binomial distribution be approximated
 * with a normal distribution?
 *
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 * @param th - threshold for n*p and n*(1-p) (default: 10)
 *
 * @return logical value indicating whether the binomial distribution can be approximated as a normal distribution
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
bool math::BinomDist::normalApprox(
          const I& n,
          const F& p,
          const F& th
       ) throw(math::StatisticsException)
{
    math::BinomDist::__private::__checkParams<F, I>(
        static_cast<F>(0), n, p);

    return ( ( static_cast<F>(n)*p >= th ) &&
              (static_cast<F>(n)*(static_cast<F>(1)-p) >= th) );
}


/**
 * Binomial distribution's probability mass function (pmf)
 *
 * @param k - number of successes
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 *
 * @return probability mass function for the given 'k'
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::BinomDist::pmf(
          const I& k,
          const I& n,
          const F& p
        ) throw(math::StatisticsException)
{
    math::BinomDist::__private::__checkParams<F, I>(
        k, n, p);

    /*
     * Binomial distribution's probability mass function is defined as:
     *
     *         /   \
     *         | n |    k        n-k
     *   pmf = |   | * p  * (1-p)
     *         | k |
     *         \   /
     */

    const F ONE = static_cast<F>(1);
    F bnm;

    try
    {
        /*
         * This function to evaluate binom is usually
         * more efficient....
         */

        bnm = static_cast<F>( math::IntCombinatorics::binom<I>(n, k) );
    }
    catch ( const math::CombinatoricsException& cex )
    {
        /*
         * ... but if it failed due to an integer overflow,
         * it is evaluated via the beta function:
         * 
         *   /   \
         *   | n |              1
         *   |   | = -----------------------
         *   | k |    (n+1) * B(n-k+1, k+1)
         *   \   /
         * 
         */

        bnm = ONE /
                ( (static_cast<F>(n) + ONE) *
                   math::SpecFun::beta<F>(static_cast<F>(n-k)+ONE, static_cast<F>(n)+ONE) );
    }

    return bnm *
           math::IntExponentiator::power<F, I>(p, k) *
           math::IntExponentiator::power<F, I>(ONE-p, n-k);

}


/**
 * Probability of a<= X<= b, where any '<=' may be replaced
 * by '<', depending on 'incLower' and 'incUpper'.
 *
 * If strictly less operator is desired at any side of the interval,
 * 'lowerInc' and/or 'upperInc' must be set to 'true'
 *
 * @note 'a' and 'b' are not interchangeable.
 *
 * @note If the actual probability should be negative, it will be
 *       "rounded" to 0.
 *
 * @param a - lower number of successes
 * @param b - upper number of successes
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 * @param incLower - should the lower limit be included into the probability (default: TRUE)
 * @param incUpper - should the upper limit be included into the probability (default: TRUE)
 *
 * @return P(a<=X<=b) or P(a<X<=b) or P(a<=X<b) or P(a<X<b), depending on 'incLower' and 'incUpper'
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::BinomDist::probInt(
          const I& a,
          const I& b,
          const I& n,
          const F& p,
          const bool incLower,
          const bool incUpper
        ) throw(math::StatisticsException)
{
    /*
     * P(a<=X<=b) is defined as:
     *
     *                  b
     *                ----- /   \
     *                \     | n |    i        n-i
     *   P(a<=X<=b) =  >    |   | * p  * (1-p)
     *                /     | i |
     *                ----- \   /
     *                 i=a
     *
     * Alternatively the probability can be evaluated as
     * the difference of both cdf's
     *
     *   P(a<=X<=b) = cdf(b, n, p, incUpper) - cdf(a, n, p, incLower)
     */

    const F PP = math::BinomDist::prob<F, I>(b, n, p, incUpper) -
                 math::BinomDist::prob<F, I>(a, n, p, !incLower);

    // "round" negative probabilities to 0
    return std::max(PP, static_cast<F>(0) );
}


/**
 * Binomial distribution's cumulative distribution function (cdf),
 * denoting the probability, that the number of successes is
 * less than (or greater than, depending on 'lowerTail' and 'incl') 'k'.
 *
 * Depending on both parameters, the following probability is returned:
 *
 * - lowerTail=true  and incl=true:    P(X<=k)
 * - lowerTail=true  and incl=false:   P(X<k)
 * - lowerTail=false and incl=true:    P(X>=k)
 * - lowerTail=false and incl=false:   P(X>k)
 *
 * @param k - number of successes
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 * @param incl - should k be included into the cdf (default: true)
 * @param lowerTail - if 'true', returns P(X<=k), otherwise P(X>k) (default: true)
 *
 * @return see above
 *
 * @throw StatisticsException if any argument is invalid or if evaluation of cdf failed
 */
template <typename F, typename I>
F math::BinomDist::prob(
          const I& k,
          const I& n,
          const F& p,
          const bool incl,
          const bool lowerTail
        ) throw(math::StatisticsException)
{
    math::BinomDist::__private::__checkParams<F, I>(
        k, n, p);

    /*
     * Binomial distribution's cumulative distribution function (cdf)
     * is defined as:
     *
     *               k
     *             ----- /   \
     *             \     | n |    i        n-i
     *   P(X<=k) =  >    |   | * p  * (1-p)
     *             /     | i |
     *             ----- \   /
     *              i=0
     *
     * However, it is more convenient to obtain the probability
     * via the regularized incomplete beta function as mentioned
     * at the page 338 of the Numerical Recipes, 3rd Edition:
     *
     * - P(X<=k) = 1 - Ip(k+1, n-k)
     * - P(X<k)  = 1 - Ip(k, n-k+1)
     * - P(X>=k) = Ip(k, n-k+1)
     * - P(X>k)  = Ip(k+1, n-k)
     */

    try
    {
        // Tolerance for the incomplete beta function
        const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                      static_cast<F>(PROBDIST_PROB_TOL_DEN);

        /*
         * A special case when n==0 and hence k==0:
         * If 'incl' equals true, 1 is returned, 0 otherwise:
         */
        if ( static_cast<I>(0) == n )
        {
            return ( true==incl ? static_cast<F>(1) : static_cast<F>(0) );
        }

        // Handling a situation when k==0
        if ( static_cast<I>(0) == k )
        {
            // P(X<=0) = pmf(0)
            if ( true==lowerTail && true==incl )
            {
                return math::BinomDist::pmf<F, I>(static_cast<I>(0), n, p);
            }

            // P(X<0) = 0
            if ( true==lowerTail && false==incl )
            {
                return static_cast<F>(0);
            }

            // P(X>=0) = 1
            if ( false==lowerTail && true==incl )
            {
                return static_cast<F>(1);
            }

            // P(X>0) = 1 - pmf(0)
            if ( false==lowerTail && false==incl )
            {
                return ( static_cast<F>(1) - math::BinomDist::pmf<F, I>(static_cast<I>(0), n, p) );
            }
        }

        // Handling a situation when k==n
        if ( k == n )
        {
            // P(X<=n) = 1
            if ( true==lowerTail && true==incl )
            {
                return static_cast<F>(1);
            }

            // P(X<n) = 1 - pmf(n)
            if ( true==lowerTail && false==incl )
            {
                return ( static_cast<F>(1) - math::BinomDist::pmf<F, I>(n, n, p) );
            }

            // P(X>=n) = pmf(n)
            if ( false==lowerTail && true==incl )
            {
                return math::BinomDist::pmf<F, I>(n, n, p);
            }

            // P(X>n) = 0
            if ( false==lowerTail && false==incl )
            {
                return static_cast<F>(0);
            }
        }


        /*
         * See the comments above to find out which incomplete beta function
         * must be obtained.
         */
        const F Ip = ( lowerTail != incl ?
                math::SpecFun::incBetaLowerReg<F>(
                     static_cast<F>(k),
                     static_cast<F>(n-k) + static_cast<F>(1),
        		     p, TOL ) :
                math::SpecFun::incBetaLowerReg<F>(
                     static_cast<F>(k) + static_cast<F>(1),
                     static_cast<F>(n-k),
                     p, TOL )   );

        // subtract the inc. beta function from 1 if necessary
        const F pr = ( true==lowerTail ?
                       static_cast<F>(1) - Ip : Ip );

        return pr;
    }
    catch ( const math::SpecFunException& sfex )
    {
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }
}


/**
 * Quantile function for the specified binomial distribution.
 *
 * If 'lowerTail' equals TRUE, and 'smallest' equals TRUE,
 * it returns the smallest 'k' that satisfies:
 *    sum(i=0, k, pmf(i, n, p)) >= prob
 *
 * If 'lowerTail' equals TRUE and 'smallest' equals FALSE,
 * it returns the largest 'k' that satisfies:
 *    sum(i=0, k, pmf(i, n, p)) <= prob
 *
 * If 'lowerTail' equals FALSE and 'smallest' equals TRUE,
 * it returns the smallest 'k' that satisfies:
 *    sum(i=k, n, pmf(i, n, p)) <= prob
 *
 * If 'lowerTail' equals FALSE and 'smallest' equals FALSE,
 * it returns the largest 'k' that satisfies:
 *    sum(i=k, n, pmf(i, n, p)) >= prob
 *
 * @param prob - probability (must be greater than 0 and less than 1)
 * @param n - number of trials
 * @param p - success probability (default: 0.5)
 * @param smallest - see above (default: true)
 * @param lowerTail - see above (default: true)
 *
 * @return quantile for the given 'prob', depending on 'smallest' and 'lowerTail'. See above for details.
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
I math::BinomDist::quant(
          const F& prob,
          const I& n,
          const F& p,
          const bool smallest,
          const bool lowerTail
        ) throw(math::StatisticsException)
{
    const I Imax = std::numeric_limits<I>::max();
    const I ONE = static_cast<I>(1);
    const I ZERO = static_cast<I>(0);

    // sanity check
    math::BinomDist::__private::__checkParams<F, I>(static_cast<F>(0), n, p);
    math::BinomDist::__private::__checkParams<F, I>(static_cast<F>(0), static_cast<F>(0), prob);

    // Handling both "corner cases" (p==0 or p==1)
    if ( prob < math::NumericUtil::getEPS<F>() )
    {
        return ( true==lowerTail ? ZERO : n );
    }

    if ( prob > (static_cast<F>(1)-math::NumericUtil::getEPS<F>()) )
    {
        return (true==lowerTail ? n : ZERO );
    }


    /*
     * Proper handling of cases when 'lowerTail' equals FALSE.
     * This situation can be processed exactly the same way
     * as with lower tail probabilities, with a few modifications.
     * Note that binomial distribution is in general not symmetric.
     * Hence the complement of 'prob' must be obtained. Additionally
     * 'smallest' must be inverted.
     */

    const F PB = ( true==lowerTail ? p : static_cast<F>(1)-p );
    const bool SV = ( true==lowerTail ? smallest : !smallest );

    /*
     * Algorithm proposed in Numerical Recipes, 3rd Edition, page 339.
     * Many integer checks have been added to improve robustness.
     */

    I k = static_cast<I>( math::BinomDist::mean<F, I>(n, PB) );
    I step = ONE;
    I kl, ku;
    F cdf = math::BinomDist::prob<F, I>(k, n, PB, false);

    if ( prob < cdf )
    {
        do
        {
            k = ( k<step ? ZERO : k-step );
            step = ( step<(Imax-step) ? step * static_cast<I>(2) : Imax );
            cdf = math::BinomDist::prob<F, I>(k, n, PB, false);
        }
        while ( prob < cdf );

        const I s2 = step / static_cast<I>(2);
        kl = k;
        ku = std::min( ( k<=(Imax-s2) ? k + s2 : Imax ), n );
    }
    else
    {
        do
        {
            k = std::min( ( k<=(Imax-step) ? k+step : Imax ), n );
            step = ( step<(Imax-step) ? step * static_cast<I>(2) : Imax );
            cdf = math::BinomDist::prob<F, I>(k, n, PB, false);
            if ( k == n )
            {
                // prevent possible infinite loop if 'cdf'
                // is still less than 'prob'
                break;  // out of do - while
            }
        }
        while ( prob > cdf );

        const I s2 = step / static_cast<I>(2);
        kl = ( k>=s2 ? k-s2 : ZERO );
        ku = k;
    }

    // Bisection between 'kl' and 'ku'
    while ( (ku-kl) > ONE )
    {
        /*
         * Equivalent to (kl + ku) / 2, but more robust
    	 * Note: when a positive number is casted to an integer,
         * it is always rounded down.
         */
        k = kl + (ku - kl) / static_cast<I>(2);

        if ( prob < math::BinomDist::prob<F, I>(k, n, PB, false) )
        {
            ku = k;
        }
        else
        {
            kl = k;
        }
    }

    // Finally adjust the result according to 'smallest'
    k = kl;
    cdf = math::BinomDist::prob<F, I>(k, n, PB, true);

    if ( true==SV && cdf>=prob )
    {
        for ( ; k>ZERO && math::BinomDist::prob<F, I>(k-ONE, n, PB, true)>prob; --k );
    }
    else if ( true==SV && cdf<prob )
    {
        for ( ; k<n && math::BinomDist::prob<F, I>(k, n, PB, true)<prob; ++k );
    }
    else if ( false==SV && cdf>=prob )
    {
        for( ; k>ZERO && math::BinomDist::prob<F, I>(k, n, PB, true)>prob  ; --k );
    }
    else  // SV==false && cdf<prob
    {
        for ( ; k<n && math::BinomDist::prob<F, I>(k+ONE, n, PB, true)<prob; ++k );
    }

    /*
     * When handling upper tail probabilities, the distribution was "flipped"
     * about the mode. For that reason, 'k' must be converted back to
     * the original distribution.
     */
    return ( true==lowerTail ? k : n-k );
}
