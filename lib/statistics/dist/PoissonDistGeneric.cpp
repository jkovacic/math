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
 * Implementation of functions within the namespace PoissonDist
 * that perform various Poisson distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability mass function, etc.
 */


/*
 * More details about the Poisson distribution:
 * https://en.wikipedia.org/wiki/Poisson_distribution
 */


// no #include "PoissonDistGeneric.h" !!!
#include <cmath>
#include <algorithm>
#include <limits>

# include "../settings/probdist_settings.h"

#include "util/NumericUtil.hpp"
#include "int_util/IntUtilGeneric.hpp"
#include "int_util/IntExponentiatorGeneric.hpp"
#include "specfun/SpecFunGeneric.hpp"

#include "exception/StatisticsException.hpp"
#include "exception/SpecFunException.hpp"


namespace math {  namespace PoissonDist {  namespace __private {

/*
 * Checks validity of given parameters.
 *
 * 'k' must be non-negative and 'lambda' must
 * be strictly greater than 0.
 * 
 * @param k - number of successes to check
 * @param lambda - mean number of successes to check
 * 
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
void __checkParams(const I& k, const F& lambda) throw(math::StatisticsException)
{
    /*
     * lambda > 0
     * k >= 0
     */

    if ( true == math::IntUtil::isNegative<I>(k) ||
         lambda < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }
}

}}}  // namespace math::PoissonDist::__private



/**
 * Poisson distribution's probability mass function (pmf)
 * 
 * @param k - number of successes
 * @param lambda - mean number of successes
 * 
 * @return probability mass function for the given 'k'
 * 
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::PoissonDist::pmf(
          const I& k,
          const F& lambda
        ) throw( math::StatisticsException)
{
    /*
     * Probability mass function is defined as:
     * 
     *                          k
     *                    lambda      -lambda
     *   pmf = P(X==k) = --------- * e
     *                      k!
     * 
     */

    // sanity check
    math::PoissonDist::__private::__checkParams<F, I>(k, lambda);

    /*
     * The provided implementation of factorial for integers can quickly blow up 
     * beyond the I's range, hence it is more convenient to apply the following
     * relation between factorial and gamma function:
     * 
     *   k! = gamma(k+1)
     */

    return math::IntExponentiator::power<F, I>(lambda, k) *
           std::exp(-lambda) /
           math::SpecFun::gamma<F>(static_cast<F>(k) + static_cast<F>(1));
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
 * @param lambda - mean number of successes
 * @param incLower - should the lower limit be included into the probability (default: TRUE)
 * @param incUpper - should the upper limit be included into the probability (default: TRUE)
 * 
 * @return P(a<=X<=b) or P(a<X<=b) or P(a<=X<b) or P(a<X<b), depending on 'incLower' and 'incUpper'
 * 
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::PoissonDist::probInt(
          const I& a,
          const I& b,
          const F& lambda,
          const bool incLower,
          const bool incUpper
        ) throw (math::StatisticsException)
{
    // sanity check
    math::PoissonDist::__private::__checkParams<F, I>(a, lambda);
    math::PoissonDist::__private::__checkParams<F, I>(b, lambda);

    /* The probability can be evaluated as the difference of both cdf's
     *
     *   P(a<=X<=b) = cdf(b, n, p, incUpper) - cdf(a, n, p, incLower)
     */

    const F PP = math::PoissonDist::prob<F, I>(b, lambda, incUpper) -
                 math::PoissonDist::prob<F, I>(a, lambda, !incLower);

    // "round" negative probabilities to 0
    return std::max(PP, static_cast<F>(0) );
}


/**
 * Poisson distribution's cumulative distribution function (cdf),
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
 * @param lambda - mean number of successes
 * @param incl - should 'k' be included into the cdf (default: TRUE)
 * @param lowerTail - if 'true', returns P(X<=k) or P(X<k), otherwise P(X>k) or P(X>=k) (default: TRUE)
 * 
 * @return see above
 * 
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
F math::PoissonDist::prob(
          const I& k,
          const F& lambda,
          const bool incl,
          const bool lowerTail
        ) throw(math::StatisticsException)
{
    // sanity check
    math::PoissonDist::__private::__checkParams<F, I>(k, lambda);

    /*
     * The cumulative distribution function (cdf) could be obtained
     * as the sum of corresponding probability mass functions.
     * However, often it is more convenient to obtain the probability
     * via the regularized upper incomplete gamma function as mentioned
     * at the pp. 336- 337 of the Numerical Recipes, 3rd Edition: 
     * 
     * - P(X<k) = Q(k, lambda)
     * - P(X<=k) = Q(k+1, lambda)
     * - P(X>k) = 1 - Q(k+1, lambda)
     * - P(X>=k) = 1 - Q(k, lambda)
     * 
     * where Q(a,b) denotes the regularized upper incomplete gamma function.
     */

    const F ONE = static_cast<F>(1);

    // Tolerance for the incomplete gamma function
    const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                  static_cast<F>(PROBDIST_PROB_TOL_DEN);

    /*
     * If k==0, the incomplete gamma function may not be defined.
     * In this case the probability is obtained from pmf.
     */
    if ( k == static_cast<I>(0) )
    {
        // P(X<=0) = pmf(0, lambda)
        if ( true==incl && true==lowerTail )
        {
            return math::PoissonDist::pmf<F>(k, lambda);
        }

        // P(X<0) = 0
        if ( false==incl && true==lowerTail )
        {
            return static_cast<F>(0);
        }

        // P(X>=0) = 1
        if ( true==incl && false==lowerTail )
        {
            return ONE;
        }

        // P(X>0) = 1 - pmf(0, lambda)
        if ( false==incl && false==lowerTail )
        {
            return ( ONE - math::PoissonDist::pmf<F>(k, lambda) );
        }
    }  // if k==0

    const F q = ( (false==incl && true==lowerTail) || (true==incl && false==lowerTail) ?
                   math::SpecFun::incGammaUpperReg<F>(static_cast<F>(k), lambda, TOL) :
                   math::SpecFun::incGammaUpperReg<F>(static_cast<F>(k)+ONE, lambda, TOL) );

    return ( true==lowerTail ? q : ONE-q );
}


/**
 * Quantile function for the specified Poisson distribution.
 *
 * If 'lowerTail' equals TRUE, and 'smallest' equals TRUE,
 * it returns the smallest 'k' that satisfies:
 *    sum(i=0, k, pmf(i, lambda)) >= prob
 *
 * If 'lowerTail' equals TRUE and 'smallest' equals FALSE,
 * it returns the largest 'k' that satisfies:
 *    sum(i=0, k, pmf(i, lambda)) <= prob
 *
 * If 'lowerTail' equals FALSE and 'smallest' equals TRUE,
 * it returns the smallest 'k' that satisfies:
 *    sum(i=k, inf, pmf(i, lambda)) <= prob
 *
 * If 'lowerTail' equals FALSE and 'smallest' equals FALSE,
 * it returns the largest 'k' that satisfies:
 *    sum(i=k, inf, pmf(i, lambda)) >= prob
 * 
 * @param prob - probability (must be greater than 0 and less than 1)
 * @param lambda - mean number of successes
 * @param smallest - see above (default: TRUE)
 * @param lowerTail - see above (default: TRUE)
 * 
 * @return quantile for the given 'prob', depending on 'smallest' and 'lowerTail'. See above for details.
 *  
 * @throw StatisticsException if any argument is invalid
 */
template <typename F, typename I>
I math::PoissonDist::quant(
          const F& prob,
          const F& lambda,
          const bool smallest,
          const bool lowerTail
        ) throw(math::StatisticsException)
{
    // sanity check
    math::PoissonDist::__private::__checkParams(static_cast<I>(0), lambda);

    if ( prob < math::NumericUtil::getEPS<F>() || 
         prob > (static_cast<F>(1) - math::NumericUtil::getEPS<F>()) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    /*
     * If lowerTail==false, the algorithm is (almost!) identical
     * except that 'prob' must be replaced by its complement.
     */
    const F PP = ( true==lowerTail ? prob : static_cast<F>(1)-prob );
    const bool cdfIncl = !lowerTail;

    const I Imax = std::numeric_limits<I>::max();

    const I ZERO = static_cast<I>(0);
    const I ONE = static_cast<I>(1);

    /*
     * Algorithm proposed in Numerical Recipes, 3rd Edition, page 337.
     * Additional integer range checks have been introduced to improve robustness.
     */

    I k, kl, ku, step;


    // Initial boundaries of the search interval:

    k = std::max(static_cast<I>(lambda), static_cast<I>(5));
    step = ONE;

    F cdf = math::PoissonDist::prob<F, I>(k, lambda, cdfIncl);
    if ( PP < cdf )
    {
        do
        {
            k = ( k>=step ? k-step : ZERO );
            step = ( step<(Imax-step) ? step * static_cast<I>(2) : Imax );
            cdf = math::PoissonDist::prob<F, I>(k, lambda, cdfIncl);
        }
        while ( PP < cdf );

        kl = k;
        const I s2 = step / static_cast<I>(2);
        ku = ( k<=(Imax-s2) ? k + s2 : Imax );
    }
    else
    {
        do
        {
            k = ( k<=(Imax-step) ? k+step : Imax );
            step = ( step<(Imax-step) ? step * static_cast<I>(2) : Imax );
            cdf = math::PoissonDist::prob<F, I>(k, lambda, cdfIncl);

            if ( Imax == k )
            {
                // prevent possible infinite loop if 'cdf'
                // is still less than 'prob'
                break;  // out of do-while
            }
        }
        while ( PP > cdf );

        const I s2 = step / static_cast<I>(2);
        kl = ( k>=s2 ? k-s2 : ZERO );
        ku = k;
    }

    // Narrow the search interval using bisection:
    while ( (ku-kl) > ONE )
    {
        /*
         * Equivalent to (kl + ku) / 2, but more robust
         * Note: when a positive number is casted to an integer,
         * it is always rounded down.
         */
        k = kl + (ku - kl) / static_cast<I>(2);

        if ( PP < math::PoissonDist::prob<F, I>(k, lambda, cdfIncl) )
        {
            ku = k;
        }
        else
        {
            kl = k;
        }
    }

    // Final adjustment of the result according to 'smallest'
    k = kl;
    cdf = math::PoissonDist::prob<F, I>(k, lambda, !cdfIncl);

    if ( true==smallest && cdf>=PP )
    {
        for ( ; k>ZERO && math::PoissonDist::prob<F, I>(k-ONE, lambda, !cdfIncl)>PP; --k );
    }
    else if ( true==smallest && cdf<PP )
    {
        for ( ; math::PoissonDist::prob<F, I>(k, lambda, !cdfIncl)<PP; ++k )
        {
            if ( Imax == k )
            {
                throw math::StatisticsException(math::StatisticsException::INTEGER_OUT_OF_RANGE);
            }
        }
    }
    else if ( false==smallest && cdf>=PP )
    {
        for ( ; k>ZERO && math::PoissonDist::prob<F, I>(k, lambda, !cdfIncl)>PP; --k );
    }
    else  // if ( false==smallest && cdf<PP ) 
    {
        for ( ; math::PoissonDist::prob<>(k+ONE, lambda, !cdfIncl)<PP; ++k )
        {
            if ( k == (Imax - ONE) )
            {
                throw math::StatisticsException(math::StatisticsException::INTEGER_OUT_OF_RANGE);
            }
        }
    }

    return k;
}
