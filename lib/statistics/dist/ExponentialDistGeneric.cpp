/*
Copyright 2017, Jernej Kovacic

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
 * Implementation of functions within the namespace ExponentialDist
 * that perform various exponential distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability distribution function, etc.
 */


/*
 * More details about the exponential distribution:
 * https://en.wikipedia.org/wiki/Exponential_distribution
 */


// no #include "ExponentialDistGeneric.h" !!!
#include <cmath>
#include <algorithm>

#include "util/NumericUtil.hpp"
#include "exception/StatisticsException.hpp"


namespace math {  namespace ExponentialDist {  namespace __private {

/*
 * Checks validity of given parameters.
 *
 * 'x' must be non-negative and 'lambda' must
 * be strictly greater than 0.
 *
 * @param x - argument to check
 * @param lambda - rate of the exponential distribution to check
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
void __checkParams(const F& x, const F& lambda) throw(math::StatisticsException)
{
    /*
     * lambda > 0
     * x >= 0
     */

    if ( x < static_cast<F>(0) ||
         lambda < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }
}

}}}  // namespace math::ExponentialDist::__private




/**
 * Exponential distribution's probability distribution function (pdf)
 *
 * @param x - value to be calculated its probability distribution function
 * @param lambda - rate of the exponential distribution
 *
 * @return probability distribution function for the given 'x'
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ExponentialDist::pdf(
          const F& x,
          const F& lambda
        ) throw( math::StatisticsException)
{
    /*
     * Probability distribution function of the exponential distribution
     * is defined as:
     *
     *                             -lambda * x
     *   pdf = P(X==x) = lambda * e
     *
     */

    // sanity check
    math::ExponentialDist::__private::__checkParams<F>(x, lambda);

    return lambda * std::exp(-lambda * x);
}


/**
 * Probability in the exponential distribution with specified
 * parameters  that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param lambda - rate of the exponential distribution
 *
 * @return P(a<X<b) or P(b<X<a)
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ExponentialDist::probInt(
          const F& a,
          const F& b,
          const F& lambda
        ) throw (math::StatisticsException)
{
    /*
     *              b
     *              /
     *   P(a<X<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *              a
     */

    // sanity check
    math::ExponentialDist::__private::__checkParams<F>(a, lambda);
    math::ExponentialDist::__private::__checkParams<F>(b, lambda);

    const F from = std::min(a, b);
    const F to = std::max(a, b);

    return ( math::ExponentialDist::prob<F>(to, lambda) -
             math::ExponentialDist::prob<F>(from, lambda) );
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified exponential distribution.
 *
 * @param x - value
 * @param lambda - rate of the exponential distribution
 * @param lowerTail - if true, returns P(t<x), otherwise P(t>x) (default: true)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ExponentialDist::prob(
      const F& x,
      const F& lambda,
      const bool lowerTail
    ) throw (math::StatisticsException)
{
    /*
     * Cumulative distribution function of the exponential distribution
     * can be obtained as:
     *
     *                          -lambda * x
     *   cdf(x) = P(X<x) = 1 - e
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::ExponentialDist::__private::__checkParams<F>(x, lambda);

    const F el = std::exp(-lambda * x);

    // the return value depends on 'lowerTail':
    return ( true==lowerTail ? static_cast<F>(1)-el : el );
}


/**
 * Quantile function for the specified exponential distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @param p - probability (must be greater or equal to 0 and less than 1)
 * @param lambda - rate of the exponential distribution
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ExponentialDist::quant(
      const F& p,
      const F& lambda,
      const bool lowerTail
    ) throw (math::StatisticsException)
{
    // sanity check
    math::ExponentialDist::__private::__checkParams<F>(static_cast<F>(1), lambda);
    if ( p < static_cast<F>(0) ||
         p >= ( static_cast<F>(1)-math::NumericUtil::getEPS<F>() ) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    /*
     * Quantile of the exponential distribution can be obtained as
     *
     *          ln(1-p)
     *   x = - ---------
     *          lambda
     */


    // The actual probability in the expressions above depends on 'lowerTail':
    const F P = (true==lowerTail ? p : static_cast<F>(1)-p);

    return -std::log(static_cast<F>(1) - P) / lambda;
}
