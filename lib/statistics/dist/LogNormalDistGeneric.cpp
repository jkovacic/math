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
 * Implementation of functions within the namespace LogNormalDist
 * that perform various log-normal distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, calculation of z statistics, probability
 * distribution function, etc.
 */


/*
 * More details about the log-normal distribution:
 * https://en.wikipedia.org/wiki/Log-normal_distribution
 */


// No #include "LogNormalDistGeneric.hpp" !!!
#include <cmath>
#include <algorithm>

#include "../settings/probdist_settings.h"

#include "util/math_constant.h"
#include "util/NumericUtil.hpp"
#include "specfun/SpecFunGeneric.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/StatisticsException.hpp"


namespace math {  namespace LogNormalDist {  namespace __private
{

/*
 * Checks the validity of the given standard deviation.
 * It must be strictly greater than zero.
 *
 * @param sigma - standard deviation to check
 *
 * @throw StatisticsException if 'sigma' is not valid.
 */
template <typename F>
void __checkSigma(const F& sigma)
{
    if ( sigma < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_STDEV);
    }
}


/*
 * Checks the validity of the given argument 'x'.
 * It must be strictly greater than zero.
 *
 * @param x - argument to check
 *
 * @throw StatisticsException if 'x' is not valid.
 */
template <typename F>
void __checkArg(const F& x)
{
    if ( x < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }
}

}}}  // namespace math::LogNormalDist::__private


/**
 * Obtains the z-statistics (a.k.a. standard score) of 'x', i.e.
 * the distance of 'ln x' from the mean, expressed in number of
 * standard deviations.
 *
 * @param x - value to be calculated its z-statistics
 * @param mu - log-normal distribution's mean value (default: 0)
 * @param sigma - log-normal distribution's standard deviation (default: 1)
 *
 * @return z statistics of 'x' for the given log-normal distribution parameters
 *
 * @throw StatisticsException if 'x' or the standard deviation is invalid
 */
template <typename F>
F math::LogNormalDist::getZ(
          const F& x,
          const F& mu,
          const F& sigma
        )
{
    /*
     *
     *         ln(x) - mu
     *   z = --------------
     *           sigma
     */

    // sanity check:
    math::LogNormalDist::__private::__checkSigma<F>(sigma);
    math::LogNormalDist::__private::__checkArg<F>(x);

    return (std::log(x) - mu) / sigma;
}


/**
 * Converts z-statistics (a.k.a standard score) back to the
 * original value for the log-normal distribution with specified
 * parameters.
 *
 * @param z - z-statistics to be converted
 * @param mu - log-normal distribution's mean value (default: 0)
 * @param sigma - log-normal distribution's standard deviation (default: 1)
 *
 * @return x value within the log-normal distribution
 *
 * @throw StatisticsException if the standard deviation is invalid
 */
template <typename F>
F math::LogNormalDist::getX(
          const F& z,
          const F& mu,
          const F& sigma
        )
{
    /*
     *
     *      mu + z * sigma
     * x = e
     *
     */

    // sanity check:
    math::LogNormalDist::__private::__checkSigma<F>(sigma);

    return std::exp(mu + z * sigma);
}


/**
 * Value of the log-normal distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @param x - value to be calculated its probability distribution function
 * @param mu - log-normal distribution's mean value (default: 0)
 * @param sigma - log-normal distribution's standard deviation (default: 1)
 *
 * @return pdf(x, mu, sigma)
 *
 * @throw StatisticsException if 'x' or the standard deviation is invalid
 */
template <typename F>
F math::LogNormalDist::pdf(
      const F& x,
      const F& mu,
      const F& sigma
    )
{
    /*
     * Probability density function of the log-normal distribution:
     *
     *                                2
     *                    (ln(x) - mu)
     *                - -----------------
     *                               2
     *                      2 * sigma
     *              e
     *   pdf =   ----------------------------
     *                            ------+
     *              x * sigma * \/ 2*pi
     *
     */

    // sanity check:
    math::LogNormalDist::__private::__checkSigma<F>(sigma);
    math::LogNormalDist::__private::__checkArg<F>(x);

    const F z = (std::log(x) - mu) / sigma;
    return static_cast<F>(MATH_CONST_SQRT_INV_2_PI) *
           std::exp( -z*z / static_cast<F>(2) ) / (sigma * x);
}



/**
 * Probability in the log-normal distribution with specified
 * parameters  that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param mu - log-normal distribution's mean value (default: 0)
 * @param sigma - log-normal distribution's standard deviation (default: 1)
 *
 * @return P(a<X<b) or P(b<X<a)
 *
 * @throw StatisticsException if either 'a', 'b' or the standard deviation is invalid
 */
template <typename F>
F math::LogNormalDist::probInt(
      const F& a,
      const F& b,
      const F& mu,
      const F& sigma
    )
{
    /*
     *              b
     *              /
     *   P(a<X<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *              a
     */

    // sanity check:
    math::LogNormalDist::__private::__checkSigma<F>(sigma);
    math::LogNormalDist::__private::__checkArg<F>(a);
    math::LogNormalDist::__private::__checkArg<F>(b);

    const F from = std::min(a, b);
    const F to = std::max(a, b);

    return math::LogNormalDist::prob<F>(to, mu, sigma, true) -
           math::LogNormalDist::prob<F>(from, mu, sigma, true);
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified log-normal distribution.
 *
 * @param x - value
 * @param mu - log-normal distribution's mean value (default: 0)
 * @param sigma - log-normal distribution's standard deviation (default: 1)
 * @param lowerTail - if true, returns P(t<x), otherwise P(t>x) (default: true)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if 'x' or the standard deviation is invalid
 */
template <typename F>
F math::LogNormalDist::prob(
      const F& x,
      const F& mu,
      const F& sigma,
      const bool lowerTail
    )
{
    /*
     *
     * As evident from:
     * https://en.wikipedia.org/wiki/Log-normal_distribution
     * cdf can be expressed as:
     *
     *                +-                             -+
     *             1  |         /    ln(x) - mu     \ |
     *   cdf(x) = --- | 1 + erf | ----------------- | | =
     *             2  |         \  sigma * sqrt(2)  / |
     *                +-                             -+
     *
     *             1         /      ln(x) - mu     \
     *          = --- * erfc | - ----------------- |
     *             2         \    sigma * sqrt(2)  /
     *
     * where 'erfc' is the so called complement error function, implemented
     * in the namespace math::SpecFun.
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::LogNormalDist::__private::__checkSigma<F>(sigma);
    math::LogNormalDist::__private::__checkArg<F>(x);

    // Tolerance for the last Taylor series term
    const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                  static_cast<F>(PROBDIST_PROB_TOL_DEN);

    /*
     *          ln(x) - mu            sqrt(2) * (ln(x) - mu)
     *   z = -----------------  =  ---------------------------
     *        sigma * sqrt(2)                2 * sigma
     */

    const F z = static_cast<F>(MATH_CONST_SQRT_INV_2) * ( std::log(x) - mu ) / sigma;

    /*
     * Calculate the return value:
     *
     *   cdf = 0.5 * (1 + erf(z)) = 0.5 * erfc(-z)
     */

    const F cdf = math::SpecFun::erfc<F>(-z, TOL) / static_cast<F>(2);

    // the return value depends on 'lowerTail':
    return ( true==lowerTail ? cdf : static_cast<F>(1)-cdf );
}


/**
 * Quantile function for the specified log-normal distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param mu - log-normal distribution's mean value (default: 0)
 * @param sigma - log-normal distribution's standard deviation (default: 1)
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if either 'sigma' or 'p' is invalid
 */
template <typename F>
F math::LogNormalDist::quant(
      const F& p,
      const F& mu,
      const F& sigma,
      const bool lowerTail
    )
{
    // sanity check
    math::LogNormalDist::__private::__checkSigma<F>(sigma);
    if ( p <= math::NumericUtil::getEPS<F>() ||
         p >= static_cast<F>(1)-math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    /*
     * Quantile is a value of 'x' that solves the nonlinear equation:
     *
     *   cdf(x) = p
     *
     * If the cumulative density function is expressed as:
     *
     *                          1         /      ln(x) - mu     \
     *   cdf(x,mu,sigma) = p = --- * erfc | - ----------------- |
     *                          2         \    sigma * sqrt(2)  /
     *
     * the following expression for 'x' can be derived quickly:
     *
     *        mu - sqrt(2) * sigma * erfcInv(2 * p)
     *   x = e
     *
     */

    try
    {
        // The actual probability in the expressions above depends on 'lowerTail':
        const F P = (true==lowerTail ? p : static_cast<F>(1)-p);

        // Tolerance for the algorithm that evaluates erfcInv():)
        const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                      static_cast<F>(PROBDIST_PROB_TOL_DEN);

        return std::exp( mu - math::SpecFun::erfcInv<F>(static_cast<F>(2) * P, TOL) *
                              sigma * static_cast<F>(MATH_CONST_SQRT_2) );
    }
    catch ( const math::SpecFunException& sfex )
    {
        /*
         * The validity of 'p' is checked beforehand.
         * Although very unlikely, it is theoretically possible that
         * the algorithm to evaluate erfcInv will not converge.
         */
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }
}
