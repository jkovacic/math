/*
Copyright 2014, Jernej Kovacic

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
 * Implementation of functions within the namespace NormalDist
 * that perform various normal distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, calculation of z statistics, probability
 * distribution function, etc.
 */


/*
 * More details about the normal distribution:
 * https://en.wikipedia.org/wiki/Normal_distribution
 */


// No #include "NormalDistGeneric.hpp" !!!
#include <cmath>
#include <algorithm>

#include "../settings/stat_settings.h"

#include "util/math_constant.h"
#include "util/NumericUtil.hpp"
#include "specfun/SpecFunGeneric.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/StatisticsException.hpp"

namespace math {  namespace NormalDist {  namespace __private
{

/*
 * Checks the validity of the given standard deviation.
 * It must be strictly greater than zero.
 *
 * @param sigma - standard deviation to check
 *
 * @throw StatisticsException if 'sigma' is not valid.
 */
template <class T>
void __checkSigma(const T& sigma) throw (math::StatisticsException)
{
    if ( sigma < NumericUtil::getEPS<T>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_STDEV);
    }
}

}}}  // namespace math::NormalDist::__private



/**
 * Obtains the z-statistics (a.k.a. standard score) of 'x', i.e.
 * the distance of 'x' from the mean, expressed in number of
 * standard deviations.
 *
 * @param x - value to be calculated its z-statistics
 * @param mu - normal distribution's mean value (default: 0)
 * @param sigma - normal distribution's standard deviation (default: 1)
 *
 * @return z statistics of 'x' for the given normal distribution parameters
 *
 * @throw StatisticsException if the standard deviation is invalid
 */
template <class T>
T math::NormalDist::getZ(
          const T& x,
          const T& mu,
          const T& sigma
        ) throw (math::StatisticsException)
{
    /*
     *
     *         x - mu
     *   z = -----------
     *         sigma
     */

    // sanity check:
    math::NormalDist::__private::__checkSigma<T>(sigma);

    return (x - mu) / sigma;
}


/**
 * Converts z-statistics (a.k.a standard score) back to the
 * original value for the normal distribution with specified
 * parameters.
 *
 * @param z - z-statistics to be converted
 * @param mu - normal distribution's mean value (default: 0)
 * @param sigma - normal distribution's standard deviation (default: 1)
 *
 * @return x value within the normal distribution
 *
 * @throw StatisticsException if the standard deviation is invalid
 */
template <class T>
T math::NormalDist::getX(
          const T& z,
          const T& mu,
          const T& sigma
        ) throw (math::StatisticsException)
{
    /*
     * x = mu + z * sigma
     */

    // sanity check:
	math::NormalDist::__private::__checkSigma<T>(sigma);

    return mu + z * sigma;
}


/**
 * Value of the normal distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @param x - value to be calculated its probability distribution function
 * @param mu - normal distribution's mean value (default: 0)
 * @param sigma - normal distribution's standard deviation (default: 1)
 *
 * @return pdf(x, mu, sigma)
 *
 * @throw StatisticsException if the standard deviation is invalid
 */
template <class T>
T math::NormalDist::pdf(
      const T& x,
      const T& mu,
      const T& sigma
    ) throw (math::StatisticsException)
{
    /*
     * Probability density function of the normal distribution:
     *
     *                            2
     *                    (x - mu)
     *                - --------------
     *                             2
     *                    2 * sigma
     *              e
     *   pdf =   -----------------------
     *                        ------+
     *              sigma * \/ 2*pi
     *
     */
    
    // sanity check:
    math::NormalDist::__private::__checkSigma<T>(sigma);

    const T z = (x - mu) / sigma;
    return static_cast<T>(MATH_CONST_SQRT_INV_2_PI) *
           std::exp( -z*z / static_cast<T>(2) ) / sigma;
}



/**
 * Probability in the normal distribution with specified
 * parameters  that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param mu - normal distribution's mean value (default: 0)
 * @param sigma - normal distribution's standard deviation (default: 1)
 *
 * @return P(a<x<b) or P(b<x<a)
 *
 * @throw StatisticsException if the standard deviation is invalid
 */
template <class T>
T math::NormalDist::probInt(
      const T& a,
      const T& b,
      const T& mu,
      const T& sigma
    ) throw (math::StatisticsException)
{
    /*
     *              b
     *              /
     *   P(a<x<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *              a
     */

    // sanity check:
    math::NormalDist::__private::__checkSigma<T>(sigma);

    const T from = std::min(a, b);
    const T to = std::max(a, b);

    return math::NormalDist::prob<T>(to, mu, sigma, true) -
           math::NormalDist::prob<T>(from, mu, sigma, true);
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified normal distribution.
 *
 * @param x - value
 * @param mu - normal distribution's mean value (default: 0)
 * @param sigma - normal distribution's standard deviation (default: 1)
 * @param lowerTail - if true, returns P(t<x), otherwise P(t>x) (default: true)
 *
 * @return P(t<x) if 'lowerTail' equals true, P(t>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if the standard deviation is invalid
 */
template <class T>
T math::NormalDist::prob(
      const T& x,
      const T& mu,
      const T& sigma,
      bool lowerTail
    ) throw (math::StatisticsException)
{
    /*
     * The cdf represents probability that a value from the normal
     * distribution is smaller than a specified value 'x' and can be
     * expressed as:
     *
     *                     x                   x
     *                     /              1    /
     *   cdf(x) = P(t<x) = | pdf(t) dt = --- + | pdf(t) dt
     *                     /              2    /
     *                   -inf                 mu
     *
     * One approach to calculate this would be numerical integration,
     * however there are more efficient methods available.
     *
     * As evident from:
     * https://en.wikipedia.org/wiki/Normal_distribution
     * cdf can be further expressed as:
     *
     *                +-                             -+
     *             1  |         /      x - mu       \ |
     *   cdf(x) = --- | 1 + erf | ----------------- | | =
     *             2  |         \  sigma * sqrt(2)  / |
     *                +-                             -+
     *
     *             1         /        x - mu       \
     *          = --- * erfc | - ----------------- |
     *             2         \    sigma * sqrt(2)  /
     * 
     * where 'erf' is the so called error function, implemented
     * in the namespace math::SpecFun.
     *
     * For the upper tail probability ( P(t>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::NormalDist::__private::__checkSigma<T>(sigma);

    // Tolerance for the last Taylor series term
    const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) /
                  static_cast<T>(STAT_DIST_PROB_TOL_DEN);

    /*
     *            x - mu            sqrt(2) * (x - mu)
     *   z = -----------------  =  --------------------
     *        sigma * sqrt(2)            2 * sigma
     */

    const T z = static_cast<T>(MATH_CONST_SQRT_INV_2) * ( x - mu ) / sigma;

    /*
     * Calculate the return value:
     *
     *   cdf = 0.5 * (1 + erf(z)) = 0.5 * erfc(-z)
     */

    const T cdf = math::SpecFun::erfc<T>(-z, TOL) / static_cast<T>(2);

    // the return value depends on 'lowerTail':
    return ( true==lowerTail ? cdf : static_cast<T>(1)-cdf );
}


/**
 * Quantile function for the specified normal distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(t<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(t>x) = p
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param mu - normal distribution's mean value (default: 0)
 * @param sigma - normal distribution's standard deviation (default: 1)
 * @param lowerTail - if true, returns q(t<x), otherwise q(t>x) (default: true)
 *
 * @return x: P(t<x) if 'lowerTail' equals true, x: P(t>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if either 'sigma' or 'p' is invalid
 */
template <class T>
T math::NormalDist::quant(
      const T& p,
      const T& mu,
      const T& sigma,
      bool lowerTail
    ) throw (math::StatisticsException)
{
    // sanity check
    math::NormalDist::__private::__checkSigma<T>(sigma);
    if ( p<=static_cast<T>(0) || p>=static_cast<T>(1) )
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
     *                          1         /        x - mu       \
     *   cdf(x,mu,sigma) = p = --- * erfc | - ----------------- |
     *                          2         \    sigma * sqrt(2)  /
     * 
     * the following expression for 'x' can be derived quickly:
     * 
     *   x = mu - sqrt(2) * sigma * erfcInv(2 * p)
     */

    try
    {
        // The actual probability in the expressions above depends on 'lowerTail':
        const T P = (true==lowerTail ? p : static_cast<T>(1)-p);

        // Tolerance for the algorithm that evaluates erfcInv():)
        const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) / 
                      static_cast<T>(STAT_DIST_PROB_TOL_DEN);

        return mu - math::SpecFun::erfcInv<T>(static_cast<T>(2) * P, TOL) * 
                    sigma * static_cast<T>(MATH_CONST_SQRT_2);
    }
    catch ( math::SpecFunException& spex )
    {
        /*
         * The validity of 'p' is checked beforehand.
         * Although very unlikely, it is theoretically possible that
         * the algorithm to evaluate erfcInv will not converge.
         */
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }
}
