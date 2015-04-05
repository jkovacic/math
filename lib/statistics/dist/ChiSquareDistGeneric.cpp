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
 * Implementation of functions within the namespace ChiSquareDist
 * that perform various Chi-square distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability distribution function, etc.
 */


/*
 * More details about the Chi-squared distribution:
 * https://en.wikipedia.org/wiki/Chi-squared_distribution
 */


// No #include "ChiSquareDistGeneric.hpp" !!!
#include <cmath>
#include <algorithm>

#include "../settings/stat_settings.h"

#include "util/NumericUtil.hpp"
#include "specfun/SpecFunGeneric.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/StatisticsException.hpp"

namespace math {  namespace ChiSquareDist {  namespace __private
{

/*
 * Checks the validity of arguments 'x' (must be greater than or equal
 * to 0) and 'df' (must be strictly greater than 0).
 *
 * @param x - value to check
 * @param df - degrees of freedom to check
 *
 * @throw StatisticsException if any argument is invalid.
 */
template <typename F>
void __checkParams(const F& x, const F& df) throw (math::StatisticsException)
{
    if ( x < static_cast<F>(0) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    if (df < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_DF);
    }
}

}}}  // namespace math::ChiSquareDist::__private



/**
 * Value of the chi-square distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @param x - value to be calculated its probability distribution function
 * @param df - degrees of freedom
 *
 * @return pdf(x, df)
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ChiSquareDist::pdf(
      const F& x,
      const F& df
    ) throw (math::StatisticsException)
{
    /*
     * Probability density function of the chi-squared distribution:
     *
     *               1          -df/2    df/2 - 1    -x/2
     *   pdf = ------------- * 2      * x         * e
     *          gamma(df/2)
     *
     */

    // sanity check:
    math::ChiSquareDist::__private::__checkParams<F>(x, df);

    const F df2 = df / static_cast<F>(2);
    return std::pow(static_cast<F>(2), -df2) * std::pow(x, df2 - static_cast<F>(1) ) *
           std::exp(-x / static_cast<F>(2)) / math::SpecFun::gamma<F>(df2);
}



/**
 * Probability in the chi-squared distribution with the specified
 * degrees of freedom that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param df - degrees of freedom
 *
 * @return P(a<X<b) or P(b<X<a)
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ChiSquareDist::probInt(
      const F& a,
      const F& b,
      const F& df
    ) throw (math::StatisticsException)
{
    /*
     *              b
     *              /
     *   P(a<X<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *              a
     */

    // Sanity check will be performed by prob()

    const F from = std::min(a, b);
    const F to = std::max(a, b);

    return math::ChiSquareDist::prob<F>(to, df, true) -
           math::ChiSquareDist::prob<F>(from, df, true);
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified chi-squared distribution.
 *
 * @param x - value
 * @param df - degrees of freedom
 * @param lowerTail - if true, returns P(X<x), otherwise P(X>x) (default: true)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ChiSquareDist::prob(
      const F& x,
      const F& df,
      const bool lowerTail
    ) throw (math::StatisticsException)
{
    /*
     * The 'cdf' represents probability that a value from the chi-squared
     * distribution is smaller than a specified value 'x' and can be
     * expressed as:
     *
     *           /  df    x  \
     *   cdf = P | ----, --- |
     *           \   2    2  /
     *
     * where P(a,x) denotes the regularized lowered incomplete gamma function,
     * implemented in the namespace math::SpecFun.
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::ChiSquareDist::__private::__checkParams<F>(x, df);

    try
    {
        // Tolerance for the algorithm that evaluates the incomplete gamma function:
        const F TOL = static_cast<F>(STAT_DIST_PROB_TOL_NUM) /
                      static_cast<F>(STAT_DIST_PROB_TOL_DEN);


        const F cdf = math::SpecFun::incGammaLowerReg<F>(
              df / static_cast<F>(2), x / static_cast<F>(2), TOL);

        // the return value depends on 'lowerTail':
        return ( true==lowerTail ? cdf : static_cast<F>(1)-cdf );
    }
    catch ( const math::SpecFunException& sfex )
    {
        throw math::StatisticsException::OPERATION_FAILED;
    }
}


/**
 * Quantile function for the chi-squared distribution with the specified
 * number of degrees of freedom.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param df - degrees of freedom
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F math::ChiSquareDist::quant(
      const F& p,
      const F& df,
      const bool lowerTail
    ) throw (math::StatisticsException)
{
    // sanity check
    math::ChiSquareDist::__private::__checkParams<F>(static_cast<F>(1), df);
    if ( p <= math::NumericUtil::getEPS<F>() ||
         p >= static_cast<F>(1) - math::NumericUtil::getEPS<F>() )
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
     *                      /  df    x  \
     *   cdf(df,x) = p =  P | ----, --- |
     *                      \   2    2  /
     *
     * the following expression for 'x' can be derived quickly:
     *
     *            -1 /  df     \
     *   x = 2 * P   | ----, p |
     *               \   2     /
     */

    try
    {
        // The actual probability in the expressions above depends on 'lowerTail':
        const F P = (true==lowerTail ? p : static_cast<F>(1)-p);

        // Tolerance for the algorithm that evaluates pInv():)
        const F TOL = static_cast<F>(STAT_DIST_PROB_TOL_NUM) /
                      static_cast<F>(STAT_DIST_PROB_TOL_DEN);

        return math::SpecFun::incGammaLowerRegInv<F>(df / static_cast<F>(2), P, TOL) *
               static_cast<F>(2);
    }
    catch ( const math::SpecFunException& sfex )
    {
        /*
         * The validity of 'p' is checked beforehand.
         * Although very unlikely, it is theoretically possible that
         * the algorithm to evaluate incGammaLowerRegInv will 
         * not converge.
         */
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }
}
