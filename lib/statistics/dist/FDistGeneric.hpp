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
 * An internal header file, it should not be included directly.
 * @headername{FDist.h}
 *
 * Declaration of functions within the namespace FDist
 * that perform various F-distribution (a.k.a.Fisher - Snedecor distribution)
 * related operations, such as calculation of upper and lower tail
 * probabilities and quantiles, probability distribution function, etc.
 */


/*
 * More details about the F distribution (a.k.a. Fisher - Snedecor distibution):
 * https://en.wikipedia.org/wiki/F-distribution
 */


#ifndef _MATH_FDISTGENERIC_HPP_
#define _MATH_FDISTGENERIC_HPP_


#include <cmath>
#include <algorithm>

#include "../settings/probdist_settings.h"

#include "util/NumericUtil.hpp"
#include "specfun/SpecFunGeneric.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various F-distribution (a.k.a.Fisher-Snedecor distribution)
 * related operations, such as calculation of upper and lower
 * tail probabilities and quantiles, probability distribution
 * function, etc.
 */
namespace FDist
{


namespace __private
{

/*
 * Checks the validity of arguments 'x' (must be greater than or equal
 * to 0) and 'd1' and 'd2' (both must be strictly greater than 0).
 *
 * @param x - value to check
 * @param d1 - F-distribution's parameter 1 (degrees of freedom) to check
 * @param d2 - F-distribution's parameter 2 (degrees of freedom) to check
 *
 * @throw StatisticsException if any argument is invalid.
 */
template <typename F>
void __checkParams(const F& x, const F& d1, const F& d2)
{
    if ( x < static_cast<F>(0) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    if ( d1 <= math::NumericUtil::getEPS<F>() ||
         d2 <= math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_DF);
    }
}

}  // namespace math::FDist::__private



/**
 * Value of the F-distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @param x - value to be calculated its probability distribution function
 * @param d1 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 * @param d2 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 *
 * @return pdf(x, d1, d2)
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F pdf(
	  const F& x,
	  const F& d1,
	  const F& d2
	)
{
    /*
     * Probability density function of the F-distribution:
     *
     *                  +-------------------+
     *                 /        d1     d2
     *         +      /   (d1*x)   * d2
     *          \    /  ------------------
     *           \  /               d1+d2
     *            \/     (d1*x + d2)
     *   pdf = -------------------------------
     *                   /  d1    d2  \
     *             x * B | ----, ---- |
     *                   \   2     2  /
     */

    // sanity check:
    math::FDist::__private::__checkParams<F>(x, d1, d2);

    F pdf = static_cast<F>(0);

    try
    {
        if ( true == math::NumericUtil::isZero<F>(x) )
        {
            /*
             * x ~= 0 :
             *
             * This situation requires special handling to prevent possible
             * division by zero:
             *
             * The pdf can be rewritten to:
             *
             *               d1/2     d2/2    d1/2 - 1
             *             d1     * d2     * x
             *   pdf = ----------------------------------
             *                                 (d1+d2)/2
             *          B(d1/2,d2/2) * (d1*x+d2)
             *
             * From this expression it is evident that the pdf is only defined when
             * 'd1' is strictly greater than 2, when pdf equals 0. When 'd1' is less than
             * or equal to 2, exponentiation of zero by a negative value (or by 0) would occur
             * which is not defined (blows up towards +infinity).
             */

            if ( d1 <= static_cast<F>(2) )
            {
                throw math::StatisticsException(math::StatisticsException::UNDEFINED);
            }
            else
            {
                pdf = static_cast<F>(0);
            }
        }
        else
        {
            pdf = std::sqrt(std::pow(d1*x, d1) * std::pow(d2, d2) / std::pow(d1*x+d2, d1+d2) ) /
                  ( x * math::SpecFun::beta(d1/static_cast<F>(2), d2/static_cast<F>(2)) );
        }
    }
    catch ( const math::SpecFunException& sfex )
    {
        // this exception should never be thrown
    }

    return pdf;
}



/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified F-distribution.
 *
 * @param x - value
 * @param d1 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 * @param d2 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 * @param lowerTail - if true, returns P(X<x), otherwise P(X>x) (default: true)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F prob(
	  const F& x,
	  const F& d1,
	  const F& d2,
	  const bool lowerTail = true
   )
{
    /*
     * The 'cdf' represents probability that a value from the
     * F distribution is smaller than a specified value 'x' and can be
     * expressed as:
     *
     *            /  d1    d2  \
     *   cdf = I  | ----, ---- |
     *          F \   2     2  /
     *
     * where I(a,b) denotes the regularized lowered incomplete beta function,
     * implemented in the namespace math::SpecFun,
     * and F denotes the following expression:
     *
     *         d1 * x
     *   F = -----------
     *        d1*x + d2
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::FDist::__private::__checkParams<F>(x, d1, d2);

    try
    {
        // Tolerance for the algorithm that evaluates the incomplete gamma function:
        const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                      static_cast<F>(PROBDIST_PROB_TOL_DEN);


        const F f = d1 * x / (d1 * x + d2);

        const F cdf = math::SpecFun::incBetaLowerReg<F>(
              d1 / static_cast<F>(2), d2 / static_cast<F>(2), f, TOL );

        // the return value depends on 'lowerTail':
        return ( true==lowerTail ? cdf : static_cast<F>(1)-cdf );
    }
    catch ( const math::SpecFunException& sfex )
    {
        throw math::StatisticsException::OPERATION_FAILED;
    }
}



/**
 * Probability in the F-distribution with the specified values of 'd1' and 'd2'
 * (degrees of freedom) that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param d1 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 * @param d2 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 *
 * @return P(a<X<b) or P(b<X<a)
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F probInt(
	  const F& a,
	  const F& b,
	  const F& d1,
	  const F& d2
	)
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

    return math::FDist::prob<F>(to, d1, d2, true) -
           math::FDist::prob<F>(from, d1, d2, true);
}



/**
 * Quantile function for the F-distribution with the specified
 * parameters 'd1' and 'd2' (degrees of freedom).
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param d1 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 * @param d2 - F-distribution's parameter 1 (degrees of freedom), must be greater than 0
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if any argument is invalid
 */
template <typename F>
F quant(
	  const F& p,
	  const F& d1,
	  const F& d2,
	  const bool lowerTail = true
	)
{
    // sanity check
    math::FDist::__private::__checkParams<F>(static_cast<F>(1), d1, d2);
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
     *                /  d1    d2  \
     *   cdf = p = I  | ----, ---- |
     *              F \   2     2  /
     *
     * the following expression for 'F' can be derived quickly:
     *
     *        -1 /  d1    d2  \
     *   F = I   | ----, ---- |
     *        p  \   2     2  /
     *
     * Then 'x' can be derived from 'F' as follows:
     *
     *           F * d2
     *   x = --------------
     *        d1 * (1 - F)
     */

    try
    {
        // The actual probability in the expressions above depends on 'lowerTail':
        const F P = (true==lowerTail ? p : static_cast<F>(1)-p);

        // Tolerance for the algorithm that evaluates iInv():
        const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                      static_cast<F>(PROBDIST_PROB_TOL_DEN);

        const F f = math::SpecFun::incBetaLowerRegInv<F>(
                                            d1 / static_cast<F>(2),
                                            d2 / static_cast<F>(2),
                                            P,
                                            TOL );

        return f * d2 / ( d1 * (static_cast<F>(1) - f) );
    }
    catch ( const math::SpecFunException& sfex )
    {
        /*
         * The validity of 'p' is checked beforehand.
         * Although very unlikely, it is theoretically possible that
         * the algorithm to evaluate incBetaLowerRegInv will
         * not converge.
         */
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }
}


}  // namespace FDist

}  // namespace math


#endif  // _MATH_FDISTGENERIC_HPP_
