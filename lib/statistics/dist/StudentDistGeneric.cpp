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
 * Implementation of functions within the namespace StudentDist
 * that perform various Student's (Gosset's) distribution related
 * operations, such as calculation of upper and lower tail
 * probabilities and quantiles, calculation of t-statistics,
 * probability distribution function, etc.
 */


/*
 * More details about the Student's distribution:
 * https://en.wikipedia.org/wiki/Student%27s_t-distribution
 */


// no #include "StudentDistGeneric.hpp" !!!
#include <cmath>
#include <algorithm>

#include "../settings/probdist_settings.h"
#include "util/math_constant.h"
#include "util/NumericUtil.hpp"
#include "specfun/SpecFunGeneric.hpp"

#include "exception/StatisticsException.hpp"
#include "exception/SpecFunException.hpp"


// Implementation of "private" functions

namespace math {  namespace StudentDist {  namespace __private
{

/*
 * Checks the validity of the given standard deviation
 * and degrees of freedom, both must be strictly greater
 * than zero.
 *
 * @param sigma - standard deviation to check
 * @param df - number of degrees of freedom to check
 *
 * @throw StatisticsException if 'sigma' or 'df' is not valid.
 */
template <typename F>
void __checkParams(const F& sigma, const F& df)
{
    if ( sigma < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_STDEV);
    }

    if ( df < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_DF);
    }
}

}}}  // namespace math::StudentDist::__private



/**
 * Obtains the t-statistics (a.k.a. standard score) of 'x', i.e.
 * the distance of 'x' from the mean, expressed in number of
 * standard deviations, divided by square root of the number
 * of elements.
 *
 * @param x - value to be calculated its t-statistics
 * @param n - number of elements in the sample (default:1)
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return t-statistics of 'x' for the given Student's distribution parameters
 *
 * @throw StatisticsException if 'n' or 'sigma' is invalid
 */
template <typename F, typename I>
F math::StudentDist::getT(
        const F& x,
        const I& n,
        const F& mu,
        const F& s
      )
{
    /*
     *
     *           x - mu         sqrt(n) * (x - mu)
     *   t  =  -----------  =  --------------------
     *          s/sqrt(n)               s
     *
     */

    // sanity check
    math::StudentDist::__private::__checkParams<F>(s, static_cast<F>(n));

    return (x-mu) * std::sqrt(static_cast<F>(n)) / s;
}


/**
 * Converts t-statistics (a.k.a standard score) back to the
 * original value for the Student's distribution with specified
 * parameters.
 *
 * @param t - t-statistics to be converted
 * @param n - number of elements in the sample (default: 1)
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return x value within the normal distribution
 *
 * @throw StatisticsException if 'n' or 'sigma' is invalid
 */
template <typename F, typename I>
F math::StudentDist::getX(
        const F& t,
        const I& n,
        const F& mu,
        const F& s
      )
{
    /*
     *
     *              t * s
     *   x = mu + ---------
     *             sqrt(n)
     *
     */

    // sanity check
    math::StudentDist::__private::__checkParams<F>(s, static_cast<F>(n));

    return mu + t * s / std::sqrt(static_cast<F>(n));
}


/**
 * Value of the Student's distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @param x - value to be calculated its probability distribution function
 * @param df - degrees of freedom
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return pdf(x, df, mu, sigma)
 *
 * @throw StatisticsException if 'df' or 'sigma' is invalid
 */
template <typename F>
F math::StudentDist::pdf(
        const F& x,
        const F& df,
        const F& mu,
        const F& sigma
      )
{
    // sanity check
    math::StudentDist::__private::__checkParams<F>(sigma, df);

    /*
     *                                                                        df+1
     *                                       +-                         -+ - ------
     *                 G((df+1)/2)           |      1     /  x - mu  \ 2 |     2
     * pdf = ------------------------------- | 1 + ---- * | -------- |   |
     *        G(df/2) * sqrt(pi*df) * sigma  |      df    \   sigma  /   |
     *                                       +-                         -+
     */

    F pdf = static_cast<F>(-1);

    try
    {
        const F t = (x - mu) / sigma;

        pdf = math::SpecFun::gamma<F>((df + static_cast<F>(1)) / static_cast<F>(2)) *
              static_cast<F>(MATH_CONST_SQRT_INV_PI) /
              ( math::SpecFun::gamma<F>(df / static_cast<F>(2)) *
                std::sqrt(df) * sigma );

        pdf *= std::pow(static_cast<F>(1) + t*t/df, -(df + static_cast<F>(1)) / static_cast<F>(2));
    }
    catch ( const math::SpecFunException& sfex )
    {
        // this exception should never be thrown.
    }

    return pdf;
}


/**
 * Probability in the Student's distribution with specified
 * parameters  that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param df - degrees of freedom
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return P(a<X<b) or P(b<X<a)
 *
 * @throw StatisticsException if 'df' or 'sigma' is invalid
 */
template <typename F>
F math::StudentDist::probInt(
        const F& a,
        const F& b,
        const F& df,
        const F& mu,
        const F& sigma
      )
{
    /*
     *              b
     *              /
     *   P(a<X<b) = | pdf(t) dt
     *              /
     *              a
     */

    // sanity check
    math::StudentDist::__private::__checkParams<F>(sigma, df);
    
    const F from = std::min(a, b);
    const F to = std::max(a, b);

    return math::StudentDist::prob<F>(to, df, true, mu, sigma) -
           math::StudentDist::prob<F>(from, df, true, mu, sigma);
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified Student's distribution.
 *
 * @param x - value
 * @param df - degrees of freedom
 * @param lowerTail - if true, returns P(X<x), otherwise P(X>x) (default: true)
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if 'df' or 'sigma' is invalid
 */
template <typename F>
F math::StudentDist::prob(
        const F& x,
        const F& df,
        const bool lowerTail,
        const F& mu,
        const F& sigma
      )
{
    /*
     * The cdf represents probability that a value from the Student's
     * distribution is smaller than a specified value 'x' and can be
     * expressed as:
     *
     *                     x                   x
     *                     /              1    /
     *   cdf(x) = P(X<x) = | pdf(t) dt = --- + | pdf(t) dt
     *                     /              2    /
     *                   -inf                 mu
     *
     * One approach to calculate this would be numerical integration,
     * however there are more efficient methods available.
     *
     * As evident from:
     * https://en.wikipedia.org/wiki/Student%27s_t-distribution#Cumulative_distribution_function
     * 't' can be defined as:
     * 
     *              df
     *  t = --------------------
     *            /  x - mu  \2
     *       df + | -------- | 
     *            \  sigma   /
     * 
     * then the cdf can be further expressed as:
     *
     *             /  1     /  df    1  \
     *             | --- I  | ----, --- |        <== x <= mu
     *             /  2   t \   2    2  /
     *   cdf(t) = {
     *             \      1     /  df    1  \
     *             | 1 - --- I  | ----, --- |    <== x > mu
     *             \      2   t \   2    2  /
     *
     * where It(a,b) is the incomplete beta function, implemented
     * in the namespace math::SpecFun.
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */
    
    // sanity check
    math::StudentDist::__private::__checkParams<F>(sigma, df);

    F cdf = static_cast<F>(-1);

    try
    {
        // Tolerance for the incomplete beta function
        const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                      static_cast<F>(PROBDIST_PROB_TOL_DEN);

        /*
         *            df
         *   t = --------------------
         *             /  x - mu  \2
         *        df + | -------- |
         *             \  sigma   /
         */

        const F ttemp = (x - mu) / sigma;

        // if 'x' is very close to 'mu', the probability is exactly 0.5
        if ( true == math::NumericUtil::isZero<F>(ttemp) )
        {
            return static_cast<F>(1) / static_cast<F>(2);
        }

        const F t = df / ( df + ttemp*ttemp );

        /*
         * Evaluate:
         * 
         *   I  (df/2, 1/2) / 2
         *    t
         */
        
        cdf = math::SpecFun::incBetaLowerReg<F>(
                df / static_cast<F>(2),
                static_cast<F>(1) / static_cast<F>(2),
                t,
                TOL );
        cdf /= static_cast<F>(2);

        // Adjust the cdf depending on 'x' (w.r.t. 'mu') and 'lowerTail':
        if ( ( x>mu && true==lowerTail ) ||
             ( x<mu && false==lowerTail ) )
        {
            cdf = static_cast<F>(1) - cdf;
        }
    }
    catch ( const math::SpecFunException& sfex )
    {
        // should only occur in a very unlikely event
        // that incomplete beta does not converge
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }

    return cdf;
}


/**
 * Quantile function for the specified Student's distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param df - degrees of freedom
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if either 'p' or 'df' or 'sigma' is invalid
 */
template <typename F>
F math::StudentDist::quant(
        const F& p,
        const F& df,
        const bool lowerTail,
        const F& mu,
        const F& sigma
      )
{
    // sanity check
    math::StudentDist::__private::__checkParams<F>(sigma, df);
    if ( p <= math::NumericUtil::getEPS<F>() ||
         p >= static_cast<F>(1)-math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    /*
     * Quantile is a value of 'x' that solves the nonlinear equation:
     *   cdf(x) = p
     *
     * If x < mu ==> p < 0.5
     * In this case 't' (as defined in comments of prob()) can be calculated as:
     * 
     *        -1 /  df    1  \
     *   t = I   | ----, --- |
     *        2p \   2    2  /
     * 
     * If p > 0.5, its complement (1-p) should be obtained and a
     * solution being greater than 'mu' should be considered.
     * 
     * If 't' is defined as:
     * 
     *                df
     *   t = --------------------
     *             /  x - mu  \2
     *        df + | -------- |
     *             \   sigma  /
     * 
     * then 'x' can be quickly derived from this equation. Note that
     * the equation has two solutions, consider the one being either smaller
     * or greater than 'mu', depending on the initial 'p'.
     * 
     *                        +---------------+
     *   x = mu +/- sigma * \/ df * ( 1/t - 1)
     * 
     * From properties of the regularized incomplete beta function
     * and its inverse it is evident that 't' will always be greater than 0
     * (except when p==0 ==> x=-inf) and less than 1 (except when p==0.5 ==> x=mu) 
     * hence the expression for 'x' is defined always except the two special cases
     * that are handled separately.
     */

    F x = static_cast<F>(0);

    // separate handling of p==1/2:
    if ( true==math::NumericUtil::isZero<F>(p-static_cast<F>(1)/static_cast<F>(2)) )
    {
        return mu;
    }

    try
    {
        const F P = (true==lowerTail ? p : static_cast<F>(1)-p);

        // consider the p<1/2:
        const F PP = std::min(P, static_cast<F>(1)-P);

        // Tolerance for the inverse incomplete beta function:
        const F TOL = static_cast<F>(PROBDIST_PROB_TOL_NUM) /
                      static_cast<F>(PROBDIST_PROB_TOL_DEN);

        // t = Iinv(df/2, 1/2, 2*p):
        const F t = math::SpecFun::incBetaLowerRegInv<F>(
                    df / static_cast<F>(2),
                    static_cast<F>(1) / static_cast<F>(2),
                    static_cast<F>(2) * PP,
                    TOL );

        // From t obtain tt = sqrt(df * (1/t - 1)):
        const F tt = sigma * std::sqrt(df * (static_cast<F>(1) - t) / t);
    
        // And finally add or subtract 'tt' to/from 'mu', depending on the initial 'p':
        x = ( P > static_cast<F>(1) / static_cast<F>(2) ? mu + tt : mu - tt );
    }
    catch ( const math::SpecFunException& sfex )
    {
        // This exception can only be thrown in an unlikely event
        // that the inverse incomplete beta function does not converge.
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }

    return x;
}
