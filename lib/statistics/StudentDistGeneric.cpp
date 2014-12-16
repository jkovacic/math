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
#include <cstddef>
#include <cmath>
#include <algorithm>

#include "../settings/stat_settings.h"
#include "util/NumericUtil.hpp"
#include "util/IFunctionGeneric.hpp"
#include "specfun/SpecFunGeneric.hpp"
#include "calculus/IntegGeneric.hpp"

#include "exception/StatisticsException.hpp"
#include "exception/FunctionException.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/CalculusException.hpp"


// Implementation of "private" functions

namespace math
{

namespace StudentDist
{

namespace __private
{

/*
 * Extension of IFunctionGeneric whose function 'func'
 * returns Student's probability distribution for the given 'x'.
 * Note that distribution's mean value, standard deviation and
 * degrees of freedom are passed via the constructor.
 *
 * @note It is assumed that the standard deviation and degrees
 *       of freedom are valid (strictly greater than 0). This check
 *       should be performed beforehand.
 */
template <class T>
class StudentDistPdf : public math::IFunctionGeneric<T>
{

private:
    const T m_mu;            // normal distribution's mean value
    const T m_sigma;         // normal distribution's standard deviation
    const T m_df;            // degrees of freedom
    T m_term;                // the multiplicand of the pdf, depends on df and sigma


public:

    /*
     * Constructor, assigns Student's distribution's mean value,
     * standard deviation and number of degrees of freedom that
     * are used for calculation of pdf.
     *
     * @note The constructor assumes that the standard deviation and
     *       degress of freedom are valid (strictly greater than 0)
     *       and does not check  this.
     *
     * @param df - degrees of freedom, must be greater than 0
     * @param mu - Student's distribution's mean value (default 0)
     * @param sigma - Student's distribution's standard deviation, must be greater than 0 (default: 1)
     */
    StudentDistPdf( const T& df, const T& mu, const T& sigma) :
        m_mu(mu), m_sigma(sigma), m_df(df)
    {
        /*
         * The term below is not very trivial to compute, hence its value
         * is calculated at instantiation of this class and stored into a variable:
         *
         *
         *                    G((df+1)/2)                            1
         * m_term = ------------------------------- = ---------------------------------
         *           G(df/2) * sqrt(pi*df) * sigma     B(1/2, df/2) * sqrt(df) * sigma
         */

        try
        {
            this->m_term = static_cast<T>(1) /
                           (math::SpecFun::beta<T>(
                                static_cast<T>(1)/static_cast<T>(2), df/2) *
                                std::sqrt(df) * sigma );
        }
        catch ( math::FunctionException& fex )
        {
            // If a valid 'df' is assumed, this exception
            // should never be thrown.
        }
    }


    /*
     * "Reimplementation" of 'func' that returns the value of the
     * Student's distribution's probability distribution function
     * (pdf) for the given x, and previously assigned number of
     * degrees of freedom, mean value and standard deviation.
     *
     * @param x - value to be calculated the probability distribution function
     *
     * @return pdf(x, df, mu, sigma)
     *
     * @throw FunctionException (never by this function)
     */
    T func(const T& x) const throw (math::FunctionException)
    {
        /*
         *                                                                        df+1
         *                                       +-                         -+ - ------
         *                 G((df+1)/2)           |      1     /  x - mu  \ 2 |     2
         * pdf = ------------------------------- | 1 + ---- * | -------- |   |
         *        G(df/2) * sqrt(pi*df) * sigma  |      df    \   sigma  /   |
         *                                       +-                         -+
         */

  	    const T t = (x - this->m_mu) / this->m_sigma;

  	    /*
  	     * If an exponent is not integer, the power can be calculated as:
  	     *
  	     *    b      b*ln(a)
  	     *   a   =  e
  	     */
        return this->m_term * std::exp((-(this->m_df+static_cast<T>(1))/static_cast<T>(2)) *
               std::log(static_cast<T>(1) + t*t/this->m_df));
    }

};  // class StudentDistPdf


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
template <class T>
void __checkParams(const T& sigma, const T& df) throw (math::StatisticsException)
{
    if ( sigma < NumericUtil::getEPS<T>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_STDEV);
    }

    if ( df < math::NumericUtil::getEPS<T>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_DF);
    }
}

}  // namespace __private
}  // namespace StudentDist
}  // namespace math



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
template <class T>
T math::StudentDist::getT(
        const T& x,
        size_t n,
        const T& mu,
        const T& s
      ) throw (math::StatisticsException)
{
    /*
     *
     *           x - mu         sqrt(n) * (x - mu)
     *   t  =  -----------  =  --------------------
     *          s/sqrt(n)               s
     *
     */

    // sanity check
    math::StudentDist::__private::__checkParams<T>(s, static_cast<T>(n));

    return (x-mu) * std::sqrt(static_cast<T>(n)) / s;
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
template <class T>
T math::StudentDist::getX(
        const T& t,
        size_t n,
        const T& mu,
        const T& s
      ) throw (math::StatisticsException)
{
    /*
     *
     *              t * s
     *   x = mu + ---------
     *             sqrt(n)
     *
     */

    // sanity check
    math::StudentDist::__private::__checkParams<T>(s, static_cast<T>(n));

    return mu + t * s / std::sqrt(static_cast<T>(n));
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
template <class T>
T math::StudentDist::pdf(
        const T& x,
        const T& df,
        const T& mu,
        const T& sigma
      ) throw (math::StatisticsException)
{
    // sanity check
    math::StudentDist::__private::__checkParams<T>(sigma, df);

    T retVal = static_cast<T>(-1);

    try
    {
        const math::StudentDist::__private::StudentDistPdf<T> pdf(df, mu, sigma);

        retVal = pdf.func(x);
    }
    catch ( const math::FunctionException &fex )
    {
        // 'pdf' is defined everywhere, hence
        // this exception should never be  thrown.
    }

    return retVal;
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
 * @return P(a<x<b) or P(b<x<a)
 *
 * @throw StatisticsException if 'df' or 'sigma' is invalid
 */
template <class T>
T math::StudentDist::probInt(
        const T& a,
        const T& b,
        const T& df,
        const T& mu,
        const T& sigma
      ) throw (math::StatisticsException)
{
    // sanity check
    math::StudentDist::__private::__checkParams<T>(sigma, df);

    /*
     *              b
     *              /
     *   P(a<x<b) = | pdf(t) dt
     *              /
     *              a
     */

    T retVal = static_cast<T>(-1);

    try
    {
        const T HSTEP = static_cast<T>(STAT_NUM_INTEG_STEP_NUM) / static_cast<T>(STAT_NUM_INTEG_STEP_NUM);

        const T from = std::min(a, b);
        const T to = std::max(a, b);
        const math::StudentDist::__private::StudentDistPdf<T> pdf(df, mu, sigma);

        retVal = math::Integ::integH<T>(pdf, from, to, HSTEP, STAT_INTEG_ALG);
    }
    catch ( const math::CalculusException &cex )
    {
        // 'pdf' is defined everywhere, the same holds for integration,
        // hence this exception should never be  thrown.
    }

    return retVal;
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified Student's distribution.
 *
 * @param x - value
 * @param df - degrees of freedom
 * @param lowerTail - if true, returns P(t<x), otherwise P(t>x) (default: true)
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return P(t<x) if 'lowerTail' equals true, P(t>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if 'df' or 'sigma' is invalid
 */
template <class T>
T math::StudentDist::prob(
        const T& x,
        const T& df,
        bool lowerTail,
        const T& mu,
        const T& sigma
      ) throw (math::StatisticsException)
{
    // sanity check
    math::StudentDist::__private::__checkParams<T>(sigma, df);

    /*
     * Lower tail probability ( P(t<x) ):
     *
     *            x                   x
     *            /              1    /
     *   P(t<x) = | pdf(t) dt = --- + | pdf(t) dt
     *            /              2    /
     *          -inf                  mu
     *
     * (also valid for x<mu).
     *
     *
     * Upper tail probability ( P(t>x) ):
     *
     *                               x                   mu
     *                          1    /              1    /
     *   P(t>x) = 1 - P(t<x) = --- - | pdf(t) dt = --- + | pdf(t) dt
     *                          2    /              2    /
     *                              mu                   x
     */

    T retVal = static_cast<T>(-1);

    try
    {
        const T HSTEP = static_cast<T>(STAT_NUM_INTEG_STEP_NUM) / static_cast<T>(STAT_NUM_INTEG_STEP_NUM);

        const T from = ( true==lowerTail ? mu : x );
        const T to = ( true==lowerTail ? x : mu );
        const math::StudentDist::__private::StudentDistPdf<T> pdf(df, mu, sigma);

        retVal = static_cast<T>(1) / static_cast<T>(2) +
                 math::Integ::integH<T>(pdf, from, to, HSTEP, STAT_INTEG_ALG);
    }
    catch ( math::CalculusException& cex )
    {
        // 'pdf' is defined everywhere, the same holds for integration,
        // hence this exception should never be  thrown.
    }

    return retVal;
}


/**
 * Quantile function for the specified Student's distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(t<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(t>x) = p
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param df - degrees of freedom
 * @param lowerTail - if true, returns q(t<x), otherwise q(t>x) (default: true)
 * @param mu - Student's distribution's mean value (default: 0)
 * @param sigma - Student's distribution's standard deviation (default: 1)
 *
 * @return x: P(t<x) if 'lowerTail' equals true, x: P(t>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if either 'p' or 'df' or 'sigma' is invalid
 */
template <class T>
T math::StudentDist::quant(
        const T& p,
        const T& df,
        bool lowerTail,
        const T& mu,
        const T& sigma
      ) throw (math::StatisticsException)
{
    // sanity check
    math::StudentDist::__private::__checkParams<T>(sigma, df);
    if ( p<=static_cast<T>(0) || p>=static_cast<T>(1) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }


    /*
     * Quantile is a value of 'x' that solves the nonlinear equation:
     *   cdf(x) = p
     *
     * Since the cumulative distribution function (cdf) is continuous, smooth and
     * monotonically increasing and its derivative (pdf) is also known and relatively
     * simple to obtain, the Newton - Raphson root finding method is chosen that is
     * guaranteed to be reasonably precise and reasonably fast.
     *
     * However, the cdf is a special function (a result of integration) that cannot
     * be expressed analytically, the root finding method, implemented by
     * math::IntegrationGeneric<T>::newton() will not be used as it would require
     * integration over the same intervals at each iteration, which may
     * be time consuming. Instead a slightly modified procedure will be applied.
     * When a new value of 'x' is calculated, the "integration" will just update
     * the previous value of cdf by the integral over the difference in x's:
     *
     *              cdf(x0) - p
     *    x = x0 - -------------
     *                pdf(x0)
     *
     *                       x
     *                       /
     *    cdf(x) = cdf(x0) + | pdf(t) dt
     *                       /
     *                      x0
     *
     * The procedure is iterated until cdf(x) gets reasonably close to
     * the desired probability 'p'. The algorithm is expected to converge
     * quite quickly, so no limit in the number of iterations is imposed.
     */

    T retVal = static_cast<T>(0);

    try
    {
        const T P = (true==lowerTail ? p : static_cast<T>(1)-p);
        const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) / static_cast<T>(STAT_DIST_PROB_TOL_DEN);
        const T HSTEP = static_cast<T>(STAT_NUM_INTEG_STEP_NUM) / static_cast<T>(STAT_NUM_INTEG_STEP_NUM);

        const math::StudentDist::__private::StudentDistPdf<T> pdf(df, mu, sigma);

        // The algorithm will start at x = mu, where the cdf is known
        // to equal exactly 1/2.
        T cdf0 = static_cast<T>(1) / static_cast<T>(2);
        T x0 = mu;

        // This initialization is necessary to perform the first loop condition check:
        T x = x0;
        T cdf = cdf0;

        // Iterate until 'f' gets close enough to 'p'
        while ( false == math::NumericUtil::isZero<T>(cdf-P, TOL) )
        {
            x = x0 - (cdf-P) / pdf.func(x0);
            cdf = cdf0 + math::Integ::integH<T>(pdf, x0, x, HSTEP, STAT_INTEG_ALG);
            x0 = x;
            cdf0 = cdf;
        }

        retVal = x;
    }
    catch ( math::FunctionException& fex )
    {
        // 'pdf' is defined everywhere, hence
        // this exception should never be  thrown.
    }
    catch ( math::CalculusException& cex )
    {
        // 'pdf' and 'cdf' are defined everywhere, consequently
      	// the integral of cdf is defined for any integration interval.
        // Hence this exception should never be  thrown.
    }

    return retVal;
}
