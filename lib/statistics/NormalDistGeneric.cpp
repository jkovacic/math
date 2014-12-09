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

#include "util/NumericUtil.hpp"
#include "util/IFunctionGeneric.hpp"
#include "exception/FunctionException.hpp"
#include "exception/StatisticsException.hpp"
#include "calculus/IntegGeneric.hpp"
#include "exception/CalculusException.hpp"


// Optimal step size (0.0001)
#define HSTEP              ( static_cast<T>(1) / static_cast<T>(10000) )
// Selected Integration algorithm for calculation of probability intervals
#define INT_ALG            math::EIntegAlg::SIMPSON_3_8

namespace math
{

namespace NormalDist
{

namespace __private
{

/*
 * Extension of IFunctionGeneric whose function 'func'
 * returns normal probability distribution for the given 'x'.
 * Note that distribution's mean value and standard deviation
 * are passed via the constructor.
 *
 * @note It is assumed that the standard deviation is valid
 *       (strictly greater than 0). This check should be performed
 *       beforehand.
 */
template <class T>
class NormDistPdf : public math::IFunctionGeneric<T>
{

private:
    const T m_mu;            // normal distribution's mean value
    const T m_sigma;         // normal distribution's standard deviation

public:
    /*
     * Constructor, assigns normal distribution's mean value
     * and standard deviation that are used for calculation of pdf.
     *
     * @note The constructor assumes that the standard deviation is
     *       valid (strictly greater than 0) and does not check  this.
     *
     * @param mu - normal distribution's mean value (default 0)
     * @param sigma - normal distribution's standard deviation, must be greater than 0 (default: 1)
     */
    NormDistPdf(
                const T& mu = static_cast<T>(0),
                const T& sigma = static_cast<T>(1) ) :
        m_mu (mu), m_sigma (sigma)
    {
        /*
         * Both properties have been assigned above,
         * so there is nothing else to do.
         */
    }

    /*
     * "Reimplementation" of 'func' to return the value of the
     * probability distribution function (pdf) for the given x,
     * and previously assigned mean value and standard deviation.
     *
     * @param x - value to be calculated the probability distribution function
     *
     * @return pdf(x, mu, sigma)
     *
     * @throw FunctionException (never by this function)
     */
    T func(const T& x) const throw (math::FunctionException)
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

        /*
         * High precision approximation for the constant 1 / sqrt(2*pi).
         *
         * According to this topic:
         * http://stackoverflow.com/questions/476212/what-is-the-precision-of-long-double-in-c
         * the precision for long double is around 18 digits.
         * In Maxima this expression is evaluated in even higher precision
         * (25 digits) and copied to the expression below. The compiler will
         * cast it and round to the precision required by the type T.
         (%i1)  fpprec : 25$
         (%i2)  bfloat(1/sqrt(2 * %pi));
         (%o2)  3.989422804014326779399461b−1
         */

    	const T z = (x - m_mu) / m_sigma;
        return static_cast<T>(0.3989422804014326779399461L) *
               std::exp( -z*z / static_cast<T>(2) ) / m_sigma;
    }
};  // class NormDistPdf


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

}  // namespace __private
}  // namespace NormalDist
}  // namespace math



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
    // sanity check:
    math::NormalDist::__private::__checkSigma<T>(sigma);

    T retVal = static_cast<T>(-1);

    try
    {
        const math::NormalDist::__private::NormDistPdf<T> pdf(mu, sigma);

        retVal = pdf.func(x);
    }
    catch ( math::FunctionException& fex )
    {
        // 'pdf' is defined everywhere, hence
        // this exception should never be  thrown.
    }

    return retVal;
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
     *   P(a<x<b) = | pdf(t) dt
     *              /
     *              a
     */

    // sanity check:
    math::NormalDist::__private::__checkSigma<T>(sigma);

    T retVal = static_cast<T>(0);

    try
    {
        const T from = std::min(a, b);
        const T to = std::max(a, b);
        const math::NormalDist::__private::NormDistPdf<T> pdf(mu, sigma);

        retVal = math::Integ::integH<T>(pdf, from, to, HSTEP, INT_ALG);
    }
    catch ( math::CalculusException& cex )
    {
        // 'pdf' is defined everywhere, the same holds for integration,
    	// hence this exception should never be  thrown.
    }

    return retVal;
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

    // sanity check
    math::NormalDist::__private::__checkSigma<T>(sigma);

    T retVal = static_cast<T>(-1);

    try
    {
        const T from = ( true==lowerTail ? mu : x );
        const T to = ( true==lowerTail ? x : mu );
        const math::NormalDist::__private::NormDistPdf<T> pdf(mu, sigma);

        retVal = static_cast<T>(1) / static_cast<T>(2) +
                 math::Integ::integH<T>(pdf, from, to, HSTEP, INT_ALG);
    }
    catch ( math::CalculusException& cex )
    {
        // 'pdf' is defined everywhere, the same holds for integration,
        // hence this exception should never be  thrown.
    }

    return retVal;
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
     *   cdf(x) = p
     *
     * Since the cumulative distribution function (cdf) is continuous, smooth and
     * monotonically increasing and its derivative (pdf) is also known and relatively
     * simple obtain, the Newton - Raphson root finding method is chosen that is
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
        const T TOL = static_cast<T>(1) / static_cast<T>(10000);
        const math::NormalDist::__private::NormDistPdf<T> pdf(mu, sigma);

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
            cdf = cdf0 + math::Integ::integH<T>(pdf, x0, x, HSTEP, INT_ALG);
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


// Undef all macros:
#undef HSTEP
#undef INT_ALG
