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
#include <cstddef>
#include <cmath>
#include <algorithm>

#include "../settings/stat_settings.h"

#include "util/math_constant.h"
#include "util/NumericUtil.hpp"
#include "util/IFunctionGeneric.hpp"
#include "exception/FunctionException.hpp"
#include "exception/StatisticsException.hpp"
#include "root_find/RootFindGeneric.hpp"
#include "exception/RootFindException.hpp"

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

        const T z = (x - this->m_mu) / this->m_sigma;
        return static_cast<T>(MATH_CONST_SQRT_INV_2_PI) *
               std::exp( -z*z / static_cast<T>(2) ) / this->m_sigma;
    }
};  // class NormDistPdf


/*
 * Extension of IFunctionGeneric whose function 'func'
 * returns normal distribution's cumulative distribution function (cdf)
 * for the given 'x'. Note that distribution's mean value and standard deviation
 * are passed via the constructor.
 *
 * @note It is assumed that the standard deviation is valid
 *       (strictly greater than 0). This check should be performed
 *       beforehand.
 *
 * @note As this class is also passed to a root finding method, 'func'
 *       actually returns the value of cdf, subtracted by the desired
 *       probability 'p'. By default 'p' is set to 0 and must be set via
 *       the setP() method.
 */
template <class T>
class NormDistCdf : public math::IFunctionGeneric<T>
{

private:
    const T m_mu;            // normal distribution's mean value
    const T m_sigma;         // normal distribution's standard deviation
    T m_p;                   // the actual value of cdf will be subtracted by this value

public:

    /*
     * Constructor, assigns normal distribution's mean value
     * and standard deviation that are used for calculation of cdf.
     * The subtrahend is automatically set to 0.
     *
     * @note The constructor assumes that the standard deviation is
     *       valid (strictly greater than 0) and does not check  this.
     *
     * @param mu - normal distribution's mean value (default 0)
     * @param sigma - normal distribution's standard deviation, must be greater than 0 (default: 1)
     */
    NormDistCdf(
                const T& mu = static_cast<T>(0),
                const T& sigma = static_cast<T>(1) ) :
        m_mu (mu), m_sigma (sigma), m_p(static_cast<T>(0))
    {
        /*
         * All properties have been assigned above,
         * so there is nothing else to do.
         */
    }


    /*
     * Sets the subtrahend to subtract from the actual cdf.
     * This is useful when the class is passed to a root finding method.
     * The constructor will set this subtrahend to 0.
     *
     * @param p - desired subtrahend
     */
    void setP(const T& p)
    {
        this->m_p = p;
    }


    /*
     * "Reimplementation" of 'func' to return the value of the
     * cumulative distribution function (cdf) for the given x,
     * and previously assigned mean value and standard deviation.
     * The cdf value is additionally subtracted by the subtrahend
     * 'p' that can be set via setP().
     *
     * @param x - value to be calculated the cumulative distribution function
     *
     * @return cdf(x, mu, sigma) - p
     *
     * @throw FunctionException (never by this function)
     */
    T func(const T& x) const throw (math::FunctionException)
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
         *   cdf(x) = --- | 1 + erf | ----------------- | |
         *             2  |         \  sigma * sqrt(2)  / |
         *                +-                             -+
         *
         * where erf is the so called error function, defined as:
         *
         *                       inf
         *                2       /  -t^2
         *   erf(x) = ----------  | e    dt
         *             sqrt(pi)   /
         *                        x
         *
         * The definite integral cannot be calculated analytically,
         * however the exponential can be expanded to a Taylor series
         * and each term can be integrated separately. As evident from
         * https://en.wikipedia.org/wiki/Error_function#Taylor_series
         * the error function can thus be expanded into the following
         * Taylor series:
         *
         *                       +-      3     5      7      9        -+
         *                2      |      z     z      z      z          |
         *   erf(z) ~ ---------- | z - --- + ---- - ---- + ----- - ... |
         *             sqrt(pi)  |      3     10     42     216        |
         *                       +-                                   -+
         *
         *                        inf               i
         *                       -----            -----     2
         *                2      \        z       |   |   -z
         *   erf(z) ~ ----------  >   --------- * |   | -------
         *             sqrt(pi)  /     2*i + 1    |   |    j
         *                       -----            |   |
         *                        i=0              j=1
         *
         * The implemented algorithm will add Taylor series terms until
         * a term's absolute value drops below the specified tolerance.
         */


        // Tolerance for the last Taylor series term
        const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) /
                      static_cast<T>(STAT_DIST_PROB_TOL_DEN);

        /*
         *            x - mu            sqrt(2) * (x - mu)
         *   z = -----------------  =  --------------------
         *        sigma * sqrt(2)            2 * sigma
         */
        const T z = static_cast<T>(MATH_CONST_SQRT_INV_2) *
                    ( x - this->m_mu ) / this->m_sigma;

        // Subterms of the Taylor series:
        T zt = z * static_cast<T>(2) * static_cast<T>(MATH_CONST_SQRT_INV_PI );
        T t = zt;

        // Initial value of the 'erf':
        T erf = t;

        // Add Taylor series terms to 'erf' until term's abs. value
        // drops below TOL:
        for ( size_t i=1; false==math::NumericUtil::isZero<T>(t, TOL); ++i )
        {
            // update the product term from the algorithm described above:
            zt *= -z*z / static_cast<T>(i);
            // new Taylor series element:
            t = zt / (static_cast<T>(2) * static_cast<T>(i) + static_cast<T>(1) );
            erf += t;
        }

        /*
         * Finally calculate the return value:
         *
         * cdf = 0.5 * (1 + erf)
         *
         * Additionally m_p is subtracted from the cdf.
         */
        return ( erf + static_cast<T>(1) ) / static_cast<T>(2) - this->m_p;
    }

};  // class NormDistCdf


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
     *   P(a<x<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *              a
     */

    // sanity check:
    math::NormalDist::__private::__checkSigma<T>(sigma);

    const T from = std::min(a, b);
    const T to = std::max(a, b);

    return math::NormalDist::prob<T>(to, mu, sigma) -
           math::NormalDist::prob<T>(from, mu, sigma);
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
     * Lower tail probability ( P(t<x) ) is calculated inside the class
     * NormDistCdf. See comments inside its implementation for more details.
     *
     * For the upper tail probability ( P(t>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::NormalDist::__private::__checkSigma<T>(sigma);

    const math::NormalDist::__private::NormDistCdf<T> cdf(mu, sigma);

    T retVal = static_cast<T>(-1);

    try
    {
        retVal = cdf.func(x);

        if ( false == lowerTail )
        {
            retVal = static_cast<T>(1) - retVal;
        }
    }
    catch ( const math::FunctionException& fex )
    {
        // 'cdf' is defined everywhere,
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
     */

    T retVal = static_cast<T>(0);

    try
    {
        const T P = (true==lowerTail ? p : static_cast<T>(1)-p);

        // Tolerance for the Newton - Raphson root finding method:
        const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) / 
                      static_cast<T>(STAT_DIST_PROB_TOL_DEN);

        const math::NormalDist::__private::NormDistPdf<T> pdf(mu, sigma);
        math::NormalDist::__private::NormDistCdf<T> cdf(mu, sigma);
        cdf.setP(P);

        retVal = math::RootFind::newton<T>(cdf, pdf, mu, TOL, static_cast<size_t>(-1));
    }
    catch ( math::RootFindException& rex )
    {
        // This exception can only be thrown in an unlikely event that the
        // root find algorithm could not converge despite the large 
        // number of allowed iterations.
        throw math::StatisticsException(math::StatisticsException::OPERATION_FAILED);
    }

    return retVal;
}
