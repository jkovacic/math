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
#include "util/math_constant.h"
#include "util/NumericUtil.hpp"
#include "util/IFunctionGeneric.hpp"
#include "specfun/SpecFunGeneric.hpp"
#include "root_find/RootFindGeneric.hpp"

#include "exception/StatisticsException.hpp"
#include "exception/FunctionException.hpp"
#include "exception/SpecFunException.hpp"
#include "exception/RootFindException.hpp"


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
            // The first version (with two gamma functions) is computationally
            // more efficient as evaluation of the beta function requires
            // evaluation of three different gamma functions.
            this->m_term = math::SpecFun::gamma<T>((df + static_cast<T>(1)) / static_cast<T>(2)) *
                           static_cast<T>(MATH_CONST_SQRT_INV_PI) /
                           ( math::SpecFun::gamma<T>(df / static_cast<T>(2)) *
                             std::sqrt(df) * sigma );
        }
        catch ( math::SpecFunException& sfex )
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

        return this->m_term * std::pow(static_cast<T>(1) + t*t/this->m_df,
                                       -(this->m_df + static_cast<T>(1)) / static_cast<T>(2));
    }

};  // class StudentDistPdf


/*
 * Extension of IFunctionGeneric whose function 'func'
 * returns Student's distribution's cumulative distribution function (cdf)
 * for the given 'x'. Note that distribution's degrees of freedom, mean value 
 * and standard deviation are passed via the constructor.
 *
 * @note It is assumed that the degrees of freedom and the standard deviation 
 *       are valid (strictly greater than 0). This check should be performed
 *       beforehand.
 *
 * @note As this class is also passed to a root finding method, 'func'
 *       actually returns the value of cdf, subtracted by the desired
 *       probability 'p'. By default 'p' is set to 0 and must be set via
 *       the setP() method.
 */
template <class T>
class StudDistCdf : public math::IFunctionGeneric<T>
{

private:
    const T m_df;            // Student's distribution degrees of freedom
    const T m_mu;            // Student's distribution's mean value
    const T m_sigma;         // Student's distribution's standard deviation
    T m_p;                   // the actual value of cdf will be subtracted by this value

public:

    /*
     * Constructor, assigns Student's distribution's degrees of freedom,
     * mean value and standard deviation that are used for calculation of cdf.
     * The subtrahend is automatically set to 0.
     *
     * @note The constructor assumes that the standard deviation is
     *       valid (strictly greater than 0) and does not check  this.
     *
     * @param df - Student's distribution degrees of freedom
     * @param mu - Student's distribution's mean value (default 0)
     * @param sigma - Student's distribution's standard deviation, must be greater than 0 (default: 1)
     */
    StudDistCdf(
                const T& df,
                const T& mu = static_cast<T>(0),
                const T& sigma = static_cast<T>(1) ) :
        m_df(df), m_mu (mu), m_sigma (sigma), m_p(static_cast<T>(0))
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
     * and previously assigned degrees of freedom, mean value and 
     * standard deviation. The cdf value is additionally subtracted
     * by the subtrahend 'm_p' that can be set via setP().
     *
     * @param x - value to be calculated the cumulative distribution function
     *
     * @return cdf(x, df, mu, sigma) - p
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
         * https://en.wikipedia.org/wiki/Student%27s_t-distribution#Cumulative_distribution_function
         * 't' can be defined as:
         * 
         *              df
         *  t = --------------------
         *            /  t - mu  \2
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
         */

        try
        {
            // Tolerance for the incomplete beta function
            const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) /
                          static_cast<T>(STAT_DIST_PROB_TOL_DEN);

            /*
             *            df
             *   t = --------------------
             *             /  x - mu  \2
             *        df + | -------- |
             *             \  sigma   /
             */

            const T ttemp = (x - this->m_mu) / this->m_sigma;

            // if 'x' is very close to 'mu', the probability is exactly 0.5
            if ( true == math::NumericUtil::isZero<T>(ttemp) )
            {
                return static_cast<T>(1) / static_cast<T>(2) - this->m_p;
            }

            const T t = this->m_df / ( this->m_df + ttemp*ttemp );

            /*
             * Evaluate:
             * 
             *   I  (df/2, 1/2) / 2
             *    t
             */
        
            T cdf = math::SpecFun::incBetaLowerReg<T>( 
                    this->m_df / static_cast<T>(2),
                    static_cast<T>(1) / static_cast<T>(2),
                    t,
                    TOL );
            cdf /= static_cast<T>(2);

            // Adjust the cdf if x>mu:
            if ( x > this->m_mu )
            {
                cdf = static_cast<T>(1) - cdf;
            }

            // Additionally m_p is subtracted from the cdf.
            return cdf - this->m_p;
        }
        catch ( const math::SpecFunException& spex )
        {
            throw math::FunctionException(math::FunctionException::UNDEFINED);
        }
    }

};  // class StudDistCdf


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
    /*
     *              b
     *              /
     *   P(a<x<b) = | pdf(t) dt
     *              /
     *              a
     */

    // sanity check
    math::StudentDist::__private::__checkParams<T>(sigma, df);
    
    const T from = std::min(a, b);
    const T to = std::max(a, b);

    return math::StudentDist::prob<T>(to, df, true, mu, sigma) -
           math::StudentDist::prob<T>(from, df, true, mu, sigma);
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
    /*
     * Lower tail probability ( P(t<x) ) is calculated inside the class
     * StudDistCdf. See comments inside its implementation for more details.
     *
     * For the upper tail probability ( P(t>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */
    
    // sanity check
    math::StudentDist::__private::__checkParams<T>(sigma, df);

    const math::StudentDist::__private::StudDistCdf<T> cdf(df, mu, sigma);

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
     * simple to evaluate, the Newton - Raphson root finding method is chosen that is
     * guaranteed to be reasonably precise and reasonably fast.
     */

    T retVal = static_cast<T>(0);

    try
    {
        const T P = (true==lowerTail ? p : static_cast<T>(1)-p);

        // Tolerance for the Newton - Raphson root finding method:
        const T TOL = static_cast<T>(STAT_DIST_PROB_TOL_NUM) / 
                      static_cast<T>(STAT_DIST_PROB_TOL_DEN);

        const math::StudentDist::__private::StudentDistPdf<T> pdf(df, mu, sigma);
        math::StudentDist::__private::StudDistCdf<T> cdf(df, mu, sigma);
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
