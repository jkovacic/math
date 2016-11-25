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
 * Implementation of functions within the namespace ContUniformDist
 * that perform various continuous uniform (a.k.a. rectangular) distribution
 * related operations, such as calculation of upper and lower tail probabilities
 * and quantiles, probability distribution function, etc.
 */


/*
 * More details about the continuous normal distribution:
 * https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
 */


// No #include "ContUniformDistGeneric.hpp" !!!
#include <algorithm>

#include "util/NumericUtil.hpp"
#include "exception/StatisticsException.hpp"

namespace math {  namespace ContUniformDist {  namespace __private
{

/*
 * Checks that both continuous uniform distribution's parameters
 * are not equal, which would cause division by zero.
 *
 * @param min - first continuous uniform distribution's parameter, typically its lower limit
 * @param max - second continuous uniform distribution's parameter, typically its upper limit
 *
 * @throw StatisticsException if the parameters are equal
 */
template <typename F>
void __checkParams(const F& min, const F& max) throw (math::StatisticsException)
{
    if ( true == math::NumericUtil::isZero<F>(min-max) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }
}

}}}  // namespace math::ContUniformDist::__private



/**
 * Value of the continuous uniform distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @note 'min' and 'max' are interchangeable.
 *
 * @param x - value to be calculated its probability distribution function
 * @param min - first continuous uniform distribution's parameter, typically its lower limit (default: 0)
 * @param max - second continuous uniform distribution's parameter, typically its upper limit (default: 1)
 *
 * @return pdf(x, min, max)
 *
 * @throw StatisticsException if the parameters are equal
 */
template <typename F>
F math::ContUniformDist::pdf(
      const F& x,
      const F& min,
      const F& max
    ) throw (math::StatisticsException)
{
    /*
     * Probability density function of the continuous uniform distribution:
     *
     *          /
     *          |      0       <==  x<min || x>max
     *          /
     *   pdf = { 
     *          \      1
     *          |  ---------   <==  min <= x <= max
     *          \   max-min
     */

    // sanity check:
    math::ContUniformDist::__private::__checkParams<F>(min, max);

    const F a = std::min(min, max);
    const F b = std::max(min, max);

    return ( x<a || x>b ? static_cast<F>(0) : static_cast<F>(1)/(b-a) );
}


/**
 * Probability in the continuous uniform distribution with specified
 * parameters  that 'x' is greater than 'a' and less  than 'b'
 * (or vice versa if b<a).
 *
 * @note 'a' and 'b' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @note 'min' and 'max' are interchangeable.
 *
 * @param a - lower limit of the interval
 * @param b - upper limit of the interval
 * @param min - first continuous uniform distribution's parameter, typically its lower limit (default: 0)
 * @param max - second continuous uniform distribution's parameter, typically its upper limit (default: 1)
 *
 * @return P(a<X<b) or P(b<X<a)
 *
 * @throw StatisticsException if the parameters are equal
 */
template <typename F>
F math::ContUniformDist::probInt(
      const F& a,
      const F& b,
      const F& min,
      const F& max
    ) throw (math::StatisticsException)
{
    /*
     *              b
     *              /
     *   P(a<X<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *              a
     */

    // sanity check will be performed by prob()

    const F from = std::min(a, b);
    const F to = std::max(a, b);

    return math::ContUniformDist::prob<F>(to, min, max, true) -
           math::ContUniformDist::prob<F>(from, min, max, true);
}


/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified continuous uniform distribution.
 *
 * @note 'min' and 'max' are interchangeable.
 *
 * @param x - value
 * @param min - first continuous uniform distribution's parameter, typically its lower limit (default: 0)
 * @param max - second continuous uniform distribution's parameter, typically its upper limit (default: 1)
 * @param lowerTail - if true, returns P(t<x), otherwise P(t>x) (default: true)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if the parameters are equal
 */
template <typename F>
F math::ContUniformDist::prob(
      const F& x,
      const F& min,
      const F& max,
      const bool lowerTail
    ) throw (math::StatisticsException)
{
    /*
     * The cdf represents probability that a value from the normal
     * distribution is smaller than a specified value 'x' and can be
     * expressed as:
     *
     *                     x
     *                     /
     *   cdf(x) = P(X<x) = | pdf(t) dt
     *                     /
     *                   -inf
     *
     * For the continuous uniform distribution cdf
     * can be calculated as:      
     *
     *          /     0      <==  x < min
     *          |
     *          /     1      <==  x > max
     *   cdf = { 
     *          \   x-min
     *          | ---------  <==  min <= x <= max
     *          \  max-min
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::ContUniformDist::__private::__checkParams<F>(min, max);

    const F a = std::min(min, max);
    const F b = std::max(min, max);

    F cdf = static_cast<F>(0);
    
    if ( x < a )
    {
        cdf = static_cast<F>(0);
    }
    else if ( x > b )
    {
        cdf = static_cast<F>(1);
    }
    else
    {
        cdf = (x - a) / (b - a);
    }

    // the return value depends on 'lowerTail':
    return ( true==lowerTail ? cdf : static_cast<F>(1)-cdf );
}


/**
 * Quantile function for the specified continuous uniform distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @note 'min' and 'max' are interchangeable.
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param min - first continuous uniform distribution's parameter, typically its lower limit (default: 0)
 * @param max - second continuous uniform distribution's parameter, typically its upper limit (default: 1)
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if the parameters are equal
 */
template <typename F>
F math::ContUniformDist::quant(
      const F& p,
      const F& min,
      const F& max,
      const bool lowerTail
    ) throw (math::StatisticsException)
{
    // sanity check
    math::ContUniformDist::__private::__checkParams<F>(min, max);
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
     * When 0 < p < 1,
     * x will be a value between 'min' and 'max':
     *
     *   x = a * (1-p) + p*b
     */


    // The actual probability in the expressions above depends on 'lowerTail':
    const F P = (true==lowerTail ? p : static_cast<F>(1)-p);

    const F a = std::min(min, max);
    const F b = std::max(min, max);

    return ( a * (static_cast<F>(1) - P) + P * b);
}
