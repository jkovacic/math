/*
Copyright 2016, Jernej Kovacic

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
 * @headername{TriangularDist.h}
 *
 * Declaration of functions within the namespace TriangularDist
 * that perform various triangular distribution related operations,
 * such as calculation of upper and lower tail probabilities
 * and quantiles, probability distribution function, etc.
 */


/*
 * More details about the triangular distribution:
 * https://en.wikipedia.org/wiki/Triangular_distribution
 */


#ifndef _MATH_TRIANGULARDISTGENERIC_HPP_
#define _MATH_TRIANGULARDISTGENERIC_HPP_


#include <algorithm>
#include <cmath>

#include "util/NumericUtil.hpp"
#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform various triangular
 * distribution related operations, such as calculation of upper and
 * lower tail probabilities and quantiles, probability distribution
 * function, etc.
 */
namespace TriangularDist
{


namespace __private
{

/*
 * Checks if all triangular distribution's parameter satisfy
 * the criteria: a < c < b
 *
 * @param a - triangular distribution's lower limit
 * @param b - triangular distribution's upper limit
 * @param c - triangular distribution's mode
 *
 * @throw StatisticsException if the parameters are invalid
 */
template <typename F>
void __checkParams(const F& a, const F& b, const F& c)
{
    if ( ( (c-a) < math::NumericUtil::getEPS<F>() ) ||
           (b-c) < math::NumericUtil::getEPS<F>() )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }
}

}  // namespace math::TriangularDist::__private



/**
 * Value of the triangular distribution's probability distribution
 * function (pdf) for the given 'x'.
 *
 * @note 'a', 'b' and 'c' must satisfy the criteria: a < c < b
 *
 * @param x - value to be calculated its probability distribution function
 * @param a - triangular distribution's lower limit
 * @param b - triangular distribution's upper limit
 * @param c - triangular distribution's mode
 *
 * @return pdf(x, a, b, c)
 *
 * @throw StatisticsException if the parameters 'a', 'b' and 'c' are invalid
 */
template <typename F>
F pdf(
	  const F& x,
	  const F& a,
	  const F& b,
	  const F& c
	)
{
    /*
     * Probability density function of the triangular distribution:
     *
     *          /
     *          |       0          <==  x <= a
     *          |
     *          |     2*(x-a)
     *          | -------------    <==  a < x <= c
     *          /  (b-a)*(c-a)
     *   pdf = {
     *          \     2*(b-x)
     *          |  -------------   <==  c < x < b
     *          |   (b-a)*(b-c)
     *          |
     *          |       0          <==  x >= b
     *          \
     */

    // sanity check:
    math::TriangularDist::__private::__checkParams<F>(a, b, c);

    F pdf = static_cast<F>(0);

    if ( a<x && x<=c )
    {
        pdf = 2*(x-a) / ((b-a)*(c-a));
    }
    else if ( c<x && x<b )
    {
        pdf = 2*(b-x) / ((b-a)*(b-c));
    }
    // else pdf = static_cast<F>(0);

    return pdf;
}



/**
 * Probability that a value is greater or smaller (depending on 'lowerTail')
 * than the given value 'x' for the specified triangular distribution.
 *
 * @note 'a', 'b' and 'c' must satisfy the criteria: a < c < b
 *
 * @param x - value
 * @param a - triangular distribution's lower limit
 * @param b - triangular distribution's upper limit
 * @param c - triangular distribution's mode
 * @param lowerTail - if true, returns P(t<x), otherwise P(t>x) (default: true)
 *
 * @return P(X<x) if 'lowerTail' equals true, P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if the parameters 'a', 'b' and 'c' are invalid
 */
template <typename F>
F prob(
	  const F& x,
	  const F& a,
	  const F& b,
	  const F& c,
	  const bool lowerTail = true
	)
{
    /*
     * The cdf represents probability that a value from the triangular
     * distribution is smaller than a specified value 'x' and can be
     * expressed as:
     *
     *                     x
     *                     /
     *   cdf(x) = P(X<x) = | pdf(t) dt
     *                     /
     *                   -inf
     *
     * For the triangular distribution, cdf can be calculated as:
     *
     *          /
     *          |        0            <==   x <= a
     *          |
     *          |           2
     *          |      (x-a)
     *          |  -------------      <==   a < x <= c
     *          |   (b-a)*(c-a)
     *          /
     *   cdf = {
     *          \
     *          |              2
     *          |         (b-x)
     *          | 1 - -------------   <==  c < x < b
     *          |      (b-a)*(b-c)
     *          |
     *          |        1            <==  x >= b
     *          \
     *
     * For the upper tail probability ( P(X>x) ) just return the complement
     * of the lower tail probability: 1 - cdf(x)
     */

    // sanity check
    math::TriangularDist::__private::__checkParams<F>(a, b, c);

    F cdf = static_cast<F>(0);

    if ( x <= a )
    {
        cdf = static_cast<F>(0);
    }
    else if ( a<x && x<=c )
    {
        const F xa = x - a;

        cdf = xa*xa / ((b-a)*(c-a));
    }
    else if ( c<x && x<b )
    {
        const F bx = b - x;

        cdf = static_cast<F>(1) - ( bx*bx / ((b-a)*(b-c)) );
    }
    else  // x >= b
    {
        cdf = static_cast<F>(1);
    }

    // the return value depends on 'lowerTail':
    return ( true==lowerTail ? cdf : static_cast<F>(1)-cdf );
}



/**
 * Probability in the triangular distribution with specified
 * parameters  that 'x' is greater than 'from' and less  than 'to'.
 *
 * @note 'from' and 'to' are interchangeable, the result's sign will
 *       not change in this case.
 *
 * @note 'a', 'b' and 'c' must satisfy the criteria: a < c < b
 *
 * @param from - lower limit of the interval
 * @param to - upper limit of the interval
 * @param a - triangular distribution's lower limit
 * @param b - triangular distribution's upper limit
 * @param c - triangular distribution's mode
 *
 * @return P(from<X<to) or P(to<X<from)
 *
 * @throw StatisticsException if the parameters 'a', 'b' and 'c' are invalid
 */
template <typename F>
F probInt(
	  const F& from,
	  const F& to,
	  const F& a,
	  const F& b,
	  const F& c
	)
{
    /*
     *              to
     *              /
     *   P(a<X<b) = | pdf(t) dt  =  cdf(b) - cdf(a)
     *              /
     *            from
     */

    // sanity check will be performed by prob()

    const F upper = std::max<F>(from, to);
    const F lower = std::min<F>(from, to);

    return math::TriangularDist::prob<F>(upper, a, b, c, true) -
           math::TriangularDist::prob<F>(lower, a, b, c, true);
}



/**
 * Quantile function for the specified triangular distribution.
 *
 * If 'lowerTail' equals true, it returns such 'x' that P(X<x) = p
 * If 'lowerTail' equals false, it returns such 'x' that P(X>x) = p
 *
 * @note The function always returns a value between 'a' and 'b'.
 *
 * @note 'a', 'b' and 'c' must satisfy the criteria: a < c < b
 *
 * @param p - probability (must be greater than 0 and less than 1)
 * @param a - triangular distribution's lower limit
 * @param b - triangular distribution's upper limit
 * @param c - triangular distribution's mode
 * @param lowerTail - if true, returns q(X<x), otherwise q(X>x) (default: true)
 *
 * @return x: P(X<x) if 'lowerTail' equals true, x: P(X>x) if 'lowerTail' equals false
 *
 * @throw StatisticsException if parameters 'a', 'b', 'c' or 'p' are invalid
 */
template <typename F>
F quant(
	  const F& p,
	  const F& a,
	  const F& b,
	  const F& c,
	  const bool lowerTail = true
	)
{
    // sanity check
    math::TriangularDist::__private::__checkParams<F>(a, b, c);
    if ( p < static_cast<F>(0) ||
         p > static_cast<F>(1) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    // The actual probability in the expressions above depends on 'lowerTail':
    const F P = (true==lowerTail ? p : static_cast<F>(1)-p);


    /*
     * Quantile is a value of 'x' that solves the nonlinear equation:
     *
     *   cdf(x) = P
     *
     * When P ~= 0, any x<=a is valid, the function will return 'a'.
     *
     * When P ~= 1, any x>=b is valid, the function will return 'b'.
     *
     * When 0 < P < 1, x will be a value between 'a' and 'b' as shown below:
     *
     * At the mode, the CDF is equal to:
     *
     *             c - a
     *   CDF(c) = -------
     *             b - a
     *
     * Further evaluation of the quantile depends whether P is less than or
     * greater than CDF(c):
     *
     *                /
     *                |             a                     <== P ~= 0
     *                |
     *                |
     *                |     -+  +-----------------+
     *                | a +   \/ P * (b-a) * (c-a)        <== 0 < P <= (c-a)/(b-a)
     *                /
     *    quant(P) = {
     *                \
     *                |     -+  +---------------------+
     *                | b -   \/ (1-P) * (b-a) * (b-c)    <== (c-a)/(b-a) < P < 1
     *                |
     *                |
     *                |             b                     <== P ~= 1
     *                \
     */

    F retVal = static_cast<F>(0);

    if ( P < math::NumericUtil::getEPS<F>() )
    {
        retVal = a;
    }
    else if ( P > (static_cast<F>(1)-math::NumericUtil::getEPS<F>()) )
    {
        retVal = b;
    }
    else if ( P <= ( (c-a) / (b-a) ) )
    {
        retVal = a + std::sqrt( P * (b-a) * (c-a) );
    }
    else  // (c-a)/(b-a) < P < 1
    {
        retVal = b - std::sqrt( (static_cast<F>(1)-P) * (b-a) * (b-c) );
    }

    return retVal;
}


}  // namespace TriangularDist

}  // namespace math

#endif  // _MATH_TRIANGULARDISTGENERIC_HPP_
