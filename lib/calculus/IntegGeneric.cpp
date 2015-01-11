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
 * Implementation of the namespace Integ with functions that
 * perform numerical integration of continuous functions.
 */


// No #include "IntegGeneric.hpp" !!!!
#include <cstddef>
#include <cmath>
#include <algorithm>

#include "util/NumericUtil.hpp"
#include "exception/FunctionException.hpp"

#include "../settings/omp_settings.h"
#include "omp/omp_header.h"
#include "omp/omp_coarse.h"


// Implementation of "private" functions

namespace math
{

namespace Integ
{

namespace __private
{

/*
 * Numerical integration using the rectangle method.
 *
 * For more info about the method, see:
 * https://en.wikipedia.org/wiki/Rectangle_method
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating steps
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __rectangle(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& b,
        size_t n
      ) throw(math::CalculusException)
{
    /*
     *
     *     b                     N-1
     *     /                    -----
     *     |                    \
     *     | f(x) dx   ~=   h *  >  f(a + i*h)
     *     |                    /
     *     /                    -----
     *     a                     i=0
     *
     * where h = (b - a) / N
     */

    try
    {
        // The algorithm requires evaluation of the function in
        // N points, the same number as integrating intervals

        const size_t& N = n;
        const F h = (b-a) / static_cast<F>(N);

        F sum = static_cast<F>(0);

        // Coarse grained parallelism
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(f, a, N) \
                    reduction(+ : sum)
        {
            const size_t thnr = omp_get_thread_num();
            const size_t nthreads  = omp_get_num_threads();
            const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
            const size_t istart = elems_per_thread * thnr;
            const size_t iend = std::min(istart + elems_per_thread, N);

            F tempSum = static_cast<F>(0);
            F xi = a + static_cast<F>(istart) * h;
            for (size_t i = istart;
                 i < iend;
                 ++i, xi += h )
            {
                tempSum += f.func( xi );
            }

            // update 'sum' in a thread safe manner
            sum += tempSum;
        }  // omp parallel

        return sum * h;
    } // try
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


/*
 * General function that implements closed Newton - Cotes formulae of all
 * degrees (e.g. trapezoidal, Simpson's, Boole's rule, etc.) to obtain proper
 * definite integrals.
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating steps
 * @param degree - degree of the method
 * @param coef - array of coefficients to multiply non-boundary points (must be of size 'degree'!)
 * @param bCoef - coefficient to multiply both boundary points
 * @param hCoef - coefficient to multiply the step size to obtain the final result
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __closedNewtonCotes(
               const math::IFunctionGeneric<F>& f,
               const F& a,
               const F& b,
               size_t n,
               size_t degree,
               const F* coef,
               const F& bCoef,
               const F& hCoef
             ) throw(math::CalculusException)
{
    /*
     * More info about Newton - Cotes formulae:
     * - https://en.wikipedia.org/wiki/Newton–Cotes_formulas
     * - http://mathworld.wolfram.com/Newton-CotesFormulas.html
     */

    try
    {
        /*
         * With N integrating intervals, the function must be evaluated
         * in N+1 points. Two points ('a' and 'b') are handled separately,
         * the remaining N-1 points are processed by a for loop.
         */

        // N must be divisible by 'degree' !
        const size_t N = n + ( 0!=n%degree ? degree-n%degree: 0 );
        const F h = (b-a) / static_cast<F>(N);

        // Do not preinitialize sum to hCoef * (f(a) + f(b)) now as it may
        // be reset back to 0 by the reduction clause
        F sum = static_cast<F>(0);

        // Coarse grained parallelism
        #pragma omp parallel num_threads(ompIdeal(N-1)) \
                    if((N-1)>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(f, a, coef, degree) \
                    reduction(+ : sum)
        {
            const size_t thnr = omp_get_thread_num();
            const size_t nthreads  = omp_get_num_threads();
            const size_t elems_per_thread = ((N-1) + nthreads - 1) / nthreads;
            const size_t istart = elems_per_thread * thnr + 1;
            const size_t iend = std::min(istart + elems_per_thread, N);

            F tempSum = static_cast<F>(0);
            F xi = a + static_cast<F>(istart) * h;
            for ( size_t i=istart;
                  i < iend;
                  ++i, xi += h )
            {
                tempSum += coef[i%degree] * f.func(xi);
            }

            // update sum in a thread safe manner
            sum += tempSum;
        }  // omp parallel

        // finally add the remaining two points (at 'a' and 'b'):
        sum += bCoef * ( f.func(a) + f.func(b) );

        return sum * hCoef *  h;
    }  // try
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


/*
 * Numerical integration using the trapezoidal rule.
 *
 * For more info about the method, see:
 * https://en.wikipedia.org/wiki/Trapezoidal_rule
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating steps
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __trapezoidal(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& b,
        size_t n
      ) throw(math::CalculusException)
{
    /*
     *
     *     b                      N-1
     *     /                     -----
     *     |               h     \   /                         \
     *     | f(x) dx   ~= ---  *  >  | f(a+i*h) + f(a+(i+1)*h) |  =
     *     |               2     /   \                         /
     *     /                     -----
     *     a                      i=0
     *
     *
     *         /                   N-2          \
     *         |                  -----         |
     *     h   |                  \             |
     *  = ---  | f(a) + f(b) + 2   >  f(a+i*h)  |  =
     *     2   |                  /             |
     *         |                  -----         |
     *         \                   i=1          /
     *
     *
     *        /                      N-2          \
     *        |                     -----         |
     *        |   f(a) + f(b)       \             |
     *  =  h  |  -------------  +    >  f(a+i*h)  |
     *        |        2            /             |
     *        |                     -----         |
     *        \                      i=1          /
     *
     *
     * where h = (b - a) / N
     */

    const F c[1] = { static_cast<F>(1) };

    return __closedNewtonCotes<F>(
            f, a, b, n, 1, c, static_cast<F>(1)/static_cast<F>(2), static_cast<F>(1)
        );
}


/*
 * Numerical integration using the Simpson's rule.
 *
 * For more info about the method, see:
 * https://en.wikipedia.org/wiki/Simpson%27s_rule
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating steps
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __simpson(
       const math::IFunctionGeneric<F>& f,
       const F& a,
       const F& b,
       size_t n
     ) throw(math::CalculusException)
{
    /*
     *   b
     *   /                  /                                                             \
     *   |               h  |                                                             |
     *   | f(x) dx  ~=  --- | f(x0) + 4*f(x1) + 2*f(x2) + 4*f(x3) + 4*f(x4) + ... + f(xN) |
     *   |               3  |                                                             |
     *   /                  \                                                             /
     *   a
     *
     *  where h = (b-1) / N,  xi = a + i *h,
     *  and N must be an even number (divisible by 2)
     */

    const F c[2] = { static_cast<F>(2),
                     static_cast<F>(4) };

    return math::Integ::__private::__closedNewtonCotes<F>(
            f, a, b, n, 2, c, static_cast<F>(1), static_cast<F>(1)/static_cast<F>(3)
        );
}


/*
 * Numerical integration using the Simpson's 3/8 rule.
 *
 * For more info about the method, see:
 * https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson.27s_3.2F8_rule
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating steps
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __simpson38(
      const math::IFunctionGeneric<F>& f,
      const F& a,
      const F& b,
      size_t n
    ) throw(math::CalculusException)
{
    /*
     *   b
     *   /                    /                                                                                 \
     *   |               3*h  |                                                                                 |
     *   | f(x) dx  ~=  ----- | f(x0) + 3*f(x1) + 3*f(x2) + 2*f(x3) + 3*f(x4) + 3*f(x5) + 2*f(x6) + ... + f(xN) |
     *   |                8   |                                                                                 |
     *   /                    \                                                                                 /
     *   a
     *
     *  where h = (b-1) / N,  xi = a + i *h,
     *  and N must be divisible by 3
     */

    const F c[3] = { static_cast<F>(2),
                     static_cast<F>(3),
                     static_cast<F>(3) };

    return math::Integ::__private::__closedNewtonCotes<F>(
            f, a, b, n, 3, c, static_cast<F>(1), static_cast<F>(3)/static_cast<F>(8)
        );
}


/*
 * Numerical integration using the Boole's rule.
 *
 * For more info about the method, see:
 * https://en.wikipedia.org/wiki/Boole%27s_rule
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating steps
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __boole(
      const math::IFunctionGeneric<F>& f,
      const F& a,
      const F& b,
      size_t n
    ) throw(math::CalculusException)
{
    /*
     *   b
     *   /                    /
     *   |               2*h  |
     *   | f(x) dx  ~=  ----- | 7*f(x0) + 32*f(x1) + 12*f(x2) + 32*f(x3) + 14*f(x4) +
     *   |               45   |
     *   /                    \
     *   a
     *
     *                                                                           \
     *                                                                           |
     *               + 32*f(x5) + 12*f(x6) + 32*f(x7) + 14*f(x8) + ... + 7*f(xN) |
     *                                                                           |
     *                                                                           /
     * 
     *  where h = (b-1) / N,  xi = a + i *h,
     *  and N must be divisible by 3
     */

    const F c[4] = { static_cast<F>(14),
                     static_cast<F>(32),
                     static_cast<F>(12),
                     static_cast<F>(32) };

    return math::Integ::__private::__closedNewtonCotes<F>(
            f, a, b, n, 4, c, static_cast<F>(7), static_cast<F>(2)/static_cast<F>(45)
        );
}


}  // namespace __private
}  // namepsace Integ
}  // namespace math



/**
 * Performs numerical integration using the selected algorithm.
 *
 * 'a' and 'b' can be interchanged, in this case the opposite value
 * will be returned.
 *
 * @note The actual number of integrating intervals may be slightly increased
 *       if required by the selected algorithm.
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param n - desired number of integrating intervals (default: 10000)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between 'a' and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integ(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& b,
        size_t n,
        math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    // sanity check:
    if ( n <= 0 )
    {
        throw math::CalculusException(math::CalculusException::NOT_ENOUGH_STEPS);
    }

    // If both boundaries are the same, the result is equal to 0
    if ( true == math::NumericUtil::isZero<F>(a-b) )
    {
        return static_cast<F>(0);
    }

    /*
     * Individual integration algorithm functions expect 'b' to be
     * greater than 'a'. If this is not the case, swap the boundaries
     * and revert the result's sign.
     */
    F retVal = static_cast<F>(1);

    const F from = std::min(a, b);
    const F to   = std::max(a, b);
    if ( a > b )
    {
        retVal = -retVal;
    }

    switch(algorithm)
    {
        case math::EIntegAlg::RECTANGLE :
        {
            retVal *= math::Integ::__private::__rectangle<F>(f, from, to, n);
            break;
        }

        case math::EIntegAlg::TRAPEZOIDAL :
        {
            retVal *= math::Integ::__private::__trapezoidal<F>(f, from, to, n);
            break;
        }

        case math::EIntegAlg::SIMPSON :
        {
            retVal *= math::Integ::__private::__simpson<F>(f, from, to , n);
            break;
        }

        case math::EIntegAlg::SIMPSON_3_8 :
        {
            retVal *= math::Integ::__private::__simpson38<F>(f, from, to, n);
            break;
        }

        case math::EIntegAlg::BOOLE :
        {
            retVal *= math::Integ::__private::__boole<F>(f, from, to, n);
            break;
        }

        default :
            throw math::CalculusException(math::CalculusException::UNSUPPORTED_ALGORITHM);
    }  // switch

    return retVal;
}



/**
 * Performs numerical integration using the selected algorithm.
 *
 * 'a' and 'b' can be interchanged, in this case the opposite value
 * will be returned.
 *
 * @note The actual size of integrating step may be slightly decreased
 *       if required by the selected algorithm.
 *
 * @param f - instance of a class with the function to integrate
 * @param a - lower bound of the integration interval
 * @param b - upper bound of the integration interval
 * @param h - desired size of an integrating step (default: 0.0001)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between 'a' and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integH(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& b,
        const F& h,
        math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    if ( h < math::NumericUtil::getEPS<F>() )
    {
        throw math::CalculusException(math::CalculusException::INVALID_STEP);
    }

    // Just obtain the (approximate) number of integrating intervals
    // from the size
    const size_t n = static_cast<size_t>(std::floor(std::abs(b-a) / h) +
    		static_cast<F>(1) / static_cast<F>(2) ) + 1;

    return math::Integ::integ<F>( f, a, b, n, algorithm );
}
