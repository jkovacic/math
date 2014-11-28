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
 * Implementation of the class IntegGeneric and related classes that
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
 * @param n - desired number of integrating intervals
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between 'a' and 'b'
 *
 * @see IFunctionGeneric
 */
template <class T>
T math::IntegGeneric<T>::integ(
        const math::IFunctionGeneric<T>& f,
        const T& a,
        const T& b,
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
    if ( true == math::NumericUtil::isZero<T>(a-b) )
    {
        return static_cast<T>(0);
    }

    /*
     * Individual integration algorithm functions expect 'b' to be
     * greater than 'a'. If this is not the case, swap the boundaries
     * and revert the result's sign.
     */
    T retVal = static_cast<T>(1);

    const T from = std::min(a, b);
    const T to   = std::max(a, b);
    if ( a > b )
    {
        retVal = -retVal;
    }

    switch(algorithm)
    {
        case math::EIntegAlg::RECTANGLE :
        {
            retVal *= __rectangle(f, from, to, n);
            break;
        }

        case math::EIntegAlg::TRAPEZOIDAL :
        {
            retVal *= __trapezoidal(f, from, to, n);
            break;
        }

        case math::EIntegAlg::SIMPSON :
        {
            retVal *= __simpson(f, from, to , n);
            break;
        }

        case math::EIntegAlg::SIMPSON_3_8 :
        {
            retVal *= __simpson38(f, from, to, n);
            break;
        }

        case math::EIntegAlg::BOOLE :
        {
            retVal *= __boole(f, from, to, n);
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
 * @param h - desired size of an integrating step
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f.func() between 'a' and 'b'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between 'a' and 'b'
 *
 * @see IFunctionGeneric
 */
template <class T>
T math::IntegGeneric<T>::integH(
        const math::IFunctionGeneric<T>& f,
        const T& a,
        const T& b,
        const T& h,
        math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    if ( h < math::NumericUtil::getEPS<T>() )
    {
        throw math::CalculusException(math::CalculusException::INVALID_STEP);
    }

    // Just obtain the (approximate) number of integrating intervals
    // from the size
    const size_t n = static_cast<size_t>(std::floor(std::abs(b-a) / h) +
    		static_cast<T>(1) / static_cast<T>(2) ) + 1;

    return integ( f, a, b, n, algorithm );
}


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
template <class T>
T math::IntegGeneric<T>::__rectangle(
        const math::IFunctionGeneric<T>& f,
        const T& a,
        const T& b,
        size_t n
      ) throw(math::CalculusException)
{
    /*
     *
     *     b                    N-1
     *     /                   -----
     *     |                   \
     *     | f(x) dx   ~   h *  >  f(a + i*h)
     *     |                   /
     *    /                    -----
     *    a                     i=0
     *
     * where h = (b - a) / N
     */

    try
    {
        // The algorithm requires evaluation of the function in
        // N points, the same number as integrating intervals

        const size_t& N = n;
        const T h = (b-a) / static_cast<T>(N);

        T sum = static_cast<T>(0);

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

            T tempSum = static_cast<T>(0);
            T xi = a + static_cast<T>(istart) * h;
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
template <class T>
T math::IntegGeneric<T>::__trapezoidal(
        const math::IFunctionGeneric<T>& f,
        const T& a,
        const T& b,
        size_t n
      ) throw(math::CalculusException)
{
    /*
     *
     *     b                     N-1
     *     /                    -----
     *     |              h     \   /                         \
     *     | f(x) dx   ~ ---  *  >  | f(a+i*h) + f(a+(i+1)*h) |  =
     *     |              2     /   \                         /
     *    /                     -----
     *    a                      i=0
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

    const T c[1] = { static_cast<T>(1) };

    return __closedNewtonCotes(
            f, a, b, n, 1, c, static_cast<T>(1)/static_cast<T>(2), static_cast<T>(1)
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
template <class T>
T math::IntegGeneric<T>::__simpson(
       const math::IFunctionGeneric<T>& f,
       const T& a,
       const T& b,
       size_t n
     ) throw(math::CalculusException)
{
    /*
     *   b
     *   /                 /                                                             \
     *   |              h  |                                                             |
     *   | f(x) dx  ~  --- | f(x0) + 4*f(x1) + 2*f(x2) + 4*f(x3) + 4*f(x4) + ... + f(xN) |
     *   |              3  |                                                             |
     *  /                  \                                                             /
     *  a
     *
     *  where h = (b-1) / N,  xi = a + i *h,
     *  and N must be an even number (divisible by 2)
     */

    const T c[2] = { static_cast<T>(2),
                     static_cast<T>(4) };

    return __closedNewtonCotes(
            f, a, b, n, 2, c, static_cast<T>(1), static_cast<T>(1)/static_cast<T>(3)
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
template <class T>
T math::IntegGeneric<T>::__simpson38(
      const math::IFunctionGeneric<T>& f,
      const T& a,
      const T& b,
      size_t n
    ) throw(math::CalculusException)
{
    /*
     *   b
     *   /                   /                                                                                 \
     *   |              3*h  |                                                                                 |
     *   | f(x) dx  ~  ----- | f(x0) + 3*f(x1) + 3*f(x2) + 2*f(x3) + 3*f(x4) + 3*f(x5) + 2*f(x6) + ... + f(xN) |
     *   |               8   |                                                                                 |
     *  /                    \                                                                                 /
     *  a
     *
     *  where h = (b-1) / N,  xi = a + i *h,
     *  and N must be divisible by 3
     */

    const T c[3] = { static_cast<T>(2),
                     static_cast<T>(3),
                     static_cast<T>(3) };

    return __closedNewtonCotes(
            f, a, b, n, 3, c, static_cast<T>(1), static_cast<T>(3)/static_cast<T>(8)
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
template <class T>
T math::IntegGeneric<T>::__boole(
      const math::IFunctionGeneric<T>& f,
      const T& a,
      const T& b,
      size_t n
    ) throw(math::CalculusException)
{
    /*
     *   b
     *   /                   /
     *   |              2*h  |
     *   | f(x) dx  ~  ----- | 7*f(x0) + 32*f(x1) + 12*f(x2) + 32*f(x3) + 14*f(x4) +
     *   |              45   |
     *  /                    \
     *  a
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

    const T c[4] = { static_cast<T>(14),
                     static_cast<T>(32),
                     static_cast<T>(12),
                     static_cast<T>(32) };

    return __closedNewtonCotes(
            f, a, b, n, 4, c, static_cast<T>(7), static_cast<T>(2)/static_cast<T>(45)
        );
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
template <class T>
T math::IntegGeneric<T>::__closedNewtonCotes(
               const math::IFunctionGeneric<T>& f,
               const T& a,
               const T& b,
               size_t n,
               size_t degree,
               const T* coef,
               const T& bCoef,
               const T& hCoef
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
        const T h = (b-a) / static_cast<T>(N);
        
        // Do not preinitialize sum to hCoef * (f(a) + f(b)) now as it may
        // be reset back to 0 by the reduction clause
        T sum = static_cast<T>(0);

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

            T tempSum = static_cast<T>(0);
            T xi = a + static_cast<T>(istart) * h;
            for ( size_t i=istart;
                  i < iend;
                  ++i,  xi += h )
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
