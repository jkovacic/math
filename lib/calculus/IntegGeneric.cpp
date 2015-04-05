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
 * Checks validity of the given integration step.
 * It must be positive and not less than EPS for the given
 * type F.
 *
 * If 'h' is not valid, CalculusException is thrown.
 *
 * @param h - step to check
 *
 * @throw CalculusException if 'h' is invalid
 */
template <typename F>
inline void __checkStep(const F& h) throw(math::CalculusException)
{
    if ( h < math::NumericUtil::getEPS<F>() )
    {
        throw math::CalculusException(math::CalculusException::INVALID_STEP);
    }
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
 * @return definite integral of f(x) between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __rectangle(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& b,
        const size_t n
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
                tempSum += f( xi );
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
 * @return definite integral of f(x) between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __closedNewtonCotes(
               const math::IFunctionGeneric<F>& f,
               const F& a,
               const F& b,
               const size_t n,
               const size_t degree,
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
                tempSum += coef[i%degree] * f(xi);
            }

            // update sum in a thread safe manner
            sum += tempSum;
        }  // omp parallel

        // finally add the remaining two points (at 'a' and 'b'):
        sum += bCoef * ( f(a) + f(b) );

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
 * @return definite integral of f(x) between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __trapezoidal(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& b,
        const size_t n
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
 * @return definite integral of f(x) between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __simpson(
       const math::IFunctionGeneric<F>& f,
       const F& a,
       const F& b,
       const size_t n
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
 * @return definite integral of f(x) between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __simpson38(
      const math::IFunctionGeneric<F>& f,
      const F& a,
      const F& b,
      const size_t n
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
 * @return definite integral of f(x) between 'a' and 'b'
 *
 * @throw CalculusException if function is not defined between 'a' and 'b'
 */
template <typename F>
F __boole(
      const math::IFunctionGeneric<F>& f,
      const F& a,
      const F& b,
      const size_t n
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


/*
 * An interface (i.e. a pure virtual class) for classes that
 * perform integration by substitution.
 *
 * The derived classes should take a function class (an instance of
 * math::IFunctionGeneric) and provide all necessary functions for
 * integration by substitution. The algorithm is derived from the
 * fundamental theorem of calculus:
 *
 *                -1
 *     b       phi  (b)
 *     /          /
 *     |f(x) dx = | f(phi(t)) * phi'(t) dt
 *     /          /
 *     a          -1
 *              phi (a)
 *
 * More details about integration by substitution:
 * https://en.wikipedia.org/wiki/Integration_by_substitution
 */
template <typename F>
class IIntegSubst : public math::IFunctionGeneric<F>
{
public:

    /*
     * Inverse of phi(x) a.k.a. phi^(-1)(x).
     * Typically used to obtain new integration limits.
     *
     * @param x
     *
     * @return phi^(-1)(x)
     *
     * @throw FunctionException if inverse of phi is invalid for the given 'x'
     */
    virtual F inv(const F& x) const throw(math::FunctionException) = 0;

    /*
     * Inverse of phi(x) for x = -inf
     * In most cases the exact return value should be 0 however
     * often -eps should be returned instead.
     *
     * @return phi^(-1)(-inf)
     *
     * @throw FunctionExceptionif phi^(-1) is undefined at x = -inf
     */
    virtual F invNegInf() const throw(math::FunctionException) = 0;

    /*
     * Inverse of phi(x) for x = +inf
     * In most cases the exact return value should be 0 however
     * often +eps should be returned instead.
     *
     * @return phi^(-1)(inf)
     *
     * @throw FunctionExceptionif phi^(-1) is undefined at x = inf
     */
    virtual F invPosInf() const throw(math::FunctionException) = 0;

    /*
     * Returns TRUE if the following condition is satisfied:
     *
     *                   -1         -1
     *     a < b  ==> phi  (a) > phi  (b)
     *
     * Typically this indicates that integration limits should be swapped
     * when integrating by substitution.
     *
     * @return true or false.
     */
    virtual bool swappedLimits() const = 0;

    // inherited from IFunctionGeneric:
    // virtual F operator()(const F& x) const throw(math::FunctionException)

    /*
     * Destructor
     */
    virtual ~IIntegSubst()
    {
        // nothing to do
    }
};


/*
 * Implementation of IIntegSubst for phi(x) = 1/x,
 * i.e reciprocal value of 'x'
 */
template <typename F>
class IntegSubstRec : public IIntegSubst<F>
{

private:
    const math::IFunctionGeneric<F>& f;
    const F EPS;

public:

    /*
     * Constructor
     *
     * @param fg - an instance of IFunctionGeneric with the function to integrate
     * @param eps - eps value for inv at both infinities (default: system EPS for the given type F)
     */
    IntegSubstRec(
        const math::IFunctionGeneric<F>& fg,
        const F& eps = math::NumericUtil::getEPS<F>() ) :
            f(fg), EPS(eps)
    {
        /*
         * TODO validity of 'eps' should be checked here.
         * However this class is private within this file and
         * it is assumed that a valid 'eps' will always be passed
         * by functions that utilize it.
         */
    }

    /*
     * Implementation of operator() that returns the value of the
     * substituted function.
     *
     * @param x - input argument
     *
     * @return value of the substituted function for the given 'x'
     *
     * @throw FunctionException if the substituted function is undefined at 'x'
     */
    F operator()(const F& x) const throw(math::FunctionException)
    {
        if ( true == math::NumericUtil::isZero<F>(x) )
        {
            throw math::FunctionException(math::FunctionException::UNDEFINED);
        }

        const F u = static_cast<F>(1) / x;

        return f(u) * u * u;
    }

    /*
     * Inverse of 1/x, in this case it also equals 1/x
     *
     * @param x - input argument
     *
     * @return 1/x
     *
     * @throw FunctionException if inv() is not defined at the given 'x' (when x==0)
     */
    F inv(const F& x) const throw(math::FunctionException)
    {
        if ( true == math::NumericUtil::isZero<F>(x) )
        {
        	throw math::FunctionException(math::FunctionException::UNDEFINED);
        }

        return static_cast<F>(1) / x;
    }

    /*
     * Inverse of -inf.
     * It returns -eps instead of 0
     *
     * @return -eps
     *
     * @throw never thrown from this class
     */
    F invNegInf() const throw(math::FunctionException)
    {
        return -this->EPS;
        //return static_cast<F>(0);
    }

    /*
     * Inverse of +inf.
     * It returns +eps instead of 0
     *
     * @return eps
     *
     * @throw never thrown from this class
     */
    F invPosInf() const throw(math::FunctionException)
    {
        return this->EPS;
        //return static_cast<F>(0);
    }

    /*
     * Indicates that
     *   a < b  ==>  1/a > 1/b
     *
     * @return always true
     */
    bool swappedLimits() const
    {
        return true;
    }
};



/*
 * Swaps two variables if necessary, i.e. if
 * 'swapped' equals FALSE.
 *
 * @param swapped - if 'false', 'v1' and'v2' will be swapped
 * @param v1 - the first variable
 * @param v2 - the second variable
 */
template <typename F>
inline void __checkSwap(const bool swapped, F& v1, F& v2)
{
    if ( false == swapped )
    {
        std::swap(v1, v2);
    }
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
 * @return definite integral of f(x) between 'a' and 'b'
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
        const size_t n,
        const math::EIntegAlg::alg algorithm
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
 * @return definite integral of f(x) between 'a' and 'b'
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
        const math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    math::Integ::__private::__checkStep<F>(h);

    // Just obtain the (approximate) number of integrating intervals
    // from the size
    const size_t n = static_cast<size_t>(std::floor(std::abs(b-a) / h) +
    		static_cast<F>(1) / static_cast<F>(2) ) + 1;

    return math::Integ::integ<F>( f, a, b, n, algorithm );
}


/**
 * Improper definite integral between -inf and a finite limit:
 *
 *      b
 *      /
 *      | f(x) dx
 *      /
 *    -inf
 *
 * If 'b' is positive, the integration is performed in two parts:
 * 'f' in the range between -inf and 'bp' (a.k.a. the improper range) is
 * integrated by substitution, in the range between 'bp' and 'b'
 * (a.k.a. the proper part) it is integrated as a proper integral.
 * If 'b' is negative, 'f' is integrated as the improper integral
 * over the entire range.
 *
 * @note the result may blow up if the finite integral of 'f' does not
 *       exist for the given range
 *
 * @param f - function to integrate
 * @param b - upper limit of the integration range
 * @param nimp - desired number of integration intervals for the "improper part" (default: 10000)
 * @param nprop - desired number of integration intervals for the "proper part" when 'b' is positive (default: 10000)
 * @param bp - breakpoint between the "proper" and "improper" part, ignored if 'b' is negative (default: -5)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f(x) between -inf and 'b'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between -inf and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integImpNegInf(
        const math::IFunctionGeneric<F>& f,
        const F& b,
        const size_t nimp,
        const size_t nprop,
        const F& bp,
        const math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    // sanity check
    if ( b >= static_cast<F>(0) && bp >= static_cast<F>(0) )
    {
        throw math::CalculusException(math::CalculusException::INVALID_BREAKPOINT);
    }

    try
    {
    	F retVal = static_cast<F>(0);

        const math::Integ::__private::IntegSubstRec<F> fsub(f);

        const bool swappped = fsub.swappedLimits();

        F from;
        F to;

        if ( b < -math::NumericUtil::getEPS<F>() )
        {
            /*
             * 'b' is negative, only the improper integral by substitution
             * is needed.
             */

            from = fsub.inv(b);
            to = fsub.invNegInf();
            math::Integ::__private::__checkSwap<F>(swappped, from, to);

            retVal = math::Integ::integ<F>(fsub, from, to, nimp, algorithm);
        }
        else
        {
            /*
             * 'b' is positive, the proper integral (between 'bp' and 'b')
             * and the improper integral by substition are needed.
             */

            from = fsub.inv(bp);
            to = fsub.invNegInf();
            math::Integ::__private::__checkSwap<F>(swappped, from, to);

            retVal =
               math::Integ::integ<F>(fsub, from, to, nimp, algorithm) +
               math::Integ::integ<F>(f, bp, b, nprop, algorithm);
        }

        return retVal;
    }
    catch ( const math::FunctionException& fex )
    {
    	throw math::CalculusException(math::CalculusException::UNDEFINED);
    }

}


/**
 * Improper definite integral between a finite limit and inf:
 *
 *     inf
 *      /
 *      | f(x) dx
 *      /
 *      a
 *
 * If 'a' is negative, the integration is performed in two parts:
 * 'f' in the range between 'a' and 'bp' (a.k.a. the proper range) is
 * integrated as a proper integral, in the range between 'bp' and inf
 * (a.k.a. the improper part) it is integrated as an improper integral
 * (by substitution). If 'a' is positive, 'f' is integrated as the improper
 * integral over the entire range.
 *
 * @note the result may blow up if the finite integral of 'f' does not
 *       exist for the given range
 *
 * @param f - function to integrate
 * @param a - lower limit of the integration range
 * @param nimp - desired number of integration intervals for the "improper part" (default: 10000)
 * @param nprop - desired number of integration intervals for the "proper part" when 'a' is negative (default: 10000)
 * @param bp - breakpoint between the "proper" and "improper" part, ignored if 'a' is positive (default: 5)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f(x) between 'a' and inf
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between -inf and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integImpPosInf(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const size_t nimp,
        const size_t nprop,
        const F& bp,
        const math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    // sanity check
    if ( a <= static_cast<F>(0) && bp <= static_cast<F>(0) )
    {
        throw math::CalculusException(math::CalculusException::INVALID_BREAKPOINT);
    }

    try
    {
        F retVal = static_cast<F>(0);

        const math::Integ::__private::IntegSubstRec<F> fsub(f);

        const bool swapped = fsub.swappedLimits();

        F from;
        F to;

        if ( a > math::NumericUtil::getEPS<F>() )
        {
            /*
             * 'a' is positive, only the improper integral by substitution
             * is needed.
             */

            from = fsub.invPosInf();
            to = fsub.inv(a);
            math::Integ::__private::__checkSwap<F>(swapped, from, to);

            retVal = math::Integ::integ<F>(fsub, from, to, nimp, algorithm);
        }
        else
        {
            /*
             * 'a' is negative, the proper integral (between 'a' and 'bp')
             * and the improper integral by substitution are needed.
             */

            from = fsub.invPosInf();
            to = fsub.inv(bp);
            math::Integ::__private::__checkSwap<F>(swapped, from, to);

            retVal =
               math::Integ::integ<F>(f, a, bp, nprop, algorithm) +
               math::Integ::integ<F>(fsub, from, to, nimp, algorithm);
        }

        return retVal;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


/**
 * Improper definite integral between -inf and +inf:
 *
 *     inf
 *      /
 *      | f(x) dx
 *      /
 *    -inf
 *
 * The integration is performed in three parts. Over both "improper" ranges
 * (between -inf and 'bpneg' and between 'bppos' and inf), 'f' is integrated
 * by substitution. Over the middle range (between 'bpneg' and 'bppos'), 'f'
 * is integrated as a proper integral.
 *
 * @note the result may blow up if the finite integral of 'f' does not
 *       exist for the given range
 *
 * @param f - function to integrate
 * @param nimp - desired number of integration intervals for each "improper part" (default: 10000)
 * @param nprop - desired number of integration intervals for the "proper part" (default: 10000)
 * @param bpneg - breakpoint between the negative "improper" and the "proper" part (default: -5)
 * @param bppos - breakpoint between the "proper" and the positive "improper" part (default: 5)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f(x) between -inf and +inf
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between -inf and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integImp(
        const math::IFunctionGeneric<F>& f,
        const size_t nimp,
        const size_t nprop,
        const F& bpneg,
        const F& bppos,
        const math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    // sanity check
    if (bpneg >= static_cast<F>(0) || bppos <= static_cast<F>(0) )
    {
        throw math::CalculusException(math::CalculusException::INVALID_BREAKPOINT);
    }

    try
    {
        const math::Integ::__private::IntegSubstRec<F> fsub(f);

        const bool swapped = fsub.swappedLimits();

        F retVal;
        F from;
        F to;

        /*
         * Between -inf and 'bpneg', 'f' is integrated as an
         * improper integral by substitution
         */

        from = fsub.inv(bpneg);
        to = fsub.invNegInf();
        math::Integ::__private::__checkSwap<F>(swapped, from, to);

        retVal = math::Integ::integ<F>(fsub, from, to, nimp, algorithm);

        /*
         * Between 'bpneg' and 'bppos', 'f' is integrated as a proper integral
         */

        retVal += math::Integ::integ<F>(f, bpneg, bppos, nprop, algorithm);

        /*
         * Between 'bppos' and +inf, 'f' is integrated as an
         * improper integral by substitution
         */

        from = fsub.invPosInf();
        to =  fsub.inv(bppos);
        math::Integ::__private::__checkSwap<F>(swapped, from, to);

        retVal += math::Integ::integ<F>(fsub, from, to, nimp, algorithm);

        return retVal;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


/**
 * Improper definite integral between -inf and a finite limit:
 *
 *      b
 *      /
 *      | f(x) dx
 *      /
 *    -inf
 *
 * If 'b' is positive, the integration is performed in two parts:
 * 'f' in the range between -inf and 'bp' (a.k.a. the improper range) is
 * integrated by substitution, in the range between 'bp' and 'b'
 * (a.k.a. the proper part) it is integrated as a proper integral.
 * If 'b' is negative, 'f' is integrated as the improper integral
 * in the entire range.
 *
 * @note the result may blow up if the finite integral of 'f' does not
 *       exist for the given range
 *
 * @note The actual size of integrating step may be slightly decreased
 *       if required by the selected algorithm.
 *
 * @param f - function to integrate
 * @param b - upper limit of the integration range
 * @param himp - desired size of an integrating step for the "improper part" (default: 10000)
 * @param hprop - desired size of an integrating step for the "proper part" when 'b' is positive (default: 10000)
 * @param bp - breakpoint between the "proper" and "improper" part, ignored if 'b' is negative (default: -5)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f(x) between -inf and 'b'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between 'a' and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integImpNegInfH(
        const math::IFunctionGeneric<F>& f,
        const F& b,
        const F& himp,
        const F& hprop,
        const F& bp,
        const math::EIntegAlg::alg algorithm
     ) throw(math::CalculusException)
{
    // sanity check
    if (  b >= static_cast<F>(0) && bp >= static_cast<F>(0) )
    {
        throw math::CalculusException(math::CalculusException::INVALID_BREAKPOINT);
    }

    math::Integ::__private::__checkStep<F>(himp);
    math::Integ::__private::__checkStep<F>(hprop);

    try
    {
        F retVal = static_cast<F>(0);

        const math::Integ::__private::IntegSubstRec<F> fsub(f);

        const bool swapped = fsub.swappedLimits();

        F from;
        F to;

        if ( b < -math::NumericUtil::getEPS<F>() )
        {
            /*
             * 'b' is negative, only the improper integral by substitution
             * is needed.
             */

            from = fsub.inv(b);
            to = fsub.invNegInf();
            math::Integ::__private::__checkSwap<F>(swapped, from, to);

            retVal = math::Integ::integH<F>(fsub, from, to, himp, algorithm);
        }
        else
        {
            /*
             * 'b' is positive, the proper integral (between 'bp' and 'b')
             * and the improper integral by substitution are needed.
             */

            from = fsub.inv(bp);
            to = fsub.invNegInf();
            math::Integ::__private::__checkSwap<F>(swapped, from, to);

            retVal =
               math::Integ::integH<F>(fsub, from, to, himp, algorithm) +
               math::Integ::integH<F>(f, bp, b, hprop, algorithm);
        }

        return retVal;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


/**
 * Improper definite integral between a finite limit and inf:
 *
 *     inf
 *      /
 *      | f(x) dx
 *      /
 *      a
 *
 * If 'a' is negative, the integration is performed in two parts:
 * 'f' in the range between 'a' and 'bp' (a.k.a. the proper range) is
 * integrated as a proper integral, in the range between 'bp' and inf
 * (a.k.a. the improper part) it is integrated as an improper integral
 * (by substitution). If 'a' is negative, 'f' is integrated as the improper
 * integral over the entire range.
 *
 * @note the result may blow up if the finite integral of 'f' does not
 *       exist for the given range
 *
 * @param f - function to integrate
 * @param a - lower limit of the integration range
 * @param himp - desired size of an integrating step for the "improper part" (default: 0.0001)
 * @param hprop - desired size of an integrating step for the "proper part" when 'a' is negative (default: 0.0001)
 * @param bp - breakpoint between the "proper" and "improper" part, ignored if 'a' is positive (default: 5)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f(x) between 'a' and inf
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between -inf and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integImpPosInfH(
        const math::IFunctionGeneric<F>& f,
        const F& a,
        const F& himp,
        const F& hprop,
        const F& bp,
        const math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    // sanity check
    if ( a <= static_cast<F>(0) && bp <= static_cast<F>(0) )
    {
	    throw math::CalculusException(math::CalculusException::INVALID_BREAKPOINT);
    }

    math::Integ::__private::__checkStep<F>(himp);
    math::Integ::__private::__checkStep<F>(hprop);

    try
    {
        F retVal = static_cast<F>(0);

        const math::Integ::__private::IntegSubstRec<F> fsub(f);

        const bool swapped = fsub.swappedLimits();

        F from;
        F to;

        if ( a > math::NumericUtil::getEPS<F>() )
        {
            /*
             * 'a' is positive, only the improper integral by substitution
             * is needed.
             */

            from = fsub.invPosInf();
            to = fsub.inv(a);
            math::Integ::__private::__checkSwap<F>(swapped, from, to);

            retVal = math::Integ::integH<F>(fsub, from, to, himp, algorithm);
        }
        else
        {
            /*
             * 'a' is negative, the proper integral (between 'a' and 'bp')
             * and the improper integral by substitution are needed.
             */

            from = fsub.invPosInf();
            to = fsub.inv(bp);
            math::Integ::__private::__checkSwap<F>(swapped, from, to);

            retVal =
               math::Integ::integH<F>(f, a, bp, hprop, algorithm) +
               math::Integ::integH<F>(fsub, from, to, himp, algorithm);
        }

        return retVal;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


/**
 * Improper definite integral between -inf and +inf:
 *
 *     inf
 *      /
 *      | f(x) dx
 *      /
 *    -inf
 *
 * The integration is performed in three parts. Over both "improper" ranges
 * (between -inf and 'bpneg' and between 'bppos' and inf), 'f' is integrated
 * by substitution. Over the middle range (between 'bpneg' and 'bppos'), 'f'
 * is integrated as a proper integral.
 *
 * @note the result may blow up if the finite integral of 'f' does not
 *       exist for the given range
 *
 * @param f - function to integrate
 * @param himp - desired size of an integrating step for both "improper" parts (default: 0.0001)
 * @param hprop - desired size of an integrating step for the "proper part" (default: 0.0001)
 * @param bpneg - breakpoint between negative "improper" and the "proper" part (default: -5)
 * @param bppos - breakpoint between the "proper" and the positive "improper" part (default: 5)
 * @param algorithm - one of the supported algorithms to obtain the definite integral (default: SIMPSON)
 *
 * @return definite integral of f(x) between -inf and +inf
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined between -inf and 'b'
 *
 * @see IFunctionGeneric
 */
template <typename F>
F math::Integ::integImpH(
        const math::IFunctionGeneric<F>& f,
        const F& himp,
        const F& hprop,
        const F& bpneg,
        const F& bppos,
        const math::EIntegAlg::alg algorithm
      ) throw(math::CalculusException)
{
    // sanity check
    if (bpneg >= static_cast<F>(0) || bppos <= static_cast<F>(0) )
    {
        throw math::CalculusException(math::CalculusException::INVALID_BREAKPOINT);
    }

    math::Integ::__private::__checkStep<F>(himp);
    math::Integ::__private::__checkStep<F>(hprop);

    try
    {
        const math::Integ::__private::IntegSubstRec<F> fsub(f);

        const bool swapped = fsub.swappedLimits();

        F retVal;
        F from;
        F to;

        /*
         * Between -inf and 'bpneg', 'f' is integrated as an
         * improper integral by substitution
         */

        from = fsub.inv(bpneg);
        to = fsub.invNegInf();
        math::Integ::__private::__checkSwap<F>(swapped, from, to);

        retVal = math::Integ::integH<F>(fsub, from, to, himp, algorithm);

        /*
         * Between 'bpneg' and 'bppos', 'f' is integrated as a proper integral
         */

        retVal += math::Integ::integH<F>(f, bpneg, bppos, hprop, algorithm);

        /*
         * Between 'bppos' and +inf, 'f' is integrated as an
         * improper integral by substitution
         */

        from = fsub.invPosInf();
        to = fsub.inv(bppos);
        math::Integ::__private::__checkSwap<F>(swapped, from, to);

        retVal += math::Integ::integH<F>(fsub, from, to, himp, algorithm);

        return retVal;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}
