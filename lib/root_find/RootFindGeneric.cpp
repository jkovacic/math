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
 * Implementation of the class RootFindGeneric with several
 * root finding algorithms.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */


// No #include "RootFindGeneric.hpp" !!!
#include <cstddef>
#include <algorithm>

#include "util/NumericUtil.hpp"
#include "exception/FunctionException.hpp"
#include "util/IFunctionGeneric.hpp"
#include "calculus/DiffGeneric.hpp"



/**
 * Bisection method to find a root of a nonlinear function between
 * 'from' and 'to'.
 *
 * The method converges to a root if the following conditions are satisfied:
 * - 'f.func' is continuous between 'from' and 'two'
 * - "f.func's" values at 'from' and 'two' have opposite signs
 * - 'f.func' has exactly one root between 'from' and 'two'
 *
 * @note The method will swap 'from' and 'to' if 'from' is greater than 'to'.
 *
 * @param f - instance of a class with the function to find its root
 * @param from - lower boundary of the search interval
 * @param to - upper boundary of the search interval
 * @param EPSX - the smallest width of the search interval (default: NumericUtil::EPS)
 * @param EPSY - tolerance (default: NumericUtil::EPS)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 */
template <class T>
T math::RootFindGeneric<T>::bisection(
           const math::IFunctionGeneric<T>& f,
           const T& from,
           const T& to,
           const T& EPSX,
           const T& EPSY
         ) throw(math::RootFindException)
{
    /*
     * The bisection method is described in detail at:
     * https://en.wikipedia.org/wiki/Bisection_method
     */

    // sanity check
    if ( EPSX < math::NumericUtil<T>::getEPS() ||
         EPSY < math::NumericUtil<T>::getEPS() )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        T a, fa;
        T b, fb;
        T c, fc;
        short int sa, sb, sc;

        a = std::min(from, to);
        b = std::max(from, to);
        fa = f.func(a);
        fb = f.func(b);
        sa = math::NumericUtil<T>::sign(fa);
        sb = math::NumericUtil<T>::sign(fb);

        // Maybe a boundary is already a root?
        if ( true == math::NumericUtil<T>::isZero(fa, EPSY) )
        {
            return a;
        }

        if ( true == math::NumericUtil<T>::isZero(fb, EPSY) )
        {
            return b;
        }

        // In order to converge, fa and fb must be of opposite signs
        if ( (sa * sb) > 0)
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        // Start of the actual bisection method.
        // Proceed until the interval narrows down to EPSX.
        while ( false == math::NumericUtil<T>::isZero(b-a, EPSX) )
        {
            c = (a+b) / static_cast<T>(2);
            fc = f.func(c);
            sc = math::NumericUtil<T>::sign(fc);

            // Maybe c is a root?
            if ( true == math::NumericUtil<T>::isZero(c, EPSY) )
            {
                // if it is, just exit the loop
                break;  // out of while (!isZero(b-a))
            }

            // 'c' is not the root (yet).
            // Set new boundaries depending on fc's sign.
            if ( sa == sc )
            {
                a = c;
                fa = fc;
                sa = sc;
            }
            else
            {
                b = c;
                fb = fc;
                sb = sc;
            }

            // 'fa' and 'fb' must still be of opposite signs!
            // If the function has more than one root,this condition may be violated.
            if ( (sa * sb) > 0 )
            {
                throw math::RootFindException(math::RootFindException::NO_CONVERGENCE);
            }
        }

        return c;
    }  // try
    catch ( const math::FunctionException& fex )
    {
        throw math::RootFindException(math::RootFindException::UNDEFINED);
    }
}


/**
 * Regula falsi (or false position) method to find a root of a
 * nonlinear function between 'from' and 'two'.
 *
 *  The function converges to a root if the following conditions
 *  are satisfied:
 * - 'f.func' is continuous between 'from' and 'to'
 * - "f.func's" values at 'from' and 'two' have opposite signs
 * - 'f.func' has exactly one root between 'from' and 'to'
 *
 * @note The method will swap 'from' and 'to' if 'from' is greater than 'to'.
 *
 * @param f - instance of a class with the function to find its root
 * @param from - lower boundary of the search interval
 * @param to - upper boundary of the search interval
 * @param EPSX - the smallest width of the search interval (default: NumericUtil::EPS)
 * @param EPSY - tolerance (default: NumericUtil::EPS)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 */
template <class T>
T math::RootFindGeneric<T>::regulaFalsi(
           const math::IFunctionGeneric<T>& f,
           const T& from,
           const T& to,
           const T& EPSX,
           const T& EPSY
         ) throw(math::RootFindException)
{
    /*
     * The regula falsi method is described in detail at:
     * https://en.wikipedia.org/wiki/False_position_method
     */

    // sanity check
    if ( EPSX < math::NumericUtil<T>::getEPS() ||
         EPSY < math::NumericUtil<T>::getEPS() )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        T a, fa;
        T b, fb;
        T c, fc;
        short int sa, sb, sc;

        a = std::min(from, to);
        b   = std::max(from, to);
        fa = f.func(a);
        fb = f.func(b);
        sa = math::NumericUtil<T>::sign(fa);
        sb = math::NumericUtil<T>::sign(fb);

        // Maybe a boundary is already a root?
        if ( true == math::NumericUtil<T>::isZero(fa, EPSY) )
        {
            return a;
        }

        if ( true == math::NumericUtil<T>::isZero(fb, EPSY) )
        {
             return b;
        }

        // In order to converge, fa and fb must be of opposite signs
        if ( (sa * sb) > 0)
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        // Start of the actual regula falsi method.
        // Proceed until the interval narrows down to EPSX.
        while ( false == math::NumericUtil<T>::isZero(b-a, EPSX) )
        {
            c = a - fa * (b-a) / (fb-fa);
            fc = f.func(c);
            sc = math::NumericUtil<T>::sign(fc);

            // Maybe c is a root?
            if ( true == math::NumericUtil<T>::isZero(c, EPSY) )
            {
                // if it is, just exit the loop
                break;  // out of while (!isZero(b-a))
            }

            // 'c' is not the root (yet).
            // Set new boundaries depending on fc's sign.
            if ( sa == sc )
            {
                a = c;
                fa = fc;
                sa = sc;
            }
            else
            {
                b = c;
                fb = fc;
                sb = sc;
            }

            // 'fa' and 'fb' must still be of opposite signs!
            // If the function has more than one root,this condition may be violated.
            if ( (sa * sb) > 0 )
            {
                throw math::RootFindException(math::RootFindException::NO_CONVERGENCE);
            }
        }

        return c;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::RootFindException(math::RootFindException::UNDEFINED);
    }
}


/**
 * Secant method to find one root of a nonlinear function.
 *
 * @note The method is not guaranteed to converge towards any
 *       root even if 'f.func' has roots.
 *
 * @note If 'f.func' has more than root, the method will converge
 *       (not guaranteed) to one root only.
 *
 * @param f - instance of a class with the function to find its root
 * @param r0 - first initial value
 * @param r1 - second initial value
 * @param EPSY - tolerance (default: NumericUtil::EPS)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 */
template <class T>
T math::RootFindGeneric<T>::secant(
           const math::IFunctionGeneric<T>& f,
           const T& r0,
           const T& r1,
           const T& EPSY,
           size_t Nmax
         ) throw(math::RootFindException)
{
    /*
     * The secant method is described in detail at:
     * https://en.wikipedia.org/wiki/Secant_method
     */

    // sanity check
    if ( EPSY < math::NumericUtil<T>::getEPS() ||
         0 == Nmax )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        T xo = r0;
        T xn = r1;
        T fo = f.func(xo);
        T fn = f.func(xn);

        // Maybe xo is already a root?
        if ( true == math::NumericUtil<T>::isZero(fo, EPSY) )
        {
            return xo;
        }

        // Start of the actual secant method
        for ( size_t iter = Nmax;
              false==math::NumericUtil<T>::isZero(fn, EPSY) && iter > 0;
              --iter )
        {
            const T df = fo - fn;
        	if ( true == math::NumericUtil<T>::isZero(df) )
        	{
                throw math::RootFindException(math::RootFindException::NO_CONVERGENCE);
        	}

            // Obtain new 'xn' (and 'fn') and update 'xo' (and 'fo')
            const T c = xo - fo * (xo-xn) / df;
            xo = xn;
            fo = fn;
            xn = c;
            fn = f.func(xn);
        }  // for iter

        // Finally check whether the algorithm has converged to any root
        if ( false == math::NumericUtil<T>::isZero(fn, EPSY) )
        {
            // apparently not
            throw math::RootFindException(math::RootFindException::NO_CONVERGENCE);
        }

        return xn;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::RootFindException(math::RootFindException::UNDEFINED);
    }
}


/**
 * Newton - Raphson (or shorter: Newton's) method to find
 * one root of a nonlinear function.
 *
 * The function's slope (derivation) must be known in advance and passed
 * as the second argument.
 *
 * @note The method is not guaranteed to converge towards any root
 *       even if 'f.func' has roots
 *
 * @note If 'f.func' has more than one root, the method will converge
 *       (not guaranteed) to one root only
 *
 * @param f - instance of a class with the function to find its root
 * @param diff - instance of a class with the derivation of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param EPSY - tolerance (default: NumericUtil::EPS)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see quasiNewton
 */
template <class T>
T math::RootFindGeneric<T>::newton(
           const math::IFunctionGeneric<T>& f,
           const math::IFunctionGeneric<T>& diff,
           const T& x0,
           const T& EPSY,
           size_t Nmax
         ) throw(math::RootFindException)
{
    /*
     * The Newton - Raphson method is described in detail at:
     * https://en.wikipedia.org/wiki/Newton%27s_method
     */

    // The algorithm is implemented in __newtonCommon:
    return __newtonCommon(f, diff, x0, EPSY, static_cast<T>(0), Nmax, true);
}


/**
 * Newton - Raphson (or shorter: Newton's) method to find
 * one root of a nonlinear function.
 *
 * Unlike 'newton', this method does not require a derivation of 'f'.
 * Instead it performs numerical differentiation using the default
 * differentiation method, i.e EDiffMethod::CENTRAL.
 *
 * @note The method is not guaranteed to converge towards any root
 *       even if 'f.func' has roots
 *
 * @note If 'f.func' has more than one root, the method will converge
 *       (not guaranteed) to one root only
 *
 * @param f - instance of a class with the function to find its root
 * @param x0 - starting point of the algorithm
 * @param EPSY - tolerance (default: NumericUtil::EPS)
 * @param h - step size for numerical derivation (default: 0.001)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see newton
 * @see EdiffMethod
 */
template <class T>
T math::RootFindGeneric<T>::quasiNewton(
           const math::IFunctionGeneric<T>& f,
           const T& x0,
           const T& EPSY,
           const T& h,
           size_t Nmax
         ) throw(math::RootFindException)
{
    /*
     * The Newton - Raphson method is described in detail at:
     * https://en.wikipedia.org/wiki/Newton%27s_method
     */

    // The algorithm is implemented in __newtonCommon:
    return __newtonCommon(f, f, x0, EPSY, h, Nmax, false);
}


/*
 * Common implementation of the Newton - Raphson algorithm.
 * This is a private method, called by both public methods.
 * Depending on 'diffFunc', it will evaluate f's slope either
 * via the provided function 'diff' or numerically via the central
 * method. In the latter case 'diff' is ignored.
 *
 * @param f - instance of a class with the function to find its root
 * @param diff - instance of a class with the derivation of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param EPSY - tolerance
 * @param h - step size for numerical derivation (ignored if 'diffFunc'==true)
 * @param Nmax - maximum number of iterations
 * @param diffFunc - if true evaluate f's slope via 'diff', otherwise numerically
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see newton
 * @see quasiNewton
 */
template <class T>
T math::RootFindGeneric<T>::__newtonCommon(
           const math::IFunctionGeneric<T>& f,
           const math::IFunctionGeneric<T>& diff,
           const T& x0,
           const T& EPSY,
           const T& h,
           size_t Nmax,
           bool diffFunc
         ) throw(math::RootFindException)
{
    /*
     * The Newton - Raphson method is described in detail at:
     * https://en.wikipedia.org/wiki/Newton%27s_method
     */

    // sanity check
    if ( EPSY < math::NumericUtil<T>::getEPS() ||
         0 == Nmax )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        T x = x0;
        T fx = f.func(x);

        // start of the actual Newton - Raphson method
        for ( size_t iter = Nmax;
              false == math::NumericUtil<T>::isZero(fx, EPSY) && iter > 0;
              --iter )
        {
            const T dx = ( true == diffFunc ?
                      diff.func(x) :
                      math::DiffGeneric<T>::diff(f, x,h, math::EDiffMethod::CENTRAL) );

            if ( true == math::NumericUtil<T>::isZero(dx) )
            {
                // cannot continue as division by zero would occur inthe next step
                throw math::RootFindException(math::RootFindException::ZERO_SLOPE);
            }

            // Update 'x':
            x -= fx / dx;
            fx = f.func(x);
        }  // for iter

        // Has the algorithm converged to any root?
        if ( false == math::NumericUtil<T>::isZero(fx, EPSY) )
        {
            // apparently not
            throw math::RootFindException(math::RootFindException::NO_CONVERGENCE);
        }

        return x;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::RootFindException(math::RootFindException::UNDEFINED);
    }
    catch ( const math::CalculusException& cex )
    {
        // cex is only thrown if derivation function is undefined near x
        throw math::RootFindException(math::RootFindException::UNDEFINED);;
    }
}
