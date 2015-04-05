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
 * Implementation of functions within the namespace RootFind with several
 * root finding algorithms.
 */


// No #include "RootFindGeneric.hpp" !!!
#include <cstddef>
#include <algorithm>

#include "../settings/rootfind_settings.h"
#include "util/NumericUtil.hpp"
#include "exception/FunctionException.hpp"
#include "util/IFunctionGeneric.hpp"
#include "calculus/DiffGeneric.hpp"


// Definition of "private" functions

namespace math
{

namespace RootFind
{

namespace __private
{

/*
 * Common implementation of the Newton - Raphson algorithm.
 * This is a private method, called by both public methods.
 * Depending on 'diffFunc', it will evaluate f's slope either
 * via the provided function 'diff' or numerically via the central
 * method. In the latter case 'diff' is ignored.
 *
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param diff - instance of a class with the derivation of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance
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
template <typename F>
F __newtonCommon(
           const math::IFunctionGeneric<F>& f,
           const math::IFunctionGeneric<F>& diff,
           const F& x0,
           const F& epsy,
           const F& h,
           const size_t Nmax,
           bool diffFunc
         ) throw(math::RootFindException)
{
    /*
     * The Newton - Raphson method is described in detail at:
     * https://en.wikipedia.org/wiki/Newton%27s_method
     */

    // sanity check
    if ( 0 == Nmax )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        const F EPSY = std::max(epsy, math::NumericUtil::getEPS<F>() );

        // sanity check of the final epsilon:
        if ( EPSY <= static_cast<F>(0) )
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        F x = x0;
        F fx = f(x);

        // start of the actual Newton - Raphson method
        for ( size_t iter = Nmax;
              false == math::NumericUtil::isZero<F>(fx, EPSY) && iter > 0;
              --iter )
        {
            const F dx = ( true == diffFunc ?
                      diff(x) :
                      math::Diff::diff<F>(f, x, h, ROOTFIND_DEFAULT_DIFF_METHOD) );

            if ( true == math::NumericUtil::isZero<F>(dx) )
            {
                // cannot continue as division by zero would occur in the next step
                throw math::RootFindException(math::RootFindException::ZERO_SLOPE);
            }

            // Update 'x':
            x -= fx / dx;
            fx = f(x);
        }  // for iter

        // Has the algorithm converged to any root?
        if ( false == math::NumericUtil::isZero<F>(fx, EPSY) )
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
        throw math::RootFindException(math::RootFindException::UNDEFINED);
    }
}


/*
 * Common implementation of the Halley's algorithm which is actually
 * extended Newton's method.
 *
 * This is a private method, called by public methods.
 *
 * Depending on 'diffFunc', it will evaluate f's derivatives either
 * via the provided functions 'diff' and 'diff2' or numerically via
 * the central method. In the latter case 'diff' and 'diff2' are ignored.
 *
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param d - instance of a class with the derivation of 'f.func'
 * @param d2 - instance of a class with the 2nd order derivative of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance
 * @param h - step size for numerical derivation (ignored if 'diffFunc'==true)
 * @param Nmax - maximum number of iterations
 * @param diffFunc - if true evaluate f's slope via 'diff', otherwise numerically
 * @param mod - if 'true', adjust the Halley's term to never be less than 0.8
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see halley
 * @see quasiHalley
 * @see halleyMod
 * @see quasiHalleyMod
 */
template <typename F>
F __halleyCommon(
        const math::IFunctionGeneric<F>& f,
        const math::IFunctionGeneric<F>& d,
        const math::IFunctionGeneric<F>& d2,
        const F& x0,
        const F& epsy,
        const F& h,
        const size_t Nmax,
        const bool diffFunc,
        const bool mod
      ) throw (math::RootFindException)
{
    /*
     * The Halley's method is described in detail at:
     * https://en.wikipedia.org/wiki/Halley%27s_method
     */

    // sanity check
    if ( 0 == Nmax )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        const F EPSY = std::max(epsy, math::NumericUtil::getEPS<F>() );

        // sanity check of the final epsilon:
        if ( EPSY <= static_cast<F>(0) )
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        // Evaluate 'f' at the initial point
        F x = x0;
        F fx = f(x);

        // start of the actual Halley's method
        for ( size_t iter = Nmax;
              false == math::NumericUtil::isZero<F>(fx, EPSY) && iter > 0;
              --iter )
        {
            // Evaluate both derivatives at 'x':
            const F dx = ( true == diffFunc ?
                     d(x) :
                     math::Diff::diff<F>(f, x, h, ROOTFIND_DEFAULT_DIFF_METHOD) );

            const F d2x = ( true== diffFunc ?
                     d2(x) :
                     math::Diff::diff2<F>(f, x, h) );

            // check the first derivative to prevent possible division by zero
            if ( true == math::NumericUtil::isZero<F>(dx) )
            {
                // cannot continue as division by zero would occur in the next step
                throw math::RootFindException(math::RootFindException::ZERO_SLOPE);
            }

            /*
             * x will be updated as follows:
             *
             *                                f(xi)
             *  x(i+1) = x(i) - ------------------------------------
             *                            /      f(xi) * f''(xi)  \
             *                   f'(xi) * | 1 - ----------------- |
             *                            \       2 * f'(xi)^2    /
             */

            // Obtain the Halley's term:
            F term = static_cast<F>(1) - fx * d2x  / (static_cast<F>(2) * d2x * d2x );

            // Theoretically the term can also equal 0...
            if ( false==mod && true==math::NumericUtil::isZero<F>(term) )
            {
                // cannot continue as division by zero would occur in the next step
                throw math::RootFindException(math::RootFindException::ZERO_SLOPE);
            }

            // Adjust the term to min. 0.8 if required
            if ( true == mod )
            {
                term = std::max( static_cast<F>(4) / static_cast<F>(5),
                                 std::min(static_cast<F>(6) / static_cast<F>(5), term) );
            }

            // Finally update x
            x -= fx / (dx * term);
            fx = f(x);
        } // for iter

        // Has the algorithm converged to any root?
        if ( false == math::NumericUtil::isZero<F>(fx, EPSY) )
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
        throw math::RootFindException(math::RootFindException::UNDEFINED);
    }
}


}  // namespace __private
}  // namespace RootFind
}  // namespace math


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
 * @note 'epsx' and/or 'epsy' will be increased to the system epsilon for the type
 *       if their values are too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param from - lower boundary of the search interval
 * @param to - upper boundary of the search interval
 * @param epsx - the smallest width of the search interval (default: 1e-4)
 * @param epsy - tolerance (default: 1e-6)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 */
template <typename F>
F math::RootFind::bisection(
           const math::IFunctionGeneric<F>& f,
           const F& from,
           const F& to,
           const F& epsx,
           const F& epsy
         ) throw(math::RootFindException)
{
    /*
     * The bisection method is described in detail at:
     * https://en.wikipedia.org/wiki/Bisection_method
     */

    try
    {
        const F EPSX = std::max(epsx, math::NumericUtil::getEPS<F>() );
        const F EPSY = std::max(epsy, math::NumericUtil::getEPS<F>() );

        // sanity check of final epsilons:
        if ( EPSX <= static_cast<F>(0) ||
             EPSY <= static_cast<F>(0) )
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        F a, fa;
        F b, fb;
        F c, fc;
        short int sa, sb, sc;

        a = std::min(from, to);
        b = std::max(from, to);
        fa = f(a);
        fb = f(b);
        sa = math::NumericUtil::sign<F>(fa);
        sb = math::NumericUtil::sign<F>(fb);

        // Maybe a boundary is already a root?
        if ( true == math::NumericUtil::isZero<F>(fa, EPSY) )
        {
            return a;
        }

        if ( true == math::NumericUtil::isZero<F>(fb, EPSY) )
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
        while ( false == math::NumericUtil::isZero<F>(b-a, EPSX) )
        {
            c = (a+b) / static_cast<F>(2);
            fc = f(c);
            sc = math::NumericUtil::sign<F>(fc);

            // Maybe c is a root?
            if ( true == math::NumericUtil::isZero<F>(c, EPSY) )
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
 * @note 'epsx' and/or 'epsy' will be increased to the system epsilon for the type
 *       if their values are too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param from - lower boundary of the search interval
 * @param to - upper boundary of the search interval
 * @param epsx - the smallest width of the search interval (default: 1e-4)
 * @param epsy - tolerance (default: 1e-6)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 */
template <typename F>
F math::RootFind::regulaFalsi(
           const math::IFunctionGeneric<F>& f,
           const F& from,
           const F& to,
           const F& epsx,
           const F& epsy
         ) throw(math::RootFindException)
{
    /*
     * The regula falsi method is described in detail at:
     * https://en.wikipedia.org/wiki/False_position_method
     */

    try
    {
        const F EPSX = std::max(epsx, math::NumericUtil::getEPS<F>() );
        const F EPSY = std::max(epsy, math::NumericUtil::getEPS<F>() );

        // sanity check of final epsilons:
        if ( EPSX <= static_cast<F>(0) ||
             EPSY <= static_cast<F>(0) )
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        F a, fa;
        F b, fb;
        F c, fc;
        short int sa, sb, sc;

        a = std::min(from, to);
        b = std::max(from, to);
        fa = f(a);
        fb = f(b);
        sa = math::NumericUtil::sign<F>(fa);
        sb = math::NumericUtil::sign<F>(fb);

        // Maybe a boundary is already a root?
        if ( true == math::NumericUtil::isZero<F>(fa, EPSY) )
        {
            return a;
        }

        if ( true == math::NumericUtil::isZero<F>(fb, EPSY) )
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
        while ( false == math::NumericUtil::isZero<F>(b-a, EPSX) )
        {
            c = a - fa * (b-a) / (fb-fa);
            fc = f(c);
            sc = math::NumericUtil::sign<F>(fc);

            // Maybe c is a root?
            if ( true == math::NumericUtil::isZero<F>(c, EPSY) )
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
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param r0 - first initial value
 * @param r1 - second initial value
 * @param epsy - tolerance (default: 1e-6)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 */
template <typename F>
F math::RootFind::secant(
           const math::IFunctionGeneric<F>& f,
           const F& r0,
           const F& r1,
           const F& epsy,
           const size_t Nmax
         ) throw(math::RootFindException)
{
    /*
     * The secant method is described in detail at:
     * https://en.wikipedia.org/wiki/Secant_method
     */

    // sanity check
    if ( 0 == Nmax )
    {
        throw math::RootFindException(math::RootFindException::INVALID_ARGS);
    }

    try
    {
        const F EPSY = std::max(epsy, math::NumericUtil::getEPS<F>() );

        // sanity check of the final epsilon:
        if ( EPSY <= static_cast<F>(0) )
        {
            throw math::RootFindException(math::RootFindException::INVALID_ARGS);
        }

        F xo = r0;
        F xn = r1;
        F fo = f(xo);
        F fn = f(xn);

        // Maybe xo is already a root?
        if ( true == math::NumericUtil::isZero<F>(fo, EPSY) )
        {
            return xo;
        }

        // Start of the actual secant method
        for ( size_t iter = Nmax;
              false==math::NumericUtil::isZero<F>(fn, EPSY) && iter > 0;
              --iter )
        {
            const F df = fo - fn;
            if ( true == math::NumericUtil::isZero<F>(df) )
            {
                throw math::RootFindException(math::RootFindException::NO_CONVERGENCE);
            }

            // Obtain new 'xn' (and 'fn') and update 'xo' (and 'fo')
            const F c = xo - fo * (xo-xn) / df;
            xo = xn;
            fo = fn;
            xn = c;
            fn = f(xn);
        }  // for iter

        // Finally check whether the algorithm has converged to any root
        if ( false == math::NumericUtil::isZero<F>(fn, EPSY) )
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
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param diff - instance of a class with the derivation of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance (default: 1e-6)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see quasiNewton
 */
template <typename F>
F math::RootFind::newton(
           const math::IFunctionGeneric<F>& f,
           const math::IFunctionGeneric<F>& diff,
           const F& x0,
           const F& epsy,
           const size_t Nmax
         ) throw(math::RootFindException)
{
    /*
     * The Newton - Raphson method is described in detail at:
     * https://en.wikipedia.org/wiki/Newton%27s_method
     */

    // The algorithm is implemented in __newtonCommon:
    return math::RootFind::__private::__newtonCommon<F>(
            f, diff, x0, epsy, static_cast<F>(0), Nmax, true );
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
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance (default: 1e-6)
 * @param h - step size for numerical derivation (default: 0.001)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see newton
 */
template <typename F>
F math::RootFind::quasiNewton(
           const math::IFunctionGeneric<F>& f,
           const F& x0,
           const F& epsy,
           const F& h,
           const size_t Nmax
         ) throw(math::RootFindException)
{
    /*
     * The Newton - Raphson method is described in detail at:
     * https://en.wikipedia.org/wiki/Newton%27s_method
     */

    // The algorithm is implemented in __newtonCommon:
    return math::RootFind::__private::__newtonCommon<F>(
            f, f, x0, epsy, h, Nmax, false );
}


/**
 * Halley's method to find one root of a nonlinear function.
 *
 * The function's 1st and 2nd order derivations must be known
 * in advance and passed as the second and third argument.
 *
 * @note The method is not guaranteed to converge towards any root
 *       even if 'f.func' has roots
 *
 * @note If 'f.func' has more than one root, the method will converge
 *       (not guaranteed) to one root only
 *
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param diff - instance of a class with the derivation of 'f.func'
 * @param diff2 - instance of a class with the 2nd order derivative of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance (default: 1e-6)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see quasiHalley
 * @see halleyMod
 * @see quasiHalleyMod
 */
template <typename F>
F math::RootFind::halley(
       const math::IFunctionGeneric<F>& f,
       const math::IFunctionGeneric<F>& diff,
       const math::IFunctionGeneric<F>& diff2,
       const F& x0,
       const F& epsy,
       const size_t Nmax
     ) throw (math::RootFindException)
{
    /*
     * The Halley's method is described in detail at:
     * https://en.wikipedia.org/wiki/Halley%27s_method
     */

    // The algorithm is implemented in __halleyCommon:
    return math::RootFind::__private::__halleyCommon<F>(
        f, diff, diff2, x0, epsy, static_cast<F>(0), Nmax, true, false);
}


/**
 * Halley's method to find one root of a nonlinear function.
 *
 * Unlike 'halley', this method does not require any derivation of 'f'.
 * Instead it performs numerical differentiation.
 *
 * @note The method is not guaranteed to converge towards any root
 *       even if 'f.func' has roots
 *
 * @note If 'f.func' has more than one root, the method will converge
 *       (not guaranteed) to one root only
 *
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance (default: 1e-6)
 * @param h - step size for numerical derivation (default: 0.001)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see halley
 * @see halleyMod
 * @see quasiHalleyMod
 */
template <typename F>
F math::RootFind::quasiHalley(
       const math::IFunctionGeneric<F>& f,
       const F& x0,
       const F& epsy,
       const F& h,
       const size_t Nmax
     ) throw (math::RootFindException)
{
    /*
     * The Halley's method is described in detail at:
     * https://en.wikipedia.org/wiki/Halley%27s_method
     */

    // The algorithm is implemented in __halleyCommon:
	return math::RootFind::__private::__halleyCommon<F>(
        f, f, f, x0, epsy, h, Nmax, false, false);
}


/**
 * Modified Halley's method to find one root of a nonlinear function.
 * The Halley's term is adjusted to be never less than 0.8
 *
 * The function's 1st and 2nd order derivations must be known
 * in advance and passed as the second and third argument.
 *
 * @note The method is not guaranteed to converge towards any root
 *       even if 'f.func' has roots
 *
 * @note If 'f.func' has more than one root, the method will converge
 *       (not guaranteed) to one root only
 *
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param diff - instance of a class with the derivation of 'f.func'
 * @param diff2 - instance of a class with the 2nd order derivative of 'f.func'
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance (default: 1e-6)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see halley
 * @see quasiHalley
 * @see quasiHalleyMod
 */
template <typename F>
F math::RootFind::halleyMod(
       const math::IFunctionGeneric<F>& f,
       const math::IFunctionGeneric<F>& diff,
       const math::IFunctionGeneric<F>& diff2,
       const F& x0,
       const F& epsy,
       const size_t Nmax
     ) throw (math::RootFindException)
{
    /*
     * The Halley's method is described in detail at:
     * https://en.wikipedia.org/wiki/Halley%27s_method
     */

    // The algorithm is implemented in __halleyCommon:
	return math::RootFind::__private::__halleyCommon<F>(
        f, diff, diff2, x0, epsy, static_cast<F>(0), Nmax, true, true);
}


/**
 * Halley's method to find one root of a nonlinear function.
 * The Halley's term is adjusted to be never less than 0.8
 *
 * Unlike 'halley', this method does not require any derivation of 'f'.
 * Instead it performs numerical differentiation.
 *
 * @note The method is not guaranteed to converge towards any root
 *       even if 'f.func' has roots
 *
 * @note If 'f.func' has more than one root, the method will converge
 *       (not guaranteed) to one root only
 *
 * @note 'epsy' will be increased to the system epsilon for the type
 *       if its value is too small.
 *
 * @param f - instance of a class with the function to find its root
 * @param x0 - starting point of the algorithm
 * @param epsy - tolerance (default: 1e-6)
 * @param h - step size for numerical derivation (default: 0.001)
 * @param Nmax - maximum number of iterations (default: 10000)
 *
 * @return root of the function 'f' if it exists
 *
 * @throw RootFindException if the root could not be found
 *
 * @see halley
 * @see halleyMod
 * @see quasiHalleyMod
 */
template <typename F>
F math::RootFind::quasiHalleyMod(
       const math::IFunctionGeneric<F>& f,
       const F& x0,
       const F& epsy,
       const F& h,
       const size_t Nmax
     ) throw (math::RootFindException)
{
    /*
     * The Halley's method is described in detail at:
     * https://en.wikipedia.org/wiki/Halley%27s_method
     */

    // The algorithm is implemented in __halleyCommon:
    return math::RootFind::__private::__halleyCommon<F>(
        f, f, f, x0, epsy, h, Nmax, false, true);
}
