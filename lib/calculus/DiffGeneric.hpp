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
 * An internal header file, it should not be included directly.
 * @headername{Diff.h}
 *
 * Declaration and implementation of functions within the namespace Diff
 * that perform numerical differentiation of continuous functions.
 */


#ifndef _MATH_DIFF_GENERIC_HPP_
#define _MATH_DIFF_GENERIC_HPP_

#include "../settings/calc_settings.h"
#include "util/NumericUtil.hpp"
#include "util/IFunctionGeneric.hpp"
#include "exception/CalculusException.hpp"
#include "exception/FunctionException.hpp"


namespace math
{


/**
 * A namespace with several implemented algorithms for
 * numerical differentiation.
 * 
 * @note Numerical differentiation is typically not very accurate.
 *       Make sure that the step 'h' is not too small to avoid
 *       large floating point rounding errors.
 */
namespace Diff
{


namespace __private
{

/*
 * Numerical differentiation of a function at the given
 * point 'x' using the forward method.
 *
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the slope
 * @param h - step size
 *
 * @return estimation of function's slope at the point 'x'
 *
 * @throw FunctionException if the function is not defined near 'x'
 */
template <typename F>
F __forwardDiff(
        const math::IFunctionGeneric<F>& f,
        const F& x,
        const F& h
      )
{
    /*
     *
     *   df(x)        f(x+h) - f(x)
     *  -------  ~=  ---------------
     *    dx                h
     *
     */

    return (f(x+h) - f(x)) / h;
}


/*
 * Numerical differentiation of a function at the given
 * point 'x' using the backward method.
 *
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the slope
 * @param h - step size
 *
 * @return estimation of function's slope at the point 'x'
 *
 * @throw FunctionException if the function is not defined near 'x'
 */
template <typename F>
F __backwardDiff(
        const math::IFunctionGeneric<F>& f,
        const F& x,
        const F& h
      )
{
    /*
     *
     *   df(x)        f(x) - f(x-h)
     *  -------  ~=  ---------------
     *    dx                h
     *
     */

    return (f(x) - f(x-h)) / h;
}


/*
 * Numerical differentiation of a function at the given
 * point 'x' using the central method.
 *
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the slope
 * @param h - step size
 *
 * @return estimation of function's slope at the point 'x'
 *
 * @throw FunctionException if the function is not defined near 'x'
 */
template <typename F>
F __centralDiff(
        const math::IFunctionGeneric<F>& f,
        const F& x,
        const F& h
      )
{
    /*
     *
     *   df(x)        f(x+h) - f(x-h)
     *  -------  ~=  -----------------
     *    dx               2 * h
     *
     */

    return (f(x+h) - f(x-h)) / ( static_cast<F>(2) * h );
}


/*
 * Numerical differentiation of a function at the given
 * point 'x' using the 5 point method.
 *
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the slope
 * @param h - step size
 *
 * @return estimation of function's slope at the point 'x'
 *
 * @throw FunctionException if the function is not defined near 'x'
 */
template <typename F>
F __5pointDiff(
        const math::IFunctionGeneric<F>& f,
        const F& x,
        const F& h
      )
{
    /*
     *
     *   df(x)        -f(x+2h) + 8*f(x+h) - 8*f(x-h) + f(x-2h)
     *  -------  ~=  ------------------------------------------
     *    dx                            12*h
     *
     */

    return ( -f(x + static_cast<F>(2) * h) +
              static_cast<F>(8) * f(x+h) -
              static_cast<F>(8) * f(x-h) +
              f(x - static_cast<F>(2) * h) ) /
           (static_cast<F>(12) * h);
}


}  // namespace __private



/**
 * @brief Supported algorithms for numerical differentiation
 */
struct EDiffMethod
{
    enum method
    {
        FORWARD,             /// Forward difference method
        BACKWARD,            /// Backward difference method
        CENTRAL,             /// Central difference method
        FIVE_POINT,          /// 5 point method
    };
};  // struct EDiffMethod



/**
 * Performs numerical differentiation of a function at the given
 * point 'x' using the selected algorithm.
 *
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the slope
 * @param h - step size (default: 0.0001)
 * @param method - one of the supported methods (default: CENTRAL)
 *
 * @return estimation of function's slope at the point 'x'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined near 'x'
 */
template <typename F>
F diff(
           const IFunctionGeneric<F>& f,
           const F& x,
           const F& h = static_cast<F>(DIFF_STEP_NUM) / static_cast<F>(DIFF_STEP_DEN),
            const EDiffMethod::method method = DIFF_DEFAULT_METHOD
          )
{
    // sanity check:
    if ( h < math::NumericUtil::getEPS<F>() )
    {
        throw math::CalculusException(math::CalculusException::INVALID_STEP);
    }

    try
    {
        F retVal;

        switch (method)
        {
            case math::Diff::EDiffMethod::FORWARD :
            {
                retVal = math::Diff::__private::__forwardDiff<F>(f, x, h);
                break;
            }

            case math::Diff::EDiffMethod::BACKWARD :
            {
                retVal = math::Diff::__private::__backwardDiff<F>(f, x, h);
                break;
            }

            case math::Diff::EDiffMethod::CENTRAL :
            {
                retVal = math::Diff::__private::__centralDiff<F>(f, x, h);
                break;
            }

            case math::Diff::EDiffMethod::FIVE_POINT :
            {
                retVal = math::Diff::__private::__5pointDiff<F>(f, x, h);
                break;
            }

            default :
                throw math::CalculusException(math::CalculusException::UNSUPPORTED_ALGORITHM);
        }  //switch

        return retVal;
    }  // try
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}



/**
 * Numerical 2nd order derivative of the given function,
 * using the 2nd order central method.
 * 
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the 2nd order derivative
 * @param h - step size (default: 0.0001)
 *
 * @return estimation of d2(x)/dx2 at the point 'x'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined near 'x'
 */
template <typename F>
F diff2(
           const IFunctionGeneric<F>& f,
           const F& x,
           const F& h = static_cast<F>(DIFF_STEP_NUM) / static_cast<F>(DIFF_STEP_DEN)
         )
{
    /*
     * 2nd order central method:
     * 
     *     2
     *    d f(x)        f(x+h) - 2*f(x) + f(x-h)
     *   --------  ~=  --------------------------
     *       2                     2
     *     dx                     h
     * 
     */

    const F h2 = h * h;

    // sanity check
    if ( h < math::NumericUtil::getEPS<F>() ||
         true == math::NumericUtil::isZero<F>(h2) )
    {
        throw math::CalculusException(math::CalculusException::INVALID_STEP);
    }
    
    try
    {
        return ( f(x+h) - static_cast<F>(2) * f(x) + f(x-h) ) / h2;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::CalculusException(math::CalculusException::UNDEFINED);
    }
}


}  // namespace Diff

}  // namespace math


#endif  // _MATH_DIFF_GENERIC_HPP_
