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
 * Implementation of functions within the namespace Diff that
 * perform numerical differentiation of continuous functions.
 */


// No "include "DiffGeneric.hpp" !!!
#include "util/NumericUtil.hpp"


// Implementation of "private" functions

namespace math
{

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
template <class T>
T __forwardDiff(
        const math::IFunctionGeneric<T>& f,
        const T& x,
        const T& h
      ) throw(math::FunctionException)
{
    /*
     *
     *   df(x)       f(x+h) - f(x)
     *  -------  ~  ---------------
     *    dx               h
     *
     */

    return (f.func(x+h) - f.func(x)) / h;
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
template <class T>
T __backwardDiff(
        const math::IFunctionGeneric<T>& f,
        const T& x,
        const T& h
      ) throw(math::FunctionException)
{
    /*
     *
     *   df(x)       f(x) - f(x-h)
     *  -------  ~  ---------------
     *    dx               h
     *
     */

    return (f.func(x) - f.func(x-h)) / h;
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
template <class T>
T __centralDiff(
        const math::IFunctionGeneric<T>& f,
        const T& x,
        const T& h
      ) throw(math::FunctionException)
{
    /*
     *
     *   df(x)       f(x+h) - f(x-h)
     *  -------  ~  -----------------
     *    dx              2 * h
     *
     */

    return (f.func(x+h) - f.func(x-h)) / ( static_cast<T>(2) * h );
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
template <class T>
T __5pointDiff(
        const math::IFunctionGeneric<T>& f,
        const T& x,
        const T& h
      ) throw(math::FunctionException)
{
    /*
     *
     *   df(x)       -f(x+2h) + 8*f(x+h) - 8*f(x-h) + f(x-2h)
     *  -------  ~  ------------------------------------------
     *    dx                           12*h
     *
     */

    return ( -f.func(x + static_cast<T>(2) * h) +
              static_cast<T>(8) * f.func(x+h) -
              static_cast<T>(8) * f.func(x-h) +
              f.func(x - static_cast<T>(2) * h) ) /
           (static_cast<T>(12) * h);
}


}  // namespace __private
}  // namespace Diff
}  // namespace math



/**
 * Performs numerical differentiation of a function at the given
 * point 'x' using the selected algorithm.
 *
 * @param f - instance of a class with the function to differentiate
 * @param x - point to estimate the slope
 * @param h - step size
 * @param method - one of the supported methods (default: CENTRAL)
 *
 * @return estimation of function's slope at the point 'x'
 *
 * @throw CalculusException if input arguments are invalid or the function is not defined near 'x'
 */
template <class T>
T math::Diff::diff(
        const math::IFunctionGeneric<T>& f,
        const T& x,
        const T& h,
        math::EDiffMethod::method method
      ) throw(math::CalculusException)
{
    // sanity check:
    if ( h < math::NumericUtil::getEPS<T>() )
    {
        throw math::CalculusException(math::CalculusException::INVALID_STEP);
    }

    try
    {
        T retVal;

        switch (method)
        {
            case math::EDiffMethod::FORWARD :
            {
                retVal = math::Diff::__private::__forwardDiff<T>(f, x, h);
                break;
            }

            case math::EDiffMethod::BACKWARD :
            {
                retVal = math::Diff::__private::__backwardDiff<T>(f, x, h);
                break;
            }

            case math::EDiffMethod::CENTRAL :
            {
                retVal = math::Diff::__private::__centralDiff<T>(f, x, h);
                break;
            }

            case math::EDiffMethod::FIVE_POINT :
            {
                retVal = math::Diff::__private::__5pointDiff<T>(f, x, h);
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
