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
 * Implementation of special functions in the namespace SpecFun.
 */


// no #include "SpecFunGeneric.hpp" !!!
#include <cmath>
#include <complex>

#include "util/NumericUtil.hpp"


// Implementation of "private" functions

namespace math
{

namespace SpecFun
{

namespace __private
{


/*
 * Checks whether the argument 'x' is located in the
 * left half plane, i.e. its real part is less or equal to 0.
 *
 * @param x - argument to check
 *
 * @return 'true' if Re(x)<=0, 'false' otherwise
 */
template <class T>
bool __leftHalfPlane(const T& x)
{
    // left half plane conveniently includes anything less than EPS
    return ( x < math::NumericUtil::getEPS<T>() );
}


/*
 * Partial "specialization" of __leftHalfPlane for complex numbers.
 *
 * Checks whether the argument 'x' is located in the
 * left half plane, i.e. its real part is less or equal to 0.
 *
 * @param x - argument to check
 *
 * @return 'true' if Re(x)<=0, 'false' otherwise
 */
template <class T>
bool __leftHalfPlane(const std::complex<T>& x)
{
    return ( std::real(x) < math::NumericUtil::getEPS<T>() );
}


/*
 * Checks whether 'x' is located in the so called "middle segment",
 * i.e. its real part is greater than 0 and less than 1.
 *
 * @param x - argument to check
 *
 * @return 'true' if 0<Re(x)<1, 'false' otherwise
 */
template <class T>
bool __midSegment(const T& x)
{
    return ( x > math::NumericUtil::getEPS<T>() &&
             x < static_cast<T>(1) );
}


/*
 * Partial "specialization" of __midSegment for complex numbers.
 *
 * Checks whether 'x' is located in the so called "middle segment",
 * i.e. its real part is greater than 0 and less than 1.
 *
 * @param x - argument to check
 *
 * @return 'true' if 0<Re(x)<1, 'false' otherwise
 */
template <class T>
bool __midSegment(const std::complex<T>& x)
{
    return ( std::real(x) > math::NumericUtil::getEPS<T>() &&
             std::real(x) < static_cast<T>(1) );
}


/*
 * A convenient function that efficiently calculates logarithm of
 * gamma(x). The function expects that argument's real part is
 * greater than 1. If this is not the case, several reflection
 * methods exist.
 *
 * @param x - input argument (its real part should be greater than 0)
 *
 * @return ln( Gamma(x) )
 */
template <class T>
T __lnGamma(const T& x)
{
    /*
     * Gamma function G(t) is defined as:
     *
     *         inf
     *          /
     *          |  t-1    -x
     *   G(t) = | x    * e   dx
     *          |
     *          /
     *          0
     *
     * There are several approaches to calculate its approximation,
     * one of the most convenient ones is the Lanczos approximation
     * https://en.wikipedia.org/wiki/Lanczos_approximation
     *
     *                                   z-0.5
     *             _  -----+    (z+g-0.5)
     *   G(z) ~=    \/ 2*pi  * ---------------- * Lg(z)
     *                              z+g-0.5
     *                             e
     *
     * where Lg(z) is defined as:
     *
     *                 N-1
     *                -----
     *                \      ck
     *   Lg(z) = c0 +  >   -------
     *                /     z+k-1
     *                -----
     *                 k=1
     *
     * Coefficients c0 .. c(N-1) depend on chosen parameters g and N.
     * They can be calculated by the Paul Godfrey's algorithm as described at:
     * http://mrob.com/pub/ries/lanczos-gamma.html
     * Using Maxima it is possible to calculate them in an arbitrary precision.
     * The code below is slightly modified code from the link above:
     *
       load("diag");

       Dc(n) := diag(makelist(2*double_factorial(2*k-1),k,0,n));

       cmatrix_element[row,col]:=
          if is(col>row) then 0
          elseif row=1 and col=1 then 1/2
          else (-1)^(row+col)*4^(col-1)*(row-1)*(row+col-3)!/(row-col)!/(2*col-2)!;

       C(n) := genmatrix(cmatrix_element, n+1);

       f(g,n):=sqrt(2)*(%e/(2*(n+g)+1))^(n+1/2);

       Dr(k) := diag(append([1],makelist(-(2*n+2)!/(2*n!*(n+1)!),n,0,k-1)));

       bmatrix_element[row,col] :=
           if row = 1 then 1
           elseif is(row > col) then 0
           else (-1)^(col-row)*binomial(col+row-3,2*row-3);

       B(k) := genmatrix(bmatrix_element,k+1);

       lanczos_coeff(g, n) :=
           block([M : (Dr(n) . B(n)) . (C(n) . Dc(n)),
                   f : transpose(matrix(makelist(f(g,k), k, 0, n)))],
               (M . f));

       tlg(g,n) := lanczos_coeff(g,n-1)*exp(g)/sqrt(2*%pi);
     *
     * When very high accuracy is not necessary, Boost recommends the following values for
     * g and N:
     *   g = 1.428456135094165802001953125
     *   N = 6
     *
     * More details at:
     * http://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/lanczos.html
     *
     * For the 25 digit precision, the Maxima code above returns the following values of
     * Lanczos coefficients:
     (%i1)  fpprec : 25$
     (%i2)  bfloat(tlg(1.428456135094165802001953125, 6));

              [  1.000000123880027543803602b0  ]
              [  9.727644985297789802780476b−1 ]
     (%o2)    [  9.79090078038154694104146b−3  ]
              [ −6.669128561042791409832734b−3 ]
              [  3.802685303772964895033826b−3 ]
              [ −1.117185821247582043028556b−3 ]
     *
     * This function returns natural logarithm of G(x)
     * Applying the properties of logarithms, it can be further
     * simplified to:
     *
     *   ln( Gzx) ) ~= 0.5*ln(2*pi) + ln(Lg(z)) + (z-0.5)*ln(z+g-0.5) - (z+g-0.5)     *
     */


    /*
     * High precision (25 digits) numerical constants
     * ln(sqrt(2*pi)),
     * obtained by Maxima:
     *
     (%i1)  fpprec : 25$
     (%i2)  bfloat(log(sqrt(2*%pi)));
     (%o2)  9.189385332046727417803297b−1
    */

	// Chosen number of Lanczos coefficients:
    #define N_LANCZOS    ( 6 )

	// Chosen parameter 'g' casted to T:
    const T g = static_cast<T>(1.428456135094165802001953125L);

     // An array with Lanczos coefficients casted to T:
    const T c[ N_LANCZOS ] =
    {
        static_cast<T>(1.000000123880027543803602L),
        static_cast<T>(0.9727644985297789802780476L),
        static_cast<T>(0.00979090078038154694104146L),
        -static_cast<T>(0.006669128561042791409832734L),
        static_cast<T>(0.003802685303772964895033826L),
        -static_cast<T>(0.001117185821247582043028556L)
    };


    // x + g - 0.5:
    const T term = x + g - static_cast<T>(1)/static_cast<T>(2);

    /*
     *              N-1
     *             -----
     *             \       c[i]
     * Lg = c[0] +  >   ----------
     *             /     x + i -1
     *             -----
     *              i=1
     */
    T lg = c[0];
    for ( size_t i=1; i<N_LANCZOS; ++i )
    {
        lg += c[i] / ( x + static_cast<T>(i) - static_cast<T>(1) );
    }

    /*
     * ln G(x) = ln sqrt(2*pi) + ln Lg + (x-0.5) * ln(x+g-0.5) - (x+g-0.5)
     */
    return static_cast<T>(0.9189385332046727417803297L) +
           std::log(lg) + (x - static_cast<T>(1)/static_cast<T>(2)) *
           std::log(term) - term;


    // The macro is not needed anymore and will be undefined
    #undef N_LANCZOS
}

}  // namespace __private
}  // namespace SpecFun
}  // namespace math



/**
 * Gamma function.
 *
 * Complex numbers are also supported.
 *
 * @param x - input argument
 *
 * @return Gamma(x)
 *
 * @throw SpecFunexception if gamma function is not defined for the given 'x' (when 'x' is a negative integer)
 */
template <class T>
T math::SpecFun::gamma(const T& x) throw (math::SpecFunException)
{
    /*
     * Gamma function G(t) is defined as:
     *
     *         inf
     *          /
     *          |  t-1    -x
     *   G(t) = | x    * e   dx
     *          |
     *          /
     *          0
     *
     */

    /*
     * High precision (25 digits) approximation for the
     * numerical constant pi, as obtained by Maxima:
     (%i1)  fpprec : 25$
     (%i2)  bfloat(%pi);
     (%o2)  3.141592653589793238462643b0
     */
    const T pi = static_cast<T>(3.141592653589793238462643L);

    // sin(pi*x), only applicable when Re(x)<0
    T st;

    /*
     * __lnGamma expects real part of 'x' to be greater than 1,
     * otherwise the algorithm may not be accurate enough or may
     * return NaN. When Re(x) is smaller than 1, a reflection method
     * is applied. In this case another value of the input argument
     * ('z') to __lnGamma must be obtained first as implemented below:
     */
    T z;
    if ( true == math::SpecFun::__private::__leftHalfPlane(x) )
    {
        /*
         * Re(x) < 0:
         *
         * The following reflection method will be applied:
         *
         *                   pi
         *   G(1-x) = ------------------
         *             G(x) * sin(pi*x)
         *
         *  In this case the real part of (1-x) will always be greater than 1.
         */

        st = std::sin(pi*x);

        // if sin(pi*x) equals 0, the gamma function is undefined:
        if ( true == math::NumericUtil::isZero<T>(st) )
        {
            throw math::SpecFunException(math::SpecFunException::UNDEFINED);
        }

        z = static_cast<T>(1) - x;

    }
    else if ( true == math::SpecFun::__private::__midSegment(x) )
    {
        /*
         * 0 < Re(x) < 1:
         *
         * This reflection method will be applied:
         *
         *   G(x+1) = x * G(x)
         *
         * In this case the real part of (x+1) will always be greater than 1.
         */

        z = x + static_cast<T>(1);
    }
    else
    {
        /*
         * Re(x) > 1:
         *
         * In this case __lnGamma(x) is always accurate enough
         * and no reflection method is necessary.
         */

        z = x;
    }

    /*
     * Now the real part of 'z' is guaranteed to be greater than 1.
     * ln G(z) can be safely calculated and exponentiated:
     */
    T g_temp = std::exp( math::SpecFun::__private::__lnGamma<T>(z) );

    /*
     * Finally obtain the actual value of the gamma function,
     * depending on the applied reflection method:
     */
    T retVal;
    if ( true == math::SpecFun::__private::__leftHalfPlane(x) )
    {
        /*
         * Re(x) < 0:
         *
         * From
         *
         *                   pi
         *   G(1-x) = ------------------
         *             G(x) * sin(pi*x)
         *
         * G(x) can be obtained as:
         *
         *                   pi
         *   G(x) = --------------------
         *           G(1-x) * sin(pi*x)
         */

        retVal = pi / (g_temp * st);
    }
    else if ( true == math::SpecFun::__private::__midSegment(x) )
    {
        /*
         * 0 < Re(x) < 1:
         *
         * From
         *
         *   G(x+1) = x * G(x)
         *
         * G(x) can be obtained as:
         *
         *           G(x+1)
         *   G(x) = --------
         *             x
         */

        retVal = g_temp / x;
    }
    else
    {
       /*
        * Re(x) > 1:
        *
        * In this case nothing else is necessary to do and
        * 'g_temp' can be returned.
        */

       retVal = g_temp;
    }

    return retVal;
}


/**
 * Beta function.
 *
 * Complex numbers are also supported.
 *
 * @param x - first input argument
 * @param y - second input argument
 *
 * @return Beta(x, y)
 *
 * @throw SpecFunException if function is not defined for any 'x' or 'y' (if any input is a negative integer)
 */
template <class T>
T math::SpecFun::beta(const T& x, const T& y) throw (math::SpecFunException)
{
    /*
     * Beta function B(x,y) is defined as:
     *
     *            1
     *            /
     *            |  x-1        y-1
     *   B(x,y) = | t    * (1-t)     dt
     *            |
     *            /
     *            0
     *
     * It can be shown that B(x,y) can be expressed with gamma
     * functions:
     *
     *             G(x) * G(y)
     *   B(x,y) = -------------
     *               G(x+y)
     */

    return math::SpecFun::gamma<T>(x) *
           math::SpecFun::gamma<T>(y) /
           math::SpecFun::gamma<T>(x+y);
}
