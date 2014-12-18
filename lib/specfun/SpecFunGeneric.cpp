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
#include <cstddef>

#include "specfun/lanczos_coef.h"
#include "util/NumericUtil.hpp"
#include "util/math_constant.h"
#include "specfun/CtdFracGeneric.hpp"


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
     *                                  z-0.5
     *               -----+    (z+g-0.5)
     *   G(z) ~=   \/ 2*pi  * ---------------- * Lg(z)
     *                             z+g-0.5
     *                            e
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
     * Lanczos coefficients c0 .. c(N-1) and other necessary parameters are 
     * specified in the header "specfun/lanczos_coef.h".
     * Please see that header file for more details how they were obtained.
     * 
     * This function returns natural logarithm of G(x)
     * Applying the properties of logarithms, it can be further
     * simplified to:
     *
     *   ln( Gzx) ) ~= 0.5*ln(2*pi) + ln(Lg(z)) + (z-0.5)*ln(z+g-0.5) - (z+g-0.5)
     */


    // Chosen parameter 'g' casted to T:
    const T g = static_cast<T>(LANCZOS_G);

    // A handy macro to cast Lancos coefficients to T:
    #define GEN_CAST_COEF(COEF)             static_cast<T>(COEF),

     // An array with Lanczos coefficients casted to T:
    const T c[ NR_LANCZOS_COEF ] =
    {
        LANCZOS_COEF_ARRAY( GEN_CAST_COEF )
    };

    // The macro is not needed anymore
    #undef GEN_CAST_COEF
    
    // x + g - 0.5:
    const T term = x + g - static_cast<T>(1)/static_cast<T>(2);

    /*
     *              N-1
     *             -----
     *             \        c[i]
     * Lg = c[0] +  >   -----------
     *             /     x + i - 1
     *             -----
     *              i=1
     */
    T lg = c[0];
    for ( size_t i=1; i<NR_LANCZOS_COEF; ++i )
    {
        lg += c[i] / ( x + static_cast<T>(i) - static_cast<T>(1) );
    }

    /*
     * ln G(x) = ln sqrt(2*pi) + ln Lg + (x-0.5) * ln(x+g-0.5) - (x+g-0.5)
     */
    return static_cast<T>(MATH_CONST_LOG_SQRT_2_PI) +
           std::log(lg) + (x - static_cast<T>(1)/static_cast<T>(2)) *
           std::log(term) - term;

}

}  // namespace __private
}  // namespace SpecFun
}  // namespace math



/**
 * Gamma function.
 *
 * Complex numbers are also supported.
 *
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
 * @param x - input argument
 *
 * @return Gamma(x)
 *
 * @throw SpecFunexception if gamma function is not defined for the given 'x' (when 'x' is a negative integer)
 */
template <class T>
T math::SpecFun::gamma(const T& x) throw (math::SpecFunException)
{

    const T PI = static_cast<T>(MATH_CONST_PI);

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

        st = std::sin(PI*x);

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
         *                  pi
         *   G(x) = --------------------
         *           G(1-x) * sin(pi*x)
         */

        retVal = PI / (g_temp * st);
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
 * Beta function B(x,y) is defined as:
 *
 *            1
 *            /
 *            |  x-1        y-1          gamma(x) * gamma(y)
 *   B(x,y) = | t    * (1-t)     dt  =  ---------------------
 *            |                               gamma(x+y)
 *            /
 *            0
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


// implementation of an auxiliary "private" function:
namespace math {  namespace SpecFun {  namespace __private {

/*
 * Evaluates an incomplete gamma function. The exact kind of the returned
 * value depends on parameters 'upper' and 'reg'.
 *
 * @note both 'a' and 'x' must be strictly greater than 0
 *
 * @param a - first input argument
 * @param x - second input argument, the integration limit
 * @param upper - selects whether the upper (if 'true') or lower incomplete gamma functionis returned
 * @param reg - if 'true', the regularized gamma function is returned, i.e. divided by gamma(a)
 * @param tol - tolerance (default: 1e-6)
 */
template <class T>
T __incGamma(
                 const T& a,
                 const T& x,
                 bool upper,
                 bool reg,
                 const T& tol
               ) throw(math::SpecFunException)
{

    /*
     * When x > (1+a), the following continued fraction will be evaluated:
     *
     *                            1*(1-a)
     *   cf = (x-a+1) - ---------------------------
     *                                  2*(2-a)
     *                    (x-a+3) - ---------------
     *                               (x-a+5) - ...
     *
     * 'b_i' and 'a_i' are defined as follows:
     *
     * * a(x,i) = - i * (i-a)
     * * b(x,i) = x - a + 1 + 2*i
     */

    // Derive a class from ICtdFracFuncGeneric that
    // properly implements 'fa' and 'fb'.
	// 'a' will be the class's internal parameter
    class CtdF : public math::CtdFrac::ICtdFracFuncGeneric<T>
    {

    private:
        const T m_a;      // parameter 'a'

    public:
        /*
         * Constructor, sets value of 'm_a'
         *
         * @param a - parameter 'a' from the definition of incomplete gamma function
         */
        CtdF(const T& a) : m_a(a)
        {
            // nothing else to do
        }

        // a(x,i) = -i * (i-a)
        T fa(const T& x, size_t i) const throw (math::FunctionException)
        {
        	(void) x;
            const T f = static_cast<T>(i);
            return -f * (f - this->m_a);
        }

        // b(x,i) = x - a + 1 + 2*i
        T fb(const T& x, size_t i) const throw (math::FunctionException)
        {
            return x - this->m_a + static_cast<T>(1) + static_cast<T>(2*i);
        }
    };  // class CtdF


    // An instance of Ctdf that implements 'fa' and 'fb':
    CtdF coef(a);

    // sanity check:
    if ( a < math::NumericUtil::getEPS<T>() ||
         x < math::NumericUtil::getEPS<T>() )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    /*
     * The algorithm for numerical approximation of the incomplete gamma function
     * as proposed in:
     *
     *   William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
     *   Numerical Recipes, The Art of Scientific Computing, 3rd Edition,
     *   Cambridge University Press, 2007
     *
     * When x > (a+1), the upper gamma function can be evaluated as
     *
     *                -x    a
     *               e   * x
     *  G(a,x) ~= --------------
     *               cf(a,x)
     *
     * where 'cf(a,x) is the continued fraction defined above, its coefficients
     * 'a_i' and 'b_i' are implented in 'coef'.
     *
     * When x < (a+1), it is more convenient to apply the following Taylor series
     * that evalutes the lower incomplete gamma function:
     *                         inf
     *                        -----
     *              -x   a    \        G(a)       i
     *   g(a,x) ~= e   *x  *   >    ---------- * x
     *                        /      G(a+i+i)
     *                        -----
     *                         i=0
     *
     * Applying the following property of the gamma function:
     *
     *   G(a+1) = a * G(a)
     *
     * The Taylor series above can be further simplified to:
     *
     *                         inf
     *                        -----              i
     *              -x   a    \                 x
     *   g(a,x) ~= e   *x  *   >    -------------------------
     *                        /      a * (a+1) * ... * (a+i)
     *                        -----
     *                         i=0
     *
     * Once either a lower or an upper incomplete gamma function is evaluted,
     * the other value may be quickly obtained by applying the following
     * property of the incomplete gamma function:
     *
     *   G(a,x) + g(a,x) = G(a)
     *
     * A value of a regularized incomplete gamma function is obtained
     * by dividing g(a,x) or G(a,x) by G(a).
     */

    // This factor is common to both algorithms described above:
    T ginc = std::exp(-x) * std::pow(x, a);

    if ( x > (a + static_cast<T>(1)) )
    {
        /*
         * x > (a + 1)
         *
         * In this case evaluate the upper gamma function as described above.
         */

    	const T G = ( true==upper && false==reg ? static_cast<T>(0) : math::SpecFun::gamma<T>(a) );

        ginc /=  math::CtdFrac::ctdFrac<T>(coef, x);

        /*
         * Apply properties of the incomplete gamma function
         * if anything else except a generalized upper incomplete
         * gamma function is desired.
         */
        if ( false == upper )
        {
            ginc = G - ginc;
        }

        if ( true == reg )
        {
            ginc /= G;
        }
    }
    else
    {
        /*
         * x < (a + 1)
         *
         * In this case evaluate the lower gamma function as described above.
         */
    	const T G = ( false==upper && false==reg ? static_cast<T>(0) : math::SpecFun::gamma<T>(a) );

        // Initial term of the Taylor series at i=0:
        ginc /= a;
        T term = ginc;

        // Proceed the Taylor series for i=1, 2, 3... until it converges:
        for ( size_t i = 1; false==math::NumericUtil::isZero<T>(term, tol); ++i )
        {
            term *= x / (static_cast<T>(i) + a);
            ginc += term;
        }

        /*
         * Apply properties of the incomplete gamma function
         * if anything else except a generalized lower incomplete
         * gamma function is desired.
         */
        if ( true == upper )
        {
            ginc = G - ginc;
        }

        if ( true == reg )
        {
            ginc /= G;
        }
    }

    return ginc;
}

}}}  // namespace math::SpecFun::__private



/**
 * Generalized upper incomplete gamma function, defined as:
 *
 *           inf
 *            /
 *            |  a-1    -t
 *   G(a,x) = | t    * e   dt
 *            |
 *            /
 *            x
 *
 * @note 'a' and 'x' must be greater than 0
 *
 * @param a - first input argument
 * @param x - second input argument, the lower integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return G(a,x)
 *
 * @throw SpecFunException if 'a' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incGammaUpper(
               const T& a,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incGamma<T>(a, x, true, false, tol);
}


/**
 * Generalized lower incomplete gamma function, defined as:
 *
 *            x
 *            /
 *            |  a-1    -t
 *   g(a,x) = | t    * e   dt
 *            |
 *            /
 *            0
 *
 * @note 'a' and 'x' must be greater than 0
 *
 * @param a - first input argument
 * @param x - second input argument, the upper integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return g(a,x)
 *
 * @throw SpecFunException if 'a' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incGammaLower(
               const T& a,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incGamma<T>(a, x, false, false, tol);
}


/**
 * Regularized upper incomplete gamma function, defined as:
 *
 *                  inf
 *                   /
 *              1    |  a-1    -t
 *   Q(a,x) = ------ | t    * e   dt
 *             G(a)  |
 *                   /
 *                   x
 *
 * @note 'a' and 'x' must be greater than 0
 *
 * @param a - first input argument
 * @param x - second input argument, the lower integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return Q(a,x)
 *
 * @throw SpecFunException if 'a' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incGammaUpperReg(
               const T& a,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incGamma<T>(a, x, true, true, tol);
}


/**
 * Regularized lower incomplete gamma function, defined as:
 *
 *                   x
 *                   /
 *              1    |  a-1    -t
 *   P(a,x) = ------ | t    * e   dt
 *             G(a)  |
 *                   /
 *                   0
 *
 * @note 'a' and 'x' must be greater than 0
 *
 * @param a - first input argument
 * @param x - second input argument, the upper integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return P(a,x)
 *
 * @throw SpecFunException if 'a' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incGammaLowerReg(
               const T& a,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incGamma<T>(a, x, false, true, tol);
}


/**
 * Error function, defined as:
 *
 *                        x
 *                2       /  -t^2
 *   erf(x) = ----------  | e    dt
 *             sqrt(pi)   /
 *                        0
 *
 * @param x - input argument
 * @param tol - tolerance (default: 1e-6)
 *
 * @return erf(x)
 */
template <class T>
T math::SpecFun::erf(const T& x, const T& tol)
{
    /*
     * The definite integral could be calculated numerically, however
     * more efficients methods exist. The exponential can be expanded
     * into a Taylor series and each term can be integrated separately.
     * As evident from
     * https://en.wikipedia.org/wiki/Error_function#Taylor_series
     * the error function can thus be expanded into the following
     * Taylor series:
     *
     *                        +-      3     5      7      9        -+
     *                 2      |      x     x      x      x          |
     *   erf(x) ~= ---------- | x - --- + ---- - ---- + ----- - ... |
     *              sqrt(pi)  |      3     10     42     216        |
     *                        +-                                   -+
     *
     *                         inf               i
     *                        -----            -----     2
     *                 2      \        x       |   |   -x
     *   erf(x) ~= ----------  >   --------- * |   | -------
     *              sqrt(pi)  /     2*i + 1    |   |    j
     *                        -----            |   |
     *                         i=0              j=1
     *
     * The error function can also be obtained via the incomplete
     * gamma function:
     *
     *   erf(x) = P(1/2, x^2)   for x>0
     *
     * If x is negative, apply the following property:
     *
     *   erf(-x) = -erf(x)
     */


    const T t = x * x;

    // If 't' is very close to 0, the incomplete gamma function
    // might throw an exception, to prevent this,
    // this is handled separately.
    if ( true == math::NumericUtil::isZero<T>(t) )
    {
        return static_cast<T>(0);
    }

    T retVal = math::SpecFun::incGammaLowerReg<T>(static_cast<T>(1)/static_cast<T>(2), t, tol);

    if ( x < static_cast<T>(0) )
    {
        retVal = -retVal;
    }

    return retVal;
}


/**
 * Complementary error function, defined as:
 *
 *                        inf
 *                 2       /  -t^2
 *   erfc(x) = ----------  | e    dt  =  1 - erf(x)
 *              sqrt(pi)   /
 *                         x
 *
 * @param x - input argument
 * @param tol - tolerance (default: 1e-6)
 *
 * @return erfc(x)
 */
template <class T>
T math::SpecFun::erfc (const T& x, const T& tol)
{
    /*
     * The complementary error function can also be obtained via the incomplete
     * gamma function:
     *
     *   erfc(x) = Q(1/2, x^2)   for x>0
     *
     * If x is negative, apply the following property:
     *
     *   erfc(-x) = 2 - erfc(x)
     */


    const T t = x * x;

    // If 't' is very close to 0, the incomplete gamma function
    // might throw an exception, to prevent this,
    // this is handled separately.
    if ( true == math::NumericUtil::isZero<T>(t) )
    {
        return static_cast<T>(1);
    }

    T retVal = math::SpecFun::incGammaUpperReg<T>(static_cast<T>(1)/static_cast<T>(2), t, tol);

    if ( x < static_cast<T>(0) )
    {
        retVal = static_cast<T>(2) - retVal;
    }

    return retVal;
}
