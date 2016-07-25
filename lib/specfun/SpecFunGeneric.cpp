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


/*
 * The following two publications are often referred in this file:
 *
 * - [Numerical Recipes] or shorter [NR]:
 *      William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
 *      Numerical Recipes, The Art of Scientific Computing, 3rd Edition,
 *      Cambridge University Press, 2007
 *      http://www.nr.com
 *
 * - [Abramowitz & Stegun] or shorter [AS]:
 *      Milton Abramowitz, Irene A. Stegun
 *      Handbook of Mathematical Functions with Formulas, Graphs and Mathematical Tables
 *      National Bureau of Standards,
 *      Applied Mathematics Series 55
 *      Issued in June 1964
 *
 *      The book is copyright free and can be downloaded from:
 *      http://www.cs.bham.ac.uk/~aps/research/projects/as/book.php
 */


// Implementation of "private" functions
namespace math {  namespace SpecFun {  namespace __private
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
inline bool __leftHalfPlane(const T& x)
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
inline bool __leftHalfPlane(const std::complex<T>& x)
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
inline bool __midSegment(const T& x)
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
inline bool __midSegment(const std::complex<T>& x)
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

    // A handy macro to cast Lanczos coefficients to T:
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
    for ( std::size_t i=1; i<NR_LANCZOS_COEF; ++i )
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

}}}  // namespace math::SpecFun::__private



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

    const T gxy = math::SpecFun::gamma<T>(x + y);

    // handle a very unlikely occassion that G(x,y) gets very close to 0:
    if ( true == math::NumericUtil::isZero<T>(gxy) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    return math::SpecFun::gamma<T>(x) *
           math::SpecFun::gamma<T>(y) /
           gxy;
}


// implementation of auxiliary "private" functions:
namespace math {  namespace SpecFun {  namespace __private {

/*
 * Extension of ICtdFracFuncGeneric, that evaluates coefficients of
 * the following continued fraction, necessary to evaluate the
 * incomplete gamma function:
 * 
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
 * 
 * 'a' and 'x' are the class' internal parameters
 */
template< class T>
class __CtdFIncGamma : public math::CtdFrac::ICtdFracFuncGeneric<T>
{

private:
    const T m_a;      // parameter 'a'
    const T m_x;      // input argument to the incomplete gamma function ('x') 

public:
    /*
     * Constructor, sets value of 'm_a'
     *
     * @param a - parameter 'a' from the definition of incomplete gamma function
     * @param x - input argument to the incomplete gamma function
     */
    __CtdFIncGamma(const T& a, const T& x) : m_a(a), m_x(x)
    {
        // nothing else to do
    }

    // a(x,i) = -i * (i-a)
    T fa(const std::size_t i) const throw (math::FunctionException)
    {
        const T f = static_cast<T>(i);
        return -f * (f - this->m_a);
    }

    // b(x,i) = x - a + 1 + 2*i
    T fb(const std::size_t i) const throw (math::FunctionException)
    {
        return this->m_x - this->m_a + static_cast<T>(1) + 
                    static_cast<T>(2) * static_cast<T>(i);
    }
};  // class __CtdFIncGamma


/*
 * Extension of ICtdFracFuncGeneric, that evaluates coefficients of
 * the following continued fraction, necessary to evaluate the
 * incomplete beta function:
 *
 *                     a1
 *   cf = 1 + ---------------------
 *                       a2
 *             1 + ---------------
 *                         a3
 *                  1 + ---------
 *                       1 + ...
 *
 * 'b_i' and 'a_i' are defined as follows:
 *
 * * b(x,i) = 1
 * *
 * *           /    (a+m) * (a+b+m) * x
 * *           | - -------------------------   <== i = 2*m+1
 * *           /    (a+2*m) * (a + 2*m + 1)
 * * a(x,i) = {
 * *           \       m * (b-m) * x
 * *           |   ---------------------       <== i = 2*m
 * *           \    (a+2*m-1) * (a+2*m)
 * 
 * 'a', 'b' and 'c' are the class' internal parameters
 */
template <class T>
class __CtdFIncBeta : public math::CtdFrac::ICtdFracFuncGeneric<T>
{
private:

    const T m_a;      // parameter 'a'
    const T m_b;      // parameter 'b'
    const T m_x;      // input argument to the incomplete beta function ('x')

public:
    /*
     * Constructor, sets values of 'm_a' and 'm_b'
     *
     * @param a - parameter 'a' from the definition of the incomplete beta function
     * @param b - parameter 'b' from the definition of the incomplete beta function
     * @param x - input argument to the incomplete beta function
     */
    __CtdFIncBeta(const T& a, const T& b, const T& x) : m_a(a), m_b(b), m_x(x)
    {
        // nothing else to do
    }

    /*
     * a(x,i) = -(a+m)*(a+b+m)*x / ( (a+2m)*(a+2m+1) )   when i=2*m+1
     * a(x,i) = m*(b-m)*x / ( (a+2m-1)*(a+2m) )          when i=2*m
     */
    T fa(const std::size_t i) const throw(math::FunctionException)
    {
        const std::size_t m = i / 2;
        T ai = static_cast<T>(0);

        if ( 1 == i%2 )
        {
            // 'i' is odd, i.e i=2*m+1
            ai = -(this->m_a + static_cast<T>(m)) * (this->m_a + this->m_b + static_cast<T>(m)) * this->m_x /
                  ( (this->m_a + static_cast<T>(i) - static_cast<T>(1) ) * (this->m_a + static_cast<T>(i) ) );
        }
        else
        {
            // 'i' is even, i.e i=2*m
            ai = static_cast<T>(m) * (this->m_b - static_cast<T>(m)) * this->m_x /
                 ( (this->m_a + static_cast<T>(i) - static_cast<T>(1)) * (this->m_a + static_cast<T>(i)) );
        }

        return ai;
    }

    // b(x,i) = 1
    T fb(const std::size_t i) const throw(math::FunctionException)
    {
        (void) i;
        return static_cast<T>(1);
    }
};  // class __CtdFIncBeta
    
    
/*
 * Evaluates an incomplete gamma function. The exact kind of the returned
 * value depends on parameters 'upper' and 'reg'.
 *
 * @note both 'a' and 'x' must be strictly greater than 0
 *
 * @param a - parameter of the incomplete gamma function
 * @param x - second input argument, the integration limit
 * @param upper - should the upper (if 'true') or the lower (if 'false) inc. gamma function be returned
 * @param reg - if 'true', the regularized gamma function is returned, i.e. divided by gamma(a)
 * @param tol - tolerance (default: 1e-6)
 */
template <class T>
T __incGamma(
                 const T& a,
                 const T& x,
                 const bool upper,
                 const bool reg,
                 const T& tol
               ) throw(math::SpecFunException)
{
    // An instance of __CtdFIncGamma that implements 'fa' and 'fb':
    const math::SpecFun::__private::__CtdFIncGamma<T> coef(a, x);

    // sanity check:
    if ( a < math::NumericUtil::getEPS<T>() ||
         x < math::NumericUtil::getEPS<T>() )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    /*
     * The algorithm for numerical approximation of the incomplete gamma function
     * as proposed by [Numerical Recipes], section 6.2:
     *
     * When x > (a+1), the upper gamma function can be evaluated as
     *
     *                 -x    a
     *                e   * x
     *   G(a,x) ~= --------------
     *                cf(a,x)
     *
     * where 'cf(a,x) is the continued fraction defined above, its coefficients
     * 'a_i' and 'b_i' are implented in 'coef'.
     *
     * When x < (a+1), it is more convenient to apply the following Taylor series
     * that evaluates the lower incomplete gamma function:
     * 
     *                          inf
     *                         -----
     *              -x    a    \        G(a)       i
     *   g(a,x) ~= e   * x  *   >    ---------- * x
     *                         /      G(a+1+i)
     *                         -----
     *                          i=0
     *
     * Applying the following property of the gamma function:
     *
     *   G(a+1) = a * G(a)
     *
     * The Taylor series above can be further simplified to:
     *
     *                          inf
     *                         -----              i
     *              -x    a    \                 x
     *   g(a,x) ~= e   * x  *   >    -------------------------
     *                         /      a * (a+1) * ... * (a+i)
     *                         -----
     *                          i=0
     *
     * Once either a lower or an upper incomplete gamma function is evaluated,
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

        const T G = ( true==upper && false==reg ? 
            static_cast<T>(0) : 
            math::SpecFun::gamma<T>(a) );

        ginc /=  math::CtdFrac::ctdFrac<T>(coef);

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
            // Note: if a>0, gamma(a) is always greater than 0
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
        const T G = ( false==upper && false==reg ? 
            static_cast<T>(0) : 
            math::SpecFun::gamma<T>(a) );

        // Initial term of the Taylor series at i=0:
        ginc /= a;
        T term = ginc;

        // Proceed the Taylor series for i = 1, 2, 3... until it converges:
        T at = a;
        std::size_t i = 1;
        for ( i = 1; 
              false==math::NumericUtil::isZero<T>(term, tol) && i<SPECFUN_MAX_ITER; 
              ++i )
        {
            at += static_cast<T>(1);
            term *= x / at;
            ginc += term;
        }

        // has the series converged?
        if ( i >= SPECFUN_MAX_ITER )
        {
            throw math::SpecFunException(math::SpecFunException::NO_CONVERGENCE);
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
            // Note: if a>0, gamma(a) is always greater than 0
            ginc /= G;
        }
    }

    return ginc;
}


/*
 * Partial "specialization" of __incGamma for complex numbers.
 * 
 * Evaluates an incomplete gamma function. The exact kind of the returned
 * value depends on parameters 'upper' and 'reg'.
 *
 * @note unlike at real numbers, incomplete gamma function is defined
 *       virtually everywhere on the complex plane except a = negative integer
 *
 * @param a - parameter of the incomplete gamma function
 * @param x - second input argument, the integration limit
 * @param upper - should the upper (if 'true') or the lower (if 'false) inc. gamma function be returned
 * @param reg - if 'true', the regularized gamma function is returned, i.e. divided by gamma(a)
 * @param tol - tolerance (default: 1e-6)
 */
template <class T>
std::complex<T> __incGamma(
                 const std::complex<T>& a,
                 const std::complex<T>& x,
                 const bool upper,
                 const bool reg,
                 const std::complex<T>& tol
               ) throw(math::SpecFunException)
{
    /*
     * The lower incomplete gamma function can be expanded into:
     *
     *                          inf
     *                         -----              i
     *              -x    a    \                 x
     *   g(a,x) ~= e   * x  *   >    -------------------------
     *                         /      a * (a+1) * ... * (a+i)
     *                         -----
     *                          i=0
     *
     * Once either a lower incomplete gamma function is evaluated,
     * the upper one may be quickly obtained by applying the following
     * property of the incomplete gamma function:
     *
     *   G(a,x) + g(a,x) = G(a)
     *
     * A value of a regularized incomplete gamma function is obtained
     * by dividing g(a,x) or G(a,x) by G(a).
     */

    // If 'a' equals zero, division by zero would occur
    if ( true == math::NumericUtil::isZero<std::complex<T> >(a) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    // G = gamma(a) or (0,0) when the value is not used
    const std::complex<T> G = ( true==reg || true==upper ? 
               math::SpecFun::gamma<std::complex<T> >(a) : 
               std::complex<T>(static_cast<T>(0)) );

    std::complex<T> at = a;

    // The first term of the series
    std::complex<T> ginc = std::pow(x, a) * std::exp(-x) / a;
    std::complex<T> term = ginc;

    // proceed the series until it converges
    std::size_t i = 0;
    for ( i = 0;
          false == math::NumericUtil::isZero<std::complex<T> >(term, tol) &&
            i < SPECFUN_MAX_ITER;
          ++i )
    {
        at += static_cast<T>(1);

        // if 'a' is a negative integer, sooner or later this exception will be thrown
        if ( true == math::NumericUtil::isZero<std::complex<T> >(at) )
        {
            throw math::SpecFunException(math::SpecFunException::UNDEFINED);
        }
        
        term *= x / at;
        ginc += term;
    }

    // has the series converged?
    if ( i >= SPECFUN_MAX_ITER )
    {
        throw math::SpecFunException(math::SpecFunException::NO_CONVERGENCE);
    }

    /*
     * Apply properties of the incomplete gamma function
     * if anything else except a non-regularized lower incomplete
     * gamma function is desired.
     */
    if ( true == upper )
    {
        ginc = G - ginc;
    }

    if ( true == reg )
    {
        // very unlikely but check it anyway
        if ( true == math::NumericUtil::isZero<std::complex<T> >(G) )
        {
            throw math::SpecFunException(math::SpecFunException::UNDEFINED);
        }

        ginc /= G;
    }

    return ginc;
}


/*
 * Evaluates an incomplete beta function. The exact kind of the returned
 * value depends on parameters 'upper' and 'reg'.
 *
 * @note both 'a', 'b' and 'x' must be strictly greater than 0,
 *       x must be greater than 0 and less than 1
 *
 * @param a - parameter 'a' of the beta function
 * @param b - parameter 'b' of the beta function
 * @param x - the integration limit
 * @param lower - should the lower (if 'true') or the upper (if 'false) inc. beta function be returned
 * @param reg - if 'true', the regularized beta function is returned, i.e. divided by beta(a, b)
 * @param tol - tolerance (default: 1e-6)
 */
template <class T>
T __incBeta(
          const T& a,
          const T& b,
          const T& x,
          const bool lower,
          const bool reg,
          const T& tol
        ) throw (math::SpecFunException)
{
    // sanity check
    if ( a < math::NumericUtil::getEPS<T>() ||
         b < math::NumericUtil::getEPS<T>() ||
         x < static_cast<T>(0) || x > static_cast<T>(1) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    /*
     * The algorithm is described in detail in [Numerical Recipes], section 6.4
     *
     * If x < (a+1)/(a+b+2), the incomplete beta function will be evaluated as:
     *
     *                  a        b
     *                 x  * (1-x)
     *   Bx(a,b) ~= ------------------
     *                a * cf(x,a,b)
     *
     * where cf(x,a,b) is the continued fraction defined above, its
     * coefficients are implemented in 'coef'.
     *
     * If x > (a+1)/(a+b+2), the algorithm will converge faster if the following
     * property of the incomplete beta function is applied:
     *
     *   B (a, b) + B   (b, a) = B(a, b)
     *    x          1-x
     *
     * and B_1-x(b,a) can be evaluated as described above.
     *
     * Once either a lower or an upper incomplete beta function is evaluated,
     * the other value may be quickly obtained by applying the following
     * property of the incomplete beta function:
     *
     *   Bx(a,b) + bx(a,b) = B(a,b)
     *
     * A value of a regularized incomplete beta function is obtained
     * by dividing Bx(a,b) or bx(a,b) by B(a,b).
     */

    // checks if x is "large", i.e. greater than the threshold defined above:
    const bool xlarge = (x > (a+static_cast<T>(1)) / (a+b+static_cast<T>(2)) );

    // If x is large, the parameters must be swapped
    const T xn = ( false==xlarge ? x :  static_cast<T>(1)-x );
    const T an = ( false==xlarge ? a : b );
    const T bn = ( false==xlarge ? b : a );

    // B = B(a,b) or 0 when the value is not used
    const T B = ( 
        true==reg || (false==xlarge && false==lower) || (true==xlarge  && true==lower) ?
        math::SpecFun::beta<T>(a, b) :
        static_cast<T>(0) );

    // An instance of __CtdFIncBeta that implements 'fa' and 'fb':
    const math::SpecFun::__private::__CtdFIncBeta<T> coef(an, bn, xn);

    T binc;

    if ( true == math::NumericUtil::isZero<T>(xn) )
    {
        /*
         * Both boundary conditions are handled separately:
         *   B0(a,b) = 0  and  B1(a,b) = B(a,b)
         */

        binc = ( true==xlarge ? B : static_cast<T>(0) );
    }
    else
    {
        // 'x' is somewhere between 0 and 1, apply the algorithm described above

        binc = std::pow(xn, an) * std::pow(static_cast<T>(1)-xn, bn) / an;
        binc /= math::CtdFrac::ctdFrac<T>(coef, tol);
    }

    /*
     * When x is "large", the algorithm actually returns the
     * upper incomplete beta function!
     *
     * Depending on the requested result ('lower') adjust
     * 'binc' if necessary
     */
    if ( (false==xlarge && false==lower) ||
         (true==xlarge  && true==lower) )
    {
        binc = B - binc;
    }

    // Finally regularize the result if requested (via 'reg')
    if ( true == reg )
    {
        // Just in case handle the very unlikely case
        if ( true == math::NumericUtil::isZero<T>(B) )
        {
            throw math::SpecFunException(math::SpecFunException::UNDEFINED);
        }

        binc /= B;
    }

    return binc;
}


/*
 * Partial "specialization" of __incBeta for complex numbers.
 * 
 * Evaluates an incomplete beta function. The exact kind of the returned
 * value depends on parameters 'upper' and 'reg'.
 *
 * @note unlike at real numbers, incomplete beta function is defined
 *       virtually everywhere except if 'a' is a negative integer, currently this
 *       function additionally requires that /x/ < 1 otherwise it may not converge.
 *
 * @param a - parameter 'a' of the beta function
 * @param b - parameter 'b' of the beta function
 * @param x - integration limit
 * @param lower - should the lower (if 'true') or the upper (if 'false) inc. beta function be returned
 * @param reg - if 'true', the regularized beta function is returned, i.e. divided by beta(a, b)
 * @param tol - tolerance (default: 1e-6)
 */
template <class T>
std::complex<T> __incBeta(
          const std::complex<T>& a,
          const std::complex<T>& b,
          const std::complex<T>& x,
          const bool lower,
          const bool reg,
          const std::complex<T>& tol
        ) throw (math::SpecFunException)
{
    /*
     * Series expansion of the incomplete beta function as explained at:
     * http://functions.wolfram.com/GammaBetaErf/Beta3/06/01/03/01/01/
     * 
     *                a    /                                   2        \
     *               x     |      a*(1-b)*x     a*(1-b)*(2-b)*x         |
     *   Bx(a,b) ~= ---- * | 1 + ----------- + ------------------ + ... | =
     *               a     |        (a+1)           2*(a+2)             |
     *                     \                                            /
     * 
     * 
     *                     /       inf                                         \
     *                a    |      -----                                     i  |
     *               x     |      \      a * (1-b) * (2-b) * ... * (i-b) * x   |
     *            = ---- * | 1 +   >    -------------------------------------- |
     *               a     |      /                  i! * (a+i)                |
     *                     |      -----                                        |
     *                     \       i=1                                         /
     * 
     */

    // Division by zero will occur if a==0
    if ( true == math::NumericUtil::isZero<std::complex<T> >(a) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    // Currently the function requires that /x/ < 1
    // when the series is guaranteed to converge
    if ( std::norm(x) >= static_cast<T>(1) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    // B = B(a,b) or (0,0) when the value is not used
    const std::complex<T> B = ( 
            true==reg || false==lower ?
            math::SpecFun::beta<std::complex<T> >(a, b) :
            std::complex<T>(static_cast<T>(0)) );

    // The first term of the series
    std::complex<T> binc = std::pow(x, a) / a;
    std::complex<T> term = binc * a;

    std::complex<T> at = a;
    std::complex<T> bt = -b;
 
    // proceed the series until it converges
    std::size_t i = 1;
    for ( i = 1;
          false == math::NumericUtil::isZero<std::complex<T> >(term, tol) &&
                i <= SPECFUN_MAX_ITER;
          ++i )
    {
        at += static_cast<T>(1);
        bt += static_cast<T>(1);
 
        // if 'a' is a negative integer, sooner or later this exception will be thrown
        if ( true == math::NumericUtil::isZero<std::complex<T> >(at) )
        {
            throw math::SpecFunException(math::SpecFunException::UNDEFINED);
        }

        term *= x * bt / static_cast<T>(i);
        
        binc += term / at;
    }

    // has the series converged?
    if ( i>= SPECFUN_MAX_ITER )
    {
        throw math::SpecFunException(math::SpecFunException::NO_CONVERGENCE);
    }

    /*
     * Apply properties of the incomplete beta function
     * if anything else except a non-regularized lower incomplete
     * beta function is desired.
     */
    if ( false == lower )
    {
        binc = B - binc;
    }
    
    if ( true == reg )
    {
        // Very unlikely but check it anyway
        if ( true == math::NumericUtil::isZero<std::complex<T> >(B) )
        {
            throw math::SpecFunException(math::SpecFunException::UNDEFINED);
        }

        binc /= B;
    }

    return binc;
}

}}}  // namespace math::SpecFun::__private



/**
 * Upper incomplete gamma function, defined as:
 *
 *           inf
 *            /
 *            |  a-1    -t
 *   G(a,x) = | t    * e   dt
 *            |
 *            /
 *            x
 *
 * @note For real arguments, 'a' and 'x' must be greater than 0.
 *       For complex arguments, the function is defined for virtually
 *       any combination of 'a' and 'x' except when 'a' is 0 or a
 *       negative integer.
 *
 * @param a - parameter of the incomplete gamma function
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
    return math::SpecFun::__private::__incGamma(a, x, true, false, tol);
}


/**
 * Lower incomplete gamma function, defined as:
 *
 *            x
 *            /
 *            |  a-1    -t
 *   g(a,x) = | t    * e   dt
 *            |
 *            /
 *            0
 *
 * @note For real arguments, 'a' and 'x' must be greater than 0.
 *       For complex arguments, the function is defined for virtually
 *       any combination of 'a' and 'x' except when 'a' is 0 or a
 *       negative integer.
 *
 * @param a - parameter of the incomplete gamma function
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
    return math::SpecFun::__private::__incGamma(a, x, false, false, tol);
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
 * @note For real arguments, 'a' and 'x' must be greater than 0.
 *       For complex arguments, the function is defined for virtually
 *       any combination of 'a' and 'x' except when 'a' is 0 or a
 *       negative integer.
 *
 * @param a - parameter of the incomplete gamma function
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
    return math::SpecFun::__private::__incGamma(a, x, true, true, tol);
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
 * @note For real arguments, 'a' and 'x' must be greater than 0.
 *       For complex arguments, the function is defined for virtually
 *       any combination of 'a' and 'x' except when 'a' is 0 or a
 *       negative integer.
 *
 * @param a - parameter of the incomplete gamma function
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
    return math::SpecFun::__private::__incGamma(a, x, false, true, tol);
}


/**
 * Lower incomplete beta function, defined as:
 *
 *             x
 *             /
 *             |  a-1        b-1
 *   Bx(a,b) = | t    * (1-t)    dt
 *             |
 *             /
 *             0
 *
 * @note For real arguments, 'a' and 'b' must be greater than 0,
 *       'x' must be greater than 0 and less than 1.
 *       For complex arguments, 'a' and 'b' can be any complex values except
 *       'a' must not be 0 or a negative integer. Currently the x is restricted
 *       to: /x/ < 1.
 *
 * @param a - parameter 'a' of the beta function
 * @param b - parameter 'b' of the beta function
 * @param x - integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return Bx(a,b)
 *
 * @throw SpecFunException if 'a', 'b' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incBetaLower(
               const T& a,
               const T& b,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incBeta(a, b, x, true, false, tol);
}


/**
 * Upper incomplete beta function, defined as:
 *
 *             1
 *             /
 *             |  a-1        b-1
 *   bx(a,b) = | t    * (1-t)    dt = B(a,b) - Bx(a,b)
 *             |
 *             /
 *             x
 *
 * @note For real arguments, 'a' and 'b' must be greater than 0,
 *       'x' must be greater than 0 and less than 1.
 *       For complex arguments, 'a' and 'b' can be any complex values except
 *       'a' must not be 0 or a negative integer. Currently the x is restricted
 *       to: /x/ < 1.
 *
 * @param a - parameter 'a' of the beta function
 * @param b - parameter 'b' of the beta function
 * @param x - integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return bx(a,b)
 *
 * @throw SpecFunException if 'a', 'b' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incBetaUpper(
               const T& a,
               const T& b,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incBeta(a, b, x, false, false, tol);
}


/**
 * Regularized lower incomplete beta function, defined as:
 *
 *                      x
 *                      /
 *                1     |  a-1        b-1
 *   Ix(a,b) = -------- | t    * (1-t)    dt
 *              B(a,b)  |
 *                      /
 *                      0
 *
 * @note For real arguments, 'a' and 'b' must be greater than 0,
 *       'x' must be greater than 0 and less than 1.
 *       For complex arguments, 'a' and 'b' can be any complex values except
 *       'a' must not be 0 or a negative integer. Currently the x is restricted
 *       to: /x/ < 1.
 *
 * @param a - parameter 'a' of the beta function
 * @param b - parameter 'b' of the beta function
 * @param x - integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return Ix(a,b)
 *
 * @throw SpecFunException if 'a', 'b' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incBetaLowerReg(
               const T& a,
               const T& b,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incBeta(a, b, x, true, true, tol);
}


/**
 * Regularized upper beta function, defined as:
 *
 *                      1
 *                      /
 *                1     |  a-1        b-1
 *   ix(a,b) = -------- | t    * (1-t)    dt = 1 - Ix(a,b)
 *              B(a,b)  |
 *                      /
 *                      x
 *
 * @note For real arguments, 'a' and 'b' must be greater than 0,
 *       'x' must be greater than 0 and less than 1.
 *       For complex arguments, 'a' and 'b' can be any complex values except
 *       'a' must not be 0 or a negative integer. Currently the x is restricted
 *       to: /x/ < 1.
 *
 * @param a - parameter 'a' of the beta function
 * @param b - parameter 'b' of the beta function
 * @param x - integration limit
 * @param tol - tolerance (default: 1e-6)
 *
 * @return ix(a,b)
 *
 * @throw SpecFunException if 'a', 'b' or 'x' is invalid
 */
template <class T>
T math::SpecFun::incBetaUpperReg(
               const T& a,
               const T& b,
               const T& x,
               const T& tol
             ) throw(math::SpecFunException)
{
    return math::SpecFun::__private::__incBeta(a, b, x, false, true, tol);
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
     * more efficient methods exist. The exponential can be expanded
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



// Implementation of "private" functions that evaluate
// inverses of incomplete gamma and beta
namespace math {  namespace SpecFun {  namespace __private {


/*
 * Implementation of the formula 26.2.22 in [Abramowitz & Stegun].
 *
 * Unlike the original formula, this function accepts 'p' being
 * greater than 0.5.
 *
 * @note This "private" function expects that p is greater than 0
 *       and less than 1.
 *
 * @param p - desired probability
 *
 * @return approximation of 'x' that satisfies Q(x) ~= p
 */
template <class T>
T __as26_2_22(const T& p)
{
    /*
     * If 0 < p < 0.5:
     *
     *             +--------+
     *      -+    /      1               +------------+
     *   t =  \  /  ln -----   =    -+  / (-2) * ln(p)
     *         \/       p^2           \/
     *
     * Then initial approximation of 'x' can be calculated as:
     *
     *                 2.30753 + 0.27061 * t
     *   x ~= t - ---------------------------------
     *             1 + t * (0.99229 + 0.04481 * t)
     *
     * If p > 0.5, the expressions above are performed on its complementary
     * value (1-p) and the sign of the final  'x' is reversed.
     */

    const T pp = ( p>= static_cast<T>(1) / static_cast<T>(2) ?
                     static_cast<T>(1) - p : p );

    const T t = std::sqrt(static_cast<T>(-2) * std::log(pp) );

    // [Abramowitz & Stegun], section 26.2.22:
    const T a0 = static_cast<T>(230753) / static_cast<T>(100000);   // 2.30753
    const T a1 = static_cast<T>(27061) / static_cast<T>(100000);    // 0.27061
    const T b0 = static_cast<T>(99229) / static_cast<T>(100000);    // 0.99229
    const T b1 = static_cast<T>(4481) / static_cast<T>(100000);     // 0.04481

    T x = t - (a0 + a1 * t) / (1.0 + t * (b0 + b1 * t));
    if ( p > static_cast<T>(1) / static_cast<T>(2) )
    {
        x = -x;
    }

    return x;
}


/*
 * Evaluates inverse of an incomplete gamma function. The exact kind
 * of the returned value depends on parameters 'upper' and 'reg'.
 *
 * @note both 'a' and 'g' must be strictly greater than 0
 *
 * @param a - parameter of the incomplete gamma function
 * @param g - desired value of the incomplete gamma function
 * @param upper -  if 'true', inverse of the upper incomplete gamma function will be evaluated
 * @param reg - if 'true', inverse of the regularized incomplete gamma function (divided by gamma(a)) will be evaluated
 * @param tol - tolerance
 */
template <class T>
T __invIncGamma(
             const T& a,
             const T& g,
             const bool upper,
             const bool reg,
             const T& tol
           ) throw (math::SpecFunException)
{
    // sanity check
    if ( a < math::NumericUtil::getEPS<T>() ||
         g < static_cast<T>(0) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    /*
     * The algorithm finds an inverse of the regularized lower
     * incomplete gamma function. If 'g' represents a result of
     * any other kind of the incomplete gamma function (specified)
     * by 'upper' and 'reg', it should be converted to a result of
     * the regularized lower incomplete gamma function by applying
     * well known properties of this function.
     *
     * The algorithm is proposed in [Numerical Recipes], section 6.2.1.
     *
     */

    // gamma(a)
    const T G = math::SpecFun::gamma<T>(a);

    T p = g;
    if ( false == reg )
    {
        p /= G;
    }
    if ( true == upper )
    {
        p = static_cast<T>(1) - p;
    }

    // if 'p' equals 0 (or is very close to it),
    // the inverse inc. gamma function will also equal 0.
    if ( true == math::NumericUtil::isZero<T>(p) )
    {
        return static_cast<T>(0);
    }

    /*
     * If 'p' is greater than 1, the incomplete gamma function
     * is not defined at all.
     * If 'p' equals 1 (or is very close to it), the result would
     * be infinity which is not supported by this algorithm.
     */
    if ( p>= (static_cast<T>(1) - math::NumericUtil::getEPS<T>()) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }


    /*
     * There is no direct series or continued fraction to evaluate
     * an inverse of the incomplete gamma function. On the other hand
     * there are known algorithms that estimate this value very well.
     * This algorithm implements approximation method proposed by the
     * [Numerical Recipes], section 6.2.1 and [Abramowitz & Stegun],
     * sections 26.2.22 and 26.4.17.
     *
     * Note: another method to estimate the initial value of the incomplete
     * gamma function is described in:
     *   Armido R. DiDonato, Alfred H. Morris
     *   Computation of the incomplete gamma function ratios and their inverse
     *   ACM Transactions on Mathematical Software 
     *   Volume 12 Issue 4, Dec. 1986, pp. 377-393
     *   http://dl.acm.org/citation.cfm?id=23109&coll=portal&dl=ACM
     * 
     * When an approximation is known, the Newton - Raphson method can be
     * applied to refine it further to the desired tolerance.
     */

    T x = static_cast<T>(1);

    if ( a <= static_cast<T>(1) )
    {
        /*
         * a <= 1:
         *
         * First P(a,1) is approximated by:
         *
         *   Pa = P(a,1) ~= a* (0.253 + 0.12 * a)
         *
         * If 'p' is less than Pa's complement value, 'x' is
         * approximated as:
         *
         *         a      +------+
         *        --+    /   p
         *   x ~=    \  / ------
         *            \/   1-Pa
         *
         * If 'p' is greater than 1-Pa, 'x' is approximated as:
         *
         *                1 - p
         *   x ~= 1 - ln -------
         *                 Pa
         *
         * Note: as long as 0<a<1, Pa will never get even close to 1,
         *       hence the first expression is always defined.
         */

        const T c1 = static_cast<T>(253) / static_cast<T>(1000);  // 0.253
        const T c2 = static_cast<T>(12) / static_cast<T>(100);    // 0.12
        const T Pa = a * (c1 + c2 * a);
        const T cpa = static_cast<T>(1) - Pa;

        x = ( p<cpa ?
                 std::pow(p / cpa, static_cast<T>(1) / a) :
                 static_cast<T>(1) - std::log((static_cast<T>(1) - p)/ Pa) );
    }
    else
    {
        /*
         * a > 1:
         *
         * The initial approximation of 'x' is obtained by the function
         * math::SpecFun::__private::__as26_2_22()
         * that implements [Abramowitz & Stegun], section 26.2.22.
         *
         * Finally the approximated x0 is evaluated as:
         *
         *             /       1          x      \ 3
         *   x0 ~= a * | 1 - ----- + ----------- |
         *             \      9*a     3*sqrt(a)  /
         *
         */

        x = -math::SpecFun::__private::__as26_2_22<T>(p);

        // [Abramowitz & Stegun], section 26.4.17:
        x = static_cast<T>(1) - static_cast<T>(1) / (static_cast<T>(9) * a) +
                x / ( static_cast<T>(3) * std::sqrt(a) );

        // x = a * x^3:
        x  *= a * x * x;
    }


    /*
     * When the initial value of 'x' is obtained, a root finding
     * method is applied to determine the exact inverse. This algorithm
     * applies the slightly modified Newton - Raphson method:
     * if 'x' tries to go negative, its half value is assigned
     * to the new value.
     * For that reason, the Newton - Raphson method must be reimplemented
     * as the general implementation (root_find/RootFindGeneric) cannot
     * be used.
     */

    /*
     * The Newton - Raphson algorithm also requires the differentiation of
     * the regularized lower incomplete gamma function:
     *
     *                   a-1    -x
     *    d P(a,x)      x    * e
     *   ---------- = --------------
     *       dx         gamma(a)
     *
     * Verified by Maxima:
     (%i1)  diff(gamma_incomplete_regularized(a, x), x);
     (%o1)  (x^a−1*%e^−x)/gamma(a)
     * 
     * See also 'scripts/lib/specfun.mac'
     */

    T xn;
    T f = math::SpecFun::incGammaLowerReg<T>(a, x, tol) - p;
    std::size_t i = 0;
    for ( i = 0;
          false==math::NumericUtil::isZero<T>(f, tol) && i<SPECFUN_MAX_ITER;
          ++i )
    {
        xn = x - f * G * std::exp(x) / std::pow(x, a-static_cast<T>(1));
        
        // x must not go negative!
        x = ( xn > math::NumericUtil::getEPS<T>() ?
              xn : x / static_cast<T>(2) );
        
        f = math::SpecFun::incGammaLowerReg<T>(a, x, tol) - p;
    }

    // Has the algorithm converged?
    if ( i>= SPECFUN_MAX_ITER )
    {
        throw math::SpecFunException(math::SpecFunException::NO_CONVERGENCE);
    }

    return x;
}


/*
 * Evaluates inverse of an incomplete beta function. The exact kind
 * of the returned value depends on parameters 'lower' and 'reg'.
 *
 * @note 'a', 'b' and 'y' must be strictly greater than 0
 *
 * @param a - parameter a of the incomplete beta function
 * @param b - parameter b of the incomplete beta function
 * @param y - desired value of the incomplete beta function
 * @param lower -  if 'true', inverse of the lower incomplete beta function will be evaluated
 * @param reg - if 'true', inverse of the regularized incomplete beta function (divided by beta(a,b)) will be evaluated
 * @param tol - tolerance
 */
template <class T>
T __invIncBeta(
             const T& a,
             const T& b,
             const T& y,
             const bool lower,
             const bool reg,
             const T& tol
           ) throw (math::SpecFunException)
{
    // sanity check
    if ( a < math::NumericUtil::getEPS<T>() ||
         b < math::NumericUtil::getEPS<T>() ||
         y < static_cast<T>(0) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    /*
     * The algorithm finds an inverse of the regularized lower
     * incomplete beta function. If 'y' represents a result of
     * any other kind of the incomplete beta function (specified)
     * by 'lower' and 'reg', it should be converted to a result of
     * the regularized lower incomplete beta function by applying
     * well known properties of this function.
     *
     * The algorithm is proposed in [Numerical Recipes], section 6.4.
     *
     */

    // beta(a,b)
    const T B = math::SpecFun::beta<T>(a, b);

    T p = y;

    if ( false == reg )
    {
        p /= B;
    }

    if (false == lower)
    {
        p = static_cast<T>(1) - p;
    }

    // if 'p' equals 0 (or is very close to it),
    // the inverse inc. beta function will also equal 0.
    if ( true == math::NumericUtil::isZero<T>(p) )
    {
        return static_cast<T>(0);
    }

    /*
     * If 'p' is greater than 1, the incomplete beta function
     * is not defined at all.
     * If 'p' equals 1 (or is very close to it), the result would
     * be infinity which is not supported by this algorithm.
     */
    if ( p>= (static_cast<T>(1) - math::NumericUtil::getEPS<T>()) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    /*
     * There is no known direct series or continued fraction to evaluate
     * an inverse of the incomplete beta function. On the other hand
     * there are known algorithms that estimate this value quite well.
     * This algorithm implements approximation method proposed by the
     * [Numerical Recipes], section 6.4 and [Abramowitz & Stegun],
     * section 26.5.22.
     *
     * When an approximation is known, the Newton - Raphson method can be
     * applied to refine it further to the desired tolerance.
     */

    T x = static_cast<T>(0);

    if (a>=static_cast<T>(1) && b>=static_cast<T>(1) )
    {
        /*
         * a >= 1  and  b >= 1:
         *
         * Initial 'x' is approximated as described in [Abramowitz & Stegun],
         * section 26.2.22 (implemented in a separate function).
         *
         * Then it is further refined as described in [Abramowitz & Stegun],
         * section 26.5.22:
         *
         *           1                   1
         *   ta = -------   and  tb = -------
         *         2*a-1               2*b-1
         *
         *              2
         *             x  - 3
         *   lambda = --------
         *               6
         *
         *           2
         *   h = ---------
         *        ta + tb
         *
         *        x * sqrt(h + lambda)                /           5      2   \
         *   w = ---------------------- - (tb - ta) * | lambda + --- - ----- |
         *                h                           \           6     3*h  /
         *
         *                   a
         *   x  ~=   ------------------
         *                       2*w
         *              a + b * e
         */

        // [Abramowitz & Stegun], section 26.2.22:
        x = math::SpecFun::__private::__as26_2_22<T>(p);

        const T lambda = (x*x - static_cast<T>(3)) / static_cast<T>(6);
        const T ta = static_cast<T>(1) / ( static_cast<T>(2) * a - static_cast<T>(1) );
        const T tb = static_cast<T>(1) / ( static_cast<T>(2) * b - static_cast<T>(1) );
        const T h = static_cast<T>(2) / (ta + tb);

        const T w = x * std::sqrt(h + lambda) / h - (tb - ta) *
               (lambda + static_cast<T>(5) / static_cast<T>(6) - static_cast<T>(2) / (static_cast<T>(3) * h) );

        x = a / (a + b * std::exp(static_cast<T>(2) * w));
    }
    else
    {
        /*
         * Either 'a' or 'b' or both are less than 1:
         *
         * In this case the initial 'x' is approximated as
         * described in [Numerical Recipes], section 6.4:
         *
         *
         *         1    /    a    \a              1    /    b    \b
         *   ta = --- * | ------- |    and  tb = --- * | ------- |
         *         a    \  a + b  /               b    \  a + b  /
         *
         *
         *   S = ta + tb
         *
         * If p < ta/S:
         *
         *         a    +----------+
         *   x ~=  -+  / p * S * a
         *           \/
         *
         * If p >= ta/S:
         *
         *        b    +----------+
         *   x ~= -+  / 1 - p*S*b
         *          \/
         */

        const T ta = std::pow(a/(a+b), a) / a;
        const T tb = std::pow(b/(a+b), b) / b;
        const T S = ta + tb;

        x = ( p < ta/S ?
                std::pow(p*S*a, static_cast<T>(1)/a) :
                std::pow(static_cast<T>(1)-p*S*b, static_cast<T>(1)/b) );
    }

    /*
     * When the initial value of 'x' is obtained, a root finding
     * method is applied to determine the exact inverse. This algorithm
     * applies the slightly modified Newton - Raphson method:
     * if 'x' tries to go negative or beyond 1, it is bisected between
     * its current value and the upper/lower bondary.
     * For that reason, the Newton - Raphson method must be reimplemented
     * as the general implementation (root_find/RootFindGeneric) cannot
     * be used.
     */

    /*
     * The Newton - Raphson algorithm also requires the differentiation of
     * the regularized lower incomplete beta function:
     *
     *                   a-1        b-1
     *    d Ix(a,b)     x    * (1-x)
     *   ----------- = ------------------
     *      dx             beta(a,b)
     *
     * Verified by Maxima:
     (%i1) diff(beta_incomplete_regularized(a, b, x), x);
     (%o1) ((1−x)^b−1*x^a−1)/beta(a,b)
     * 
     * See also 'scripts/lib/specfun.mac'
     */

    T xn;
    T f = math::SpecFun::incBetaLowerReg<T>(a, b, x, tol) - p;
    std::size_t i = 0;
    for ( i = 0;
          false==math::NumericUtil::isZero<T>(f, tol) && i<SPECFUN_MAX_ITER;
          ++i )
    {
        xn = x - f * B / (std::pow(x, a-static_cast<T>(1)) * std::pow(static_cast<T>(1)-x, b-static_cast<T>(1)) );

        // x must not go negative or beyond 1!
        if ( xn < math::NumericUtil::getEPS<T>() )
        {
            x /= static_cast<T>(2);
        }
        else if ( xn > (static_cast<T>(1)-math::NumericUtil::getEPS<T>()) )
        {
            x = (static_cast<T>(1) + x) / static_cast<T>(2);
        }
        else
        {
            x = xn;
        }

        f = math::SpecFun::incBetaLowerReg<T>(a, b, x, tol) - p;
    }

    // Has the algorithm converged?
    if ( i>= SPECFUN_MAX_ITER )
    {
        throw math::SpecFunException(math::SpecFunException::NO_CONVERGENCE);
    }

    return x;
}


}}}  // namespace math::SpecFun::__private



/**
 * Inverse of the lower incomplete gamma function,
 * i.e. it returns such 'x' that satisfies:
 *
 *    x
 *    /
 *    |  a-1    -t
 *    | t    * e   dt  =  g
 *    |
 *    /
 *    0
 *
 * @note 'a' must be strictly greater than 0, 'g' must be greater or equal
 *       to zero and less than gamma(a).
 *
 * @param a - parameter of the incomplete gamma function
 * @param g - desired value of the lower incomplete gamma function
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of g(a,x)
 *
 * @throw SpecFunException if 'a' or 'g' is invalid
 */
template <class T>
T math::SpecFun::incGammaLowerInv(
           const T& a,
           const T& g,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncGamma<T>(a, g, false, false, tol);
}


/**
 * Inverse of the upper incomplete gamma function,
 * i.e. it returns such 'x' that satisfies:
 *
 *   inf
 *    /
 *    |  a-1    -t
 *    | t    * e   dt  =  g
 *    |
 *    /
 *    x
 *
 * @note 'a' must be strictly greater than 0, 'g' must be greater than
 *       to zero and less or equal to gamma(a).
 *
 * @param a - parameter of the incomplete gamma function
 * @param g - desired value of the upper incomplete gamma function
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of G(a,x)
 *
 * @throw SpecFunException if 'a' or 'g' is invalid
 */
template <class T>
T math::SpecFun::incGammaUpperInv(
           const T& a,
           const T& g,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncGamma<T>(a, g, true, false, tol);
}


/**
 * Inverse of the regularized lower incomplete gamma function,
 * i.e. it returns such 'x' that satisfies:
 *
 *          x
 *          /
 *     1    |  a-1    -t
 *   ------ | t    * e   dt  =  g
 *    G(a)  |
 *          /
 *          0
 *
 * @note 'a' must be strictly greater than 0, 'g' must be greater or equal
 *       to zero and less than 1.
 *
 * @param a - parameter of the incomplete gamma function
 * @param g - desired value of the regularized lower incomplete gamma function
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of P(a,x)
 *
 * @throw SpecFunException if 'a' or 'g' is invalid
 */
template <class T>
T math::SpecFun::incGammaLowerRegInv(
           const T& a,
           const T& g,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncGamma<T>(a, g, false, true, tol);
}


/**
 * Inverse of the upper incomplete gamma function,
 * i.e. it returns such 'x' that satisfies:
 *
 *         inf
 *          /
 *     1    |  a-1    -t
 *   ------ | t    * e   dt  =  g
 *    G(a)  |
 *          /
 *          x
 *
 * @note 'a' must be strictly greater than 0, 'g' must be greater than
 *       to zero and less or equal to 1.
 *
 * @param a - parameter of the incomplete gamma function
 * @param g - desired value of the regularized upper incomplete gamma function
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of Q(a,x)
 *
 * @throw SpecFunException if 'a' or 'g' is invalid
 */
template <class T>
T math::SpecFun::incGammaUpperRegInv(
           const T& a,
           const T& g,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncGamma<T>(a, g, true, true, tol);
}


/**
 * Inverse of the lower incomplete beta function,
 * i.e. it returns such 'x' that satisfies:
 *
 *             x
 *             /
 *             |  a-1        b-1
 *   Bx(a,b) = | t    * (1-t)    dt = y
 *             |
 *             /
 *             0
 *
 * @note 'a' and 'b' must be strictly greater than 0, 'y' must be greater or equal
 *       to zero and less than beta(a,b).
 *
 * @param a - parameter 'a' of the incomplete beta function
 * @param b - parameter 'b' of the incomplete beta function
 * @param y - desired value of the lower incomplete beta function (Bx(a,b))
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of Bx(a,b)
 *
 * @throw SpecFunException if 'a', b' or 'y' is invalid
 */
template <class T>
T math::SpecFun::incBetaLowerInv(
           const T& a,
           const T& b,
           const T& y,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncBeta<T>(a, b, y, true, false, tol);
}


/**
 * Inverse of the upper incomplete beta function,
 * i.e. it returns such 'x' that satisfies:
 *
 *             1
 *             /
 *             |  a-1        b-1
 *   bx(a,b) = | t    * (1-t)    dt = y
 *             |
 *             /
 *             x
 *
 * @note 'a' and 'b' must be strictly greater than 0, 'y' must be greater or equal
 *       to zero and less than beta(a,b).
 *
 * @param a - parameter 'a' of the incomplete beta function
 * @param b - parameter 'b' of the incomplete beta function
 * @param y - desired value of the upper incomplete beta function (bx(a,b))
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of bx(a,b)
 *
 * @throw SpecFunException if 'a', b' or 'y' is invalid
 */
template <class T>
T math::SpecFun::incBetaUpperInv(
           const T& a,
           const T& b,
           const T& y,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncBeta<T>(a, b, y, false, false, tol);
}


/**
 * Inverse of the regularized lower incomplete beta function,
 * i.e. it returns such 'x' that satisfies:
 *
 *                       x
 *                       /
 *                1      |  a-1        b-1
 *   Ix(a,b) = --------  | t    * (1-t)    dt = y
 *              B(a,b)   |
 *                       /
 *                       0
 *
 * @note 'a' and 'b' must be strictly greater than 0, 'y' must be greater or equal
 *       to zero and less than 1.
 *
 * @param a - parameter 'a' of the incomplete beta function
 * @param b - parameter 'b' of the incomplete beta function
 * @param y - desired value of the regularized lower incomplete beta function (Ix(a,b))
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of Ix(a,b)
 *
 * @throw SpecFunException if 'a', b' or 'y' is invalid
 */
template <class T>
T math::SpecFun::incBetaLowerRegInv(
           const T& a,
           const T& b,
           const T& y,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncBeta<T>(a, b, y, true, true, tol);
}


/**
 * Inverse of the regularized upper incomplete beta function,
 * i.e. it returns such 'x' that satisfies:
 *
 *                      1
 *                      /
 *                1     |  a-1        b-1
 *   ix(a,b) = -------  | t    * (1-t)    dt = y
 *              B(a,b)  |
 *                      /
 *                      x
 *
 * @note 'a' and 'b' must be strictly greater than 0, 'y' must be greater or equal
 *       to zero and less than 1.
 *
 * @param a - parameter 'a' of the incomplete beta function
 * @param b - parameter 'b' of the incomplete beta function
 * @param y - desired value of the regularized upper incomplete beta function (ix(a,b))
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of ix(a,b)
 *
 * @throw SpecFunException if 'a', b' or 'y' is invalid
 */
template <class T>
T math::SpecFun::incBetaUpperRegInv(
           const T& a,
           const T& b,
           const T& y,
           const T& tol
         ) throw (math::SpecFunException)
{
    return math::SpecFun::__private::__invIncBeta<T>(a, b, y, false, true, tol);
}


/**
 * Inverse of the error function, i.e it returns such 'x'
 * that satisfies:
 *
 *   erf(x) = e
 *
 * @note 'e' must be greater than -1 and less than 1
 *
 * @param e - desired value of the error function
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of erf(x)
 *
 * @throw SpecFunException if 'e' is out of the definition range
 */
template <class T>
T math::SpecFun::erfInv(
           const T& e,
           const T& tol
         ) throw (math::SpecFunException)
{
    /*
     * The following properties of the error function will be applied:
     *
     *   x >= 0  ==>  erf(x) = P(1/2, x^2)  and  0 <= erf(x) <= 1
     *   x <= 0  ==>  erf(x) = -erf(x)  and  -1 <= erf(x) <= 0
     *
     * 'x' can be obtained as follows:
     *
     * If e>=0:
     *   x = sqrt(invP(1/2, e))
     *
     * If e<0:
     *   x = -sqrt(invP(1/2, -e)
     */

    // sanity check
    if ( e < (static_cast<T>(-1) + math::NumericUtil::getEPS<T>()) ||
         e > (static_cast<T>(1) - math::NumericUtil::getEPS<T>()) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    // handle e==0 as a special case:
    if ( true == math::NumericUtil::isZero<T>(e) )
    {
        return static_cast<T>(0);
    }

    const T arg = ( e<static_cast<T>(0) ? -e : e );
    const T x = std::sqrt(math::SpecFun::incGammaLowerRegInv<T>(
             static_cast<T>(1) / static_cast<T>(2),
             arg,
             tol ) );

    return ( e>static_cast<T>(0) ? x : -x );

}


/**
 * Inverse of the complementary error function, i.e it returns
 * such 'x' that satisfies:
 *
 *   erfc(x) = e
 *
 * @note 'e' must be greater than 0 and less than 2
 *
 * @param e - desired value of the complementary error function
 * @param tol - tolerance (default: 1e-6)
 *
 * @return inverse of erfc(x)
 *
 * @throw SpecFunException if 'e' is out of the definition range
 */
template <class T>
T math::SpecFun::erfcInv(
           const T& e,
           const T& tol
         ) throw (math::SpecFunException)
{
    /*
     * The following properties of the complementary error function will be applied:
     *
     *   x >= 0  ==>  erfc(x) = Q(1/2, x^2)  and  0 <= erfc(x) <= 1
     *   x <= 0  ==>  erfc(x) = 2 - erfc(-x)  and  1 <= erfc(x) <= 2
     *
     * 'x' can be obtained as follows:
     *
     * If e<=1:
     *   x = sqrt(invQ(1/2, e))
     *
     * If e>1:
     *   x = -sqrt(invQ(1/2, 2-e))
     */

    // sanity check
    if ( e < math::NumericUtil::getEPS<T>() ||
         e > (static_cast<T>(2) - math::NumericUtil::getEPS<T>()) )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }

    // handle e==1 as a special case:
    if ( true == math::NumericUtil::isZero<T>(e-static_cast<T>(1)) )
    {
        return static_cast<T>(0);
    }

    const T arg = ( e<static_cast<T>(1) ? e : static_cast<T>(2)-e );
    const T x = std::sqrt(math::SpecFun::incGammaUpperRegInv<T>(
              static_cast<T>(1) / static_cast<T>(2),
              arg,
              tol ) );

    return ( e<static_cast<T>(1) ? x : -x );
}
