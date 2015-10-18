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
 * Implementation of auxiliary classes and functions that efficiently
 * evaluate continued fractions.
 */


// no #include "CtdFracGeneric.hpp" !!!!
#include <cstddef>
#include <cmath>

#include "../settings/specfun_settings.h"
#include "util/NumericUtil.hpp"

#include "exception/FunctionException.hpp"
#include "exception/SpecFunException.hpp"


/*
 * ICtdFracFuncGeneric's destructor, "implemented" as an empty function
 */
template <class T>
math::CtdFrac::ICtdFracFuncGeneric<T>::~ICtdFracFuncGeneric()
{
    // empty destructor
}



/*
 * Evaluates the continued fraction:
 *
 *                       a1
 *   f = b0 + -------------------------
 *                          a2
 *              b1 + -----------------
 *                             a3
 *                    b2 + ----------
 *                          b3 + ...
 *
 * where ai and bi are functions of 'i'.
 *
 * @param ctdf - instance of a class ICtdFracFuncGeneric that returns values of 'a_i' and 'b_i'
 * @param tol - tolerance (default: 1e-6)
 *
 * @return the value of the continued fraction, specified by terms 'ai' and 'bi'
 *
 * @throw SpecFunException if 'ctdf.fa' or 'ctdf.fb' is undefined for any 'i'
 */
template <class T>
T math::CtdFrac::ctdFrac(
           const math::CtdFrac::ICtdFracFuncGeneric<T>& ctdf,
           const T& tol
         ) throw (math::SpecFunException)
{
    /*
     * The Lentz's algorithm (modified by I. J. Thompson and A. R. Barnett)
     * is applied to evaluate the continued fraction. The algorithm is
     * presented in detail in:
     *
     *   William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
     *   Numerical Recipes, The Art of Scientific Computing, 3rd Edition,
     *   Cambridge University Press, 2007
     *
     *   https://books.google.com/books?id=1aAOdzK3FegC&lpg=PA207&ots=3jNoK9Crpj&pg=PA208#v=onepage&f=false
     *
     * The procedure of the algorithm is as follows:
     *
     * - f0 = b0, if b0==0 then f0 = eps
     * - C0 = f0
     * - D0 = 0
     * - for j = 1, 2, 3, ...
     *   -- Dj = bj + aj * D_j-1, if Dj==0 then Dj = eps
     *   -- Cj = bj + aj / C_j-1, if Cj==0 then Cj = eps
     *   -- Dj = 1 / Dj
     *   -- Delta_j = Cj * Dj
     *   -- fj = f_j-1 * Delta_j
     *   -- if abs(Delta_j-1) < TOL then exit for loop
     * - return fj
     */

    try
    {
        // f0 = b0
        T f = ctdf.fb(0);

        // adjust f0 to eps if necessary
        if ( true == math::NumericUtil::isZero<T>(f) )
        {
            f = math::NumericUtil::getEPS<T>();
        }

        // c0 = f0,  d0 = 0
        T c = f;
        T d = static_cast<T>(0);

        // Initially Delta should not be equal to 1
        T Delta = static_cast<T>(0);

        size_t j = 1;
        for ( 
              j=1; 
              false == math::NumericUtil::isZero<T>(Delta-static_cast<T>(1), tol) &&
                     j <= SPECFUN_MAX_ITER ; 
              ++j )
        {
            // obtain 'aj' and 'bj'
            const T a = ctdf.fa(j);
            const T b = ctdf.fb(j);

            // dj = bj + aj * d_j-1
            d = b + a * d;
            // adjust dj to eps if necessary
            if ( true == math::NumericUtil::isZero<T>(d) )
            {
                d = math::NumericUtil::getEPS<T>();
            }

            // cj = bj + aj/c_j-1
            c = b + a /c;
            // adjust cj to eps if necessary
            if ( true == math::NumericUtil::isZero(c) )
            {
                c = math::NumericUtil::getEPS<T>();
            }

            // dj = 1 / dj
            d = static_cast<T>(1) / d;

            // Delta_j = cj * dj
            Delta = c * d;

            // fj = f_j-1 * Delta_j
            f *= Delta;

            // for loop's condition will check, if abs(Delta_j-1)
            // is less than the tolerance
        }

        // check if the algorithm has converged:
        if ( j >= SPECFUN_MAX_ITER )
        {
            throw math::SpecFunException(math::SpecFunException::NO_CONVERGENCE);
        }

        // ... if yes, return the fj
        return f;
    }
    catch ( const math::FunctionException& fex )
    {
        throw math::SpecFunException(math::SpecFunException::UNDEFINED);
    }
}
