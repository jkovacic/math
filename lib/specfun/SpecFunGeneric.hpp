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
 * @headername{SpecFunGeneric.h}
 *
 * Declaration of special functions in the namespace SpecFun.
 */

#ifndef _MATH_SPECFUNGENERIC_HPP_
#define _MATH_SPECFUNGENERIC_HPP_

#include "../settings/specfun_settings.h"

#include "exception/SpecFunException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that evaluate several
 *        special functions.
 */
namespace SpecFun
{

    template <class T>
    T gamma(const T& x) throw (SpecFunException);

    template <class T>
    T beta(const T& x, const T& y) throw (SpecFunException);


    template <class T>
    T incGammaUpper(
               const T& a,
               const T& x,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaLower(
               const T& a,
               const T& x,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaUpperReg(
               const T& a,
               const T& x,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaLowerReg(
               const T& a,
               const T& x,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaLowerInv(
               const T& a,
               const T& g,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaUpperInv(
               const T& a,
               const T& g,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaLowerRegInv(
               const T& a,
               const T& g,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incGammaUpperRegInv(
               const T& a,
               const T& g,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaLower(
               const T& x,
               const T& a,
               const T& b,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaUpper(
               const T& x,
               const T& a,
               const T& b,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaLowerReg(
               const T& x,
               const T& a,
               const T& b,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaUpperReg(
               const T& x,
               const T& a,
               const T& b,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaLowerInv(
               const T& a,
               const T& b,
               const T& y,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaUpperInv(
               const T& a,
               const T& b,
               const T& y,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaLowerRegInv(
               const T& a,
               const T& b,
               const T& y,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T incBetaUpperRegInv(
               const T& a,
               const T& b,
               const T& y,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T erf(
               const T& x,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             );


    template <class T>
    T erfc (
               const T& x,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             );


    template <class T>
    T erfInv(
               const T& e,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


    template <class T>
    T erfcInv(
               const T& e,
               const T& tol = static_cast<T>(SPECFUN_TOL_NUM)/static_cast<T>(SPECFUN_TOL_DEN)
             ) throw (SpecFunException);


}  // namespace SpecFun

}  // namespace math

// DEFINITION:
#include "specfun/SpecFunGeneric.cpp"

#endif  // _MATH_SPECFUNGENERIC_HPP_
