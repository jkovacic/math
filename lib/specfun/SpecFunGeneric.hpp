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

#include "exception/SpecFunException.hpp"


namespace math
{

namespace SpecFun
{

    template <class T>
    T gamma(const T& x) throw (SpecFunException);

    template <class T>
    T beta(const T& x, const T& y) throw (SpecFunException);

}  // namespace SpecFun

}  // namespace math

// DEFINITION:
#include "specfun/SpecFunGeneric.cpp"

#endif  // _MATH_SPECFUNGENERIC_HPP_
