/*
Copyright 2016, Jernej Kovacic

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
 * @headername{IntFunctionGeneric.h}
 *
 * Declaration of the namespace IntFunction with
 * some elementary functions, "specialized" for integer numbers.
 */

#ifndef _MATH_INTFUNCTIONGENERIC_HPP_
#define _MATH_INTFUNCTIONGENERIC_HPP_

#include "exception/IntFunctionException.hpp"


namespace math
{

/**
 * @brief A separate namespace with some elementary functions,
 *        "specialized" for integer numbers.
 *
 * @note 'I' is expected to be an integral integer type, preferably unsigned. In
 *       any other case the behaviour of functions from this namespace is
 *       unpredictable.
 */

namespace IntFunction
{

    // floor/ceil of sqrt(n)
    template <typename I>
    I intSqrt(const I& n, const bool ceil=false) throw(IntFunctionException);


    // floor/ceil of log2(n)
    template <typename I>
    I intLog2(const I& n, const bool ceil=false) throw(IntFunctionException);

}  // namespace IntFunction

}  // namespace math

// DEFINITION
#include "int_util/IntFunctionGeneric.cpp"

#endif  // _MATH_INTFUNCTIONGENERIC_HPP_
