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
 * @headername{IntExponentiator.h}
 *
 * Declaration of the namespace IntExponentiator with
 * functions for efficient calculation of powers of
 * positive integer exponents.
 */


#ifndef _MATH_INTEXPONENTIATORGENERIC_HPP_
#define _MATH_INTEXPONENTIATORGENERIC_HPP_

#include "exception/IntFactorizationException.hpp"


namespace math
{

/**
 * @brief A separate namespace with some functions that efficiently
 * perform exponentiation of positive integer exponents.
 */
namespace IntExponentiator
{
    // Efficient calculation of base^n:
    template <class T, typename I>
    T power(const T& base, const I& n) throw (IntFactorizationException);

}  // namespace IntExponentiator

}  // namespace math

// DEFINITION
#include "int_util/IntExponentiatorGeneric.cpp"

#endif  // _MATH_INTEXPONENTIATORGENERIC_HPP_
