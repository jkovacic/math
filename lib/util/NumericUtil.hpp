/*
Copyright 2011, Jernej Kovacic

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
 * @headername{NumericUtil.h}
 *
 * Declaration of the namespace NumericUtil with a collection of
 * some useful numerical utilities.
 */

#ifndef _MATH_NUMERICUTIL_HPP_
#define _MATH_NUMERICUTIL_HPP_


#include <cstddef>


namespace math
{

/**
 * @brief A namepace with some frequently used functionality, e.g.
 *        determination whether a numeric value is close enough to
 *        zero, etc.
 */
namespace NumericUtil
{

    // Does the value equal zero? Note that in numerical mathematics
    // (where mostly float or double values are in use), "equals" has a
    // different meaning than in discrete mathematics (int etc.)
    template <class T>
    bool isZero(const T& value);

    template <class T>
    bool isZero(const T& value, const T& eps);

    // Get the system specific epsilon for the specified type
    template <class T>
    T getEPS();
    
    // num's sign:
    template <class T>
    short int sign(const T& num);

    // round from a floating type T to an integer type I 
    template <class T, class I>
    I intRound(const T& n);

    template <class T, class I>
    I intFloor(const T& n);

    template <class T, class I>
    I intCeil(const T& n);

}  // namespace NumericUtil 

}  // namespace math

// DEFINITION
#include "util/NumericUtil.cpp"

#endif	// _MATH_NUMERICUTIL_HPP_
