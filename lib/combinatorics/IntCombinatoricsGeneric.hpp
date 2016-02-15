/*
Copyright 2013, Jernej Kovacic

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
 * @headername{IntCombinatoricsGeneric.h}
 *
 * Declaration of the namespace IntCombinatorics with functions 
 * for calculation of factorials, binomial coefficients, etc.
 */

#ifndef _MATH_INTCOMBINATORICSGENERIC_HPP_
#define _MATH_INTCOMBINATORICSGENERIC_HPP_

#include "exception/CombinatoricsException.hpp"

namespace math
{

/**
 * @brief A namespace with some combinatorics related functions for calculation of
 *        various factorials, binomial coefficients, etc.
 * 
 * @note 'I' is expected to be an integral integer type, preferably unsigned. In
 *       any other case the behaviour of functions from this namespace is
 *       unpredictable.
 */
namespace IntCombinatorics
{

    // factorial: N! or from * (from+1) * ... * N:
    template <typename I>
    I factorial(
            const I& N, 
            const I& from = static_cast<I>(1) ) 
          throw(CombinatoricsException);
    
    // falling factorial power: N * (N-1) * ... * (N-K+1):
    template <typename I>
    I fallingFactorial(
            const I& N, 
            const I& K ) 
          throw(CombinatoricsException);
    
    // rising factorial power: N * (N+1) * ... * (N+K-1):
    template <typename I>
    I risingFactorial(
            const I& N, 
            const I& K ) 
          throw(CombinatoricsException);
    
    // multifactorial: N * (N-K) * (N-2*K) * ... :
    template <typename I>
    I multiFactorial(
            const I& N, 
            const I& K ) 
          throw(CombinatoricsException);
    
    // double factorial: N * (N-2) * (N-4) * ... :
    template <typename I>
    I doubleFactorial( const I& N ) throw(CombinatoricsException);
    
    // binomial coefficient:
    template <typename I>
    I binom(
            const I& N, 
            const I& K ) 
          throw(CombinatoricsException);

}  // namespace IntCombinatorics

}  // namespace math


// DEFINITION
#include "combinatorics/IntCombinatoricsGeneric.cpp"

#endif	// _MATH_INTCOMBINATORICSGENERIC_HPP_
