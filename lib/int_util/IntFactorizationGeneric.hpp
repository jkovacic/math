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
 * @headername{IntFactorizationGeneric.h}
 *
 * Declaration of the namespace IntFactorization with
 * functions for factorization of integers and some 
 * prime number utilities.
 */

#ifndef _MATH_INTFACTORIZATIONGENERIC_HPP_
#define	_MATH_INTFACTORIZATIONGENERIC_HPP_

#include "exception/IntFactorizationException.hpp"

#include <map>
#include <set>

namespace math
{

/**
 * @brief A separate namespace with some integer factorization functions for prime
 * factorization, calculation of the least common multiple, greatest common divisor, etc.
 * 
 * @note 'I' is expected to be an integral integer type, preferably unsigned. In
 *       any other case the behaviour of functions from this namespace is
 *       unpredictable.
 */

namespace IntFactorization
{    

    template <typename I>
    bool isPrime(const I& n);


    // greatest common divisor
    template <typename I>
    I greatestCommonDivisor(
                        const I& first, 
                        const I& second ) 
                    throw(IntFactorizationException);


    // least common multiple
    template <typename I>
    I leastCommonMultiple(
                        const I& first, 
                        const I& second ) 
                    throw(IntFactorizationException);


    template <typename I>
    I nextPrime(const I& n) 
                    throw(IntFactorizationException);


    // the largest integer not exceeding sqrt(n)
    template <typename I>
    I intSqrt(const I& n) throw(IntFactorizationException);
 

    // prime factorization
    template <typename I>
    void factor(
                const I& n, 
                std::map<I, I>& fac ) 
            throw(IntFactorizationException);


    // list of all divisors of 'n'
    template <typename I>
    void divisors(
                const I& n,
                std::set<I>& div )
            throw(IntFactorizationException);

}  // namespace IntFactorization

}  // namespace math

// DEFINITION
#include "int_util/IntFactorizationGeneric.cpp"

#endif	// _MATH_INTFACTORIZATIONGENERIC_HPP_
