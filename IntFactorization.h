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
 * @file IntFactorization.h
 * 
 * Declaration of the class IntFactorization with static
 * functions for factorization of integers and some 
 * prime number utilities.
 */

#ifndef _MATH_INTFACTORIZATION_H_
#define	_MATH_INTFACTORIZATION_H_

#include "IntFactorizationException.h"

#include <map>
#include <set>

namespace math
{
 
struct IntFactorization
{
public:
    static bool isPrime(unsigned long int n);
    static unsigned long int gcd(unsigned long int first, unsigned long int second) throw(IntFactorizationException);
    static unsigned long int lcm(unsigned long int first, unsigned long int second) throw(IntFactorizationException);
    static unsigned long int nextPrime(unsigned long int n) throw(IntFactorizationException);

    static unsigned long int intSqrt(unsigned long int n);
    static std::map<unsigned long int, unsigned int> factor(unsigned long int n) throw(IntFactorizationException);
    static std::set<unsigned long int> divisors(unsigned long int n) throw(IntFactorizationException);
};

}  // namespace math

#endif	// _MATH_INTFACTORIZATION_H_
