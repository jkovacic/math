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
 * @file IntCombinatorics.h
 * 
 * Declaration of the class IntCombinatorics with static
 * functions for calculation of factorials, binomial coefficients, etc.
 * 
 * @author Jernej Kovacic
 */

#ifndef _MATH_INTCOMBINATORICS_H_
#define	_MATH_INTCOMBINATORICS_H_

#include "CombinatoricsException.h"

namespace math
{

struct IntCombinatorics
{
public:
    // factorial: N! or from * (from+1) * ... * N:
    static unsigned long long int factorial(unsigned long long int N, unsigned long long int from=1LL) throw(CombinatoricsException);
    
    // falling factorial power: N * (N-1) * ... * (N-K+1):
    static unsigned long long int fallingFactorial(unsigned long long int N, unsigned long long int K) throw(CombinatoricsException);
    
    // rising factorial power: N * (N+1) * ... * (N+K-1):
    static unsigned long long int risingFactorial(unsigned long long int N, unsigned long long int K) throw(CombinatoricsException);
    
    // multifactorial: N * (N-K) * (N-2*K) * ... :
    static unsigned long long int multiFactorial(unsigned long long int N, unsigned int K) throw(CombinatoricsException);
    
    // double factorial: N * (N-2) * (N-4) * ... :
    static unsigned long long int doubleFactorial(unsigned long long int N) throw(CombinatoricsException);
    
    // binomial coefficient:
    static unsigned long long int binom(unsigned long long int N, unsigned long long int K) throw(CombinatoricsException);
};
    
}  // namespace math


#endif	// _MATH_INTCOMBINATORICS_H_
