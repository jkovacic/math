/*
Copyright 2015, Jernej Kovacic

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
 *
 * Declaration of various pseudo functions within namespace PseudoFunction
 * that are monotonically increasing (or decreasing) function of their
 * official "counterparts", but more efficient to compute.
 * 
 * @note This header should only be used internally!
 */

#ifndef _MATH_PSEUDOFUNCTIONGENERIC_HPP_
#define _MATH_PSEUDOFUNCTIONGENERIC_HPP_

#include <complex>


namespace math
{

namespace PseudoFunction
{

template <class T>
inline T pabs(const T& x);

template <class T>
inline std::complex<T> pabs(const std::complex<T>& x);


template <class T>
inline bool absgt(const T& a, const T& b);

template <class T>
inline bool absgt(const std::complex<T>& a, const std::complex<T>& b);


template <class T>
inline T pabs2abs(const T& x);

template <class T>
inline std::complex<T> pabs2abs(const std::complex<T>& x);

}  // namespace PseudoFunction

}  // namespace math


// DEFINITION
#include "util/PseudoFunctionGeneric.cpp"

#endif  // _MATH_PSEUDOFUNCTIONGENERIC_HPP_
