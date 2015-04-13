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
 * Declaration of functions that perform swapping of
 * vector's elements.
 *
 * Since vectors are contiguous containers, it is possible to
 * distribute large amounts of swaps among multiple threads
 * and thus accelerate the task.
 *
 * @note All thrown exceptions (e.g. std:bad_alloc, std::out_of_range, etc.)
 *       are transferred to caller functions.
 */


#ifndef _MATH_MTSWAP_HPP_
#define _MATH_MTSWAP_HPP_

#include <vector>

namespace math
{


template <class T>
void mtswap(
        const typename std::vector<T>::iterator& first,
        const typename std::vector<T>::iterator& last,
        const typename std::vector<T>::iterator& destit );

}  // namespace math


// DEFINITION
#include "util/mtswap.cpp"


#endif  // _MATH_MTSWAP_HPP_
