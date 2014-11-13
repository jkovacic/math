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
 * @headername{mtcopy.h}
 * 
 * Declaration of functions that copy data from vectors or arrays to
 * vectors.
 *
 * Since both arrays and vectors are contiguous containers, it is
 * possible to distribute copying large amounts of data among
 * multiple threads and thus accelerate the task.
 *
 * @note All thrown exceptions (e.g. std:bad_alloc, std::out_of_range, etc.)
 *       are transferred to caller functions.
 */


#ifndef _MATH_MTCOPY_HPP_
#define _MATH_MTCOPY_HPP_

#include <cstddef>
#include <vector>

namespace math
{

template <class T>
void mtcopy(const T* first, const T* last, std::vector<T>& dest);


template <class T>
void mtcopy(const T* first, size_t len, std::vector<T>& dest);


template <class T>
void mtcopy(const std::vector<T>& src, std::vector<T>& dest);


template <class T>
void mtcopy(const std::vector<T>&src,
            size_t first,
            size_t len,
            std::vector<T>& dest);


template <class T>
void mtcopy(const typename std::vector<T>::const_iterator& first,
            const typename std::vector<T>::const_iterator& last,
            std::vector<T>& dest);

}   // namspace math


// Definition of templated functions must follow their declaration
// and cannot be compiled separately
#include "util/mtcopy.cpp"

#endif  // _MATH_MTCOPY_HPP_
