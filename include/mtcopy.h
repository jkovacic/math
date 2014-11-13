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

#ifndef _MATH_MTCOPY_H_
#define _MATH_MTCOPY_H_

#include "util/mtcopy.hpp"

#endif  /*  _MATH_MTCOPY_H_  */
