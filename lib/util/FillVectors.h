/*
Copyright 2016, Jernej Kovacic

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
 * @name
 * @author Jernej Kovacic
 * 
 * An internal header file, it should not be included directly.
 * 
 * Declaration of auxiliary functions that fill integer vector(s) with
 * elements' consecutive indices.
 * 
 * Since vectors are contiguous containers, it is possible to
 * distribute large amounts of assignments among multiple threads
 * and thus accelerate the task.
 * 
 * @note All thrown exceptions (in particular std:bad_alloc) are transferred
 *       to caller functions.
 */


#ifndef _MATH_FILLVECTORS_H_
#define	_MATH_FILLVECTORS_H_

#include <cstddef>
#include <vector>


namespace math
{

namespace util
{


void fillVectorsWithConsecutiveIndices(
        const std::size_t N,
        std::vector<std::size_t>& v1,
        std::vector<std::size_t>& v2 );

void fillVectorWithConsecutiveIndices(
        const std::size_t N,
        std::vector<std::size_t>& v );

}}  // namespace math::util

#endif	// _MATH_FILLVECTORS_H_
