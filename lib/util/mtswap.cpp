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
 * Implementation of functions that perform swapping of
 * vector's elements.
 */


#include <cstddef>
#include <vector>
#include <algorithm>

#include "../settings/omp_settings.h"
#include "omp/omp_header.h"
#include "omp/omp_coarse.h"



/**
 * Swaps all elements of two equally sized vector's ranges.
 * The regions can belong either to two different vectors
 * (of the same templated parameter) or to two ranges of
 * the same vector.
 *
 * @note The regions must not overlap!
 *
 * @param first - iterator to the initial element of the selected range
 * @param last - iterator to the final element of the selected range (will not be swapped)
 * @param destit - iterator to the first element of the second range
 */
template <class T>
void math::mtswap(
        const typename std::vector<T>::iterator& first,
        const typename std::vector<T>::iterator& last,
        const typename std::vector<T>::iterator& destit )
{
    // Nothing to do if 'last' is located before 'first'
    if ( last <= first )
    {
        return;
    }

    // Number of elements in the range
    const size_t N = last - first;

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(destit, first)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        // Iterator to the final element of the source block
        const typename std::vector<T>::const_iterator final = first + iend;

        // iterator to the current element of the source block
        typename std::vector<T>::iterator it = first + istart;

        typename std::vector<T>::iterator dest = destit + istart;
        for ( ;
	          it != final;
	          ++it, ++dest )
        {
	        std::swap<T>(*it, *dest);
        }
    }
}
