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
 * Implementation of functions that copy data from vectors or arrays
 * to vectors.
 */

#include <vector>
#include <cstddef>
#include <algorithm>

#include "../settings/omp_settings.h"
#include "omp/omp_header.h"
#include "omp/omp_coarse.h"


/**
 * Copies the selected range of an array to a vector.
 *
 * @note The function changes the actual content of 'dest'!
 *
 * @param first - pointer to the initial position of the array range
 * @param last - pointer to the final position of the range (will not be copied)
 * @param dest - reference to a vector where the selected range will be copied to
 */
template <class T>
void math::mtcopy(const T* const first, const T* const last, std::vector<T>& dest)
{
    // Number of elements in the range
    const size_t N = last - first;

    // Preallocate the destination vector
    dest.resize(N);

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(dest)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<T>::iterator it = dest.begin() + istart;
        for ( size_t idx = istart;
              idx<iend && it!=dest.end();
              ++it, ++idx )
        {
            *it = *(first+idx);
        }
    }
}


/**
 * Copies the selected range of an array to a vector.
 *
 * @note The function changes the actual content of 'dest'!
 *
 * @param first - pointer to the initial position of the array range
 * @param len - number of elements in the array range
 * @param dest - reference to a vector where the selected range will be copied to
 */
template <class T>
void math::mtcopy(const T* const first, const size_t len, std::vector<T>& dest)
{
    math::mtcopy<T>(first, first+len, dest);
}


/**
 * Copies a vector to another vector.
 *
 * @note The function changes the actual content of 'dest'!
 *
 * @param src - vector to be copied
 * @param dest - reference to a vector where the vector will be copied to
 */
template <class T>
void math::mtcopy(const std::vector<T>& src, std::vector<T>& dest)
{
    math::mtcopy<T>(src.begin(), src.end(), dest);
}


/**
 * Copies the selected range of a vector to another vector.
 *
 * @note The function changes the actual content of 'dest'!
 *
 * @param src - vector to be copied
 * @param first - position of the initial element of the selected range
 * @param len - maximum number of elements in the selected range (may be decreased if exceeds src's size)
 * @param dest - reference to a vector where the selected range will be copied to
 */
template <class T>
void math::mtcopy(const std::vector<T>&src,
                  const size_t first,
                  const size_t len,
                  std::vector<T>& dest)
{
    // First make sure that "first+len" does not exceed size_t's range
    const size_t lenlim = std::min<size_t>(len, static_cast<size_t>(-1)-first);

    // Then ensure that the range does not exceed
    // the actual src's range:
    math::mtcopy<T>(src.begin() + first,
                    src.begin() + std::min<size_t>( first+lenlim, src.size() ),
                    dest);
}


/**
 * Copies the selected range between vector's iterators to another vector.
 *
 * @note The function changes the actual content of 'dest'!
 *
 * @note Nothing will be done if 'first' is located after 'last'.
 * 
 * @param first - iterator to the initial element of the selected range
 * @param last - iterator to the final position of the range (will not be copied)
 * @param dest - reference to a vector where the selected range will be copied to
 */
template <class T>
void math::mtcopy(const typename std::vector<T>::const_iterator& first,
                  const typename std::vector<T>::const_iterator& last,
                  std::vector<T>& dest)
{
    /*
     * 'first' and 'last' are instances of a Random Access Iterator
     * with numerous properly defined operators
     */

    /*
     * Sanity check.
     * If 'first' is located after 'last', nothing will be done.
     * 
     * Note that inequality comparison operators exist
     * for random access iterators, incl. std::vector::const_iterator.
     */
    if ( last <= first )
    {
        return;
    }

    // Number of elements in the range
    const size_t N = last - first;

    // Preallocate the dest. vector
    dest.resize(N);

    // Coarse grained parallelism if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(dest, first)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        // Iterator to the final element of the source block
        const typename std::vector<T>::const_iterator final = first + iend;

        // iterator to the current element of 'src'
        typename std::vector<T>::const_iterator idx = first + istart;

        typename std::vector<T>::iterator it = dest.begin() + istart;
        for ( ;
              idx != final && it != dest.end();
              ++it, ++idx )
        {
            *it = *idx;
        }
    }
}
