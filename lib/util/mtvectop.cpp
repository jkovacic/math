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
 * Implementation of functions that perform basic operations on
 * vector elements, such as addition/subtraction of corresponding
 * elements or multiplication of all elements by the same scalar value.
 */


// no #include "mtvectop.hpp" !!!
#include <vector>
#include <cstddef>
#include <algorithm>

#include "../settings/omp_settings.h"
#include "omp/omp_header.h"
#include "omp/omp_coarse.h"


/**
 * Adds or subtracts (depending on 'add') the corresponding elements
 * of two vectors.
 *
 * If 'v1' and 'v2' are of different sizes, the operation will only
 * perform on the number of elements (N) that corresponds to the size
 * of the smaller vector.
 *
 * If dest's preallocated size is smaller than N, it will be resized to N.
 * If its size is larger than or same as N, it will not be resized.
 *
 * @param v1 - the first vector
 * @param v2 - the second vector
 * @param dest - reference of a vector where sums/differences will be written to
 * @param add - if 'true', addition will be performed, subtraction if 'false' (default: true)
 */
template <class T>
void math::mtvectadd(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<T>& dest, bool add=true)
{
    const size_t N = std::min( v1.size(), v2.size() );

    // Only resize if 'dest' is too small
    if ( dest.size() < N )
    {
        dest.resize(N);
    }

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest, v1, v2, add)
    {
        // Depending on the number of available threads,
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will copy.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thrnr;
        const size_t iend = std::min(istart + elems_per_thread, N);

        typename std::vector<T>::const_iterator it1 = v1.begin() + istart;
        typename std::vector<T>::const_iterator it2 = v2.begin() + istart;
        typename std::vector<T>::iterator it = dest.begin() + istart;
        for ( size_t i = istart;
              i<iend && it1!=v1.end() && it2!=v2.end() && it!=dest.end();
              ++it1, ++it2, ++it, ++i )
        {
            if ( &v1 != &dest )
            {
                *it = (true==add ? *it1 + *it2 : *it1 - *it2 );
            }
            else
            {
                // if 'v1' is actually the same vector as 'dest',
                // rather perform += or -= operations
                if ( true == add )
                {
                    *it += *it2;
                }
                else
                {
                    *it -= *it2;
                }
            }
        }  // for
    } // omp parallel
}


/**
 * Multiplies each vector's ('v1') element by a the same scalar value
 * and stores the products into vector 'dest'.
 *
 * If dest's preallocated size is larger than v1's, it will be resized
 * to v1's size. If it is larger than or of the same size as v1,
 * it will not be resized.
 *
 * @param v1 - vector whose elements will be multiplied by the scalar
 * @param scalar - scalar to multiply each vector's element
 * @param dest - reference of a vector where products will be written to
 */
template <class T>
void math::mtvectmult(const std::vector<T>& v1, const T& scalar, std::vector<T>& dest)
{
    const size_t N = v1.size();

    // Only resize if 'dest' is too small
    if ( dest.size() < N )
    {
        dest.resize(N);
    }

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest, v1, scalar)
    {
        // Depending on the number of available threads,
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will copy.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thrnr;
        const size_t iend = std::min(istart + elems_per_thread, N);

        typename std::vector<T>::const_iterator srcit = v1.begin() + istart;
        typename std::vector<T>::iterator destit = dest.begin() + istart;
        for ( size_t i = istart;
              i<iend && srcit!=v1.end() && destit!=dest.end();
              ++srcit, ++destit, ++i)
        {
            if ( &v1 != &dest )
            {
                *destit = *srcit * scalar;
            }
            else
            {
                // if 'v1' is actually the same vector as 'dest',
                // rather perform *= operation
            	*destit *= scalar;
            }
        }  // for
    }  // omp parallel
}


/**
 * Adds or subtracts (depending on 'add') each vector's element by a scalar
 * and stores the results into 'dest'.
 *
 * If dest's preallocated size is smaller than N, it will be resized to N.
 * If its size is larger than or same as N, it will not be resized.
 *
 * Combinations of 'add' and 'vectFirst':
 *
 * - add == true, vectFirst == true:
 *        dest[i] = v1[i] + scalar
 * - add == false, vectFirst == true:
 *        dest[i] = v1[i] - scalar
 *
 * - add == true, vectFirst == false:
 *        dest[i] = scalar + v1[i]
 *
 * - add == false, vectFirst == false:
 *        dest[i] = scalar - v1[i]
 *
 * @param v1 - vector whose elements will be added/subtracted to/from scalar
 * @param scalar - scalar to be added/subtracted to/from the vector
 * @param dest - reference of a vector where sums/differences will be written to
 * @param add - if 'true', addition will be performed, subtraction if 'false' (default: true)
 * @param vectFirst - if 'true', the first element of the desired operation is the vector's element,otherwise the scalar (default: true)
 */
template <class T>
void math::mtvectscalaradd(
        const std::vector<T>& v1,
        const T& scalar,
        std::vector<T>& dest,
        bool add,
        bool vectFirst )
{
    const size_t N = v1.size();

    // Only resize if 'dest' is too small
    if ( dest.size() < N )
    {
        dest.resize(N);
    }

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest, v1, scalar, add, vectFirst)
    {
        // Depending on the number of available threads,
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will copy.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thrnr;
        const size_t iend = std::min(istart + elems_per_thread, N);

        typename std::vector<T>::const_iterator srcit = v1.begin() + istart;
        typename std::vector<T>::iterator destit = dest.begin() + istart;
        for ( size_t i = istart;
              i<iend && srcit!=v1.end() && destit!=dest.end();
              ++srcit, ++destit, ++i)
        {
            if ( &v1 != &dest || false == vectFirst )
            {
                if ( true == vectFirst )
                {
                	if ( true == add )
                    {
                        *destit = *srcit + scalar;
                    }
                    else
                    {
                        *destit = *srcit - scalar;
                	}
                }
                else
                {
                    // vectFirst == false
                    if ( true == add )
                    {
                        *destit = scalar + *srcit;
                    }
                    else
                    {
                        *destit = scalar - *srcit;
                    }
                }
            }
            else
            {
                // if 'v1' is actually the same vector as 'dest',
                // rather use compound assignment operators.
                // Note that 'vectFirst' can only equal true.
               	if ( true == add )
               	{
               	    *destit += scalar;
               	}
               	else
               	{
               	    *destit -= scalar;
               	}
            }
        }  // for
    }  // omp parallel
}
