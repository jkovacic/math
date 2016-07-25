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
void math::mtvectadd(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<T>& dest, const bool add=true)
{
    const std::size_t N = std::min<std::size_t>( v1.size(), v2.size() );

    // Only resize if 'dest' is too small
    if ( dest.size() < N )
    {
        dest.resize(N);
    }

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest, v1, v2)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<T>::const_iterator it1 = v1.begin() + istart;
        typename std::vector<T>::const_iterator it2 = v2.begin() + istart;
        typename std::vector<T>::iterator it = dest.begin() + istart;
        for ( std::size_t i = istart;
              i<iend && it1!=v1.end() && it2!=v2.end() && it!=dest.end();
              ++it1, ++it2, ++it, ++i )
        {
            T& destCurr = *it;
            const T& v1cur = *it1;
            const T& v2cur = *it2;

            if ( &v1 != &dest )
            {
                destCurr = (true==add ? v1cur + v2cur : v1cur - v2cur );
            }
            else
            {
                // if 'v1' is actually the same vector as 'dest',
                // rather perform += or -= operations
                if ( true == add )
                {
                    destCurr += v2cur;
                }
                else
                {
                    destCurr -= v2cur;
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
    const std::size_t N = v1.size();

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
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<T>::const_iterator srcit = v1.begin() + istart;
        typename std::vector<T>::iterator destit = dest.begin() + istart;
        for ( std::size_t i = istart;
              i<iend && srcit!=v1.end() && destit!=dest.end();
              ++srcit, ++destit, ++i)
        {
            T& currDest = *destit;
            const T& currSrc = *srcit;

            if ( &v1 != &dest )
            {
                currDest = currSrc * scalar;
            }
            else
            {
                // if 'v1' is actually the same vector as 'dest',
                // rather perform *= operation
                currDest *= scalar;
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
        const bool add,
        const bool vectFirst )
{
    const std::size_t N = v1.size();

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
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<T>::const_iterator srcit = v1.begin() + istart;
        typename std::vector<T>::iterator destit = dest.begin() + istart;
        for ( std::size_t i = istart;
              i<iend && srcit!=v1.end() && destit!=dest.end();
              ++srcit, ++destit, ++i)
        {
            const T& currSrc = *srcit;
            T& currDest = *destit;

            if ( &v1 != &dest || false == vectFirst )
            {
                if ( true == vectFirst )
                {
                    if ( true == add )
                    {
                        currDest = currSrc + scalar;
                    }
                    else
                    {
                        currDest = currSrc - scalar;
                    }
                }
                else
                {
                    // vectFirst == false
                    if ( true == add )
                    {
                        currDest = scalar + currSrc;
                    }
                    else
                    {
                        currDest = scalar - currSrc;
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
                    currDest += scalar;
                }
                else
                {
                    currDest -= scalar;
                }
            }
        }  // for
    }  // omp parallel
}


/**
 * Multiplies or divides (depending on 'mult') the corresponding elements
 * of two vectors.
 *
 * If 'v1' and 'v2' are of different sizes, the operation will only
 * perform on the number of elements (N) that corresponds to the size
 * of the smaller vector.
 *
 * If dest's preallocated size is smaller than N, it will be resized to N.
 * If its size is larger than or same as N, it will not be resized.
 *
 * If any division by zero is attempted, it will be skipped and no
 * exception will be thrown. Additionally the function will return false.
 *
 * @param v1 - the first vector
 * @param v2 - the second vector
 * @param dest - reference of a vector where products/quotients will be written to
 * @param mult - if 'true', multiplication will be performed, division if 'false' (default: true)
 *
 * @return 'true' if all operations were successful, 'false' if any division by 0 was attempted
 */
template <class T>
bool math::mtvectewmult(
        const std::vector<T>& v1,
        const std::vector<T>& v2,
        std::vector<T>& dest,
        const bool mult )
{
    bool retVal = true;

    const std::size_t N = std::min<std::size_t>( v1.size(), v2.size() );

    // Only resize if 'dest' is too small
    if ( dest.size() < N )
    {
        dest.resize(N);
    }

    // Coarse grained parallelism, if OpenMP is enabled
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest, v1, v2, retVal)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<T>::const_iterator it1 = v1.begin() + istart;
        typename std::vector<T>::const_iterator it2 = v2.begin() + istart;
        typename std::vector<T>::iterator it = dest.begin() + istart;
        for ( std::size_t i = istart;
              i<iend && it1!=v1.end() && it2!=v2.end() && it!=dest.end();
              ++it1, ++it2, ++it, ++i )
        {
            const T& v1cur = *it1;
            const T& v2cur = *it2;
            T& currDest = *it;

            // Prevent possible division by zero
            if ( false == mult &&
                 true==math::NumericUtil::isZero<T>(v2cur) )
            {
                retVal = false;
                // TODO what to do if division by 0 is attempted???
                continue;  // for i
            }

            if ( &v1 != &dest )
            {
                currDest = (true==mult ? v1cur * v2cur : v1cur / v2cur );
            }
            else
            {
                // if 'v1' is actually the same vector as 'dest',
                // rather perform *= or /= operations
                if ( true == mult )
                {
                    currDest *= v2cur;
                }
                else
                {
                    currDest /= v2cur;
                }
            }
        }  // for
    } // omp parallel

    return retVal;
}
