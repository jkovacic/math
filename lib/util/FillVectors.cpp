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
 * Implementation of auxiliary functions that fill integer vector(s) with
 * elements' consecutive indices.
 */


#include <cstddef>
#include <vector>

#include "../settings/omp_settings.h"
#include "omp/omp_header.h"
#include "omp/omp_coarse.h"

#include "util/FillVectors.h"


// An anonymous namespace for "private" functions:
namespace math {  namespace util {  namespace
{
    

/*
 * Preallocates and fills vectors 'v1' and optionally 'v2' with consecutive
 * values between 0 and N-1.
 * 
 * @note The function attempts to resize 'v1' and optionally 'v2', however it does
 *       not handle the exception std::bad_alloc (thrown when not enough memory is
 *       available). The exception is transferred to the caller. It is recommended
 *       that 'v1' and optionally 'v2' are preallocated prior to calling this function.
 * 
 * @param N - size of vector(s)
 * @param v1 - first vector
 * @param v2 - second vector (ignored if 'bothVectors' is false)
 * @param bothVectors - if false, only 'v1' is preallocated and filled with initial positions
 */
void fillVectorsWithInitialConsecutiveIndices(
        const std::size_t N,
        std::vector<std::size_t>& v1,
        std::vector<std::size_t>& v2,
        const bool bothVectors
      )
{
    v1.resize(N);

    if ( true == bothVectors )
    {
        v2.resize(N);
    }
        
    /*
     * Fill 'v1' and optionally 'v2' with initial positions.
     */
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(v1, v2)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<std::size_t>::iterator it1 = v1.begin() + istart;
        typename std::vector<std::size_t>::iterator it2 = v2.begin() + istart;
        std::size_t i = istart;
        for ( std::size_t cntr=0;
              cntr<elems_per_thread && it1!=v1.end();
              ++cntr, ++i, ++it1 )
        {
            std::size_t& v1cur = *it1;
            std::size_t& v2cur = *it2;

            v1cur = i;
            if ( true == bothVectors )
            {
                v2cur = i;
                ++it2;
            }
        }

        (void) iend;
    }  // omp parallel

}

}}} // anonymous namespace



/**
 * Preallocates and fills vectors 'v1' and 'v2' with consecutive values
 * between 0 and N-1.
 * 
 * @note The function attempts to resize 'v1' and 'v2', however it does not handle
 *       the exception std::bad_alloc (thrown when not enough memory is available).
 *       The exception is transferred to the caller. It is recommended that 'v1'
 *       and 'v2' are preallocated prior to calling this function.
 * 
 * @param N - size of vector
 * @param v1 - first vector to fill
 * @param v2 - second vector to fill
 */
void math::util::fillVectorsWithConsecutiveIndices(
        const std::size_t N, 
        std::vector<std::size_t>& v1,
        std::vector<std::size_t>& v2 )
{
    math::util::fillVectorsWithInitialConsecutiveIndices(N, v1, v2, true);
}


/**
 * Preallocates and fills vector 'v' with consecutive values between 0 and N-1.
 * 
 * @note The function attempts to resize 'v', however it does not handle
 *       the exception std::bad_alloc (thrown when not enough memory is available).
 *       The exception is transferred to the caller. It is recommended that 'v'
 *       is preallocated prior to calling this function.
 * 
 * @param N - size of vector
 * @param v - vector to fill
 */
void math::util::fillVectorWithConsecutiveIndices(
        const std::size_t N, 
        std::vector<std::size_t>& v )
{
    math::util::fillVectorsWithInitialConsecutiveIndices(N, v, v, false);
}
