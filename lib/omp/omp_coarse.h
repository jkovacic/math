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
 * Useful inline functions that facilitate
 * application of coarse-grained parallelism.
 */


#ifndef _OMP_COARSE_H_
#define _OMP_COARSE_H_


#include <cstddef>

#include "../settings/omp_settings.h"


/**
 * Ideal number of threads (if permitted by the number of
 * available CPU cores) to perform a coarse-grained
 * parallelized task
 *
 * @param n - number of all elements
 * @param per_thread - ideal number of elements to be processed by a single thread
 *
 * @return ideal number of threads
 */
inline std::size_t ompIdealNrThreads( const std::size_t n, const std::size_t per_thread )
{
    return (n / per_thread) + ( 0 == (n % per_thread) ? 0 : 1 );
}


/**
 * Same as 'ompIdealNrThreads', assumes the application default value
 * for the ideal number of elements to be processed by a single thread
 * (OMP_CHUNKS_PER_THREAD, defined in settings/omp_settings.h).
 *
 * @param n - number of all elements
 *
 * @return ideal number of threads for the default number of elements per thread (OMP_CHUNKS_PER_THREAD)
 *
 * @see ompIdealNrThreads
 */
inline std::size_t ompIdeal( const std::size_t n )
{
    return ompIdealNrThreads(n, OMP_CHUNKS_PER_THREAD);
}



/**
 * A convenience macro (unfortunately C and C++ do not allow a more elegant
 * way to accomplish this) that initializes necessary variables
 * used at coarse level parallelization.
 *
 * It declares (as "const std::size_t") and initializes the following variables
 * that can be further used by parallelized algorithms:
 * - thrnr: number of the current thread
 * - nthreads: number of all allocated threads
 * - elems_per_thread: maximum number of elements processed by a thread
 * - istart: index to the first element of the range, processed by the current thread
 * - iend: index to the final element of the range, processed by the current thread
 *
 * The macro can only be used within the "#pragma omp parallel" blocks
 * that actually implement coarse grained parallelization!
 *
 * The following headers are assumed to be included beforehand:
 * - #include <cstddef> (for the definition of size_t)
 * - #include <algorithm> (for the definition of std::min)
 * - #include "omp/omp_header.h" (for the definition of omp_get_thread_num and omp_get_num_threads)
 *
 * @param N - number of all elements to be processed by threads (must be a size_t value)
 */
#define OMP_COARSE_GRAINED_PAR_INIT_VARS( N )                                         \
    const std::size_t thrnr = omp_get_thread_num();                                   \
    const std::size_t nthreads = omp_get_num_threads();                               \
                                                                                      \
    const std::size_t elems_per_thread =                                              \
          ( std::min<std::size_t>( (N), static_cast<std::size_t>(-1) - nthreads + 1 ) \
            + nthreads - 1) / nthreads;                                               \
                                                                                      \
    const std::size_t istart = elems_per_thread * thrnr;                              \
    const std::size_t iend =                                                          \
          std::min<std::size_t>( ( N ),                                               \
              ( elems_per_thread <= (static_cast<std::size_t>(-1)-istart) ?           \
                istart + elems_per_thread :                                           \
                static_cast<std::size_t>(-1) ) );


#endif  /*  _OMP_COARSE_H_  */
