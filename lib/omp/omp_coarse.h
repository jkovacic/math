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
 *
 * @note The header can also be included into C code.
 */


#ifndef _OMP_COARSE_H_
#define _OMP_COARSE_H_

#include <stddef.h>

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
inline size_t ompIdealNrThreads( size_t n, size_t per_thread )
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
inline unsigned long int ompIdeal( unsigned long int n )
{
    return ompIdealNrThreads(n, OMP_CHUNKS_PER_THREAD);
}

#endif  /*  _OMP_COARSE_H_  */
