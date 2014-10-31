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
 * @file omp-header.h
 *
 * A "dummy" header that declares some OpenMP functions regardless whether
 * OpenMP is enabled or not. If it is enabled, "omp.h" is included.
 * Otherwise some OpenMP functions will be defined as short inline functions
 * that return appropriate values for single threaded environments.
 *
 * The purpose of this header is to avoid using numerous #ifdef's within
 * applications that make the code less readable.
 *
 * @note The header can also be included into C code.
 *
 * @author Jernej Kovacic
 */


#ifndef _OMP_HEADER_H_
#define _OMP_HEADER_H_


/*
 * When a compiler is invoked with the correct OpenMP flag,
 * this macro should be predefined automatically. If this is not
 * the case, edit Makefile appropriately.
 */
#ifdef _OPENMP

    /*
     * If the macro is defined, just include the original OpenMP header.
     *
     * Note: typically the original header already does take care of
     * including C code into C++ code, so no 'extern "C"' is necessary.
     */
    #include <omp.h>

#else

    /*
     * Otherwise define necessary OpenMP commands as inline functions that
     * return typical values for single threaded situations, e.g. 1 thread,
     * its thread number equaling 0, etc.
     */


    /*
     * Returns the number of threads in the current team. In a sequential section of
     * the program omp_get_num_threads returns 1.
     */
    inline int omp_get_num_threads(void)
    {
        return 1;
    }

    /*
     * Returns a unique thread identification number within the current team. In a
     * sequential parts of the program, omp_get_thread_num always returns 0. In parallel
     * regions the return value varies from 0 to omp_get_num_threads-1 inclusive. The
     * return value of the master thread of a team is always 0.
     */
    inline int omp_get_thread_num(void)
    {
        return 0;
    }

    /*
     * Return the maximum number of threads used for the current parallel region that does
     * not use the clause num_threads.
     */
    inline int omp_get_max_threads(void)
    {
        return 1;
    }

    /*
     * Specifies the number of threads used by default in subsequent parallel sections,
     * if those do not specify a num_threads clause. The argument of omp_set_num_threads
     * shall be a positive integer.
     */
    inline void omp_set_num_threads(int num_threads)
    {
        /* nothing to do in single threaded mode */
        (void) num_threads;
    }

#endif  /* _OPENMP */

#endif  /* OMP_HEADER_H_ */
