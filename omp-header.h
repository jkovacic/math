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
 * Otherwise some OpenMP functions will be "declared" as macros that "return"
 * appropriate values for single threaded environments.
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
     */
    #include <omp.h>

#else

   /*
    * Otherwise "declare" necessary OpenMP commands as macros that "return"
    * typical values for single threaded situations, e.g. 1 thread,
    * its thread number equaling 0, etc.
    */
   #define omp_get_num_threads()                 ( 1 )
   #define omp_get_thread_num()                  ( 0 )

#endif  /* _OPENMP */

#endif  /* OMP_HEADER_H_ */
