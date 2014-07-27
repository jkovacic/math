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
 * @file omp_settings.h
 *
 * Some reasonable default settings for OpenMP.
 *
 * @author Jernej Kovacic
 */


#ifndef _MATH_OMP_SETTINGS_H_
#define _MATH_OMP_SETTINGS_H_


/**
 * If a series of simple operations (such as matrix addition)
 * is attempted, parallelization might not always be the most
 * optimal solution as manipulation of threads also costs something.
 * This macro defines the "threshold" number of operations when
 * parallelization of simple tasks might be sensible.
 */
#define OMP_CHUNKS_PER_THREAD                 ( 100 )


/**
 * In quaternion arithmetics, it rarely makes any sense to parallelize
 * operations into 4 threads (one to calculate each quaternion's element).
 * Typically this is only sensible when T (type of elements) is a class with
 * complex (slow) basic operations, such as '+', '-' and '*'.
 *
 * If this macro is set to 0, parallelization of quaternion's arithmetics
 * is disabled. To enable parallelization, set this macro to any non-zero value.
 */
#define OMP_QUAT_PARALLELIZE                  ( 0 )

#endif  // _MATH_OMP_SETTINGS_H_
