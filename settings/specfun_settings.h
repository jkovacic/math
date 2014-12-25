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
 * Default tolerance for numerical calculations in algorithms for
 * calculation of special functions.
 */


#ifndef _MATH_SPECFUN_SETTINGS_H_
#define	_MATH_SPECFUN_SETTINGS_H_


/**
 * Default tolerance for numerical algorithms that evaluate
 * special functions.
 *
 * The tolerance is specified by its numerator and denominator. This ensures
 * more flexibility for generic (templated) functions.
 */
#define SPECFUN_TOL_NUM                       ( 1 )
#define SPECFUN_TOL_DEN                       ( 1000000 )


/**
 * Maximum number of iterations
 */
#define SPECFUN_MAX_ITER                      ( 10000 )


#endif  /* _MATH_SPECFUN_SETTINGS_H_ */
