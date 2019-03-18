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
 * Default settings for implemented root finding algorithms.
 */

#ifndef _MATH_ROOTFIND_SETTINGS_H_
#define _MATH_ROOTFIND_SETTINGS_H_


/**
 * Default x and y tolerances for root finding algorithms.
 * 
 * For better flexibility with generic (templated) functions
 * the step size is given by its numerator and denominator.
 */
#define ROOTFIND_EPSX_NUM                  ( 1 )
#define ROOTFIND_EPSX_DEN                  ( 10000 )
#define ROOTFIND_EPSY_NUM                  ( 1 )
#define ROOTFIND_EPSY_DEN                  ( 1000000 )


/**
 * Default maximum number of iterations
 */
#define ROOTFIND_MAX_ITER                  ( 10000 )


/**
 * Default step size for root finding algorithms that
 * require numerical differentiation.
 * 
 * For better flexibility with generic (templated) functions
 * the step size is given by its numerator and denominator.
 */
#define ROOTFIND_DIFF_STEP_NUM             ( 1 )
#define ROOTFIND_DIFF_STEP_DEN             ( 1000 )


/**
 * Default 1st order differentiation method for algorithms
 * that perform numerical differentiation.
 */
#define ROOTFIND_DEFAULT_DIFF_METHOD       math::Diff::EDiffMethod::CENTRAL



#endif  /* _MATH_ROOTFIND_SETTINGS_H_ */
