/*
Copyright 2015, Jernej Kovacic

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
 * Default settings for implemented algorithms for
 * solving systems of linear equations.
 */

#ifndef _MATH_LINEQ_SETTINGS_H_
#define _MATH_LINEQ_SETTINGS_H_


/**
 * Default maximum number of iterations
 */
#define LINEQ_MAX_ITER                ( 10000 )


/**
 * Default convergence tolerance for iterative algorithms that
 * solve systems of linear equations.
 * 
 * For better flexibility with generic (templated) functions
 * the convergence tolerance is given by its numerator and denominator.
 */
#define LINEQ_TOL_CONV_NUM            ( 1 )
#define LINEQ_TOL_CONV_DEN            ( 1000000 )


/**
 * Default setting, specifying whether the Gauss - Jordan method
 * performs full pivoting by default?
 */
#define LINEQ_GAUSS_JORDAN_FULL_PIVOT       ( true )


/**
 * Default setting, specifying whether the internal pivoting
 * implementation for the Gauss - Jordan method should internally
 * physically swap elements of the coefficient matrix when pivoting
 * rows and/or columns. When set to TRUE, the algorithm requires
 * less additional storage for book keeping of swaps.
 */
#define LINEQ_PHYSSWAP_COEF           ( false )


#endif 	/* _MATH_LINEQ_SETTINGS_H_ */
