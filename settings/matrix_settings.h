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
 * Default settings for implemented algorithms that compute
 * various matrix properties, such as inverse, determinat or rank.
 */

#ifndef _MATH_MATRIX_SETTINGS_H_
#define _MATH_MATRIX_SETTINGS_H_


/*
 * Default settings, specifying whether internal pivoting
 * implementation for the Gauss - Jordan algorithm based operations
 * performs full pivoting. Note that full pivoting is more complex
 * but numerically more stable.
 */
#define MATRIX_INVERSE_FULL_PIVOT             ( true )
#define MATRIX_DET_FULL_PIVOT                 ( true )

/**
 * Default setting, specifying whether the internal pivoting
 * implementation for the Gauss - Jordan algorithm based operations
 * (such as matrix inversion, determinant and rank) internally physically
 * swap elements of the coefficient matrix when pivoting rows
 * and/or columns. When set to TRUE, the algorithm requires less
 * additional storage for book keeping of swaps.
 */
#define MATRIX_PHYSSWAP_COEF                  ( false )


#endif  /* _MATH_MATRIX_SETTINGS_H_ */
