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
 * Settings for numerical calculations in statistics algorithms,
 * e.g. tolerances for calculation of distributions' probabilities
 * and quantiles, etc.
 */


#ifndef _MATH_STAT_SETTINGS_H_
#define	_MATH_STAT_SETTINGS_H_


/**
 * Tolerance for numerical algorithms (e.g. for root finding methods
 * and terms of Taylor series).
 * 
 * Tolerance is specified by its numerator and denominator. This ensures
 * more flexibility for generic (templated) statistical functions.
 * 
 * Note: typically statistical functions do not require very high accuracy.
 */
#define STAT_DIST_PROB_TOL_NUM                ( 1 )
#define STAT_DIST_PROB_TOL_DEN                ( 100000 )


/**
 * Numerical integration and step size for integration algorithms
 * where applicable, e.g. Student's distribution.
 *
 * Step size is specified by its numerator and denominator.
 */
#define STAT_INTEG_ALG               math::EIntegAlg::SIMPSON_3_8
#define STAT_NUM_INTEG_STEP_NUM               ( 1 )
#define STAT_NUM_INTEG_STEP_DEN               ( 10000 )


#endif  /* _MATH_STAT_SETTINGS_H_ */
