/*
Copyright 2016, Jernej Kovacic

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
 * Settings for numerical calculations in probability distribution algorithms,
 * e.g. tolerances for calculation of distributions' probabilities and quantiles, etc.
 */


#ifndef _MATH_PROBDIST_SETTINGS_H_
#define _MATH_PROBDIST_SETTINGS_H_


/**
 * Tolerance for numerical algorithms that evaluate probability distributions'
 * probabilities, quantiles, etc.
 * 
 * Tolerance is specified by its numerator and denominator. This ensures
 * more flexibility for generic (templated) statistical functions.
 * 
 * Note: typically statistical functions do not require very high accuracy.
 */
#define PROBDIST_PROB_TOL_NUM                ( 1 )
#define PROBDIST_PROB_TOL_DEN                ( 1000000 )


/**
 * Default threshold for products n*p and n*(1-p)
 * to determine whether a binomial distribution can be
 * approximated as nearly normal.
 *
 * The parameter is specified by its numerator and denominator. This ensures
 * more flexibility for generic (templated) statistical functions.
 */
#define PROBDIST_BINOM_DIST_NORMAL_APPROX_NUM     ( 10 )
#define PROBDIST_BINOM_DIST_NORMAL_APPROX_DEN     ( 1 )


#endif  /* _MATH_PROBDIST_SETTINGS_H_ */
