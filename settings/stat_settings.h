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
#define _MATH_STAT_SETTINGS_H_


/**
 * Default method to estimate a quantile
 */
#define STAT_DEFAULT_QUANTILE_ALG             math::EQntlType::R7


/**
 * Default number of interquartile ranges (IQRs) below/above
 * Q1 and Q3 to determine whether an element is an outlier.
 * 
 * The parameter is specified by its numerator and denominator. This ensures
 * more flexibility for generic (templated) statistical functions.
 */
#define STAT_OUTLIER_IQRS_NUM                 ( 3 )
#define STAT_OUTLIER_IQRS_DEN                 ( 2 )


/**
 * Default setting, specifying whether indices of some statistical functions
 * are zero or one based.
 */
#define STAT_DEFAULT_ZERO_BASE                true


#endif  /* _MATH_STAT_SETTINGS_H_ */
