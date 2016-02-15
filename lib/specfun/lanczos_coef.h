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
 * An internal header file, it should not be included directly.
 * 
 * A header with Lanczos coefficients, necessary to evaluate
 * a gamma function. 
 */


#ifndef _MATH_SPECFUN_LANCZOS_COEF_H_
#define _MATH_SPECFUN_LANCZOS_COEF_H_


/*
 * Lanczos coefficients depend on chosen parameters g and N.
 * They can be calculated by the Paul Godfrey's algorithm as described at:
 * http://mrob.com/pub/ries/lanczos-gamma.html
 * Using Maxima it is possible to calculate the coefficients in an arbitrary precision.
 *
 * When very high accuracy is not necessary, Boost recommends the following values for
 * g and N:
 *   g = 1.428456135094165802001953125
 *   N = 6
 *
 * More details at:
 * http://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/lanczos.html
 *
 * The Maxima script is implemented in 'scripts/lib/lanczos_coef.mac'.
 * For the 25 digit precision, it returns the values below:
 */


/* Lanczos G parameter: */
#define LANCZOS_G                ( 1.428456135094165802001953125L )

/* Number of Lanczos coefficients: */
#define NR_LANCZOS_COEF          ( 6 )


/*
 * Elements in the array LANCZOS_COEF_ARRAY are provided as high precision 
 * (25 digits) long double constants. The source file will typically
 * cast them to the appropriate type. This will be facilititated using
 * so called "for-each macros". Within a macro definition, all values are 
 * enumerated as arguments to another macro, e.g. CAST. Macros that are 
 * replaced by "CAST" are then defined in source files when they are actually
 * casted. For more details about this trick, see:
 * http://cakoose.com/wiki/c_preprocessor_abuse
 */

#define LANCZOS_COEF_ARRAY(CAST) \
    CAST( 1.000000123880027543803602L ) \
    CAST( 0.9727644985297789802780476L ) \
    CAST( 0.00979090078038154694104146L ) \
    CAST( -0.006669128561042791409832734L ) \
    CAST( 0.003802685303772964895033826L ) \
    CAST( -0.001117185821247582043028556L )


#endif  /* _MATH_SPECFUN_LANCZOS_COEF_H_ */
