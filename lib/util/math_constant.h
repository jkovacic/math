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
 * The header should only be used internally.
 *
 * High precision (38 digits) values of frequently used
 * mathematical constants. The constants are defined as
 * long double values and are expected to be statically
 * casted to the desired type (e.g. float or double).
 *
 * @note The header can also be included into C code.
 */


#ifndef _MATH_MATH_CONSTANT_H_
#define _MATH_MATH_CONSTANT_H_


/* pi */
#define MATH_CONST_PI               ( 3.1415926535897932384626433832795028842L )

/* 1/pi */
#define MATH_CONST_INV_PI           ( 0.31830988618379067153776752674502872407L )

/* sqrt(pi) */
#define MATH_CONST_SQRT_PI          ( 1.7724538509055160272981674833411451828L )

/* sqrt(1/pi) */
#define MATH_CONST_SQRT_INV_PI      ( 0.56418958354775628694807945156077258585L )

/* sqrt(2/pi) */
#define MATH_CONST_SQRT_2_DIV_PI    ( 0.79788456080286535587989211986876373695L )

/* sqrt( 1/(2*pi) ) */
#define MATH_CONST_SQRT_INV_2_PI    ( 0.39894228040143267793994605993438186848L )

/* ln( 2*pi ) */
#define MATH_CONST_LOG_SQRT_2_PI    ( 0.91893853320467274178032973640561763986L )




/* sqrt(2) */
#define MATH_CONST_SQRT_2           ( 1.4142135623730950488016887242096980786L )

/* sqrt(2)/2 = sqrt(1/2) */
#define MATH_CONST_SQRT_INV_2       ( 0.70710678118654752440084436210484903929L )



/*
 * All constants were obtained by PARI/GP that uses 38 characters
 * as its default precision:
 *
 ? Pi
 %1 = 3.1415926535897932384626433832795028842
 ? 1/Pi
 %2 = 0.31830988618379067153776752674502872407
 ? sqrt(Pi)
 %3 = 1.7724538509055160272981674833411451828
 ? sqrt(1/Pi)
 %4 = 0.56418958354775628694807945156077258585
 ? sqrt(2/Pi)
 %5 = 0.79788456080286535587989211986876373695
 ? sqrt(1/(2*Pi))
 %6 = 0.39894228040143267793994605993438186848
 ? log(sqrt(2*Pi))
 %7 = 0.91893853320467274178032973640561763986
 ? sqrt(2)
 %8 = 1.4142135623730950488016887242096980786
 ? sqrt(2)/2
 %9 = 0.70710678118654752440084436210484903929

 */


#endif  /*  _MATH_MATH_CONSTANT_H_  */
