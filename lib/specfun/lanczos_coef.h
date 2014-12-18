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
#define	_MATH_SPECFUN_LANCZOS_COEF_H_


/*
 * Lanczos coefficients depend on chosen parameters g and N.
 * They can be calculated by the Paul Godfrey's algorithm as described at:
 * http://mrob.com/pub/ries/lanczos-gamma.html
 * Using Maxima it is possible to calculate them in an arbitrary precision.
 * The code below is slightly modified code from the link above:
 *
   load("diag");

   Dc(n) := diag(makelist(2*double_factorial(2*k-1),k,0,n));

   cmatrix_element[row,col]:=
      if is(col>row) then 0
      elseif row=1 and col=1 then 1/2
      else (-1)^(row+col)*4^(col-1)*(row-1)*(row+col-3)!/(row-col)!/(2*col-2)!;

   C(n) := genmatrix(cmatrix_element, n+1);

   f(g,n):=sqrt(2)*(%e/(2*(n+g)+1))^(n+1/2);

   Dr(k) := diag(append([1],makelist(-(2*n+2)!/(2*n!*(n+1)!),n,0,k-1)));

   bmatrix_element[row,col] :=
       if row = 1 then 1
       elseif is(row > col) then 0
       else (-1)^(col-row)*binomial(col+row-3,2*row-3);

   B(k) := genmatrix(bmatrix_element,k+1);

   lanczos_coeff(g, n) :=
       block([M : (Dr(n) . B(n)) . (C(n) . Dc(n)),
               f : transpose(matrix(makelist(f(g,k), k, 0, n)))],
           (M . f));

   tlg(g,n) := lanczos_coeff(g,n-1)*exp(g)/sqrt(2*%pi);
 *
 * When very high accuracy is not necessary, Boost recommends the following values for
 * g and N:
 *   g = 1.428456135094165802001953125
 *   N = 6
 *
 * More details at:
 * http://www.boost.org/doc/libs/1_57_0/libs/math/doc/html/math_toolkit/lanczos.html
 *
 * For the 25 digit precision, the Maxima code above returns the following values of
 * Lanczos coefficients:
 (%i1)  fpprec : 25$
 (%i2)  bfloat(tlg(1.428456135094165802001953125, 6));

          [  1.000000123880027543803602b0  ]
          [  9.727644985297789802780476b-1 ]
 (%o2)    [  9.79090078038154694104146b-3  ]
          [ -6.669128561042791409832734b-3 ]
          [  3.802685303772964895033826b-3 ]
          [ -1.117185821247582043028556b-3 ]
 *
 */


/* Lanczos G parameter: */
#define LANCZOS_G                ( 1.428456135094165802001953125L )

/* Number of Lanczos coefficients: */
#define NR_LANCZOS_COEF          ( 6 )


/*
 * Elements in the array LANCZOS_COEF_ARRAY are provided as high precision 
 * (25 digits)long double constants. The source file will tzpically
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
