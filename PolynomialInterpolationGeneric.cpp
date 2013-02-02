/*
Copyright 2013, Jernej Kovacic

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
 * @file PolynomialInterpolationGeneric.cpp
 *
 * Implementation of the class PolynomialInterpolationGeneric. The class
 * calculates an interpolation polynomial that goes exactly through
 * entered points.
 * 
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 *
 * @author Jernej Kovacic
 */

// Deliberately there is no #include "PolynomialInterpolationGeneric.h"
#include "MatrixGeneric.h"
#include "PolynomialGeneric.h"


/**
 * Constructor.
 * 
 * It initializes all its internal structures.
 */
template<class T>
math::PolynomialInterpolationGeneric<T>::PolynomialInterpolationGeneric()
{
    math::CurveFittingGenericAb<T>::init();
}

/**
 * Generates an interpolation polynomial that goes exactly 
 * through all entered points. The degree of the polynomial is
 * equal to nr. of points - 1.
 * 
 * @param degree - ignored as it does not make any sense at this class.
 * 
 * @throw CurveFittingException if generation of the curve failed for any reason.
 */
template<class T>
void math::PolynomialInterpolationGeneric<T>::generateCurve(unsigned int degree) throw (math::CurveFittingException)
{
    // performs necessary checks
    this->curveGenerationCheck();
    // TODO check value of degree?

    /*
         One possible algorithm would be construction of the so called
         Vandermonde matrix. However, this method is known as computationally
         inefficient and prone to significant computation errors when solving a
         system of linear equations. 
     
         Instead, the Newton's interpolation algorithm will be applied, described at: 
         http://en.wikipedia.org/wiki/Newton_polynomial
     
         Procedure of the algorithm:
         - input: points: (x0,y0), (x1,y1), ..., (xn,yn) 
         - create a (n+1)x(n+1) matrix a, fill it with zeros
         - fill the first column of the matrix with y's:
             a(i,0) = y(i)   for i = 0 .. n
         - calculate lower diagonal elements as differential quotients:
     
                       a(r,c-1) - a(r-1,c-1)     
             a(r,c) = ------------------------
                           y(r) - y(r-c)
       
              for c = 1..n, r=c..n
      
         - when this procedure is completed, the interpolation polynomial
           is expressed as:
     
             p(x) = y(0) + a(1,1)*(x-y(0)) + a(2,2)*(x-y(0))*(x-y(1)) + ... +
       
                          n-1
                         +---+     
                         |   |
              + a(n,n) * |   | (x-y(k))      
                         |   |
                          k=0
     */

    try
    {
        // number of points
        const unsigned int N = this->points.size();

        // create matrices as described above:
        math::MatrixGeneric<T> a(N, N);
        math::MatrixGeneric<T> x(N, 1);
        
        // hiding idx from the rest of the function
        {
            // As iterators are the fastest way to access linked list elements,
            // traverse the list only once and populate appropriate elements of a and b
            unsigned int idx = 0;
            for ( 
              typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator it=this->points.begin();
                        it!=this->points.end(); it++, idx++ )
                {
                    x.set(idx, 0, it->p_x);
                    a.set(idx, 0, it->p_y);
                }  // for it
        }  // hiding idx from the rest of the code
        
        // Populate lower diagonal elements of matrix a as described above
        for ( unsigned int c=1; c<N; c++ )
        {
            for ( unsigned int r=c; r<N; r++)
            {
                a.at(r, c) = (a.at(r, c-1)-a.at(r-1, c-1)) / (x.at(r, 0)-x.at(r-c, 0));
            }
        }

        // and create the polynomial as described above:
        
        // (x-y0)*(x-y1)*...*(x-yi)
        math::PolynomialGeneric<T> temp(1);
        // partial polynomial sum (sum of temp*a(i,i))
        math::PolynomialGeneric<T> sum(1);
        // 1st degree term (x-y) that will be multiplied by temp
        math::PolynomialGeneric<T> term(2);

        // Initialize the polynomials:
        
        // term will always be of form (x-yi), hence term(1)
        // is always equal to 1, while term(0) will be set depending on i
        term.set(1, math::NumericUtil<T>::ONE);

        // initial value of temp: 1 
        temp.set(0, math::NumericUtil<T>::ONE);
        // initial value of sum: a(0,0)
        sum.set(0, a.get(0, 0));
        
        // finally implement the last section of the description above
        for (unsigned int i=1; i<N; i++ )
        {
            term.set(0, -x.at(i-1, 0));
            temp *= term;
            sum += a.at(i, i) * temp;
        }  // for i

        // The algorithm is finished, assign the interpolation polynomial to poly:
        this->poly = sum;
        
        // the curve can be marked as generated
        this->curveGenerated = true;
    }  // try
    catch ( const math::MatrixException& mex )
    {
        // the only possible MatrixException is out of memory
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }
    catch ( const math::PolynomialException& pex )
    {
        // the only possible PolynomialException is out of memory
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }
}
