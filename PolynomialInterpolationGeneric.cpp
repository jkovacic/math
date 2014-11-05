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
 * @file
 * @author Jernej Kovacic
 *
 * Implementation of the class PolynomialInterpolationGeneric. The class
 * calculates an interpolation polynomial that goes exactly through
 * entered points.
 * 
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */

// Deliberately there is no #include "PolynomialInterpolationGeneric.h"
#include "PolynomialGeneric.h"
#include "CurveFittingException.h"
#include "omp_settings.h"

#include <vector>
#include <new>
#include <cstddef>


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
void math::PolynomialInterpolationGeneric<T>::generateCurve(size_t degree) throw (math::CurveFittingException)
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
                           x(r) - x(r-c)
       
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
        // As it turns out, the algorithm above can be further optimized.
        // It is sufficient to use a 1D vector (storing just the current "column") 
        // instead of a 2D matrix. It is also possible to update polynomials 
        // concurrently when updating "columns". 
        
        // number of points
        const size_t N = this->points.size();

        // create vectors to store temporary results:
        std::vector<T> a;
        std::vector<T> x;
        
        if ( N > a.max_size() )
        {
            throw math::CurveFittingException(math::CurveFittingException::CURVE_GENERATION_FAILED);
        }
        
        a.resize(N);
        x.resize(N);
        
        // hiding idx from the rest of the function
        {
            // As iterators are the fastest way to access linked list elements,
            // traverse the list only once and populate appropriate elements of a and b
            size_t idx = 0;
            for ( 
              typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator it=this->points.begin();
                        it!=this->points.end(); ++it, ++idx )
                {
                    x.at(idx) = it->p_x;
                    a.at(idx) = it->p_y;
                }  // for it
        }  // hiding idx from the rest of the code
        
        // Polynomials:
        // (x-y0)*(x-y1)*...*(x-yi)
        math::PolynomialGeneric<T> temp(true, 1);
        // partial polynomial sum (sum of temp*a(i,i))
        math::PolynomialGeneric<T> sum(true, 1);
        // 1st degree term (x-y) that will be multiplied by temp
        math::PolynomialGeneric<T> term(true, 2);

        // Initialize the polynomials:

        // term will always be of form (x-xi), hence term(1)
        // is always equal to 1, while term(0) will be set depending on c
        term.set(1, math::NumericUtil<T>::ONE);

        // initial value of temp: 1 
        temp.set(0, math::NumericUtil<T>::ONE);
        // initial value of sum: a(0,0)
        sum.set(0, a.at(0));
        
        // recalculate vector's element as differential quotients:
        // a(i) = (a(i+1)-a(i))/ (appropriate difference of x)

        /*
         * This for loop cannot be parallelized because 'a' at each iteration
         * depends on the same vector at the previous iteration
         */
        for ( size_t c=0; c<(N-1); ++c )
        {
            /*
             * On the other hand, it is possible to parallelize the inner
             * for loop. However it must be ensured that a.at(i) is updated
             * before a.at(i+1). This can be achieved using the OpenMP
             * "ordered" clause.
             */
            #pragma omp parallel for ordered default(none) shared(a, x, c)
            for ( size_t i=0; i<(N-1-c); ++i )
            {
                T el = (a.at(i+1) - a.at(i)) / (x.at(i+c+1) - x.at(i));

                // It is important that this update is performed in the right order!
                #pragma omp ordered
                a.at(i) = el;
            }  // for i
  
            
            // finally update the polynomials as described in the last section above
            term.set(0, -x.at(c));
            temp *= term;
            sum += a.at(0) * temp;
        } // for c
        
        // vectors are not needed anymore
        x.clear();
        a.clear();
                
        // The algorithm is finished, assign the interpolation polynomial to poly:
        this->poly = sum;
        
        // the curve can be marked as generated
        this->curveGenerated = true;
    }  // try
    catch ( const math::PolynomialException& pex )
    {
        // the only possible PolynomialException is out of memory
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }

    (void) degree;
}
