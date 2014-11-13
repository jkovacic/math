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
 * Implementation of the class PolynomialRegressionGeneric. The class
 * calculates a regression polynomial using the method of least squares.
 * 
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */

// Deliberately there is no #include "PolynomialRegressionGeneric.hpp" !
#include "matrix/MatrixGeneric.hpp"
#include "matrix/SqMatrixGeneric.hpp"
#include "lineq/LinearEquationSolverGeneric.hpp"
#include "util/NumericUtil.hpp"
#include "exception/CurveFittingException.hpp"

#include <cstddef>
#include <limits>


/**
 * Constructor.
 * 
 * It initializes all its internal structures.
 */
template<class T>
math::PolynomialRegressionGeneric<T>::PolynomialRegressionGeneric()
{
    // all necessary functionality is implemented by the base class's init()
    math::CurveFittingGenericAb<T>::_init();
}

/**
 * Generates a regression polynomial satisfying the least squares criteria.
 * 
 * @param degree - desired degree of the regression polynomial
 * 
 * @throw CurveFittingException if generation of the curve failed for any reason.
 */
template<class T>
void math::PolynomialRegressionGeneric<T>::generateCurve(size_t degree) throw (math::CurveFittingException)
{
    // performs necessary checks
    this->_curveGenerationCheck();
    // TODO check value of degree?

    /*
        Coefficients of the regression polynomial can be calculated as
        a solution of a system of linear equations:  A*x=b

        A is a square NxN matrix, its coefficients are calculated as follows:

                      n-1      
                    -------     
                    \        (r+c)
           A(r,c) =  |   x(i)
                    /
                    -------
                      i=0

        and

                      n-1
                    -------
                    \            r
           B(r)  =   |  y(i)*x(i)
                    /
                    -------
                      i=0

        where n is the number of points. 
     */
    
    try
    {
        if ( degree>=(std::numeric_limits<size_t>::max()/2-1) )
        {
            throw math::CurveFittingException(math::CurveFittingException::CURVE_GENERATION_FAILED);    
        }
        
        // number of polynomial's coefficients
        const size_t N = degree + 1;

        math::SqMatrixGeneric<T> a(N);
        math::MatrixGeneric<T> b(N, 1);
        
        /*
             Instead of calculating all powers of x, all points will be
             traversed once and the appropriate terms will be added to 
             matrices' elements
         */
        for ( typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator it=this->points.begin(); it!=this->points.end(); ++it )
        {
            // initial values of summands
            T aterm = math::NumericUtil<T>::ONE;
            T bterm = it->p_y;

            // i actually determines a position inside both matrices: i = r+c             
            for ( size_t i=0; i<(2*N-1); ++i )
            {
                // find all possible rows satisfying the condition above
                const size_t Rmax = ( i<=degree ? i : degree );
                for ( size_t r=0; r<=Rmax; ++r )
                {
                    // do not update anything if any element 
                    // is out of matrix's range
                    if ( r>i || i>=(N+r) )
                    {
                        continue;  // for r
                    }
                    a.at(r, i-r) += aterm;
                }  // for r
                // and update 'aterm' for the next iteration of i
                aterm *=it->p_x;

                // do not do anything if i is out of b's range 
                if ( i>degree )
                {
                    continue;  // for i
                }

                // update the b's i^th element 
                b.at(i,0) += bterm;
                // and update 'bterm' for the next iteration of i
                bterm *= it->p_x;
            }  // for i
        }  // for it

        // Matrices are filled, solve the system of linear equations
        math::LinearEquationSolverGeneric<T> leq(a, b);
        math::MatrixGeneric<T> x;
        leq.solve(x);

        // And finally fill the regression polynomial
        for ( size_t i=0; i<N; ++i )
        {
            this->poly.set(i, x.at(i, 0));
        }

        // the curve can be marked as generated
        this->curveGenerated = true;
    }  // try
    catch ( const math::MatrixException& mex )
    {
        // the only possible MatrixException is out of memory
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }
    catch ( const math::LinearEquationSolverException& leqex )
    {
        throw math::CurveFittingException(math::CurveFittingException::CURVE_GENERATION_FAILED);
    }
    catch ( const math::PolynomialException& pex )
    {
        // the only possible PolynomialException is out of memory
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }

}
