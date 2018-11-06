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
 */

// no #include "PolynomialRegressionGeneric.hpp" !!!
#include "matrix/MatrixGeneric.hpp"
#include "matrix/LinearEquationSolverGeneric.hpp"
#include "exception/MatrixException.hpp"
#include "exception/CurveFittingException.hpp"

#include <cstddef>
#include <limits>


/**
 * Constructor.
 * 
 * It initializes all its internal structures.
 */
template <typename F>
math::PolynomialRegressionGeneric<F>::PolynomialRegressionGeneric()
{
    // all necessary functionality is implemented by the base class's init()
    math::CurveFittingGenericAb<F>::_init();
}


/**
 * Generates a regression polynomial satisfying the least squares criteria.
 * 
 * @param degree - desired degree of the regression polynomial
 * 
 * @throw CurveFittingException if generation of the curve failed for any reason.
 */
template <typename F>
void math::PolynomialRegressionGeneric<F>::generateCurve(const std::size_t degree)
{
    // performs necessary checks
    this->_curveGenerationCheck();
    // TODO check value of degree?

    /*
     * Coefficients of the regression polynomial can be calculated as
     * a solution of a system of linear equations:  A*x=b
     *
     * A is a square NxN matrix, its coefficients are calculated as follows:
     *
     *               n-1      
     *              -----     
     *              \        (r+c)
     *     A(r,c) =  >   x(i)
     *              /
     *              -----
     *               i=0
     *
     * and
     *
     *               n-1
     *              -----
     *              \            r
     *     B(r)  =   >  y(i)*x(i)
     *              /
     *              -----
     *               i=0
     *
     * where 'n' denotes the number of points. 
     */
    
    try
    {
        if ( degree>=(std::numeric_limits<std::size_t>::max()/2-1) )
        {
            throw math::CurveFittingException(math::CurveFittingException::CURVE_GENERATION_FAILED);    
        }
        
        // number of polynomial's coefficients
        const std::size_t N = degree + 1;

        math::MatrixGeneric<F> a(N);
        math::MatrixGeneric<F> b(N, 1);
        
        /*
         * Instead of calculating all powers of x, all points will be
         * traversed once and the appropriate terms will be added to 
         * matrices' elements
         */
        for ( 
            typename std::list<typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator it=this->m_points.begin(); 
            it!=this->m_points.end(); 
            ++it )
        {
            // initial values of summands
            F aterm = static_cast<F>(1);
            F bterm = it->m_y;

            // i actually determines a position inside both matrices: i = r+c             
            for ( std::size_t i=0; i<(2*N-1); ++i )
            {
                // find all possible rows satisfying the condition above
                const std::size_t Rmax = ( i<=degree ? i : degree );
                for ( std::size_t r=0; r<=Rmax; ++r )
                {
                    // do not update anything if any element 
                    // is out of matrix's range
                    if ( r>i || i>=(N+r) )
                    {
                        continue;  // for r
                    }
                    a(r, i-r) += aterm;
                }  // for r
                // and update 'aterm' for the next iteration of i
                aterm *=it->m_x;

                // do not do anything if i is out of b's range 
                if ( i>degree )
                {
                    continue;  // for i
                }

                // update the b's i^th element 
                b(i, 0) += bterm;
                // and update 'bterm' for the next iteration of i
                bterm *= it->m_x;
            }  // for i
        }  // for it

        /*
         * Matrices are filled, solve the system of linear equations.
         * Note that 'degree' will usually be rather small, hence
         * the Gauss - Jordan method will be sufficient in majority of cases.
         */
        math::MatrixGeneric<F> x(b);
        const bool succ = math::LinearEquationSolver::solveGaussJordan<F>(a, b, x);

        // System of linear equations successfully solved?
        if ( false == succ )
        {
            // Obviously not
            throw math::CurveFittingException(math::CurveFittingException::CURVE_GENERATION_FAILED);
        }

        // And finally fill the regression polynomial
        for ( std::size_t i=0; i<N; ++i )
        {
            this->m_poly.set(i, x(i, 0));
        }

        // the curve can be marked as generated
        this->m_curveGenerated = true;
    }  // try
    catch ( const math::MatrixException& mex )
    {
        if ( math::MatrixException::OUT_OF_MEMORY == mex.error )
        {
            throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
        }
        else
        {
            throw math::CurveFittingException(math::CurveFittingException::CURVE_GENERATION_FAILED);
        }
    }
    catch ( const math::PolynomialException& pex )
    {
        // the only possible PolynomialException is out of memory
        throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
    }

}
