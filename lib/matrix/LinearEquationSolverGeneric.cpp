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
 * Implementation of functions within namespace LinearEquationSolver.
 */


// no #include "LinearEquationSolverGeneric.hpp" !!!
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <complex>

#include "exception/MatrixException.hpp"
#include "util/NumericUtil.hpp"
#include "matrix/MatrixGeneric.hpp"
#include "../settings/omp_settings.h"



// A namespace with "private" functions:
namespace math {  namespace LinearEquationSolver {  namespace __private
{


/*
 * A "pseudo absolute value" of 'x'.
 * 
 * The function is used in situations when values must be ordered by
 * their absolute values and sometimes it is more convenient/more efficient to
 * obtain values that are a monotonically increasing function of
 * actual absolute values.
 * 
 * The general implementation returns the actual absolute value which
 * can be obtained efficiently for most scalar types.
 * 
 * @param x - value whose absolute value is returned
 * 
 * @return absolute value of 'x' 
 */
template <class T>
inline T pabs(const T& x)
{
    return std::abs(x);
}


/*
 * Partial "specialization" of 'pabs' for complex numbers.
 * 
 * This function returns a square of the actual absolute value and is
 * as such more efficient because no additional calculation of square root
 * (not a fast operation) is necessary.
 * 
 * @param x - a complex value
 * 
 * @return square of the absolute value of 'x', returned as a complex value with imag. part equal to 0
 */
template <class T>
inline std::complex<T> pabs(const std::complex<T>& x)
{
    return std::complex<T>( std::norm(x), static_cast<T>(0) );
}


/*
 * A convenience function that compares two (absolute) values.
 * 
 * @param a - first value
 * @param b - second value
 * 
 * @return true if a>b, false otherwise
 */
template <class T>
inline bool absgt(const T& a, const T& b)
{
    return (a > b);
}


/*
 * Partial "specialization" of 'absgt' for complex numbers.
 * The function compares real parts of 'a' and 'b'.
 * 
 * @param a - first complex value
 * @param b - second complex value
 * 
 * @return true if re(a)>re(b), false otherwise
 */
template <class T>
inline bool absgt(const std::complex<T>& a, const std::complex<T>& b)
{
    return ( std::real(a) > std::real(b) );
}


/*
 * Finds the row with the highest value (by absolute value) of the
 * p.th column.
 *
 * 'p' is returned immediately if it is out of a's range.
 *
 * @param a - a matrix of coefficients
 * @param p - index of the desired column
 *
 * @return row number with the highest absolute value of the p.th element
 */
template <class T>
size_t findPivot(
        const math::MatrixGeneric<T>& a,
        const size_t p )
{
    const size_t N = a.nrRows();

    if ( p >= N || p >= a.nrColumns() )
    {
        return p;
    }

    size_t r = p;
    T maxPiv = static_cast<T>(0);

    for ( size_t i=p; i<N; ++i )
    {
        const T elabs =
            math::LinearEquationSolver::__private::pabs( a(i, p) );

        // actually equivalent to:
        // if ( elabs > maxPiv )
        if ( true == math::LinearEquationSolver::__private::absgt(elabs, maxPiv ) )
        {
            maxPiv = elabs;
            r = i;
        }
    }

    return r;
}

}}}  // namespace math::LinearEquationsolver::__private



/**
 * Solves the system of linear equations using the Gauss - Jordan elimination
 * method and returns its unique solution if it exists.
 * 
 * Number of coef's columns must be equal to the number of term's rows.
 * 
 * If unique solution does not exist (i.e. determinant of 'coef' is 0), an
 * exception will be thrown.
 * 
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' are invalid or internal allocation of memory failed
 */
template <class T>
bool math::LinearEquationSolver::solveGaussJordan(
          const math::MatrixGeneric<T>& coef,
          const math::MatrixGeneric<T>& term,
          math::MatrixGeneric<T>& sol
        ) throw (math::MatrixException)
{
    // Partial pivoting is implemented

    // Sanity check
    if ( false == coef.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    /*
     * The Gaussian elimination algorithm is implemented:
     * multiples of coef's and term's lines are added to other lines until
     * "coef" appears as a unit matrix. In this case the modified "term" is a
     * unique solution of a system of linear equations. More details about
     * the algorithm at: http://en.wikipedia.org/wiki/Gaussian_elimination
     */
    const size_t N = coef.nrColumns();          // Nr. of unknowns
    const size_t NT = term.nrColumns();         // Nr. of terms' columns

    // Check the dimensions
    if ( N != term.nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    math::MatrixGeneric<T> temp(coef);
    sol = term;

    /*
     * Try to convert the 'temp' into an identity matrix
     * by appropriately adding multiples of other lines to each line
     * (incl. lines of 'retVal')
     */

    /*
     * The first part of the algorithm will implement partial pivoting
     * of rows. Among the remaining (N-i) rows it will find the one with
     * highest absolute value of the row's i.th element and swap the lines
     * (equivalent to rearranging the equations in a different order).
     * Additionally it will subtract multiples of the i.th row from all
     * subsequent rows (r>i) so their i.th column will be equal to 0. This
     * requires plenty of additions/subtractions of individual rows so
     * parallelization of this for loop is not possible due to race conditions.
     * It will be possible to parallelize certain parts of this loop, though.
     */
    for ( size_t i=0; i<N; ++i )
    {
        // Find the most appropriate row to pivot with this one
        const size_t pr = math::LinearEquationSolver::__private::findPivot(temp, i);

        // If even the highest absolute value equals 0, there is
        // no unique solution of the system of linear equations
        if ( true == math::NumericUtil::isZero<T>(temp(pr, i)) )
        {
            return false;
        }

        // Swap the rows of 'temp' and 'sol' if necessary
        if ( pr != i )
        {
            temp.swapRows(i, pr);
            sol.swapRows(i, pr);
        }

        // Set the i.th column of all other rows (r>i) to 0 by
        // adding the appropriate multiple of the i.th row
        #pragma omp parallel for default(none) shared(temp, sol, i)
        for ( size_t r=i+1; r<N; r++ )
        {
            // Nothing to do if temp(r,i) is already 0.
            if ( false == math::NumericUtil::isZero<T>(temp(r, i)) )
            {
                // Subtract a multiple of the i^th row.
                const T el = temp.get(r, i) / temp.get(i, i);

                // temp(r,:) = temp(r,:)-el*temp(i,:)
                for ( size_t c=i; c<N; ++c )
                {
                    temp(r, c) -= el*temp(i, c);
                }

                // term(r,:) = term(r,:)-el*term(i,:)
                for ( size_t c=0; c<NT; ++c )
                {
                    sol(r, c) -= el * sol(i, c);
                }
            }  // if (temp(r,i) != 0
        }  // for r
    }  // for i

    /*
     * Set the diag elements of 'temp' and 'retVal' to 1 by dividing the
     * whole row by temp(r,r). Columns smaller than 'r' are already equal to 0.
     */

    /*
     * Normalizing of each row is independent from other rows so it is
     * perfectly safe to parallelize the task by rows.
     */
    #pragma omp parallel for default(none) shared(temp, sol)
    for ( size_t r=0; r<N; ++r )
    {
        const T el = temp.get(r, r);

        for ( size_t c=r; c<N; ++c )
        {
            temp(r, c) /= el;
        }

        for ( size_t c=0; c<NT; ++c )
        {
            sol(r, c) /= el;
        }
    }

    /*
     * Now the lower triangle (below diag excl.) is 0, the diagonal consists of 1,
     * The upper triangle (above the diag) must be set to 0 as well.
     *
     * Column 'c' of each row (for r<c) will be set to 0
     * by adding the appropriate multiple of c.th row
     *
     * Parallelization of this for loop is not possible due to race conditions.
     */
    for ( size_t c=1; c<N; ++c )
    {
        // The current row to apply the operation described above

        /*
         * It is possible to parallelize this for loop because a row 'c' (not included
         * into the for loop) will be added independently to each "parallelized" row
         * 'r' (between 0 and c-1 incl.), thus no race condition is possible.
         */
        #pragma omp parallel for \
                    default(none) \
                    shared(temp, sol, c) \
                    schedule(dynamic)
        for ( size_t r=0; r<c; ++r )
        {
            // Nothing to do if temp(r,c) already equals 0
            if ( false == math::NumericUtil::isZero<T>(temp(r, c)) )
            {
                /*
                 * To set temp(r,c) to 0 it is a good idea to add the c.th row to it.
                 * temp(c,i); i<c are already 0 (i.e. will not affect anything left of temp(i,c)
                 * and temp(c,c) is already 1.
                 */

                const T el = temp.get(r, c);

                // temp(r,:) = temp(r,:) - el*temp(c,:)
                for ( size_t i=c; i<N; ++i )
                {
                    temp(r, i) -= el * temp(c, i);
                }

                // term(r,:) = term(r,:) - el*term(c,:)
                for ( size_t i=0; i<NT; ++i )
                {
                    sol(r, i) -= el * sol(c, i);
                }
            }  // if temp(r,c) != 0
        }  // for r
    }  // for c

    // A unique solution has been successfully found
    return true;

}
