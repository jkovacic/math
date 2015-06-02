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
#include <vector>
#include <new>

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
 * Finds row and (optionally) column indices of the matrix' pivot element
 * (i.e. the one with the highest absolute value) past the p.th row and (depending
 * on 'fullp').
 *
 * 'p' is returned immediately if it is out of a's range.
 *
 * @param a - matrix of coefficients
 * @param p - index of the desired row/column to find a pivot for
 * @param row - reference to a variable to write the pivot's row number to
 * @param col - reference to a variable to write the pivot's column number to (not modified if 'fullp' is false)
 * @param fullp - if true, also consider elements at columns greater than 'p'
 */
template <class T>
void findPivot(
        const math::MatrixGeneric<T>& a,
        const size_t p,
        size_t& row,
        size_t& col,
        const bool fullp )
{
    const size_t N = a.nrRows();
    const size_t Ncol = a.nrColumns();

    // dimension check
    if ( p >= N || p >= Ncol )
    {
        row = p;

        if ( true == fullp )
        {
            col = p;
        }

        return;
    }

    // initial "pivot"
    size_t r = p;
    size_t c = p;
    T maxPiv = static_cast<T>(0);
    
    // the highest column number to be considered:
    const size_t CMAX = ( true==fullp ? Ncol : p+1 );

    for ( size_t i=p; i<N; ++i )
    {
        for (size_t j=p; j<CMAX; ++j )
        {
            const T elabs =
                math::LinearEquationSolver::__private::pabs( a(i, j) );

            // actually equivalent to:
            // if ( elabs > maxPiv )
            if ( true == math::LinearEquationSolver::__private::absgt(elabs, maxPiv ) )
            {
                maxPiv = elabs;
                r = i;
                c = j;
            }
        }  // for j
    }  // for i

    // finally assign 'row' and optionally 'col'
    row = r;

    if ( true == fullp )
    {
        col = c;
    }
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
 * The method performs either partial or full pivoting. While partial pivoting
 * should be sufficient in most cases, full pivoting is usually numerically
 * more stable, however it introduces additional overhead.
 * 
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * @param fullp - should the method perform full pivoting (default: TRUE)
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' are invalid or internal allocation of memory failed
 */
template <class T>
bool math::LinearEquationSolver::solveGaussJordan(
          const math::MatrixGeneric<T>& coef,
          const math::MatrixGeneric<T>& term,
          math::MatrixGeneric<T>& sol,
          const bool fullp
        ) throw (math::MatrixException)
{
    // Sanity check
    if ( false == coef.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    /*
     * The Gaussian elimination algorithm is implemented:
     * multiples of coef's and term's lines are added to other lines until
     * 'coef' appears as a unit matrix. In this case the modified 'term' is a
     * unique solution of a system of linear equations. More details about
     * the algorithm at: https://en.wikipedia.org/wiki/Gaussian_elimination
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
     * Summary of the entire algorithm:
     * 
     * Try to convert the 'temp' into an identity matrix
     * by appropriately adding multiples of other lines to each line
     * (incl. lines of 'sol'). Apply pivoting for better numerical
     * stability. If 'temp' is successfully converted to an identity
     * matrix, 'sol' will become a unique solution of the system
     * of linear equations.
     */

    /*
     * Notes about partial and full pivoting:
     * 
     * A set of linear equations A*x=b can be rewritten as:
     * 
     *          +       +        +       +               +       +    +     +
     *          | a_1,1 |        | a_1,2 |               | a_1,N |    | b_1 |
     *          | a_2,1 |        | a_2,2 |               | a_2,N |    | b_2 |
     *     x_1  |  ...  | + x_2  |  ...  | +  ... + x_N  |  ...  | =  | ... |
     *          | a_N,1 |        | a_N,2 |               | a_N,N |    | b_N |
     *          +       +        +       +               +       +    +     +
     * 
     * If any two rows of 'A' are swapped (i.e. partial pivoting), this is
     * equivalent to rearranging the equations in a different order. It is known
     * that solution of a system of linear equations is independent of the order
     * of equations. If any rearrangement of rows is performed in 'A', the same
     * rearrangement of b's rows must be performed immediately and no
     * rearrangement of the final matrix/vector 'sol' is necessary.
     * 
     * If any two columns of 'A' are swapped (i.e. full pivoting), this is
     * equivalent to rearrangement of terms in each equation, however the order
     * of equations remains unmodified. In this case rows of 'b' are not
     * rearranged immediately. However, all information about rearrangements of
     * A's columns must be noted. At the end of the algorithm (when 'b' "becomes"
     * the solution 'x'), sol's rows must be rearranged into the original order
     * of A's columns.
     */

    /*
     * Full pivoting requires that all swaps of columns are kept
     * in a separate vector. Initial values of the vector's elements
     * are equal to their indices. If any two columns of 'temp' are swapped,
     * the corresponding elements of 'colidx' will be swapped as well.
     */
    std::vector<size_t> colidx;

    if ( true == fullp )
    {
        try
        {
            colidx.resize(N);
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
        }

        for ( size_t i=0; i<N; ++i )
        {
            colidx.at(i) = i;
        }
    }  // if fullp


    /*
     * The first part of the algorithm will perform pivoting.
     * Additionally it will subtract multiples of the i.th row from all
     * subsequent rows (r>i) so their i.th column will be equal to 0. This
     * requires plenty of additions/subtractions of individual rows so
     * parallelization of this for loop is not possible due to race conditions.
     * It will be possible to parallelize certain parts of this loop, though.
     */

    for ( size_t i=0; i<N; ++i )
    {
        // Find the most appropriate row and column to pivot with this one
        size_t pr = i;
        size_t pc = i;
        math::LinearEquationSolver::__private::findPivot(temp, i, pr, pc, fullp);

        // If even the highest absolute value equals 0, there is
        // no unique solution of the system of linear equations
        if ( true == math::NumericUtil::isZero<T>(temp(pr, pc)) )
        {
            return false;
        }

        // Swap the rows of 'temp' and 'sol' if necessary
        if ( pr != i )
        {
            temp.swapRows(i, pr);
            sol.swapRows(i, pr);
        }

        // swap columns of 'temp' if necessary
        if ( true==fullp && pc!=i )
        {
            temp.swapColumns(i, pc);
            std::swap(colidx.at(i), colidx.at(pc));
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
                const T el = temp(r, i) / temp(i, i);

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
        const T el = temp(r, r);

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

    /*
     * If full pivoting was performed, it is likely that some columns
     * of 'temp' have been rearranged. In this case, sol's rows must be
     * rearranged back to the original columns of temp's columns.
     * 
     * The algorithm below will check if colidx(i) equqls 'i'.
     * If it does not, it will find the index of colidx's element that does equal
     * 'i' and swap the rows of 'sol' accordingly. Additionally this swap will
     * be reflected by swapping of corresponding elements of 'colidx'.
     */
    if ( true == fullp )
    {
        for ( size_t i=0; i<N; ++i )
        {
            // find such 'j' that colidx(j) == i
            size_t j = i;
            for ( j=i; colidx.at(j)!=i; ++j );

            // and swap sol's rows and colidx's elements
            sol.swapRows(i, j);
            std::swap( colidx.at(i), colidx.at(j) );
        }  // for i
    }  // if fullp

    // A unique solution has been successfully found
    return true;

}
