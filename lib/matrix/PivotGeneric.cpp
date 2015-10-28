/*
Copyright 2015, Jernej Kovacic

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
 * Implementation of auxiliary functions within namespace Pivot that
 * perform partial or full pivoting.
 */


// no #include "PivotGeneric.hpp" !!!
#include <cstddef>
#include <algorithm>
#include <complex>
#include <cmath>
#include <vector>
#include <new>

#include "matrix/MatrixGeneric.hpp"
#include "util/NumericUtil.hpp"

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/MatrixException.hpp"


// A namespace with "private" functions:
namespace math {  namespace Pivot {  namespace __private
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
            col = std::min(p, Ncol);
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
                math::Pivot::__private::pabs( a(i, j) );

            // actually equivalent to:
            // if ( elabs > maxPiv )
            if ( true == math::Pivot::__private::absgt(elabs, maxPiv ) )
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


/*
 * Converts the matrix 'A' into an upper diagonal or identity matrix (depending
 * on 'fullm') by adding appropriate multiples of other rows to each rows. The
 * same applies for appropriate rows of the matrix 'pB' if provided. To increase
 * numerical stability of the algorithm, either partial or full pivoting is
 * performed (depending on 'fullp'). Additionally it returns A's rank and determinant
 * (both as side products of the algorithm) if 'pRank' and/or 'pDet' are provided.
 * 
 * @note The function modifies matrices 'A' and 'pB' (if provided) and
 *       variables 'pRank' and 'pDet' if provided.
 * 
 * @param A - matrix to be converted into an upper diagonal or identity matrix
 * @param fullp - logical value indicating whether full pivoting should be performed
 * @param fullm - logical value indicating whether 'A' should be converted into identity matrix
 * @param pB - pointer to a matrix of terms (may be NULL if not necessary)
 * @param pRank - pointer to a variable to assign A's rank (may be NULL if not required) 
 * @param pDet - pointer to a variable to assign A's determinant (may be NULL if not required) 
 * 
 * @throw MatrixException if any argument is invalid or allocation of a column bookkeeping vector fails
 */
template <class T>
void __pivot(
        math::MatrixGeneric<T>& A,
        const bool fullp,
        const bool fullm,
        math::MatrixGeneric<T>* const pB,
        size_t* const pRank,
        T* const pDet
      ) throw(math::MatrixException)
{
    // sanity check
    // if B is provided it must have the same nr. of rows as A
    if ( NULL!=pB && A.nrRows()!=pB->nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const size_t NR = A.nrRows();
    const size_t NC = A.nrColumns();
    const size_t N = std::min(NR, NC);  // the actual nr. of rows and columns to process
    const size_t NT = ( NULL!=pB ? pB->nrColumns() : 0 );

    /*
     * Bookkeeping of swapped columns is only necessary when 'pB' is provided
     * and full pivoting is requested. More details follow later.
     */
    const bool needColidx = (NULL != pB) && (true==fullp);


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
     * rearrangement of the final matrix/vector 'b' is necessary.
     * 
     * If any two columns of 'A' are swapped (i.e. full pivoting), this is
     * equivalent to rearrangement of terms in each equation, however the order
     * of equations remains unmodified. In this case rows of 'b' are not
     * rearranged immediately. However, all information about rearrangements of
     * A's columns must be noted. At the end of the algorithm (when 'b' "becomes"
     * the solution 'x'), b's rows must be rearranged into the original order
     * of A's columns.
     */

    /*
     * In addition to the notes above, the following properties of determinant
     * will be utilized:
     *
     * - addition of one row/column with a multiple of another row/column does
     *   not affect the determinant
     * - swapping of two rows or columns negates the determinant
     * - if one row or column is multiplied by a nonzero value, the determinant
     *   will be multiplied by the same value
     * - if a matrix is "converted" into an upper or lower triangle
     *   (all elements above/below the diagonal are 0), the determinant
     *   is a product of all diagonal elements
     */


    /*
     * Full pivoting requires that all swaps of columns are kept
     * in a separate vector. Initial values of the vector's elements
     * are equal to their indices. If any two columns of 'temp' are swapped,
     * the corresponding elements of 'colidx' will be swapped as well.
     */
    std::vector<size_t> colidx;

    if ( true == needColidx )
    {
        try
        {
            colidx.resize(NC);
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
        }

        /*
         * Assign 'colidx' with initial positions.
         */
        #pragma omp parallel num_threads(ompIdeal(NC)) \
                    if(NC>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(colidx)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(NC);

            typename std::vector<size_t>::iterator it = colidx.begin() + istart;
            size_t i = istart;
            for ( size_t cntr=0;
                  cntr<elems_per_thread && it!=colidx.end();
                  ++cntr, ++i, ++it )
            {
                *it = i;
            }

            (void) iend;
        }  // omp parallel

    }  // if fullp


    /*
     * The first part of the algorithm will perform pivoting.
     * Additionally it will subtract multiples of the i.th row from all
     * subsequent rows (r>i) so their i.th column will be equal to 0. This
     * requires plenty of additions/subtractions of individual rows so
     * parallelization of this for loop is not possible due to race conditions.
     * It will be possible to parallelize certain parts of this loop, though.
     */

    // Is the total number of row/column swaps odd or even?
    bool oddSwaps = false;
    if ( NULL != pDet )
    {
        // initial value of the determinant
        *pDet = static_cast<T>(1);
    }

    for ( size_t i=0; i<N; ++i )
    {
        // Find the most appropriate row and column to pivot with this one
        size_t pr = i;
        size_t pc = i;
        math::Pivot::__private::findPivot(A, i, pr, pc, fullp);

        // If even the highest absolute value equals 0, there is
        // no unique solution of the system of linear equations
        if ( true == math::NumericUtil::isZero<T>( A(pr, pc)) )
        {
            // in this case , the rank will be equal to the current 'i'...
            if ( NULL != pRank )
            {
                *pRank = i;
            }

            // ... and the determinant will be zero
            if ( NULL != pDet )
            {
                *pDet = static_cast<T>(0);
            }

            // nothing else to do in this function
            return;
        }

        // Swap the rows of 'A' and 'b' if necessary
        if ( pr != i )
        {
            // additionally update the "counter" of row/column swaps
            oddSwaps = !oddSwaps;

            A.swapRows_(i, pr);
            if ( NULL != pB )
            {
                pB->swapRows_(i, pr);
            }
        }

        // swap columns of 'A' if necessary
        if ( true==fullp && pc!=i )
        {
            oddSwaps = !oddSwaps;
            A.swapColumns_(i, pc);

            // and note the column swap in 'colidx'
            if ( true == needColidx )
            {
                std::swap(colidx.at(i), colidx.at(pc));
            }
        }

        /*
         * The entire row will be divided by A(i,i) later.
         * According to one of properties of the determinant (see notes above),
         * the current 'det' will be multiplied by A(i,i) now.
         */
        if ( NULL != pDet )
        {
            *pDet *= A(i, i);
        }

        // Set the i.th column of all other rows (r>i) to 0 by
        // adding the appropriate multiple of the i.th row
        #pragma omp parallel for default(none) shared(A, pB, i)
        for ( size_t r=i+1; r<NR; r++ )
        {
            // Nothing to do if temp(r,i) is already 0.
            if ( false == math::NumericUtil::isZero<T>( A(r, i)) )
            {
                // Subtract a multiple of the i^th row.
                const T el = A(r, i) / A(i, i);

                // A(r,:) = A(r,:) - el*A(i,:)
                for ( size_t c=i; c<NC; ++c )
                {
                    A(r, c) -= el * A(i, c);
                }

                // b(r,:) = b(r,:) - el*b(i,:)
                if ( NULL != pB )
                {
                    math::MatrixGeneric<T>& B = *pB;

                    for ( size_t c=0; c<NT; ++c )
                    {
                        B(r, c) -= el * B(i, c);
                    }
                }
            }  // if A(r,i) != 0
        }  // for r
    }  // for i

    /*
     * Set the diag elements of 'A' and 'b' to 1 by dividing the
     * whole row by temp(r,r). Columns smaller than 'r' are already equal to 0.
     */

    /*
     * Normalizing of each row is independent from other rows so it is
     * perfectly safe to parallelize the task by rows.
     */
    #pragma omp parallel for default(none) shared(A, pB)
    for ( size_t r=0; r<NR; ++r )
    {
        const T el = A(r, r);

        for ( size_t c=r; c<NC; ++c )
        {
            A(r, c) /= el;
        }

        if ( NULL != pB )
        {
            math::MatrixGeneric<T>& B = * pB;

            for ( size_t c=0; c<NT; ++c )
            {
                B(r, c) /= el;
            }
        }
    }

    /*
     * 'A' is now an upper diagonal matrix with ones on the diagonal.
     * No additional processing of 'A' and 'b' is necessary if 'fullm' is false.
     */

    // If this point is reached, A is a full rank matrix (rank=min(NR, NC)):
    if ( NULL != pRank )
    {
        *pRank = N;
    }

    /*
     * Calculation of the determinant is almost complete.
     * Additionally it must be negated if the total number of
     * swaps of rows and columns is odd. For more details,
     * see properties of the determinant in the notes above.
     */
    if ( NULL!=pDet && true==oddSwaps)
    {
        *pDet = -(*pDet);
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
    if ( true == fullm )
    {
        for ( size_t c=1; c<NC; ++c )
        {
            // The current row to apply the operation described above

            /*
             * It is possible to parallelize this for loop because a row 'c' (not included
             * into the for loop) will be added independently to each "parallelized" row
             * 'r' (between 0 and c-1 incl.), thus no race condition is possible.
             */
            #pragma omp parallel for \
                        default(none) \
                        shared(A, pB, c) \
                        schedule(dynamic)
            for ( size_t r=0; r<c; ++r )
            {
                // Nothing to do if A(r,c) already equals 0
                if ( false == math::NumericUtil::isZero<T>( A(r, c)) )
                {
                    /*
                     * To set A(r,c) to 0 it is a good idea to add the c.th row to it.
                     * A(c,i); i<c are already 0 (i.e. will not affect anything left of A(i,c)
                     * and A(c,c) is already 1.
                     */

                    const T el = A(r, c);

                    // A(r,:) = A(r,:) - el*A(c,:)
                    for ( size_t i=c; i<NC; ++i )
                    {
                        A(r, i) -= el * A(c, i);
                    }

                    // A(r,:) = A(r,:) - el*A(c,:)
                    if ( NULL != pB )
                    {
                        math::MatrixGeneric<T>& B = *pB;

                        for ( size_t i=0; i<NT; ++i )
                        {
                            B(r, i) -= el * B(c, i);
                        }
                    }
                }  // if A(r,c) != 0
            }  // for r
        }  // for c
    }  // if !fullm

    /*
     * If full pivoting was performed, it is likely that some columns
     * of 'temp' have been rearranged. In this case, b's rows must be
     * rearranged back to the original columns of A's columns.
     * 
     * The algorithm below will check if colidx(i) equals 'i'.
     * If it does not, it will find the index of colidx's element that does equal
     * 'i' and swap the rows of 'sol' accordingly. Additionally this swap will
     * be reflected by swapping of corresponding elements of 'colidx'.
     */

    if ( true==needColidx && NR==NC )
    {
        math::MatrixGeneric<T>& B = *pB;

        for ( size_t i=0; i<NR; ++i )
        {
            // find such 'j' that colidx(j) == i
            size_t j = i;
            for ( j=i; colidx.at(j)!=i; ++j );

            // and swap sol's rows and colidx's elements
            B.swapRows_(i, j);
            std::swap( colidx.at(i), colidx.at(j) );
        }  // for i
    }  // if needColidx
}

}}}  // namespace math::Pivot::__private



/*
 * Solves the system of linear equations using the Gauss - Jordan elimination
 * method and returns its unique solution if it exists.
 * 
 * 'A' must be a square matrix and 'b' must have the same number of
 * rows as 'b'.
 * 
 * If unique solution does not exist (i.e. determinant of 'A' is 0), an
 * exception will be thrown.
 * 
 * The method performs either partial or full pivoting. While partial pivoting
 * should be sufficient in most cases, full pivoting is usually numerically
 * more stable, however it introduces additional overhead.
 * 
 * @param A - a square matrix with coefficients of the system of linear equations
 * @param b - a matrix with constant terms of the system of linear equations
 * @param x - a reference to a matrix to be assigned the solution of equations
 * @param fullp - should the algorithm perform full pivoting?
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'A' and 'b' are invalid or internal allocation of memory failed
 */
template <class T>
bool math::Pivot::solveGaussJordan(
        const math::MatrixGeneric<T>& A,
        const math::MatrixGeneric<T>& b,
        math::MatrixGeneric<T>& x,
        const bool fullp
      ) throw(math::MatrixException)
{
    // sanity check
    if ( false == A.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const size_t N = A.nrColumns();          // Nr. of unknowns

    // Check the dimensions
    if ( N != b.nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    math::MatrixGeneric<T> temp(A);
    x = b;
    
    size_t rank;

    math::Pivot::__private::__pivot<T>(temp, fullp, true, &x, &rank, NULL);

    // if rank(A) equals N, the system of lin. equations has a unique solution

    return ( N == rank );
}


/**
 * Calculates determinant of a square matrix.
 * 
 * @param A - a matrix to calculate its determinant
 * @param fullp - should the algorithm perform full pivoting?
 * 
 * @return determinant of 'A'
 * 
 * @throw MatrixException if 'A' is not a square matrix or if internal allocation of memory failed
 */
template <class T>
T math::Pivot::getDeterminant(
        const math::MatrixGeneric<T>& A,
        const bool fullp
      ) throw(math::MatrixException)
{
    // sanity check
    if ( false == A.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
    }

    math::MatrixGeneric<T> temp(A);
    T det;

    math::Pivot::__private::__pivot<T>(temp, fullp, false, NULL, NULL, &det);

    return det;
}
