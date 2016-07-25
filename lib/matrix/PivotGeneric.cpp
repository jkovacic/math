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
#include <vector>
#include <new>

#include "matrix/MatrixGeneric.hpp"
#include "util/NumericUtil.hpp"
#include "util/PseudoFunctionGeneric.hpp"

#include "util/FillVectors.h"

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/MatrixException.hpp"


// A namespace with "private" functions:
namespace math {  namespace Pivot {  namespace __private
{

/*
 * Returns either 'i' or the i.th element of 'v' if 'v' is provided,
 * i.e. it is not NULL.
 * 
 * @note It is assumed that 'i' will alway be less than size of 'v' (if provided)
 * 
 * @param v - pointer to a vector of permuted indices (must be NULL if not relevant)
 * @param i - position of the desired element in 'v'
 * 
 * @return 'i' or v[i] if 'v' is provided
 */
inline std::size_t __index(const std::vector<std::size_t>* const v, const std::size_t i)
{
    return ( NULL==v ? i : v->at(i) );
}


/*
 * Finds row and (optionally) column indices of the matrix' pivot element
 * (i.e. the one with the highest absolute value) past the p.th row and (depending
 * on 'fullp').
 * 
 * @note If 'prows' and/or 'pcols' is provided, returned 'row' and/or 'col'
 *       will satisfy: abs( a(prows(row, pcols(col) ) = max
 *       given that row>=p and col>=p 
 *
 * @note 'p' is returned immediately if it is out of a's range.
 *
 * @param a - matrix of coefficients
 * @param p - index of the desired row/column to find a pivot for
 * @param prows - pointer to the vector of row indices (may be NULL)
 * @param pcols - pointer to the vector of column indices (may be NULL) 
 * @param row - reference to a variable to write the pivot's row number to
 * @param col - reference to a variable to write the pivot's column number to (not modified if 'fullp' is false)
 * @param fullp - if true, also consider elements at columns greater than 'p'
 */
template <class T>
void findPivot(
        const math::MatrixGeneric<T>& a,
        const std::size_t p,
        const std::vector<std::size_t>* const prows,
        const std::vector<std::size_t>* const pcols,
        std::size_t& row,
        std::size_t& col,
        const bool fullp )
{
    const std::size_t Nrow = a.nrRows();
    const std::size_t Ncol = a.nrColumns();

    // check of dimensions
    if ( p >= Nrow || p >= Ncol )
    {
        row = std::min(p, Nrow);

        if ( true == fullp )
        {
            col = std::min(p, Ncol);
        }

        return;
    }

    // the number of rows and columns to be considered:
    const std::size_t ROWS = Nrow - p;
    const std::size_t COLS = ( true==fullp ? Ncol-p: 1 );

    // initial absolute value of "pivot"
    T globMax = math::PseudoFunction::pabs( 
        a( math::Pivot::__private::__index(prows, p),
           math::Pivot::__private::__index(pcols, p)) );

    // as we are searching within a subset of a matrix,
    // N should never be out of size_t's range
    const std::size_t N = ROWS * COLS;
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(a, globMax, row, col)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        // indices of the highest absolute value within the assigned block
        std::size_t r = p + istart / COLS;
        std::size_t c = p + istart % COLS;
        T localMax = math::PseudoFunction::pabs( 
            a( math::Pivot::__private::__index(prows, r),
               math::Pivot::__private::__index(pcols, c) ) );

        for ( std::size_t it=istart; it<iend; ++it )
        {
            // current indices
            const std::size_t i = p + it / COLS;
            const std::size_t j = p + it % COLS;

            const T elabs = math::PseudoFunction::pabs( 
              a( math::Pivot::__private::__index(prows, i),
                 math::Pivot::__private::__index(pcols, j) ) );

            if ( true == math::PseudoFunction::absgt(elabs, localMax) )
            {
                localMax = elabs;
                r = i;
                c = j;
            }
        }  // for it

        // prevent possible race condition when updating 'row' and 'col'
        #pragma omp critical(pivotgeneric_findpivot)
        {
            // Check if the local (i.e. within the assigned block ) highest
            // absolute value is greater than the global one
            if ( true == math::PseudoFunction::absgt(localMax, globMax) )
            {
                globMax = localMax;
                row = r;

                // 'col' can only be modified when 'fullp' equals true
                if ( true == fullp )
                {
                    col = c;
                }
            }
        }  // omp critical
    }  // omp parallel

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
 * @param physSwap - should the internal algorithm perform physical swapping of matrix elements
 * 
 * @throw MatrixException if any argument is invalid or allocation of a column bookkeeping vector fails
 */
template <class T>
void pivot(
        math::MatrixGeneric<T>& A,
        const bool fullp,
        const bool fullm,
        math::MatrixGeneric<T>* const pB,
        std::size_t* const pRank,
        T* const pDet,
        const bool physSwap
      ) throw(math::MatrixException)
{
    // sanity check
    // if B is provided it must have the same nr. of rows as A
    if ( NULL!=pB && A.nrRows()!=pB->nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const std::size_t NR = A.nrRows();
    const std::size_t NC = A.nrColumns();
    const std::size_t N = std::min(NR, NC);  // the actual nr. of rows and columns to process
    const std::size_t NT = ( NULL!=pB ? pB->nrColumns() : 0 );

    const bool detRequired = ( NULL != pDet );

    /*
     * Bookkeeping of swapped columns is necessary regardless of 'physSwap' when
     * 'pB' is provided and full pivoting is requested. More details follow later.
     */
    const bool fullPB = (NULL != pB) && (true==fullp);


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
    std::vector<std::size_t> colidx;
    std::vector<std::size_t> rowidx;

    try
    {
        if ( false == physSwap )
        {
            math::util::fillVectorWithConsecutiveIndices(NR, rowidx);
        }

        if ( false == physSwap || true == fullPB )
        {
            math::util::fillVectorWithConsecutiveIndices(NC, colidx);
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // pointers to vectors above, required by __index()
    std::vector<std::size_t>* const pRows = ( false==physSwap ? &rowidx : NULL );
    std::vector<std::size_t>* const pCols = ( false==physSwap ? &colidx : NULL );

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

    for ( std::size_t i=0; i<N; ++i )
    {
        // Find the most appropriate row and column to pivot with this one
        std::size_t pr = i;
        std::size_t pc = i;
        math::Pivot::__private::findPivot(A, i, pRows, pCols, pr, pc, fullp);

        // If even the highest absolute value equals 0, there is
        // no unique solution of the system of linear equations
        if ( true == math::NumericUtil::isZero<T>(
             A( math::Pivot::__private::__index(pRows, pr),
                math::Pivot::__private::__index(pCols, pc) ) ) )
        {
            // in this case , the rank will be equal to the current 'i'...
            if ( NULL != pRank )
            {
                *pRank = i;
            }

            // ... and the determinant will be zero
            if ( true == detRequired )
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

            if ( false == physSwap )
            {
                std::swap( rowidx.at(i), rowidx.at(pr) );
            }
            else
            {
                A.swapRows_(i, pr);
            }

            if ( NULL != pB )
            {
                pB->swapRows_(i, pr);
            }
        }

        // swap columns of 'A' if necessary
        if ( true==fullp && pc!=i )
        {
            oddSwaps = !oddSwaps;

            if ( true == physSwap )
            {
                A.swapColumns_(i, pc);
            }

            if ( false==physSwap || true==fullPB )
            {
                std::swap( colidx.at(i), colidx.at(pc) );
            }
        }

        // Set the i.th column of all other rows (r>i) to 0 by
        // adding the appropriate multiple of the i.th row
        #pragma omp parallel for default(none) shared(A, i)
        for ( std::size_t r=i+1; r<NR; ++r )
        {
            // Nothing to do if temp(r,i) is already 0.
            const T& Ari = A(math::Pivot::__private::__index(pRows, r), math::Pivot::__private::__index(pCols, i));
            if ( false == math::NumericUtil::isZero<T>(Ari) )
            {
                // Subtract a multiple of the i.th row.
                const T& Aii = A(math::Pivot::__private::__index(pRows, i), math::Pivot::__private::__index(pCols, i));
                const T el = Ari / Aii;

                // A(r,:) = A(r,:) - el*A(i,:)
                for ( std::size_t c=i; c<NC; ++c )
                {
                    T& Arc = A(math::Pivot::__private::__index(pRows, r), math::Pivot::__private::__index(pCols, c));
                    const T& Aic = A(math::Pivot::__private::__index(pRows, i), math::Pivot::__private::__index(pCols, c));
                    Arc -= el * Aic;
                }

                // b(r,:) = b(r,:) - el*b(i,:)
                if ( NULL != pB )
                {
                    math::MatrixGeneric<T>& B = *pB;

                    for ( std::size_t c=0; c<NT; ++c )
                    {
                        B(r, c) -= el * B(i, c);
                    }
                }
            }  // if A(r,i) != 0
        }  // for r
    }  // for i

    // If this point is reached, A is a full rank matrix (rank=min(NR, NC)):
    if ( NULL != pRank )
    {
        *pRank = N;
    }

    /*
     * Initialize the determinant to 1 if applicable.
     * Later it will be multiplied by diagonal elements of 'A'.
     * If the number of all row and column swaps is odd,
     * the sign of the determinant will be reversed, hence
     * the initial determinant will be set to -1 in this case.
     */
    if ( true == detRequired )
    {
        *pDet = ( false==oddSwaps ? static_cast<T>(1) : static_cast<T>(-1) );
    }


    /*
     * Set the diagonal elements of 'A' to 1 by dividing the
     * whole row by A(r,r). Columns smaller than 'r' are already equal to 0.
     * If applicable, the determinant is also calculated by this loop.
     * This operation only makes sense when a system of linear equations is
     * being solved, i.e. when 'b' is provided or a determinant is required.
     */

    if ( NULL!=pB || true==detRequired )
    {
        // If pB is NULL, B will be assigned to A and will be never modified
        math::MatrixGeneric<T>& B = ( NULL!=pB ? *pB : A );

        // The ideal number of threads depends on the actual task
        // that can be deducted from availability of 'B'
        const std::size_t NTHR = ( NULL!=pB ? N : ompIdeal(N) );

        /*
         * Normalizing of each row is independent from other rows so it is
         * perfectly safe to parallelize the task by rows.
         */
        #pragma omp parallel num_threads(NTHR) default(none) shared(A, B)
        {
            // Product of diagonal elements of rows assigned to the thread
            T diagProd = static_cast<T>(1);

            #pragma omp for
            for ( std::size_t r=0; r<N; ++r )
            {
                const T el = A(math::Pivot::__private::__index(pRows, r), math::Pivot::__private::__index(pCols, r));

                // Normalization of rows is actually only necessary
                // when 'B' is provided
                if ( NULL != pB )
                {
                    for ( std::size_t c=r; c<NC; ++c )
                    {
                        T& Arc = A(math::Pivot::__private::__index(pRows, r), math::Pivot::__private::__index(pCols, c));
                        Arc /= el;
                    }

                    for ( std::size_t c=0; c<NT; ++c )
                    {
                        B(r, c) /= el;
                    }
                }

                // (partial) product of diagonal elements when applicable
                if ( NULL != pDet )
                {
                    diagProd *= el;
                }
            }  // for r

            // If determinant is required multiply *pDet by all diagProds.
            // Note that this operation must be reduced to prevent race conditions.
            if ( true == detRequired )
            {
                #pragma omp critical(pivotgeneric_detmult)
                {
                    *pDet *= diagProd;
                }
            }

        }  // omp parallel

        (void) NTHR;
    }  // if pb!=NULL || detRequired

    /*
     * 'A' is now an upper diagonal matrix with ones on the diagonal.
     * No additional processing of 'A' and 'b' is necessary if 'fullm' is false.
     */


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
        for ( std::size_t c=1; c<NC; ++c )
        {
            // The current row to apply the operation described above

            /*
             * It is possible to parallelize this for loop because a row 'c' (not included
             * into the for loop) will be added independently to each "parallelized" row
             * 'r' (between 0 and c-1 incl.), thus no race condition is possible.
             */
            #pragma omp parallel for \
                        default(none) \
                        shared(A, c) \
                        schedule(dynamic)
            for ( std::size_t r=0; r<c; ++r )
            {
                // Nothing to do if A(r,c) already equals 0
                const T& Arc = A(math::Pivot::__private::__index(pRows, r), math::Pivot::__private::__index(pCols, c));
                if ( false == math::NumericUtil::isZero<T>(Arc) )
                {
                    /*
                     * To set A(r,c) to 0 it is a good idea to add the c.th row to it.
                     * A(c,i); i<c are already 0 (i.e. will not affect anything left of A(i,c)
                     * and A(c,c) is already 1.
                     */

                    const T el = Arc;

                    // A(r,:) = A(r,:) - el*A(c,:)
                    for ( std::size_t i=c; i<NC; ++i )
                    {
                        T& Ari = A(math::Pivot::__private::__index(pRows, r), math::Pivot::__private::__index(pCols, i));
                        const T& Aci = A(math::Pivot::__private::__index(pRows, c), math::Pivot::__private::__index(pCols, i));
                        Ari -= el * Aci;
                    }

                    // b(r,:) = b(r,:) - el*b(c,:)
                    if ( NULL != pB )
                    {
                        math::MatrixGeneric<T>& B = *pB;

                        for ( std::size_t i=0; i<NT; ++i )
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

    if ( true==fullPB && NR==NC )
    {
        math::Pivot::rearrangeMatrixRows(*pB, colidx);
    }  // if fullPB
}

}}}  // namespace math::Pivot::__private



/*
 * Tries to permute rows and columns of the matrix 'A' to obtain
 * a diagonally dominant matrix, i.e. with diagonal elements as
 * large as possible (by absolute values). Permutations of rows
 * and columns are written into 'rows' and 'cols'.
 * 
 * If 'pB' is not NULL, its rows will be permuted according to
 * the permutation vector 'rows'.
 * 
 * @note It is not necessary to initialize 'rows' and 'columns' as
 *       the function takes care of it.
 * 
 * @param A - a matrix to be "converted" to diagonally dominant matrix
 * @param pB - pointer to a matrix to permute its rows according to 'rows' (may be NULL if not required)
 * @param rows - reference to vector to write permutations of A's rows to
 * @param cols - reference to vector to write permutations of A's columns to
 * 
 * @return logical value indicating whether all "diagonal" elements are non-zero
 * 
 * @throw MatrixException if 'A' is not a square matrix or nr. of rows of 'B' is invalid
 */
template <class T>
bool math::Pivot::getDiagonallyDominantMatrix(
        const math::MatrixGeneric<T>& A,
        math::MatrixGeneric<T>* const pB,
        std::vector<std::size_t>& rows,
        std::vector<std::size_t>& cols
      ) throw (math::MatrixException)
{
    // sanity check
    if ( false == A.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
    }

    const std::size_t N = A.nrRows();

    if ( NULL!=pB && pB->nrRows()!=N )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    try
    {
        rows.resize(N);
        cols.resize(N);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }
    
    math::util::fillVectorsWithConsecutiveIndices(N, rows, cols);

    for ( std::size_t i=0; i<(N-1); ++i )
    {
        std::size_t pr;
        std::size_t pc;

        math::Pivot::__private::findPivot(A, i, &rows, &cols, pr, pc, true);

        // Is the pivot by absolute value greater than 0?
        if ( true == math::NumericUtil::isZero( A(rows.at(pr), cols.at(pc)) ) )
        {
            return false;
        }

        if ( i != pr )
        {
            std::swap( rows.at(i), rows.at(pr) );
            if ( NULL != pB )
            {
                pB->swapRows_(i, pr);
            }
        }

        if ( i != pc )
        {
            std::swap( cols.at(i), cols.at(pc) );
        }
    }

    return true;
}


/*
 * Rearranges the matrix' rows to the order, required
 * by full pivoting.
 * 
 * @note The function also rearranges elements of 'cols'.
 * 
 * @param x - the matrix to rearrange its rows
 * @param cols - permutation vector of columns
 */
template <class T>
void math::Pivot::rearrangeMatrixRows(
        math::MatrixGeneric<T>& x,
        std::vector<std::size_t>& cols )
{
    /*
     * The algorithm below will check if cols(i) equals 'i'.
     * If it does not, it will find the index of cols' element that does equal
     * 'i' and swap the rows of 'x' accordingly. Additionally this swap will
     * be reflected by swapping of corresponding elements of 'cols'.
     */

    const std::size_t N = x.nrRows();

    for ( std::size_t idx=0; idx<N; ++idx )
    {
        // find such 'colIdx' that satisfies: cols(colIdx) == idx
        // TODO: does it make any sense to parallelize this simple operation?
        std::size_t colIdx;
        for ( colIdx=idx; cols.at(colIdx)!=idx; ++colIdx );

#if 0
        /*
         * Experimental code that parallelizes searching within 'cols'.
         *
         * It replaces one simple line (above) and any potential
         * benefits of this "complication" have not been researched yet.
         * Until then the code will be "commented out" and will remain within
         * the #if 0 block for possible future reconsideration.
         */

        // Do not parallelize if only a handful of candidates remain
        if ( (N-idx) < OMP_CHUNKS_PER_THREAD )
        {
            for ( colIdx=idx; cols.at(colIdx)!=idx; ++colIdx );
        }
        else
        {
            /*
             * Note that this flag will be updated exactly once
             * hence no "synchronization" (e.g. critical section)
             * is necessary.
             */
            volatile bool foundFlag = false;
 
            #pragma omp parallel num_threads(ompIdeal(N-idx)) \
                    default(none) shared(cols, foundFlag, idx, colIdx)
            {
                OMP_COARSE_GRAINED_PAR_INIT_VARS(N-idx);

                const std::size_t idxStart = istart + idx;
                const std::size_t idxEnd = iend + idx;
                typename std::vector<std::size_t>::const_iterator it = cols.begin() + idxStart;
                for ( std::size_t currIdx=idxStart;
                      false==foundFlag && currIdx<idxEnd && it!=cols.end();
                      ++it, ++currIdx )
                {
                    const std::size_t& currCol = *it;
    
                    if ( currCol == idx )
                    {
                        colIdx = currIdx;
                        foundFlag = true;
                        #pragma omp flush(foundFlag)

                        break;  // out of for it
                    }

                    // Maybe another thread has updated the flag...
                    #pragma omp flush(foundFlag)
                }
            }  // omp parallel
        }
#endif

        // and swap x's rows and cols' elements if necessary
        if ( colIdx != idx )
        {
            x.swapRows_(idx, colIdx);
            std::swap(cols.at(idx), cols.at(colIdx));
        }
    }  // for idx
}


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
 * @param physSwap - should the internal algorithm perform physical swapping of matrix elements
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
        const bool fullp,
        const bool physSwap
      ) throw(math::MatrixException)
{
    // sanity check
    if ( false == A.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const std::size_t N = A.nrColumns();          // Nr. of unknowns

    // Check the dimensions
    if ( N != b.nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    math::MatrixGeneric<T> temp(A);
    x = b;
    
    std::size_t rank;

    math::Pivot::__private::pivot<T>(temp, fullp, true, &x, &rank, NULL, physSwap);

    // if rank(A) equals N, the system of lin. equations has a unique solution

    return ( N == rank );
}


/*
 * Calculates determinant of a square matrix.
 * 
 * @param A - a matrix to calculate its determinant
 * @param fullp - should the algorithm perform full pivoting?
 * @param physSwap - should the internal algorithm perform physical swapping of matrix elements
 * 
 * @return determinant of 'A'
 * 
 * @throw MatrixException if 'A' is not a square matrix or if internal allocation of memory failed
 */
template <class T>
T math::Pivot::getDeterminant(
        const math::MatrixGeneric<T>& A,
        const bool fullp,
        const bool physSwap
      ) throw(math::MatrixException)
{
    // sanity check
    if ( false == A.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
    }

    math::MatrixGeneric<T> temp(A);
    T det;

    math::Pivot::__private::pivot<T>(temp, fullp, false, NULL, NULL, &det, physSwap);

    return det;
}


/*
 * Calculates rank of a matrix.
 *
 * @param A - a matrix to calculate its rank
 * @param physSwap - should the internal algorithm perform physical swapping of matrix elements
 *
 * @return rank of 'A'
 *
 * @throw MatrixException if internal allocation of memory failed
 */
template <class T>
std::size_t math::Pivot::getRank(
        const math::MatrixGeneric<T>& A,
        const bool physSwap
      ) throw(math::MatrixException)
{
    math::MatrixGeneric<T> temp(A);
    std::size_t rank;

    math::Pivot::__private::pivot<T>(temp, true, false, NULL, &rank, NULL, physSwap);

    return rank;
}
