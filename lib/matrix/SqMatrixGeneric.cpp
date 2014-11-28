/*
Copyright 2011, Jernej Kovacic

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
 * Implementation of the class SqMatriGeneric.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */

// Deliberately there is no #include "SqMatrixGeneric.hpp" !
#include "util/mtcopy.hpp"
#include "util/NumericUtil.hpp"
#include "exception/MatrixException.hpp"
#include "matrix/MatrixGeneric.hpp"
#include "lineq/LinearEquationSolverGeneric.hpp"
#include "exception/LinearEquationSolverException.hpp"
#include "../settings/omp_settings.h"
#include "omp/omp_header.h"
#include "omp/omp_coarse.h"

#include <cstddef>
#include <stdexcept>
#include <vector>
#include <algorithm>


/**
 * Constructor.
 * Creates an instance of a square matrix with the specified number of rows and columns.
 *
 * @param dim  - number of rows and columns (default: 1)
 *
 * @throw MatrixException when allocation of memory fails or in case of incorrect
 *        input parameters (dim must be at least 1)
 */
template<class T>
math::SqMatrixGeneric<T>::SqMatrixGeneric(size_t dim) throw (math::MatrixException) :
    math::MatrixGeneric<T>(dim, dim)
{
    // Square matrices have the same number of rows and columns.
    // This is specified by the initializer list.
    // The class has no extra members, so nothing else to do
}

/**
 * Copy constructor.
 * Creates an instance of a matrix with the same dimensions as 'orig'
 * and copies its elements. 'orig' must have the same number of rows and columns,
 * otherwise an exception is thrown.
 *
 * @param orig - original matrix to be copied into this one
 *
 * @throw MatrixException if 'orig' has different number of rows and columns
 */
template<class T>
math::SqMatrixGeneric<T>::SqMatrixGeneric(const math::MatrixGeneric<T>& orig) throw(math::MatrixException) :
    math::MatrixGeneric<T>(orig)
{
    // check of input parameters
    if ( orig.nrRows() != orig.nrColumns() )
    {
        throw MatrixException(MatrixException::INVALID_DIMENSION);
    }
}

/**
 * Assignment operator (=)
 *
 * @param orig - a matrix to be copied into this one
 *
 * @return reference to this
 *
 * @throw MatrixException if memory allocation fails or attempting to assign a non-square matrix
 */
template<class T>
math::SqMatrixGeneric<T>& math::SqMatrixGeneric<T>::operator= (const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &m )
    {
        return *this;
    }

    // m must be a "square" matrix, i.e. with equal number of rows and columns
    if ( m.nrRows() != m.nrColumns() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    this->_copyElems(m);

    return *this;
}

/**
 * Modifies the matrix into a diagonal matrix.
 * Values of all diagonal elements are set to scalar, the other ones to zero.
 *
 * @param scalar - value of diagonal elements
 *
 * @return reference to itself
 *
 * @throw MatrixException
 */
template<class T>
math::SqMatrixGeneric<T>& math::SqMatrixGeneric<T>::setDiag(const T& scalar) throw(math::MatrixException)
{
    // A double for loop will traverse the matrix, its diagonal elements
    // (row == column) will be set to the scalar, others to 0

    const size_t N = this->rows;
    const size_t N2 = N * N;

    // Coarse grained parallelism:
    std::vector<T>& els = this->elems;

    #pragma omp parallel num_threads(ompIdeal(N2)) \
                if(N2>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(els, scalar)
    {
        const size_t thnr = omp_get_thread_num();
        const size_t nthreads  = omp_get_num_threads();
        const size_t elems_per_thread = (N2 + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thnr;

        typename std::vector<T>::iterator it = els.begin() + istart;
        for ( size_t idx = istart;
              idx<(istart+elems_per_thread) && it!=els.end();
              ++it, ++idx )
        {
            const size_t r = idx / N;
            const size_t c = idx % N;

            *it = ( r==c ? scalar : static_cast<T>(0) );
        }
    }  // omp parallel

    return *this;
}

/**
 * Modifies the matrix into a unit matrix (a diagonal matrix with ones on the diagonal)
 *
 * @return reference to itself
 *
 * @throw MatrixException
 */
template<class T>
math::SqMatrixGeneric<T>& math::SqMatrixGeneric<T>::setUnit() throw(math::MatrixException)
{
    // Actually this is a diagonal matrix with units (ones)
    // on its diagonal
    setDiag(static_cast<T>(1));

    return *this;
}

/**
 * Calculates matrix's determinant.
 * This operation makes sense if T is float, double, Rational, Complex.
 * The result may be wrong if T is any implementation of int !!!!!!
 *
 * @return determinant
 *
 * @throw MatrixException if allocation of memory for auxiliary variables fails
 */
template<class T>
T math::SqMatrixGeneric<T>::determinant() const throw(math::MatrixException)
{
    /*
     * The following properties of determinant will be utilized:
     *
     * - addition of one line with a multiple of another line does not affect
     *   the determinant
     * - swapping of two lines negates the determinant
     * - if a matrix is "converted" into an upper or lower triangle
     *   (all elements above/below the diagonal are 0), the determinant
     *   is a product of all diagonal elements
     */

    // Initial value. It will be negated each time two lines need to be swapped.
    // At the end of the algorithm it will be multiplied by all diagonal elements
    T retVal = static_cast<T>(1);

    // We do not want to modify the matrix, therefore its vector of
    // elements will be copied into a temporary one where any modifications
    // are permitted
    std::vector<T> temp;
    math::mtcopy(this->elems, temp);
    const size_t N = this->rows;  // number of rows (and columns)

    /*
     * First part of the algorithm just finds the first occurrence of a
     * row where A(i,i) does not equal 0. As such, this part is not
     * suitable for parallelization.
     */
    for ( size_t i=0; i<N-1; ++i )
    {
        // if temp(i,i) equals zero, swap the i^th line with
        // another one (r; r>i) satisfying temp(r,i)!=0
        // Each swap multiplies the determinant by -1

        if ( true == math::NumericUtil::isZero<T>(temp.at(this->_pos(i, i))) )
        {
            size_t r;

            // Line swap will be necessary.
            // Find the r^th line meeting the criteria above
            // If not found, the determinant will be 0.

            for ( r=i+1; r<N; ++r )
            {
                if ( false == NumericUtil::isZero<T>(temp.at(this->_pos(r, i))) )
                {
                    // Found, no need to search further,
                    // so end the for (r) loop
                    break;  // out of for (r)
                }
            }  // for r

            // if no appropriate row r was found, the determinant will be 0
            // and the method is finished
            if ( N == r )
            {
                return static_cast<T>(0);
            }

            // otherwise swap the lines one element by one

            // swap i.th and r.th line by replacing elements one by one

            // However the swapping part might conditionally be suitable for parallelization
            #pragma omp parallel for \
                        if((N-i)>OMP_CHUNKS_PER_THREAD) \
                        default(none) shared(temp, r, i)
            for ( size_t c=this->_pos(i,i); c<this->_pos(i+1,0); ++c )
            {
                std::swap(temp.at(c), temp.at(this->_pos(r, c)));
            }

            // finally, if two lines are swapped, det = -det
            retVal = -retVal;
        } // if temp(i,i) equals 0

        // Now temp(i,i) definitely does not equal 0
        // all temp(r,i) will be "set" to 0 where r>i.

        // Note that determinant is not changed if
        // a multiplier of one line (i in this algorithm)
        // is added to another line (r; r>i)

        /*
         * Main part of the algorithm. An appropriate multiplier of the i.th row will be
         * added to each row 'r' (r>i) so that temp(r,i) will be equal to zero.
         *
         * Note: only the outer for loop (for r) will be parallelized.
         */
        #pragma omp parallel for default(none) shared(temp, i)
        for ( size_t r=i+1; r<N; ++r )
        {
            // temp(r,i) will be calculated to 0 immediately.
            // However, its initial value is necessary to properly
            // calculate all other elements of the r.th row
            const T ri = temp.at(this->_pos(r, i));

            for ( size_t c=i; c<N; ++c)
            {
                // temp(r,c) = temp(r,c) - temp(i,c) * temp(r,i) / temp(i,i)
                temp.at(this->_pos(r, c)) -= temp.at(this->_pos(i, c)) * ri / temp.at(this->_pos(i, i));
            }  // for c
        }  // for r
    }  // for i


    /*
     * Now 'temp' is an upper triangular matrix so all its diagonal
     * elements can be multiplied.
     */

    // Coarse grained parallelism
    T prod = static_cast<T>(1);

    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(temp, prod)
    {
        const size_t thnr = omp_get_thread_num();
        const size_t nthreads  = omp_get_num_threads();
        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thnr;
        const size_t iend = std::min(istart+elems_per_thread, N);

        T tempProd = static_cast<T>(1);
        for ( size_t i = istart; i<iend; ++i )
        {
            tempProd *= temp.at(this->_pos(i, i));
        }

        // Multiply in a thread safe manner:
        #pragma omp critical(sqmatrix_determinant)
        {
            prod *= tempProd;
        }
    }  // omp parallel

    retVal *= prod;

    // temp not needed anymore, clean it
    temp.clear();

    return retVal;
}

/**
 * Matrix inversion.
 * R = A^(-1) if R*A = A*R = I
 * If matrix's determinant equals 0, the matrix is not invertible
 *
 * @return inverse matrix
 *
 * @throw MatrixException if the matrix is not invertible
 */
template<class T>
math::SqMatrixGeneric<T> math::SqMatrixGeneric<T>::inverse() const throw(math::MatrixException)
{
    /*
     * Elements of the inverse matrix can be calculated as quotients of so called
     * "co-determinants" and determinant of the main matrix. However, this
     * would be too complex, so a numerical Gauss - Jordan elimination algorithm
     * will be applied:
     * - a unit matrix is appended to the right of the matrix: [A|I].
     * - multiples of lines are added to lines (similar as solving a set of linear
     *   equations) until a unit matrix appears in the left half: [I|B]
     * - B is inverse matrix of A: B = A^(-1)
     * 
     * This functionality is already implemented by the class LinearEquationSolverGeneric.
     */

    try
    {
        // prepare an identity matrix NxN...
        SqMatrixGeneric<T> id(this->rows);
        id.setUnit();
        // ... and instantiate LinearEquationSolverGeneric
        LinearEquationSolverGeneric<T> leq(*this, id);

        // inverse matrix is a solution (if it exists) of the equation:
        // this * inv = id
        SqMatrixGeneric<T> retVal;
        leq.solve(retVal);

        return retVal;
    }
    catch ( const LinearEquationSolverException& leqex )
    {
        // is *this an uninvertible matrix? (determinant()=0):
        if ( math::LinearEquationSolverException::NO_UNIQUE_SOLUTION == leqex.error )
        {
            throw math::MatrixException(math::MatrixException::NON_INVERTIBLE_MATRIX);
        }
        else
        {
            // other than a case of non-invertible matrix, exception can only
            // be thrown if allocation of memory failed.
            throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
        }
    }
}

/**
 * Transpose the matrix and write all changes into it
 * (disregard the original matrix)
 *
 * @return reference to this
 *
 * @throw MatrixException (never thrown at this implementation)
 */
template<class T>
math::SqMatrixGeneric<T>& math::SqMatrixGeneric<T>::transposed() throw(math::MatrixException)
{
    // TODO: find and implement a better algorithm

    const size_t N = this->rows;       // number of rows (and columns)
    const size_t Ntr = N * (N-1) / 2;  // number of all elements to be transposed

    // Traverse the upper diagonal part of the matrix,
    // no need to reach the final row and
    // no need to transpose elements on the diagonal:

    /*
     * Notes about parallelization:
     * As the inner loop (for c) also depends on outer loop's
     * iterator (r), it is not possible to parallelize both loops
     * in an elegant way. Hence only the outer loop is parallelized.
     * As the threads' load varies, dynamic scheduling is applied.
     */
    std::vector<T>& els = this->elems;
    #pragma omp parallel for \
                if(Ntr>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(els) \
                schedule(dynamic)
    for ( size_t r=0; r<N-1; ++r )
    {
        for ( size_t c=r+1; c<N; ++c )
        {
            std::swap(els.at(this->_pos(r, c)), els.at(this->_pos(c, r)));
        }  // for c
    }  // for r

    return *this;

    // just to suppress a warning when OpenMP is not enabled
    (void) Ntr;
}

/**
 * "Reimplementation" of operator*= that multiplies a matrix by this one
 * and assigns the product to itself.
 *
 * The product must remain a square matrix, for that reason 'm' must
 * be a square matrix (it can be declared as an "ordinary" matrix, though)
 * with the same dimensions as 'this'
 *
 * @param m - matrix to be multiplied by this one
 *
 * @return reference to itself
 *
 * @throw MatrixException if dimensions do not match or if allocation of memory fails
 */
template<class T>
math::SqMatrixGeneric<T>& math::SqMatrixGeneric<T>::operator*= (const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    // Check of additional condition ('m' must have the same dimensions):
    if ( m.nrRows()!=this->rows || m.nrColumns()!=this->rows )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // If m's dimensions match, dimensions of this will be preserved.
    // Just apply the functionality of the parent class.
    // Also cast the returning reference to SqMatrixGeneric.
    return dynamic_cast<SqMatrixGeneric<T>&>(( MatrixGeneric<T>::operator*=(m) ));
}

/*
 * The following functions are inherited from the base class, however they change
 * matrix dimensions so they should not be permitted in a square matrix.
 * (insert or removal of a row/column of a square matrix results in a non-square matrix!)
 * Therefore they will automatically throw an exception if called.
 */
template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::removeRow(size_t rowNr) throw (math::MatrixException)
{
    throw math::MatrixException(MatrixException::FORBIDDEN);
    (void) rowNr;
}

template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::removeColumn(size_t colNr) throw (math::MatrixException)
{
    throw math::MatrixException(math::MatrixException::FORBIDDEN);
    (void) colNr;
}

template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::insertRow(size_t rowNr, const T& el) throw (math::MatrixException)
{
    throw math::MatrixException(math::MatrixException::FORBIDDEN);
    (void) rowNr;
    (void) el;
}

template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::insertColumn(size_t colNr, const T& el) throw (math::MatrixException)
{
    throw math::MatrixException(math::MatrixException::FORBIDDEN);
    (void) colNr;
    (void) el;
}
