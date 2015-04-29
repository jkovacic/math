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
 * Implementation of the class MatriGeneric.
 */


#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <ostream>
#include <complex>
#include <algorithm>

// no #include "MatrixGeneric.hpp" !!!
#include "util/NumericUtil.hpp"
#include "exception/MatrixException.hpp"
#include "matrix/LinearEquationSolverGeneric.hpp"
#include "util/mtcopy.hpp"
#include "util/mtvectop.hpp"
#include "util/mtswap.hpp"

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"



/*
 * Initializes a new matrix
 * 
 * @param rows    - number of rows
 * @param columns - number of columns
 *
 * @throw MatrixException when allocation of memory fails or in case of incorrect
 *        input arguments (both must be at least 1)
 */
template <class T>
void math::MatrixGeneric<T>::__init(
    const size_t rows, const size_t cols) throw (math::MatrixException)
{
    // Matrix must contain at least 1 row and at least 1 column
    if ( rows < 1 || cols < 1 )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    // even theoretically the number of vector's elements is limited
    if ( cols > this->m_elems.max_size()/rows )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        this->m_rows = rows;
        this->m_cols = cols;
        // make sure, the vector will be empty
        this->m_elems.clear();
        // allocate memory for required number of elements, initialize each of them
        this->m_elems.resize(rows * cols, static_cast<T>(0));
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }
}


/**
 * Constructor.
 * Creates an instance of a matrix with the specified number of rows and columns.
 * Elements of the matrix are set to the default value (zero).
 *
 * @param rows    - number of rows
 * @param columns - number of columns
 *
 * @throw MatrixException when allocation of memory fails or in case of incorrect
 *        input arguments (both must be at least 1)
 */
template <class T>
math::MatrixGeneric<T>::MatrixGeneric(const size_t rows, const size_t columns) throw(math::MatrixException)
{
    this->__init(rows, columns);
}

/**
 * Constructor.
 * Creates an instance of a square matrix with the specified number of rows and columns.
 *
 * @param n  - number of rows and columns
 *
 * @throw MatrixException when allocation of memory fails or in case of incorrect
 *        input arguments ('n' must be at least 1)
 */
template <class T>
math::MatrixGeneric<T>::MatrixGeneric(const size_t n) throw (math::MatrixException)
{
    this->__init(n, n);
}


/**
 * Copy constructor.
 * Creates an instance of a matrix with the same dimensions as 'orig'
 * and copies its elements.
 *
 * @param orig - original matrix to be copied into this one
 *
 * @throw MatrixException when allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T>::MatrixGeneric(const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    this->__copyElems(orig);
}


// Copy elements from one matrix into another. Used at copy constructors,
// assignment operators etc.
template <class T>
void math::MatrixGeneric<T>::__copyElems(const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    // Check if the original matrix is OK (if its m_rows and m_cols are set correctly)
    if ( orig.m_rows * orig.m_cols != orig.m_elems.size() )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    try
    {
        this->m_rows = orig.m_rows;
        this->m_cols = orig.m_cols;
        math::mtcopy<T>(orig.m_elems, this->m_elems);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Check if the correct number of elements was copied
    if ( this->m_rows * this->m_cols != this->m_elems.size() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }
}


/**
 * @return number of rows
 */
template <class T>
size_t math::MatrixGeneric<T>::nrRows() const
{
    return this->m_rows;
}


/**
 * @return number of columns
 */
template <class T>
size_t math::MatrixGeneric<T>::nrColumns() const
{
    return this->m_cols;
}


/**
 * Returns an element at the specified location (row and column)
 *
 * @param row number (starting with 0)
 * @param column number (starting with 0)
 *
 * @return element at the location
 *
 * @throw MatrixException if input parameters are out of range
 */
template <class T>
T math::MatrixGeneric<T>::get(const size_t row, const size_t column) const throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= this->m_rows || column >= this->m_cols )
    {
        // At least one input parameter is out of range
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->m_elems.at(this->__pos(row, column));
}


/**
 * Returns a read-write reference to the desired element of the matrix.
 * 
 * @note The reference should be used immediately after this call before
 *       the matrix is destroyed or its internal storage is reallocated.
 * 
 * @note This function is deprecated.
 *
 * @param row - row number (starting with 0)
 * @param column - column number (starting with 0)
 *
 * @return read-write reference to the element at the desired position
 *
 * @throw MatrixException if 'row' and/or 'column' are out of range
 * 
 * @deprecated
 */
template <class T>
T& math::MatrixGeneric<T>::at(const size_t row, const size_t column) throw (math::MatrixException)
{
    // Check if input arguments are within the matrix's range
    if ( row >= this->m_rows || column >= this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->m_elems.at(this->__pos(row, column));
}


/**
 * Returns a read-only (const) reference to the desired element of the matrix.
 * 
 * @note The reference should be used immediately after this call before
 *       the matrix is destroyed or its internal storage is reallocated.
 * 
 * @note This function is deprecated
 *
 * @param row - row number (starting with 0)
 * @param column - column number (starting with 0)
 *
 * @return read-only reference to the element at the desired position
 *
 * @throw MatrixException if 'row' and/or 'column' are out of range
 * 
 * @deprecated
 */
template <class T>
const T& math::MatrixGeneric<T>::at(const size_t row, const size_t column) const throw (math::MatrixException)
{
    /*
     * Implementation is actually the same as implementation of another at()
     * with non-const signature. A single macro could be used for both
     * implementations, however both functions are short, simple and unlikely
     * to change often.
     */

    // Check if input arguments are within the matrix's range
    if ( row >= this->m_rows || column >= this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->m_elems.at(this->__pos(row, column));
}


/**
 * Returns a read-write reference to the desired element of the matrix.
 * @note The reference should be used immediately after this call before
 *       the matrix is destroyed or its internal storage is reallocated.
 *
 * @param row - row number (starting with 0)
 * @param column - column number (starting with 0)
 *
 * @return read-write reference to the element at the desired position
 *
 * @throw MatrixException if 'row' and/or 'column' are out of range
 */
template <class T>
T& math::MatrixGeneric<T>::operator()(const size_t row, const size_t column) throw(math::MatrixException)
{
    // Check if input arguments are within the matrix's range
    if ( row >= this->m_rows || column >= this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->m_elems.at(this->__pos(row, column));
}


/**
 * Returns a read-only (const) reference to the desired element of the matrix.
 * @note The reference should be used immediately after this call before
 *       the matrix is destroyed or its internal storage is reallocated.
 *
 * @param row - row number (starting with 0)
 * @param column - column number (starting with 0)
 *
 * @return read-only reference to the element at the desired position
 *
 * @throw MatrixException if 'row' and/or 'column' are out of range
 */
template <class T>
const T& math::MatrixGeneric<T>::operator()(const size_t row, const size_t column) const throw(math::MatrixException)
{
    /*
     * Implementation is actually the same as implementation of another operator()
     * with non-const signature. A single macro could be used for both
     * implementations, however both functions are short, simple and unlikely
     * to change often.
     */

    // Check if input arguments are within the matrix's range
    if ( row >= this->m_rows || column >= this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->m_elems.at(this->__pos(row, column));
}


/**
 * Read-write reference to the element at the desired position,
 * denoted by a linear index.
 * 
 * 'idx' denotes the element of a long column vector, composed by
 *  consecutive column vectors.
 * 
 * If R denotes the number of rows, the function returns
 *   this(idx/R, idx%R)
 * 
 * @param idx - linear index of the desired element
 * 
 * @return read-write reference to the element at the linear index 'idx'
 * 
 * @throw MatrixException if 'idx' is out of the matrix's range
 */
template <class T>
T& math::MatrixGeneric<T>::operator()(const size_t idx) throw(math::MatrixException)
{
    // check if 'idx' is within elems' range
    if ( idx >= this->m_elems.size() )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    const size_t col = idx / this->m_rows;
    const size_t row = idx % this->m_rows;

    return this->m_elems.at(this->__pos(row, col));
}


/**
 * Read-only reference to the element at the desired position,
 * denoted by a linear index.
 * 
 * 'idx' denotes the element of a long column vector, composed by
 *  consecutive column vectors.
 * 
 * If R denotes the number of rows, the function returns
 *   this(idx/R, idx%R)
 * 
 * @param idx - linear index of the desired element
 * 
 * @return read-only reference to the element at the linear index 'idx'
 * 
 * @throw MatrixException if 'idx' is out of the matrix's range
 */
template <class T>
const T& math::MatrixGeneric<T>::operator()(const size_t idx) const throw(math::MatrixException)
{
    // check if 'idx' is within elems' range
    if ( idx >= this->m_elems.size() )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    const size_t col = idx / this->m_rows;
    const size_t row = idx % this->m_rows;

    return this->m_elems.at(this->__pos(row, col));
}


/**
 * Assigns value at the requested location.
 * This is the only allowed method to modify values of matrix's elements
 * (except of class's internal functions)
 *
 * @param row - row of the element to modify
 * @param column - column of the element to modify
 * @param element - value to be assigned at the requested location
 *
 * @return reference to itself
 *
 * @throw MatrixException if input parameters are out of range
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::set(const size_t row, const size_t column, const T& element) throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= this->m_rows || column >= this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    this->m_elems.at(this->__pos(row, column)) = element;

    return *this;
}


/**
 * Display the matrix to stdout
 *
 * @param str - output stream where the matrix will be displayed (default: cout)
 *
 * @throw MatrixException if attempting to access an element out of allocated range
 */
template <class T>
void math::MatrixGeneric<T>::display(std::ostream& str) const throw (math::MatrixException)
{

    /*
     * Elements of each row are separated by tabs.
     * This does not guarantee that elements of the same column will be
     * displayed under each other and in case of a long row (its output
     * is longer than terminal's line width), the overall output may look rather
     * confusing. However, this is more or less "just" an auxiliary function, mainly used for
     * testing purposes, and the effort was focused to other functionalities.
     * Anyway, it would be nice to improve it in future.
     */

    const size_t tabsPerRow = this->m_cols - 1;

    for ( size_t r=0; r<(this->m_rows); ++r )
    {
        // display elements of the row r, separated by tabs
        for ( size_t c=0; c<(this->m_cols); ++c )
        {
            str << this->m_elems.at(this->__pos(r, c));
            if ( c < tabsPerRow )
            {
                str << "\t";
            }
        }
        // the line must be terminated by a newline
        str << std::endl;
    } // for r

}


/**
 * Assignment operator (=)
 *
 * @param orig - a matrix to be copied into this one
 *
 * @return reference to this
 *
 * @throw MatrixException if memory allocation fails
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator= (const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &orig )
    {
        return *this;
    }

    this->__copyElems(orig);

    return *this;
}


/**
 * Addition operator (+=) that adds a matrix to this and assigns the sum to itself.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param m - matrix to be added to this one
 *
 * @return reference to this
 *
 * @throw MatrixException if dimensions do not match
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator+= (const math::MatrixGeneric<T>& m) throw(math::MatrixException)
{
    // Check if dimensions of both matrices match
    if ( this->m_rows != m.m_rows || this->m_cols != m.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // For a definition of matrix addition, see operator+
    math::mtvectadd<T>(this->m_elems, m.m_elems, this->m_elems, true);

    return *this;
}


/**
 * Subtraction operator (-=) that subtracts a matrix from this and assigns the difference to itself.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param matrix to be subtracted from this one
 *
 * @return reference to this
 *
 * @throw MatrixException if dimensions do not match
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator-= (const math::MatrixGeneric<T>& matrix) throw(math::MatrixException)
{
    // Check if dimensions of both matrices match
    if ( this->m_rows != matrix.m_rows || this->m_cols != matrix.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // For a definition of matrix subtraction, see operator-
    math::mtvectadd<T>(this->m_elems, matrix.m_elems, this->m_elems, false);

    return *this;
}


/**
 * Compound addition operator (+=) that adds each matrix's element
 * by a scalar and assigns the sum to itself.
 *
 * @param scalar
 *
 * @return reference to itself
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator+=(const T& scalar)
{
    math::mtvectscalaradd<T>(this->m_elems, scalar, this->m_elems, true, true);

    return *this;
}


/**
 * Compound subtraction operator (+=) that subtracts each matrix's element
 * by a scalar and assigns the difference to itself.
 *
 * @param scalar
 *
 * @return reference to itself
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator-=(const T& scalar)
{
    math::mtvectscalaradd<T>(this->m_elems, scalar, this->m_elems, false, true);

    return *this;
}


/**
 * Multiplication operator (*=) that multiplies a matrix by this one and assigns
 * the product to itself.
 * Number of columns of 'this' must be the same as number of rows of 'm',
 * otherwise multiplication is not possible.
 *
 * @param m - matrix to be multiplied by this one
 *
 * @return reference to itself
 *
 * @throw MatrixException if dimensions do not match or if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator*= (const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    // for a definition of matrix multiplication, see operator*

    /*
     * operator* will perform checking of numbers of rows/columns etc.
     * and throw an appropriate exception if necessary
     */
    math::MatrixGeneric<T> temp = *this * m;

    try
    {
        this->__copyElems(temp);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    return *this;
}


/**
 * Multiplication operator (*=) that multiplies a matrix by a scalar
 * and assigns the product to itself.
 *
 * @param scalar
 *
 * @return reference to itself
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator*=(const T& scalar)
{
    // Multiply each element by the 'scalar'
    math::mtvectmult<T>(this->m_elems, scalar, this->m_elems);

    return *this;
}


/**
 * Compound division operator (/=) that divides each element
 * of the matrix by 'scalar'and assigns the quotient to itself.
 *
 * @param scalar - divisor that divides each element of 'this'
 *
 * @return reference to itself
 *
 * @throw MatrixException if attempting to divide by 0 or allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator/=(const T& scalar) throw (math::MatrixException)
{
    // Division by zero is not permitted:
    if ( true == math::NumericUtil::isZero<T>(scalar) )
    {
        throw math::MatrixException(math::MatrixException::FORBIDDEN);
    }

    const T f = static_cast<T>(1) / scalar;

    // Multiply each element of 'this' by 1/scalar:
    math::mtvectmult(this->m_elems, f, this->m_elems);

    return *this;
}


/**
 * Element wise multiplication of 'this' by the given matrix.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param m - matrix to be element wise multiplied by this one
 *
 * @return reference to this
 *
 * @throw MatrixException if dimensions do not match
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::ewMult(const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    // Check if dimensions of both matrices match
    if ( this->m_rows != m.m_rows || this->m_cols != m.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    math::mtvectewmult<T>(this->m_elems, m.m_elems, this->m_elems, true);

    return *this;
}


/**
 * Element wise division of 'this' by the given matrix.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param m - matrix to be element wise multiplied by this one
 *
 * @return reference to this
 *
 * @throw MatrixException if dimensions do not match or any m's element equals 0
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::ewDiv(const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    // Check if dimensions of both matrices match
    if ( this->m_rows != m.m_rows || this->m_cols != m.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const bool succ = math::mtvectewmult<T>(this->m_elems, m.m_elems, this->m_elems, false);

    if ( false == succ )
    {
        throw math::MatrixException(math::MatrixException::FORBIDDEN);
    }

    return *this;
}


/*
 * Extracts the specified triangular part(s) of the matrix and
 * assigns the result to 'dest'.
 *
 * @param dest - a reference to a matrix to assign the result
 * @param upper - extract the upper triangular part (without diagonal)
 * @param lower - extract the lower triangular part (without diagonal)
 * @param diag - extract the diagonal
 */
template <class T>
void math::MatrixGeneric<T>::__triangPart(
    math::MatrixGeneric<T>& dest,
    const bool upper,
    const bool lower,
    const bool diag) const
{
    // it is assumed that dest's 'm_elems' is properly allocated

    const size_t N = this->m_elems.size();

    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        size_t r;
        size_t c;

        // iterator to the final element of the block:
        const typename std::vector<T>::const_iterator final = this->m_elems.begin() + iend;
        // iterator to the first/current element of this->m_elems:
        typename std::vector<T>::const_iterator idx = this->m_elems.begin() + istart;
        // iterator to the current element of dest.m_elems:
        typename std::vector<T>::iterator it = dest.m_elems.begin() + istart;

        for ( size_t i = istart;
                idx != final;
                ++i, ++it, ++idx )
        {
            // row and column of the current index 'i':
            r = i / this->m_cols;
            c = i % this->m_cols;

            // Either copy the current element or set it to 0,
            // depending on input arguments.
            *it = 
                ( (r==c && true==diag) ||
                   (r>c && true==lower) ||
                   (r<c && true==upper) ? 
                      *idx : static_cast<T>(0) );
        }
    }  // omp parallel
}


/**
 * Extracts the upper triangular part of the matrix and optionally
 * the diagonal.
 * 
 * @param inclDiag - should the diagonal be included as well (default: TRUE)
 * 
 * @return upper triangular part of the matrix
 * 
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::upperTriangularPart(
        const bool inclDiag
      ) const throw(math::MatrixException)
{
    math::MatrixGeneric<T> retVal(this->m_rows, this->m_cols);
    this->__triangPart(retVal, true, false, inclDiag);
    return retVal;
}


/**
 * Extracts the lower triangular part of the matrix and optionally
 * the diagonal.
 * 
 * @param inclDiag - should the diagonal be included as well (default: TRUE)
 * 
 * @return lower triangular part of the matrix
 * 
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::lowerTriangularPart(
        const bool inclDiag
      ) const throw(math::MatrixException)
{
    math::MatrixGeneric<T> retVal(this->m_rows, this->m_cols);
    this->__triangPart(retVal, false, true, inclDiag);
    return retVal;
}


/**
 * @return the diagonal part of the matrix
 * 
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::diagPart() const throw(math::MatrixException)
{
    math::MatrixGeneric<T> retVal(this->m_rows, this->m_cols);
    this->__triangPart(retVal, false, false, true);
    return retVal;
}


/**
 * @return a logical value indicating whether the matrix is square
 */
template <class T>
bool math::MatrixGeneric<T>::isSquare() const
{
    return ( this->m_rows == this->m_cols );
}


/**
 * Modifies the matrix into a diagonal matrix.
 * Values of all diagonal elements are set to scalar, the other ones to zero.
 * 
 * @note The method is only supported for square matrices
 *
 * @param scalar - value of diagonal elements
 *
 * @return reference to itself
 *
 * @throw MatrixException if the matrix is not square
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::setDiag(const T& scalar) throw(math::MatrixException)
{
    // Sanity check
    if ( false == this->isSquare() )
    {
        throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
    }

    // A double for loop will traverse the matrix, its diagonal elements
    // (row == column) will be set to the scalar, others to 0

    const size_t& N = this->m_rows;
    const size_t N2 = N * N;

    // Coarse grained parallelism:
    std::vector<T>& els = this->m_elems;

    #pragma omp parallel num_threads(ompIdeal(N2)) \
                if(N2>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(N, els, scalar)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N2);

        typename std::vector<T>::iterator it = els.begin() + istart;
        for ( size_t idx = istart;
              idx<(istart+elems_per_thread) && it!=els.end();
              ++it, ++idx )
        {
            const size_t r = idx / N;
            const size_t c = idx % N;

            *it = ( r==c ? scalar : static_cast<T>(0) );
        }

        (void) iend;
    }  // omp parallel

    return *this;
}


/**
 * Modifies the matrix into a unit matrix (a diagonal matrix with ones on the diagonal)
 * 
 * @note The method is only supported for square matrices
 *
 * @return reference to itself
 *
 * @throw MatrixException if the matrix is not square
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::setUnit() throw(math::MatrixException)
{
    // Actually this is a diagonal matrix with units (ones)
    // on its diagonal
    this->setDiag(static_cast<T>(1));

    return *this;
}


/**
 * Calculates matrix's determinant.
 * This operation makes sense if T is float, double, Rational, Complex.
 * The result may be wrong if T is any implementation of int !!!!!!
 *
 * @note The method is only supported for square matrices
 * 
 * @return determinant of the matrix
 *
 * @throw MatrixException if the matrix is not square or if allocation of memory for auxiliary variables fails
 */
template <class T>
T math::MatrixGeneric<T>::determinant() const throw(math::MatrixException)
{
    // Sanity check
    if ( false == this->isSquare() )
    {
        throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
    }

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

    /*
     * We do not want to modify the matrix, therefore its vector of
     * elements will be copied into a temporary one where any modifications
     * are permitted
     */
    std::vector<T> temp;
    math::mtcopy(this->m_elems, temp);
    const size_t& N = this->m_rows;  // number of rows (and columns)

    /*
     * First part of the algorithm just finds the first occurrence of a
     * row where A(i,i) does not equal 0. As such, this part is not
     * suitable for parallelization.
     */
    for ( size_t i=0; i<N-1; ++i )
    {
        /*
         * if temp(i,i) equals zero, swap the i^th line with
         * another one (r; r>i) satisfying temp(r,i)!=0
         * Each swap multiplies the determinant by -1
         */

        if ( true == math::NumericUtil::isZero<T>(temp.at(this->__pos(i, i))) )
        {
            size_t r;

            /*
             * Line swap will be necessary.
             * Find the r^th line meeting the criteria above
             * If not found, the determinant will be 0.
             */

            for ( r=i+1; r<N; ++r )
            {
                if ( false == NumericUtil::isZero<T>(temp.at(this->__pos(r, i))) )
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
            for ( size_t c=this->__pos(i,i); c<this->__pos(i+1, 0); ++c )
            {
                std::swap(temp.at(c), temp.at(this->__pos(r, c)));
            }

            // finally, if two lines are swapped, det = -det
            retVal = -retVal;
        } // if temp(i,i) equals 0

        /*
         * Now temp(i,i) definitely does not equal 0
         * all temp(r,i) will be "set" to 0 where r>i.
         *
         * Note that determinant is not changed if
         * a multiplier of one line (i in this algorithm)
         * is added to another line (r; r>i)
         */

        /*
         * Main part of the algorithm. An appropriate multiplier of the i.th row will be
         * added to each row 'r' (r>i) so that temp(r,i) will be equal to zero.
         *
         * Note: only the outer for loop (for r) will be parallelized.
         */
        #pragma omp parallel for default(none) shared(N, temp, i)
        for ( size_t r=i+1; r<N; ++r )
        {
            /*
             * temp(r,i) will be calculated to 0 immediately.
             * However, its initial value is necessary to properly
             * calculate all other elements of the r.th row
             */
            const T ri = temp.at(this->__pos(r, i));

            for ( size_t c=i; c<N; ++c)
            {
                // temp(r,c) = temp(r,c) - temp(i,c) * temp(r,i) / temp(i,i)
                temp.at(this->__pos(r, c)) -= 
                    temp.at(this->__pos(i, c)) * ri / temp.at(this->__pos(i, i));
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
                default(none) shared(N, temp, prod)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        T tempProd = static_cast<T>(1);
        for ( size_t i = istart; i<iend; ++i )
        {
            tempProd *= temp.at(this->__pos(i, i));
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
 * @note The method is only supported for square matrices
 * 
 * @return inverse matrix
 *
 * @throw MatrixException if the matrix is not square or not invertible
 */
template <class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::inverse() const throw(math::MatrixException)
{
    // Sanity check
    if ( false == this->isSquare() )
    {
        throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
    }

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

    // prepare an identity matrix NxN...
    math::MatrixGeneric<T> id(this->m_rows);
    id.setUnit();

    // Inverse matrix is a solution (if it exists) of the equation:
    // this * inv = id
    math::MatrixGeneric<T> retVal(id);

    const bool succ = math::LinearEquationSolver::solveGaussJordan<T>(*this, id, retVal);

    // is *this an uninvertible matrix? (determinant()=0):
    if ( false == succ )
    {
        throw math::MatrixException(math::MatrixException::NON_INVERTIBLE_MATRIX);
    }

    return retVal;

}


/**
 * Matrix transpose operation
 *
 * @return this^T
 *
 * @throw MatrixException if matrix does not contain enough elements
 */
template <class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::transpose() const throw (math::MatrixException)
{
    // If dimension of this is (n,m), dimension of its transposed matrix is (m,n)
    // T(r,c) = this(c,r)

    // Create an instance with swapped dimensions
    math::MatrixGeneric<T> retVal(this->m_cols, this->m_rows);

    /*
     * If any dimension equals 1 it is sufficient just to copy
     * the vector 'm_elems' and swap the dimensions.
     */
    if ( 1==this->m_rows || 1==this->m_cols )
    {
        math::mtcopy<T>(this->m_elems, retVal.m_elems);
        return retVal;
    }
    
    const size_t& tcols = this->m_cols;

    const size_t N = this->m_rows * this->m_cols;

    // Coarse grained parallelism
    const std::vector<T>& els = this->m_elems;

    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(retVal, els, tcols)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        for ( size_t idx=istart; idx<iend; ++idx )
        {
            const size_t r = idx / tcols;
            const size_t c = idx % tcols;

            retVal.m_elems.at(retVal.__pos(c, r)) = els.at(this->__pos(r, c));
        }
    }  // omp parallel

    return retVal;
}


/**
 * Transpose the matrix and write all changes into it
 * (disregard the original matrix)
 *
 * @return reference to this
 *
 * @throw MatrixException if not enough memory to perform the operation
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::transposed() throw (math::MatrixException)
{
    /*
     * TODO: find a memory efficient method for a general matrix!!
     */

    // If dimension of this is (n,m), dimension of its transposed matrix is (m,n)
    // T(r,c) = this(c,r)

    if ( 1==this->m_rows || 1==this->m_cols )
    {
        /*
         * If any dimension equals 1, no reallocation of elements is
         * necessary. Instead, only the dimensions must be swapped.
         */

        std::swap(this->m_rows, this->m_cols);
    }
    else if ( true == this->isSquare() )
    {
        // A specialized algorithm for square matrices
        // TODO: find and implement a better algorithm

        const size_t& N = this->m_rows;    // number of rows (and columns)
        const size_t Ntr = N * (N-1) / 2;  // number of all elements to be transposed

        /*
         * Traverse the upper diagonal part of the matrix,
         * no need to reach the final row and
         * no need to transpose elements on the diagonal:
         */

        /*
         * Notes about parallelization:
         * As the inner loop (for c) also depends on outer loop's
         * iterator (r), it is not possible to parallelize both loops
         * in an elegant way. Hence only the outer loop is parallelized.
         * As the threads' load varies, dynamic scheduling is applied.
         */
        std::vector<T>& els = this->m_elems;
        #pragma omp parallel for \
                    if(Ntr>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(els) \
                    schedule(dynamic) shared(N)
        for ( size_t r=0; r<N-1; ++r )
        {
            for ( size_t c=r+1; c<N; ++c )
            {
                std::swap(els.at(this->__pos(r, c)), els.at(this->__pos(c, r)));
            }  // for c
        }  // for r

        // just to suppress a warning when OpenMP is not enabled
        (void) Ntr;
    }
    else
    {
        // Until a better algorithm is implemented
        // just use the general transpose method
        math::MatrixGeneric<T> temp = this->transpose();

        // update the vector of elements:
        math::mtcopy<T>(temp.m_elems, this->m_elems);

        // and swap matrix's dimensions:
        std::swap( this->m_rows, this->m_cols );
    }
    
    return *this;
}


/**
 * Conjugation of  all matrix's elements.
 * For complex types, elements' imaginary parts are reversed
 * their signs. For all other types, a copy of *this
 * is returned.
 * 
 * @return *this conjugated
 * 
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::conj() const throw (math::MatrixException)
{
    math::MatrixGeneric<T> retVal(*this);
    math::__matrixprivate::__matconj(*this, retVal);

    return retVal;
}


/**
 * Removes the specified row number from the matrix.
 * It also decreases the number of rows.
 *
 * @param rowNr - the row number to remove
 *
 * @return reference to itself
 *
 * @throw MatrixException if attempting to remove the nonexistent row
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeRow(const size_t rowNr) throw (math::MatrixException)
{
    /*
     * Check of input arguments.
     * The matrix must contain at least two rows (as the updated matrix must still
     * contain at least one row), rowNr must be between 0 and rows-1.
     */
    if ( rowNr >= this->m_rows || this->m_rows <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    /*
     * Elements of the r^th row are located between r*m_cols and (r+1)*m_cols-1.
     * The first element of the next row is located at (r+1)*m_cols.
     * Row elements are contiguous and can be removed with one vector.erase() call.
     * It will also relocate remaining elements if applicable.
     */
    this->m_elems.erase(
                this->m_elems.begin() + rowNr * this->m_cols,
    		    this->m_elems.begin() + (rowNr+1) * this->m_cols );

    // Elements have been removed, update the number of rows
    --(this->m_rows);

    return *this;
}


/**
 * Removes the specified column number from the matrix.
 * It also decreases the number of columns.
 *
 * @param colNr - the column number to remove
 *
 * @return reference to itself
 *
 * @throw MatrixException if attempting to remove the nonexistent column
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeColumn(const size_t colNr) throw (math::MatrixException)
{
    /*
     * Checking of input arguments. The matrix must contain at least 2 columns
     * as the result must still contain at least one. colNr must be
     * between 0 and m_cols-1
     */
    if ( colNr >= this->m_cols || this->m_cols <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    /*
     * Elements of the column are not contiguous so their removal is a bit tricky.
     * It is best to remove them from the last (the highest row number) till
     * the first one (row=0). This way the position of the element to be removed
     * is (m_rows-i)*m_cols+colNr, m_cols is not updated yet. vector.erase() will
     * move remaining elements appropriately.
     *
     * Note: vector.erase() is by no means thread safe, so the for loop
     * should not be parallelized!
     */
    for ( size_t i=1; i<=this->m_rows; ++i )
    {
        this->m_elems.erase(
                this->m_elems.begin() + (this->m_rows - i) * this->m_cols + colNr );
    }

    // All required elements have been removed, now update the number of columns
    --(this->m_cols);

    return *this;
}


/**
  * Inserts a row in front of the rowNr^th row.
  * It also increases the number of rows.
  *
  * @param rowNr - a row will be inserted in front of this row. Valid values between 0 and rows
  * @param el - value to be assigned to all inserted elements (default: 0)
  *
  * @return reference to itself
  *
  * @throw MatrixException if invalid rowNr or if reallocation fails
  */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertRow(const size_t rowNr, const T& el) throw (math::MatrixException)
{
    // a valid rowNr value is between 0 and m_rows (incl.)
    if ( rowNr > this->m_rows )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // more elements than max. possible size of a vector.
    if ( this->m_rows > (this->m_elems.max_size() / this->m_cols - 1) )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        this->m_elems.reserve( (this->m_rows + 1) * this->m_cols );
        // a contiguous block of m_cols elements will be inserted
        // the position of rowNr*m_cols element.
        this->m_elems.insert(
                this->m_elems.begin() + rowNr * this->m_cols,
                this->m_cols, el );
    }
    catch ( const std::bad_alloc& ba )
    {
        // it is possible that reallocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion was successful, update the number of rows
    ++(this->m_rows);

    return *this;
}


/**
  * Inserts a column in front of the colNr^th column.
  * It also increases the number of columns.
  *
  * @param colNr - a column will be inserted in front of this column. Valid values between 0 and cols
  * @param el - value to be assigned to all inserted elements (default: 0)
  * 
  * @return reference to itself
  *
  * @throw MatrixException if invalid colNr or if reallocation fails
  */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertColumn(const size_t colNr, const T& el) throw (math::MatrixException)
{
    // A valid colNr is between 0 and m_cols (incl.)
    if ( colNr > this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // no more than MAX_SIZE_T elements
    if ( this->m_cols > (this->m_elems.max_size() / this->m_rows - 1) )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        /*
         * row elements will be inserted, however they are not contiguous.
         * First reserve enough memory for smoother reallocations
         */
        this->m_elems.reserve( this->m_rows * (this->m_cols + 1) );

        /*
         * Elements will be inserted step by step, with ascending row coordinate.
         * The position of each such element can be calculated as r*(m_cols+1)+colNr.
         *
         * Note: vector.insert() is by no means thread safe, so the for loop
         * should not be parallelized!
         */
        for ( size_t r = 0; r < this->m_rows; ++r )
        {
            this->m_elems.insert(
                    this->m_elems.begin() + r * (this->m_cols + 1) + colNr, el);
        }  // for r
    }  // try
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion successful, update the number of columns
    ++(this->m_cols);

    return *this;
}


/**
 * Swaps rows in the matrix.
 *
 * @param r1 - first row's number
 * @param r2 - second row's number
 *
 * @return a reference to itself
 *
 * @throw MatrixException if any input argument is out of range
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::swapRows(
        const size_t r1,
        const size_t r2
      ) throw(math::MatrixException)
{
    // Sanity check
    if ( r1>=this->m_rows || r2>=this->m_rows )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nothing to do if both parameters are equal
    if ( r1 == r2 )
    {
        return *this;
    }

    math::mtswap<T>(
        this->m_elems.begin() + r1 * this->m_cols,
        this->m_elems.begin() + (r1+1) * this->m_cols,
        this->m_elems.begin() + r2 * this->m_cols );


    return *this;
}


/**
 * Swaps columns in the matrix.
 *
 * @param c1 - first column's number
 * @param c2 - second column's number
 *
 * @return a reference to itself
 *
 * @throw MatrixException if any input argument is out of range
 */
template <class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::swapColumns(
        const size_t c1,
        const size_t c2
      ) throw(math::MatrixException)
{
    // Sanity check
    if ( c1>=this->m_cols || c2>=this->m_cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nothing to do if both parameters are equal
    if ( c1 == c2 )
    {
        return *this;
    }

    const size_t N = this->m_rows;

    // Coarse grained parallelization
    #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        for ( size_t r=istart; r<iend; ++r )
        {
            std::swap(
                this->at(r, c1),
                this->at(r, c2) );
        }
    }  // omp parallel

    return *this;
}


/**
 * Destructor
 */
template <class T>
math::MatrixGeneric<T>::~MatrixGeneric()
{
    // Vector's destructors would probably clean up this automatically.
    // Anyway let us clear the vector, just to be aware of allocated resources.
    this->m_elems.clear();

    // Other dynamically allocated memory (via malloc or new) should be freed here.
    // There are no other resources to release.
}



/**
 * Unary operator '+', returns a copy of the input argument 'f'.
 * 
 * @note Usage of this operator should be avoided
 * 
 * @param m - matrix to be copied
 * 
 * @return copy of 'm'
 * 
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator+(const math::MatrixGeneric<T>& m) throw(math::MatrixException)
{
    return m;
}


/**
 * Unary negation operator (-)
 *
 * @param m - matrix to be negated
 * 
 * @return -m
 *
 * @throw MatrixException if memory allocation fails
 */
template <class T>
math::MatrixGeneric<T> math::operator-(const math::MatrixGeneric<T>& m) throw(math::MatrixException)
{
    // There are no requirements about dimensions and no check of input arguments is necessary

    /*
     * Each element of the resulting matrix is a negated value of the element
     * at the same position:
     * N(r,c) = -m(r,c)
     */
    math::MatrixGeneric<T> temp(m.m_rows, m.m_cols);

    math::mtvectmult<T>(m.m_elems, static_cast<T>(-1), temp.m_elems);

    return temp;
}


/**
 * Addition operator (+) of two matrices.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param m1 - augend
 * @param m2 - addend
 *
 * @return m1 + m2
 *
 * @throw MatrixException if dimensions do not match
 */
template <class T>
math::MatrixGeneric<T> math::operator+(const math::MatrixGeneric<T>& m1, const math::MatrixGeneric<T>& m2) throw (math::MatrixException)
{
    // Check of dimensions. Numbers of rows and columns must match
    // otherwise addition is not possible
    if ( m1.m_rows != m2.m_rows || m1.m_cols != m2.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Each element of the sum matrix is a sum of elements at the same position
    // S(r,c) = m1(r,c) + m2(r,c)
    math::MatrixGeneric<T> temp(m1.m_rows, m2.m_cols);

    math::mtvectadd<T>(m1.m_elems, m2.m_elems, temp.m_elems, true);

    return temp;
}


/**
 * Subtraction operator (-) of two matrices.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param m1 - minuend
 * @param m2 - subtrahend
 *
 * @return m1 - m2
 *
 * @throw MatrixException if dimensions do not match
 */
template <class T>
math::MatrixGeneric<T> math::operator-(const math::MatrixGeneric<T>& m1, const math::MatrixGeneric<T>& m2) throw (math::MatrixException)
{
    // Check dimensions of both matrices. They must have the same number of rows and columns
    if ( m1.m_rows != m2.m_rows || m1.m_cols != m2.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Each element of the difference matrix is a difference of elements at the same position
    // D(r,c) = m1(r,c) - m2(r,c)
    math::MatrixGeneric<T> temp(m1.m_rows, m2.m_cols);

    math::mtvectadd<T>(m1.m_elems, m2.m_elems, temp.m_elems, false);

    return temp;
}


/**
 * Multiplication operator (*) of two matrices.
 * Number of columns of 'this' must be the same as number of rows of matrix,
 * otherwise multiplication is not possible. The resulting matrix will
 * have this.rows rows and matrix.cols columns.
 *
 * Note that matrix multiplication is not commutative (A*B != B*A).
 *
 * @param m1 - multiplicand
 * @param m2 - multiplier
 *
 * @return m1 * m2
 *
 * @throw MatrixException if dimensions do not match
 */
template <class T>
math::MatrixGeneric<T> math::operator*(const math::MatrixGeneric<T>& m1, const math::MatrixGeneric<T>& m2) throw (math::MatrixException)
{
    // Check if dimensions match (this.cols must be equal to matrix.rows)
    if ( m1.m_cols != m2.m_rows )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    /*
     * Multiplication modifies dimensions, so make sure the product will not
     * contain more elements than allowed by vector:
     */
    if ( m2.m_cols > m1.m_elems.max_size() / m1.m_rows )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    /*
     * If dimension of 'm1' is (m,n) and dimension of 'm2' is (n,o), the
     * dimension of the product will be (m,o).
     *
     *             n-1
     *            -----
     *            \
     *   P(r,c) =  >  ( m1(r,i) * m2(i,c) )
     *            /
     *            -----
     *             i=0
     *
     */

    math::MatrixGeneric<T> temp(m1.m_rows, m2.m_cols);

    #pragma omp parallel for collapse(2) \
                default(none) shared(m1, m2, temp)
    for ( size_t r=0; r<m1.m_rows; ++r )
    {
        for ( size_t c=0; c<m2.m_cols; ++c)
        {
            T sum = static_cast<T>(0);
            for ( size_t i=0; i<m1.m_cols; ++i )
            {
                sum += m1.m_elems.at(m1.__pos(r, i)) * m2.m_elems.at(m2.__pos(i, c));
            }

            temp.m_elems.at(temp.__pos(r, c)) = sum;
        }  // for c
    }  // for r

    return temp;
}


/**
 * Multiplication operator (*) for multiplication of a matrix and a scalar.
 * There are no restrictions about matrix dimensions.
 *
 * @param m - multiplicand (a matrix)
 * @param sc - multiplier (a scalar)
 *
 * @return m * sc
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator*(const math::MatrixGeneric<T>& m, const T& sc) throw (math::MatrixException)
{
    // The operation is possible on matrices with any dimension
    // so no need to check input arguments

    // Just multiply each element by the scalar:
    // P(r,c) = scalar * m(r,c)
    math::MatrixGeneric<T> retVal(m.m_rows, m.m_cols);

    math::mtvectmult<T>(m.m_elems, sc, retVal.m_elems);

    return retVal;
}


/**
 * Multiplication operator (*) of a scalar and a matrix.
 * This operation is commutative and does the same as operator*(scalar).
 *
 * @param sc - multiplicand (a scalar)
 * @param m - multiplier (a matrix)
 *
 * @return scalar * matrix
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator*(const T& sc, const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    MatrixGeneric<T> retVal = m * sc;
    return retVal;
}


/**
 * Adds each element of the matrix 'm' by the scalar 'sc'
 *
 * @param m - augend (a matrix)
 * @param sc - addend (a scalar)
 *
 * @return matrix 'res' where res(i,j) = m(i.j) + sc for each valid 'i' and 'j'
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator+(
            const math::MatrixGeneric<T>& m,
            const T& sc)
        throw (math::MatrixException)
{
    math::MatrixGeneric<T> retVal(m.m_rows, m.m_cols);

    math::mtvectscalaradd<T>(m.m_elems, sc, retVal.m_elems, true, true);

    return retVal;
}


/**
 * Adds each element of the matrix 'm' by the scalar 'sc'
 *
 * @param sc - addend (a scalar)
 * @param m - augend (a matrix)
 *
 * @return matrix 'res' where res(i,j) = sc + m(i.j) for each valid 'i' and 'j'
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator+(
            const T& sc,
            const math::MatrixGeneric<T>& m)
        throw (math::MatrixException)
{
    math::MatrixGeneric<T> retVal(m.m_rows, m.m_cols);

    math::mtvectscalaradd<T>(m.m_elems, sc, retVal.m_elems, true, false);

    return retVal;
}


/**
 * Subtracts each element of the matrix 'm' by the scalar 'sc'
 *
 * @param m - minuend (a matrix)
 * @param sc - subtrahend (a scalar)
 *
 * @return matrix 'res' where res(i,j) = m(i.j) - sc for each valid 'i' and 'j'
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator-(
            const math::MatrixGeneric<T>& m,
            const T& sc)
        throw (math::MatrixException)
{
    math::MatrixGeneric<T> retVal(m.m_rows, m.m_cols);

    math::mtvectscalaradd<T>(m.m_elems, sc, retVal.m_elems, false, true);

    return retVal;
}


/**
 * Subtracts the scalar 'sc' by each element of the matrix 'm'
 *
 * @param sc - minuend (a scalar)
 * @param m - subtrahend (a matrix)
 *
 * @return matrix 'res' where res(i,j) = sc - m(i.j) for each valid 'i' and 'j'
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator-(
            const T& sc,
            const math::MatrixGeneric<T>& m)
        throw (math::MatrixException)
{
    math::MatrixGeneric<T> retVal(m.m_rows, m.m_cols);

    math::mtvectscalaradd<T>(m.m_elems, sc, retVal.m_elems, false, false);

    return retVal;
}


/**
 * Divides each matrix's element by a scalar.
 *
 * @param m - dividend (a matrix)
 * @param sc - divisor (a scalar)
 *
 * @return m / sc
 *
 * @throw MatrixException if attempting to divide by 0 or allocation of memory fails
 */
template <class T>
math::MatrixGeneric<T> math::operator/(
            const math::MatrixGeneric<T>& m,
            const T& sc)
throw (math::MatrixException)
{
    // Division by zero is not permitted:
    if ( true == math::NumericUtil::isZero<T>(sc) )
    {
        throw math::MatrixException(math::MatrixException::FORBIDDEN);
    }

    const T f = static_cast<T>(1) / sc;
    math::MatrixGeneric<T> retVal(m.m_rows, m.m_cols);

    // Multiply all m's elements by 1/sc:
    math::mtvectmult(m.m_elems, f, retVal.m_elems);

    return retVal;
}


/**
 * Matrix element wise multiplication, equivalent to
 * Matlab's operator .*
 *
 * @param m1 - multiplicand matrix
 * @param m2 - multiplier matrix
 *
 * @return m1 .* m2
 *
 * @throw MatrixException if matrices' dimensions are not the same
 */
template <class T>
math::MatrixGeneric<T> math::matEwMult(
            const math::MatrixGeneric<T>& m1,
            const math::MatrixGeneric<T>& m2)
        throw (math::MatrixException)
{
    // Check of dimensions. Numbers of rows and columns must match
    // otherwise element wise multiplication is not possible
    if ( m1.m_rows != m2.m_rows || m1.m_cols != m2.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    MatrixGeneric<T> retVal(m1.m_rows, m1.m_cols);

    math::mtvectewmult<T>(m1.m_elems, m2.m_elems, retVal.m_elems, true);

    return retVal;
}


/**
 * Matrix element wise division, equivalent to
 * Matlab's operator ./
 *
 * @param m1 - dividend matrix
 * @param m2 - divisor matrix
 *
 * @return m1 ./ m2
 *
 * @throw MatrixException if matrices' dimensions are not the same or any divisor's element equals 0
 */
template <class T>
math::MatrixGeneric<T> math::matEwDiv(
            const math::MatrixGeneric<T>& m1,
            const math::MatrixGeneric<T>& m2)
        throw (math::MatrixException)
{
    // Check of dimensions. Numbers of rows and columns must match
    // otherwise element wise division is not possible
    if ( m1.m_rows != m2.m_rows || m1.m_cols != m2.m_cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    MatrixGeneric<T> retVal(m1.m_rows, m1.m_cols);

    const bool succ = math::mtvectewmult<T>(m1.m_elems, m2.m_elems, retVal.m_elems, false);

    if ( false == succ )
    {
        throw math::MatrixException(math::MatrixException::FORBIDDEN);
    }

    return retVal;
}


/*
 * Conjugation of all matrix elements.
 * The general version (for non complex types)
 * just copies the input matrix.
 *
 * @param m - matrix to conjugate
 * @param dest - a reference to a matrix to fill with elements of 'm' conjugated
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
void math::__matrixprivate::__matconj(
           const math::MatrixGeneric<T>& m,
           math::MatrixGeneric<T>& dest
        ) throw(math::MatrixException)
{
    dest = m;
}


/*
 * Overloading (partial "specialization") of _matconj
 * for complex numbers. It actually conjugates each element.
 *
 * @param m - matrix to conjugate
 * @param dest - a reference to a matrix to fill with elements of 'm' conjugated
 *
 * @throw MatrixException if allocation of memory fails
 */
template <class T>
void math::__matrixprivate::__matconj(
           const math::MatrixGeneric<std::complex<T> >& m,
           math::MatrixGeneric<std::complex<T> >& dest
        ) throw(math::MatrixException)
{
    /*
     * Specialization into another templated class implemented
     * as suggested here:
     * http://www.cplusplus.com/forum/general/68298/
     */
    dest = m;
    const size_t N = dest.nrRows() * dest.nrColumns();

    //Coarse grained parallelization
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(dest)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<std::complex<T> >::iterator it = dest.m_elems.begin() + istart;
        for ( size_t i=istart;
              i<iend && it!=dest.m_elems.end(); ++it, ++i )
        {
            *it = std::conj(*it);
        }
    }  // omp parallel

}
