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
#include "util/mtcopy.hpp"
#include "util/mtvectop.hpp"
#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"
#include "MatrixGeneric.hpp"


/**
 * Constructor.
 * Creates an instance of a matrix with the specified number of rows and columns.
 * Elements of the matrix are set to the default value (zero).
 *
 * @param rows    - number of rows (default: 1)
 * @param columns - number of columns (default: 1)
 *
 * @throw MatrixException when allocation of memory fails or in case of incorrect
 *        input parameters (both must be at least 1)
 */
template <class T>
math::MatrixGeneric<T>::MatrixGeneric(size_t rows, size_t columns) throw(math::MatrixException)
{
    // Matrix must contain at least 1 row and at least 1 column
    if ( rows < 1 || columns < 1 )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    // even theoretically the number of vector's elements is limited
    if ( columns > this->elems.max_size()/rows )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        this->rows = rows;
        this->cols = columns;
        // make sure, the vector will be empty
        this->elems.clear();
        // allocate memory for required number of elements, initialize each of them
        this->elems.resize(rows*columns, static_cast<T>(0));
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }
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
    _copyElems(orig);
}


// Copy elements from one matrix into another. Used at copy constructors,
// assignment operators etc.
template <class T>
void math::MatrixGeneric<T>::_copyElems(const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    // Check if the original matrix is OK (if its rows and cols are set correctly)
    if ( orig.rows * orig.cols != orig.elems.size() )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    try
    {
        this->rows = orig.rows;
        this->cols = orig.cols;
        math::mtcopy<T>(orig.elems, this->elems);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Check if the correct number of elements was copied
    if ( this->rows * this->cols != this->elems.size() )
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
    return this->rows;
}


/**
 * @return number of columns
 */
template <class T>
size_t math::MatrixGeneric<T>::nrColumns() const
{
    return this->cols;
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
T math::MatrixGeneric<T>::get(size_t row, size_t column) const throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= this->rows || column >= this->cols )
    {
        // At least one input parameter is out of range
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->elems.at(_pos(row, column));
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
T& math::MatrixGeneric<T>::at(size_t row, size_t column) throw (math::MatrixException)
{
    // Check if input arguments are within the matrix's range
    if ( row >= this->rows || column >= this->cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->elems.at(_pos(row, column));
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
const T& math::MatrixGeneric<T>::at(size_t row, size_t column) const throw (math::MatrixException)
{
    /*
     * Implementation is actually the same as implementation of another at()
     * with non-const signature. A single macro could be used for both
     * implementations, however both functions are short, simple and unlikely
     * to change often.
     */

    // Check if input arguments are within the matrix's range
    if ( row >= this->rows || column >= this->cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->elems.at(_pos(row, column));
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
T& math::MatrixGeneric<T>::operator()(size_t row, size_t column) throw(math::MatrixException)
{
    // Check if input arguments are within the matrix's range
    if ( row >= this->rows || column >= this->cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->elems.at(_pos(row, column));
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
const T& math::MatrixGeneric<T>::operator()(size_t row, size_t column) const throw(math::MatrixException)
{
    /*
     * Implementation is actually the same as implementation of another operator()
     * with non-const signature. A single macro could be used for both
     * implementations, however both functions are short, simple and unlikely
     * to change often.
     */

    // Check if input arguments are within the matrix's range
    if ( row >= this->rows || column >= this->cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return this->elems.at(_pos(row, column));
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
T& math::MatrixGeneric<T>::operator()(size_t idx) throw(math::MatrixException)
{
    // check if 'idx' is within elems' range
    if ( idx >= this->elems.size() )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    const size_t col = idx / this->rows;
    const size_t row = idx % this->rows;

    return this->elems.at(_pos(row, col));
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
const T& math::MatrixGeneric<T>::operator()(size_t idx) const throw(math::MatrixException)
{
    // check if 'idx' is within elems' range
    if ( idx >= this->elems.size() )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    const size_t col = idx / this->rows;
    const size_t row = idx % this->rows;

    return this->elems.at(_pos(row, col));
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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::set(size_t row, size_t column, const T& element) throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= this->rows || column >= this->cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    this->elems.at(_pos(row, column)) = element;

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

    for ( size_t r=0; r<(this->rows); ++r )
    {
        // display elements of the row r, separated by tabs
        for ( size_t c=0; c<(this->cols); ++c )
        {
            str << this->elems.at(_pos(r, c)) << "\t";
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

    this->_copyElems(orig);

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
    if ( this->rows != m.rows || this->cols != m.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // For a definition of matrix addition, see operator+
    math::mtvectadd<T>(this->elems, m.elems, this->elems, true);

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
    if ( this->rows != matrix.rows || this->cols != matrix.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // For a definition of matrix subtraction, see operator-
    math::mtvectadd<T>(this->elems, matrix.elems, this->elems, false);

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
    math::mtvectscalaradd<T>(this->elems, scalar, this->elems, true, true);

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
    math::mtvectscalaradd<T>(this->elems, scalar, this->elems, false, true);

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
        this->_copyElems(temp);
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
    math::mtvectmult<T>(this->elems, scalar, this->elems);

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
    math::mtvectmult(this->elems, f, this->elems);

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
    if ( this->rows != m.rows || this->cols != m.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    math::mtvectewmult<T>(this->elems, m.elems, this->elems, true);

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
    if ( this->rows != m.rows || this->cols != m.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const bool succ = math::mtvectewmult<T>(this->elems, m.elems, this->elems, false);

    if ( false == succ )
    {
        throw math::MatrixException(math::MatrixException::FORBIDDEN);
    }

    return *this;
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
    math::MatrixGeneric<T> retVal(this->cols, this->rows);

    const size_t& tcols = this->cols;

    const size_t N = this->rows * this->cols;

    // Coarse grained parallelism
    const std::vector<T>& els = this->elems;

    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(retVal, els, tcols)
    {
        const size_t thnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();
        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thnr;
        const size_t iend = std::min(istart + elems_per_thread, N);

        for ( size_t idx=istart; idx<iend; ++idx )
        {
            const size_t r = idx / tcols;
            const size_t c = idx % tcols;

            retVal.elems.at(retVal._pos(c, r)) = els.at(this->_pos(r, c));
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
     *
     * Note: as the method is not memory efficient at the moment, the function
     * is declared as virtual, allowing a more efficient implementation at
     * SqMatrixGeneric.
     */

    // If dimension of this is (n,m), dimension of its transposed matrix is (m,n)
    // T(r,c) = this(c,r)

    // Until a better algorithm is implemented
    // just use the general transpose method
    math::MatrixGeneric<T> temp = this->transpose();

    // update the vector of elements:
    math::mtcopy<T>(temp.elems, this->elems);

    // and swap matrix's dimensions:
    std::swap( this->rows, this->cols );

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
    math::MatrixGeneric<T> retVal;
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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeRow(size_t rowNr) throw (math::MatrixException)
{
    /*
     * Check of input arguments.
     * The matrix must contain at least two rows (as the updated matrix must still
     * contain at least one row), rowNr must be between 0 and rows-1.
     */
    if ( rowNr >= this->rows || this->rows <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    /*
     * Elements of the r^th row are located between r*cols and (r+1)*cols-1.
     * The first element of the next row is located at (r+1)*cols.
     * Row elements are contiguous and can be removed with one vector.erase() call.
     * It will also relocate remaining elements if applicable.
     */
    this->elems.erase(this->elems.begin()+rowNr*this->cols, this->elems.begin()+(rowNr+1)*this->cols);

    // Elements have been removed, update the number of rows
    --(this->rows);

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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeColumn(size_t colNr) throw (math::MatrixException)
{
    /*
     * Checking of input argumentss. The matrix must contain at least 2 columns
     * as the result must still contain at least one. colNr must be
     * between 0 and cols-1
     */
    if ( colNr >= this->cols || this->cols <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    /*
     * Elements of the column are not contiguous so their removal is a bit tricky.
     * It is best to remove them from the last (the highest row number) till
     * the first one (row=0). This way the position of the element to be removed
     * is (rows-i)*cols+colNr, cols is not updated yet. vector.erase() will
     * move remaining elements appropriately.
     *
     * Note: vector.erase() is by no means thread safe, so the for loop
     * should not be parallelized!
     */
    for ( size_t i=1; i<=this->rows; ++i )
    {
        this->elems.erase(this->elems.begin()+(this->rows-i)*this->cols+colNr);
    }

    // All required elements have been removed, now update the number of columns
    --(this->cols);

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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertRow(size_t rowNr, const T& el) throw (math::MatrixException)
{
    // a valid rowNr value is between 0 and rows (incl.)
    if ( rowNr > this->rows )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // more elements than max. possible size of a vector.
    if ( this->rows > (this->elems.max_size()/this->cols - 1) )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        this->elems.reserve( (this->rows+1) * this->cols );
        // a contiguous block of cols elements will be inserted
        // the position of rowNr*cols element.
        this->elems.insert(this->elems.begin()+rowNr*this->cols, this->cols, el);
    }
    catch ( const std::bad_alloc& ba )
    {
        // it is possible that reallocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion was successful, update the number of rows
    ++(this->rows);

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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertColumn(size_t colNr, const T& el) throw (math::MatrixException)
{
    // A valid colNr is between 0 and cols (incl.)
    if ( colNr > this->cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // no more than MAX_SIZE_T elements
    if ( this->cols > (this->elems.max_size()/this->rows - 1) )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        /*
         * row elements will be inserted, however they are not contiguous.
         * First reserve enough memory for smoother reallocations
         */
        this->elems.reserve( this->rows * (this->cols+1) );

        /*
         * Elements will be inserted step by step, with ascending row coordinate.
         * The position of each such element can be calculated as r*(cols+1)+colNr.
         *
         * Note: vector.insert() is by no means thread safe, so the for loop
         * should not be parallelized!
         */
        for ( size_t r = 0; r < this->rows; ++r )
        {
            this->elems.insert(this->elems.begin()+r*(this->cols+1)+colNr, el);
        }  // for r
    }  // try
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion successful, update the number of columns
    ++(this->cols);

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
    this->elems.clear();

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
    math::MatrixGeneric<T> temp(m.rows, m.cols);

    math::mtvectmult<T>(m.elems, static_cast<T>(-1), temp.elems);

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
    if ( m1.rows != m2.rows || m1.cols != m2.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Each element of the sum matrix is a sum of elements at the same position
    // S(r,c) = m1(r,c) + m2(r,c)
    math::MatrixGeneric<T> temp(m1.rows, m2.cols);

    math::mtvectadd<T>(m1.elems, m2.elems, temp.elems, true);

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
    if ( m1.rows != m2.rows || m1.cols != m2.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Each element of the difference matrix is a difference of elements at the same position
    // D(r,c) = m1(r,c) - m2(r,c)
    math::MatrixGeneric<T> temp(m1.rows, m2.cols);

    math::mtvectadd<T>(m1.elems, m2.elems, temp.elems, false);

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
    if ( m1.cols != m2.rows )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    /*
     * Multiplication modifies dimensions, so make sure the product will not
     * contain more elements than allowed by vector:
     */
    if ( m2.cols > m1.elems.max_size() / m1.rows )
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

    math::MatrixGeneric<T> temp(m1.rows, m2.cols);

    #pragma omp parallel for collapse(2) \
                default(none) shared(m1, m2, temp)
    for ( size_t r=0; r<m1.rows; ++r )
    {
        for ( size_t c=0; c<m2.cols; ++c)
        {
            T sum = static_cast<T>(0);
            for ( size_t i=0; i<m1.cols; ++i )
            {
                sum += m1.elems.at(m1._pos(r, i)) * m2.elems.at(m2._pos(i, c));
            }

            temp.elems.at(temp._pos(r, c)) = sum;
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
    math::MatrixGeneric<T> retVal(m.rows, m.cols);

    math::mtvectmult<T>(m.elems, sc, retVal.elems);

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
    math::MatrixGeneric<T> retVal(m.rows, m.cols);

    math::mtvectscalaradd<T>(m.elems, sc, retVal.elems, true, true);

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
    math::MatrixGeneric<T> retVal(m.rows, m.cols);

    math::mtvectscalaradd<T>(m.elems, sc, retVal.elems, true, false);

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
    math::MatrixGeneric<T> retVal(m.rows, m.cols);

    math::mtvectscalaradd<T>(m.elems, sc, retVal.elems, false, true);

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
    math::MatrixGeneric<T> retVal(m.rows, m.cols);

    math::mtvectscalaradd<T>(m.elems, sc, retVal.elems, false, false);

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
    math::MatrixGeneric<T> retVal(m.rows, m.cols);

    // Multiply all m's elements by 1/sc:
    math::mtvectmult(m.elems, f, retVal.elems);

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
    if ( m1.rows != m2.rows || m1.cols != m2.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    MatrixGeneric<T> retVal(m1.rows, m1.cols);

    math::mtvectewmult<T>(m1.elems, m2.elems, retVal.elems, true);

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
    if ( m1.rows != m2.rows || m1.cols != m2.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    MatrixGeneric<T> retVal(m1.rows, m1.cols);

    const bool succ = math::mtvectewmult<T>(m1.elems, m2.elems, retVal.elems, false);

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
 *
 * @return m "conjugated"
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
 *
 * @return m conjugated
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
        const size_t thnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();
        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thnr;
        const size_t iend = std::min(istart + elems_per_thread, N);

        typename std::vector<std::complex<T> >::iterator it = dest.elems.begin() + istart;
        for ( size_t i=istart;
              i<iend && it!=dest.elems.end(); ++it, ++i )
        {
            *it = std::conj(*it);
        }
    }  // omp parallel

}
