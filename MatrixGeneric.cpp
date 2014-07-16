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
* @file MatrixGeneric.cpp
*
* Implementation of the class MatriGeneric.
*
* As the class is templated, this file must not be compiled.
* Instead it must be included after the class declaration in the .h file
*
* @author Jernej Kovacic
*/

#include <new>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#    include <omp.h>
#endif

// Deliberately there is no #include "MatrixGeneric.h" !
#include "MatrixException.h"
#include "NumericUtil.h"


// A zero constant has already been defined in the class NumericUtil.
// It can only be accessed as math::NumericUtil<T>::ZERO
// As this notation is a bit long, this convenience macro is defined:
#define ZERO math::NumericUtil<T>::ZERO


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
template<class T>
math::MatrixGeneric<T>::MatrixGeneric(size_t rows, size_t columns) throw(math::MatrixException)
{
    // Matrix must contain at least 1 row and at least 1 column
    if ( rows < 1 || columns < 1 )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    // even theoretically the number of vector's elements is limited
    if ( columns > elems.max_size()/rows )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        this->rows = rows;
        this->cols = columns;
        // make sure, the vector will be empty
        elems.clear();
        // allocate memory for required number of elements, initialize each of them
        elems.resize(rows*columns, ZERO);
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }
} // MatrixGeneric::MatrixGeneric(int, int)

/**
 * Copy constructor.
 * Creates an instance of a matrix with the same dimensions as orig
 * and copies its elements.
 *
 * @param orig - original matrix to be copied into this one
 *
 * @throw MatrixException when allocation of memory fails
 */
template<class T>
math::MatrixGeneric<T>::MatrixGeneric(const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    copyElems(orig);
}

/*
 * A utility function that returns the position of element's "coordinates"
 * within the matrix's internal vector (r*cols+c). As this functionality
 * is used often, the purpose of this function is to define it only once
 * and to eliminate possibilities of typing errors.
 *
 * As the function is simple (short) and called frequently, it is declared
 * as an inline function to reduce overhead.
 *
 * Note: no checks are performed inside the function. The results of the
 * function are usually passed directly to std::vector.at() which throws
 * an exception if 'pos' is out of the vector's range.
 */
template <class T>
inline size_t math::MatrixGeneric<T>::pos(size_t row, size_t column) const
{
    return ( row * this->cols + column );
}

// Copy elements from one matrix into another. Used at copy constructors,
// assignment operators etc.
template <class T>
void math::MatrixGeneric<T>::copyElems(const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    // Check if the original matrix is OK (if its rows and cols are set correctly)
    if ( orig.rows * orig.cols != orig.elems.size() )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    // release elements
    elems.clear();

    try
    {
        this->rows = orig.rows;
        this->cols = orig.cols;
        // STL vector's assignment operator (=) will allocate the appropriate
        // size to elems and copy all its elements (instantiate them with copy
        // constructors, if applicable)
        elems = orig.elems;
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
}  // MatrixGeneric::copyElems(const MatrixGeneric&)

/**
 * @return number of rows
 */
template<class T>
size_t math::MatrixGeneric<T>::nrRows() const
{
    return this->rows;
}

/**
 * @return number of columns
 */
template<class T>
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
template<class T>
T math::MatrixGeneric<T>::get(size_t row, size_t column) const throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= rows || column >= cols )
    {
        // At least one input parameter is out of range
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Access the element at the requested location
    T retVal;
    try
    {
        retVal = elems.at(pos(row, column));
    }
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return retVal;
}  // MatrixGeneric::get

/**
 * Returns a read-write reference to the desired element od the matrix.
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
template<class T>
T& math::MatrixGeneric<T>::at(size_t row, size_t column) throw (math::MatrixException)
{
    // Check if input parameters are within the matrix's range
    if ( row >= rows || column >= cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return elems.at(pos(row, column));
}

/**
 * Returns a read-only (const) reference to the desired element od the matrix.
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
template<class T>
const T& math::MatrixGeneric<T>::at(size_t row, size_t column) const throw (math::MatrixException)
{
    /*
        Implementation is actually the same as implementation of another at()
        with non-const signature. A single macro could be used for both
        implementations, however both functions are short, simple and unlikely
        to change often.
     */

    // Check if input parameters are within the matrix's range
    if ( row >= rows || column >= cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return elems.at(pos(row, column));
}

/**
 * Assigns value at the requested location.
 * This is the only allowed method to modify values of matrix's elements
 * (except of class's internal functions)
 *
 * @param row     row of the element to modify
 * @param column  column of the element to modify
 * @param element value to be assigned at the requested location
 *
 * @return reference to itself
 *
 * @throw MatrixException if input parameters are out of range
 */
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::set(size_t row, size_t column, const T& element) throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= rows || column >= cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Attempt to modify the element
    try
    {
        elems.at(pos(row, column)) = element;
    }
    catch ( std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return *this;
}  // MatrixGeneric::set

/**
 * Display the matrix to stdout
 *
 * @param str - output stream where the matrix will be displayed (default: cout)
 *
 * @throw MatrixException if attempting to access an element out of allocated range
 */
template<class T>
void math::MatrixGeneric<T>::display(std::ostream& str) const throw (math::MatrixException)
{
    try
    {
        // Elements of each row are separated by tabs.
        // This does not guarantee that elements of the same column will be
        // displayed under each other and in case of a long row (its output
        // is longer than terminal's line width), the overall output may look rather
        // confusing. However, this is more or less "just" an auxiliary function, mainly used for
        // testing purposes, and the effort was focused to other functionalities.
        // Anyway, it would be nice to improve it in future.
        for ( size_t r=0; r<rows; ++r )
        {
            // display elements of the row r, separated by tabs
            for ( size_t c=0; c<cols; ++c )
            {
                str << elems.at(pos(r, c)) << "\t";
            }
            // the line must be terminated by a newline
            str << std::endl;
        } // for r
    } // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }
}  // MatrixGeneric::display

/**
 * Assignment operator (=)
 *
 * @param orig - a matrix to be copied into this one
 *
 * @return reference to this
 *
 * @throw MatrixException if memory allocation fails
 */
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator= (const math::MatrixGeneric<T>& orig) throw (math::MatrixException)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &orig )
    {
        return *this;
    }

    copyElems(orig);

    return *this;
}

/**
 * Addition operator (+) of two matrices.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param matrix to be added to this one
 *
 * @return *this + matrix
 *
 * @throw MatrixException if dimensions do not match
 */
template<class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::operator+ (const math::MatrixGeneric<T>& matrix) const throw (math::MatrixException)
{
    // Check of dimensions. Numbers of rows and columns must match
    // otherwise addition is not possible
    if ( rows != matrix.rows || cols != matrix.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Each element of the sum matrix is a sum of elements at the same position
    // S(r,c) = this(r,c) + matrix(r,c)
    math::MatrixGeneric<T> temp(rows, cols);

    try
    {
        // Matrices have the same number of elements, just traverse
        // them linearly and perform addition of elements at the same position
        for ( size_t i=0; i<rows*cols; ++i )
        {
            temp.elems.at(i) = elems.at(i) + matrix.elems.at(i);
        }
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return temp;
}

/**
 * Addition operator (+=) that adds a matrix to this and assigns the sum to itself.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param matrix to be added to this one
 *
 * @return reference to this
 *
 * @throw MatrixException if dimensions do not match
 */
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator+= (const math::MatrixGeneric<T>& m) throw(math::MatrixException)
{
    // Check if dimensions of both matrices match
    if ( rows != m.rows || cols != m.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // For a definition of matrix addition, see operator+
    try
    {
        for ( size_t i=0; i<rows*cols; ++i )
        {
            elems.at(i) += m.elems.at(i);
        }
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return *this;
}

/**
 * Subtraction operator (-) of two matrices.
 * Both matrices must have the same dimension (equal number of rows and columns)
 *
 * @param matrix to be subtracted from this one
 *
 * @return *this - matrix
 *
 * @throw MatrixException if dimensions do not match
 */
template<class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::operator- (const math::MatrixGeneric<T>& m) const throw (math::MatrixException)
{
    // Check dimensions of both matrices. They must have the same number of rows and columns
    if ( rows != m.rows || cols != m.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Each element of the difference matrix is a difference of elements at the same position
    // D(r,c) = this(r,c) - matrix(r,c)
    math::MatrixGeneric<T> temp(rows, cols);

    try
    {
        for ( size_t i=0; i<rows*cols; ++i )
        {
            temp.elems.at(i) = elems.at(i) - m.elems.at(i);
        }
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return temp;
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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator-= (const math::MatrixGeneric<T>& matrix) throw(math::MatrixException)
{
    // Check if dimensions of both matrices match
    if ( rows != matrix.rows || cols != matrix.cols )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // For a definition of matrix subtraction, see operator-

    try
    {
        for ( size_t i=0; i<rows*cols; ++i )
        {
            elems.at(i) -= matrix.elems.at(i);
        }
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return *this;
}

/**
 * Unary negation operator (-)
 *
 * @return -this
 *
 * @throw MatrixException if memory allocation fails
 */
template<class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::operator-() const throw(math::MatrixException)
{
    // There are no requirements about dimensions and no input check is necessary

    // Each element of the resulting matrix is a negated value of the element
    // at the same position:
    // N(r,c) = -this(r,c)
    math::MatrixGeneric<T> temp(rows, cols);
    try
    {
        for ( size_t i=0; i<rows*cols; ++i )
        {
            temp.elems.at(i) = -elems.at(i);
        }
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

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
 * @param matrix
 *
 * @return this * matrix
 *
 * @throw MatrixException if dimensions do not match
 */
template<class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::operator* (const math::MatrixGeneric<T>& matrix) const throw (math::MatrixException)
{
    // Check if dimensions match (this.cols must be equal to matrix.rows)
    if ( cols != matrix.rows )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // Multiplication modifies dimensions, so make sure the product will contain
    // more elements than allowed by vector:
    if ( matrix.cols > elems.max_size()/this->rows )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    // if dimension of this is (m,n) and dimension of matrix is (n,o), the
    // dimension of the product will be (m,o).
    // P(r,c) = sum( i=0, n, this(r,i)*matrix(i,c) )
    math::MatrixGeneric<T> temp(rows, matrix.cols);

    try
    {
        for ( size_t r=0; r<rows; ++r )
        {
            for ( size_t c=0; c<matrix.cols; ++c )
            {
                T sum = ZERO;
                for ( size_t i=0; i<cols; ++i )
                {
                    sum += elems.at(pos(r, i)) * matrix.elems.at(matrix.pos(i, c));
                }

                temp.elems.at(matrix.pos(r, c)) = sum;
            }  // for c
        }  // for r
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return temp;
}

/**
 * Subtraction operator (*=) that multiplies a matrix by this one and assigns
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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator*= (const math::MatrixGeneric<T>& m) throw (math::MatrixException)
{
    // for a definition of matrix multiplication, see operator*

    // operator* will perform checking of numbers of rows/columns etc.
    // and throw an appropriate exception if necessary
    math::MatrixGeneric<T> temp = *this * m;

    try
    {
        this->rows = temp.rows;   // not really necessary to do this...
        this->cols = temp.cols;
        this->elems.clear();
        this->elems = temp.elems;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    return *this;
}

/**
 * Multiplication operator (*) for multiplication of a matrix and a scalar.
 * There are no restrictions about matrix dimensions.
 *
 * @param scalar
 *
 * @return this * scalar
 *
 * @throw MatrixException if matrix does not contain enough elements
 */
template<class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::operator* (const T& scalar) const throw (math::MatrixException)
{
    // The operation is possible on matrices with any dimension
    // so no need to check input arguments

    // Just multiply each element by the scalar:
    // P(r,c) = scalar * this(r,c)
    math::MatrixGeneric<T> retVal(rows, cols);

    try
    {
        const size_t N = rows*cols;
        for ( size_t i=0; i<N; ++i )
        {
            retVal.elems.at(i) = elems.at(i) * scalar;
        }
    }  //try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return retVal;
}

/**
 * Multiplication operator (*=) that multiplies a matrix by a scalar
 * and assigns the product to itself.
 *
 * @param scalar
 *
 * @return reference to itself
 */
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::operator*=(const T& scalar)
{
    const size_t N = rows*cols;

    // Multiply each element by the 'scalar'
    for ( size_t i=0; i<N; ++i )
    {
        elems.at(i) *= scalar;
    }

    return *this;
}

 /**
  * Multiplication operator (*) of a scalar and a matrix.
  * This operation is commutative and does the same as operator*(scalar).
  * Since the first operand is not a matrix, it must be implemented as
  * a friend function.
  *
  * @param scalar
  * @param matrix
  *
  * @return scalar * matrix
  *
  * @throw MatrixException if matrix does not contain enough elements
  */
template<class T>
math::MatrixGeneric<T> math::operator* (const T& scalar, const math::MatrixGeneric<T>& matrix) throw (math::MatrixException)
{
    MatrixGeneric<T> retVal = matrix * scalar;
    return retVal;
}

/**
 * Matrix transpose operation
 *
 * @return this^T
 *
 * @throw MatrixException if matrix does not contain enough elements
 */
template<class T>
math::MatrixGeneric<T> math::MatrixGeneric<T>::transpose() const throw (math::MatrixException)
{
    // If dimension of this is (n,m), dimension of its transposed matrix is (m,n)
    // T(r,c) = this(c,r)

    // Create an instance with swapped dimensions
    math::MatrixGeneric<T> retVal(cols, rows);

    try
    {
        // "collect" all elements of this
        for ( size_t r=0; r<rows; ++r )
        {
            for ( size_t c=0; c<cols; ++c )
            {
                // and swap their "coordinates"
                retVal.elems.at(retVal.pos(c, r)) = elems.at(pos(r, c));
            }
        }
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::transposed() throw (math::MatrixException)
{
    // TODO: find a memory efficient method for a general matrix!!
    // Note: as the method is not memory efficient at the moment, the function
    // is declared as virtual, allowing a more efficient implementation at
    // SqMatrixGeneric.

    // If dimension of this is (n,m), dimension of its transposed matrix is (m,n)
    // T(r,c) = this(c,r)

    std::vector<T> tempElems;
    // reserve enough space for the temporary vector:
    try
    {
        tempElems.resize(rows * cols);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    for ( size_t r=0; r<rows; ++r )
    {
        for ( size_t c=0; c<cols; ++c )
        {
            tempElems.at(c*rows + r) = this->elems.at(pos(r, c));
        }  // for c
    }  // for r

    // update the vector of elements:
    this->elems = tempElems;

    // and swap matrix's dimensions:
    size_t c = this->cols;
    this->cols = this->rows;
    this->rows = c;

    // tempElems not needed anymore, clean it
    tempElems.clear();

    return *this;
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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeRow(size_t rowNr) throw (math::MatrixException)
{
    // Check of input parameters.
    // The matrix must contain at least two rows (as the updated matrix must still
    // contain at least one row), rowNr must be between 0 and rows-1.
    if ( rowNr >= rows || rows <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    // Elements of the r^th row are located between r*cols and (r+1)*cols-1.
    // The first element of the next row is located at (r+1)*cols.
    // Row elements are contiguous and can be removed with one vector.erase() call.
    // It will also relocate remaining elements if applicable.
    elems.erase(elems.begin()+rowNr*cols, elems.begin()+(rowNr+1)*cols);

    // Elements have been removed, update the number of rows
    --rows;

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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeColumn(size_t colNr) throw (math::MatrixException)
{
    // Checking of input parameters. The matrix must contain at least 2 columns
    // as the result must still contain at least one. colNr must be
    // between 0 and cols-1
    if ( colNr >= cols || cols <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    // Elements of the column are not contiguous so their removal is a bit tricky.
    // It is best to remove them from the last (the highest row number) till
    // the first one (row=0). This way the position of the element to be removed
    // is (rows-i)*cols+colNr, cols is not updated yet. vector.erase() will
    // move remaining elements appropriately
    for ( size_t i=1; i<=rows; ++i )
    {
        elems.erase(elems.begin()+(rows-i)*cols+colNr);
    }

    // All required elements have been removed, now update the number of columns
    --cols;

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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertRow(size_t rowNr, const T& el) throw (math::MatrixException)
{
    // a valid rowNr value is between 0 and rows (incl.)
    if ( rowNr > rows )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // more elements than max. possible size of a vector.
    if ( rows > (elems.max_size()/cols - 1) )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        elems.reserve( (rows+1)*cols );
        // a contiguous block of cols elements will be inserted
        // the position of rowNr*cols element.
        elems.insert(elems.begin()+rowNr*cols, cols, el);
    }
    catch ( const std::bad_alloc& ba )
    {
        // it is possible that reallocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion was successful, update the number of rows
    ++rows;

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
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertColumn(size_t colNr, const T& el) throw (math::MatrixException)
{
    // A valid colNr is between 0 and cols (incl.)
    if ( colNr > cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // no more than MAX_SIZE_T elements
    if ( cols > (elems.max_size()/rows - 1) )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        // rows elements will be inserted, however they are not contiguous.
        // First reserve enough memory for smoother reallocations
        elems.reserve( rows*(cols+1) );

        // Elements will be inserted step by step, with ascending row coordinate.
        // The position of each such element can be calculated as r*(cols+1)+colNr.
        for ( size_t r = 0; r < rows; ++r )
        {
            elems.insert(elems.begin()+r*(cols+1)+colNr, el);
        }  // for r
    }  // try
    catch ( const std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion successful, update the number of columns
    ++cols;

    return *this;
}

/**
 * Destructor
 */
template<class T>
math::MatrixGeneric<T>::~MatrixGeneric()
{
    // Vector's destructors would probably clean up this automatically.
    // Anyway let us clear the vector, just to be aware of allocated resources.
    elems.clear();

    // Other dynamically allocated memory (via malloc or new) should be freed here.
    // There are no other resources to release.
}

// This macro was defined especially for this file. To prevent any possible
// conflicts, it will be #undef'ed
#undef ZERO
