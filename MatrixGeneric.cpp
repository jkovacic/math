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
#include <limits>
#include <vector>

// Deliberately there is no #include "MatrixGeneric.h" !
#include "MatrixException.h"
#include "NumericUtil.h"


// Matrix's element of the row r and column c is accessed as "r*cols+c"
// (cols is a property of the class). This expression would be used
// frequently inside this file. To improve readability, to simplify
// code reviews and to eliminate possibilities of typing errors, this
// macro has been defined. It should only be used to access elements of
// *this. To access elements of other matrices (this occurs less frequently),
// the appropriate "version" of the above expression must be used
#define POS(r,c)    ( (r) * cols + (c) )

// A zero constant has already been defined in the class NumericUtil.
// It can only be accessed as math::NumericUtil<T>::ZERO
// As this notation is a bit long, this convenience macro is defined:
#define ZERO math::NumericUtil<T>::ZERO


/**
 * Constructor.
 * Creates an instance of a matrix with the specified number of rows and columns.
 *
 * @param rows    - number of rows (default: 1)
 * @param columns - number of columns (default: 1)
 *
 * @throw MatrixException when allocation ofmemory fails or in case of incorrect
 *        input parameters (both must be at least 1)
 */
template<class T>
math::MatrixGeneric<T>::MatrixGeneric(unsigned int rows, unsigned int columns) throw(math::MatrixException)
{
    // Matrix must contain at least 1 row and at least 1 column
    if ( rows < 1 || columns < 1 )
    {
        throw math::MatrixException(MatrixException::INVALID_DIMENSION);
    }

    // Max. INT_MAX (typically cca. 2e+9) elements will be allowed to prevent
    // problems with counter overflows. This is much much more than required
    // in most real life applications.
    const unsigned long int matSize = rows * columns;
    if ( matSize > static_cast<unsigned long int>(std::numeric_limits<int>::max() ) )
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
    catch ( std::bad_alloc &ba )
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
    catch (std::bad_alloc &ba )
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
unsigned int math::MatrixGeneric<T>::nrRows() const
{
    return this->rows;
}

/**
 * @return number of columns
 */
template<class T>
unsigned int math::MatrixGeneric<T>::nrColumns() const
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
T math::MatrixGeneric<T>::get(unsigned int row, unsigned int column) const throw (math::MatrixException)
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
        retVal = elems.at(POS(row, column));
    }
    catch ( std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return retVal;
}  // MatrixGeneric::get

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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::set(unsigned int row, unsigned int column, const T& element) throw (math::MatrixException)
{
    // Check of input parameters
    if ( row >= rows || column >= cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Attempt to modify the element
    try
    {
        elems.at(POS(row, column)) = element;
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
        unsigned int r;
        unsigned int c;

        // Elements of each row are separated by tabs.
        // This does not guarantee that elements of the same column will be
        // displayed under each other and in case of a long row (its output
        // is longer than terminal's line width), the overall output may look rather
        // confusing. However, this is more or less "just" an auxiliary function, mainly used for
        // testing purposes, and the effort was focused to other functionalities.
        // Anyway, it would be nice to improve it in future.
        for ( r=0; r<rows; r++ )
        {
            // display elements of the row r, separated by tabs
            for ( c=0; c<cols; c++ )
            {
                str << elems.at(POS(r, c)) << "\t";
            }
            // the line must be terminated by a newline
            str << std::endl;
        } // for r
    } // try
    catch ( std::out_of_range& oor )
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
        unsigned int i;
        for ( i=0; i<rows*cols; i++ )
        {
            temp.elems.at(i) = elems.at(i) + matrix.elems.at(i);
        }
    }  // try
    catch ( std::out_of_range& oor )
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
        unsigned int i;
        for ( i=0; i<rows*cols; i++ )
        {
            elems.at(i) += m.elems.at(i);
        }
    }  // try
    catch ( std::out_of_range& oor )
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
        unsigned int i;
        for ( i=0; i<rows*cols; i++ )
        {
            temp.elems.at(i) = elems.at(i) - m.elems.at(i);
        }
    }  // try
    catch ( std::out_of_range& oor )
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
        unsigned int i;
        for ( i=0; i<rows*cols; i++ )
        {
            elems.at(i) -= matrix.elems.at(i);
        }
    }  // try
    catch ( std::out_of_range& oor )
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
        unsigned int i;

        for ( i=0; i<rows*cols; i++ )
        {
            temp.elems.at(i) = -elems.at(i);
        }
    }  // try
    catch ( std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return temp;
}

/**
 * Multiplication operator (*) of two matrices.
 * Number of coumns of this must be the same as number of rows of matrix,
 * otherwise multiplication is not possible. The resulting matrix will
 * have this.rows rows and matrix.cols columns.
 *
 * Note that matrix multiplication is not comutative (A*B != B*A).
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
    // no more than INT_MAX elements
    const unsigned long int matSize = this->rows * matrix.cols;
    if ( matSize > (unsigned long int) std::numeric_limits<int>::max() )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    // if dimension of this is (m,n) and dimension of matrix is (n,o), the
    // dimension of the product will be (m,o).
    // P(r,c) = sum( i=0, n, this(r,i)*matrix(i,c) )
    math::MatrixGeneric<T> temp(rows, matrix.cols);

    try
    {
        unsigned int r;
        unsigned int c;

        for ( r=0; r<rows; r++ )
        {
            for ( c=0; c<matrix.cols; c++ )
            {
                T sum = ZERO;
                unsigned int i;
                for ( i=0; i<cols; i++ )
                {
                    sum += elems.at(POS(r, i))*matrix.elems.at(i*matrix.cols + c);
                }

                temp.elems.at(r*matrix.cols + c) = sum;
            }  // for c
        }  // for r
    }  // try
    catch ( std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return temp;
}

/**
 * Multiplication operator (*) for multiplication of a matrix and scalar.
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
    // so no need to check input parameters

    // Just multiply each element by the scalar:
    // P(r,c) = scalar * this(r,c)
    math::MatrixGeneric<T> retVal(rows, cols);

    try
    {
        unsigned int i;

        for ( i=0; i<rows*cols; i++ )
        {
            retVal.elems.at(i) = elems.at(i) * scalar;
        }
    }  //try
    catch ( std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    return retVal;
}
 /**
  * Multiplication operator (*) of a scalar and a matrix.
  * This operation is comutative and does the same as operator*(scalar).
  * Since the first operand is not a matrix, it must be implemented as
  * a friend function.
  *
  * @param scalar
  * @param matrix
  *
  * @return scalar * matrix
  *
  * @ŧhrow MatrixException if matrix does not contain enough elements
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
        unsigned int r;
        unsigned int c;

        for ( r=0; r<rows; r++ )
        {
            for ( c=0; c<cols; c++ )
            {
                // and swap their "coordinates"
                retVal.elems.at(c*rows + r) = elems.at(POS(r, c));
            }
        }
    }  // try
    catch ( std::out_of_range& oor )
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
    catch (std::bad_alloc &ba)
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    unsigned int r;
    unsigned int c;
    for ( r=0; r<rows; r++ )
    {
        for ( c=0; c<cols; c++ )
        {
            tempElems.at(c*rows + r) = this->elems.at(POS(r, c));
        }  // for c
    }  // for r

    // update the vector of elements:
    this->elems = tempElems;

    // and swap matrix's dimensions:
    c = this->cols;
    this->cols = this->rows;
    this->rows = c;

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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeRow(unsigned int rowNr) throw (math::MatrixException)
{
    // Check of input parameters.
    // The matrix must contain at least two rows (as the updated matrix must still
    // contain at least one row), rowNr must be between 0 and rows-1.
    if ( rowNr >= rows || rows <= 1 )
    {
        throw math::MatrixException(MatrixException::OUT_OF_RANGE);
    }

    // Elements od the r^th row are located between r*cols and (r+1)*cols-1.
    // The first element of the next row is located at (r+1)*cols.
    // Row elements are contiguous and can be removed with one vector.erase() call.
    // It will also relocate remainng elements if applicable.
    elems.erase(elems.begin()+rowNr*cols, elems.begin()+(rowNr+1)*cols);

    // Elements have been removed, update the number of rows
    rows--;

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
math::MatrixGeneric<T>& math::MatrixGeneric<T>::removeColumn(unsigned int colNr) throw (math::MatrixException)
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
    unsigned int i;
    for ( i=1; i<=rows; i++ )
    {
        elems.erase(elems.begin()+(rows-i)*cols+colNr);
    }

    // All required elements have been removed, now update the number of columns
    cols--;

    return *this;
}

/**
  * Inserts a row in front of the rowNr^th row.
  * It also increases the number of rows.
  * Elements of the new row will be set to a default value (typically 0).
  *
  * @param rowNr - a row will be inserted in front of this row. Valid values between 0 and rows
  *
  * @return reference to itself
  *
  * @throw MatrixException if invalid rowNr or if reallocation fails
  */
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertRow(unsigned int rowNr) throw (math::MatrixException)
{
    // a valid rowNr value is between 0 and rows (incl.)
    if ( rowNr > rows )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // no more than INT_MAX elements
    const unsigned long int matSize = (rows + 1) * cols;
    if ( matSize > (unsigned long int) std::numeric_limits<int>::max() )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        elems.reserve( (rows+1)*cols );
        // a contiguous block of cols elements will be inserted
        // the position of rowNr*cols element.
        elems.insert(elems.begin()+rowNr*cols, cols, ZERO);
    }
    catch ( std::bad_alloc& ba )
    {
        // it is possible that reallocation failed
        throw math::MatrixException(math::MatrixException::OUT_OF_MEMORY);
    }

    // Insertion was successful, update the number of rows
    rows++;

    return *this;
}

/**
  * Inserts a column in front of the colNr^th column.
  * It also increases the number of columns.
  * Elements of the new column will be set to a default value (typically 0).
  *
  * @param colNr - a column will be inserted in front of this column. Valid values between 0 and cols
  *
  * @return reference to itself
  *
  * @throw MatrixException if invalid colNr or if reallocation fails
  */
template<class T>
math::MatrixGeneric<T>& math::MatrixGeneric<T>::insertColumn(unsigned int colNr) throw (math::MatrixException)
{
    // A valid colNr is between 0 and cols (incl.)
    if ( colNr > cols )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Nr. of elements will increase, so make sure the product will contain
    // no more than INT_MAX elements
    const unsigned long int matSize = rows * (cols + 1);
    if ( matSize > (unsigned long int) std::numeric_limits<int>::max() )
    {
        throw math::MatrixException(math::MatrixException::TOO_LARGE);
    }

    try
    {
        // rows elements will be inserted, however they are not contiguous.
        // First reserve enough memory for smoother reallocations
        elems.reserve( rows*(cols+1) );

        // Elements will be inserted step by step, with ascending row cordinate.
        // The position of each such element can be calculated as r*(cols+1)+colNr.
        unsigned int r;

        for ( r = 0; r < rows; r++ )
        {
            elems.insert(elems.begin()+r*(cols+1)+colNr, ZERO);
        }  // for r
    }  // try
    catch ( std::bad_alloc& ba )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

    // Insertion successful, update the number of columns
    cols++;

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

// These macros were defined especially for this file. To prevent any possible
// conflicts, they will be #undef'ed
#undef POS
#undef ZERO
