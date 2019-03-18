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
 * An internal header file, it should not be included directly.
 * @headername{Matrix.h}
 *
 * Declaration and implementation of the class MatrixGeneric,
 * representing generic matrices.
 */

#ifndef _MATH_MATRIXGENERIC_HPP_
#define _MATH_MATRIXGENERIC_HPP_


#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <ostream>
#include <complex>
#include <algorithm>
#include <cmath>

#include "../settings/matrix_settings.h"
#include "util/NumericUtil.hpp"
#include "util/mtcopy.hpp"
#include "util/mtvectop.hpp"
#include "util/mtswap.hpp"
#include "util/PseudoFunctionGeneric.hpp"

#include "util/FillVectors.h"

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/MatrixException.hpp"


namespace math
{

/**
 * Templates are used to allow several types (T) of matrix elements.
 * Types must have implemented basic arithmetic operators (+, -, *, /),
 * otherwise build will fail (which is desired). In practice, types as float,
 * double, Rational, Complex make most sense. Integer types may be conditionally
 * acceptable (e.g. if the class is used to represent mathematical graphs, etc.)
 * but arithmetic operations (e.g. inversion of a square matrix) may return
 * incorrect results. This is true for unsigned types as well.
 */

// Forward declaration of the class is necessary...
template <class T> class MatrixGeneric;

// to declare the class's friend functions:
template <class T>
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m);

template <class T>
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m);

template <class T>
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);

template <class T>
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);

template <class T>
MatrixGeneric<T> operator*(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);

template <class T>
MatrixGeneric<T> operator*(const MatrixGeneric<T>& m, const T& sc);

template <class T>
MatrixGeneric<T> operator*(const T& sc, const MatrixGeneric<T>& m);

template <class T>
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m, const T& sc);

template <class T>
MatrixGeneric<T> operator+(const T& sc, const MatrixGeneric<T>& m);

template <class T>
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m, const T& sc);

template <class T>
MatrixGeneric<T> operator-(const T& sc, const MatrixGeneric<T>& m);

template <class T>
MatrixGeneric<T> operator/(const MatrixGeneric<T>& m, const T& sc);

template <class T>
MatrixGeneric<T> matEwMult(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);

template <class T>
MatrixGeneric<T> matEwDiv(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);



namespace __matrixprivate
{
    
template <class T>
void __matconj(MatrixGeneric<T>& m);

// overloaded "specialization" complex<T>:
template <class T>
void __matconj(MatrixGeneric<std::complex<T> >& m);

}  // namespace __matrixprivate



/**
 * @brief A class representing general matrices and their basic operators.
 */
template <class T>
class MatrixGeneric
{

private:

    std::size_t m_rows;      /// Number of rows
    std::size_t m_cols;      /// Number of columns

    /*
     * STL Vector has several advantages over arrays, allocated by new[],
     * e.g. elements can be accessed via at() which checks for range and throws
     * an exception when attempting to access sth. outside the allocated range
     * (e.g. in case of a typing error) which may result in a segmentation fault (crash).
     * Operations such as inserting or removing of elements are simplified as well.
     */
    std::vector<T> m_elems;   /// Elements of the matrix



    /*
     * Initializes a new matrix
     * 
     * @param rows    - number of rows
     * @param columns - number of columns
     *
     * @throw MatrixException when allocation of memory fails or in case of incorrect
     *        input arguments (both must be at least 1)
     */
    void __init(const std::size_t rows, const std::size_t cols)
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



    /*
     * Copy elements from one matrix into another. Used at copy constructors,
     * assignment operators etc.
     * 
     * @param orig - the original matrix
     * 
     * @throw MatrixException id+f dimensions do not match or allocation of memory fails
     */
    void __copyElems(const MatrixGeneric<T>& orig)
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
     * 
     * @param row - row number
     * @param column - column number
     * 
     * @return element's position witin the matrix's internal vector
     */
    inline std::size_t __pos(const std::size_t row, const std::size_t column) const
    {
        return ( row * this->m_cols + column );
    }



    /*
     * Extracts the specified triangular part(s) of the matrix and
     * assigns the result to itself.
     *
     * @param upper - extract the upper triangular part (without diagonal)
     * @param lower - extract the lower triangular part (without diagonal)
     * @param diag - extract the diagonal
     */
    void __triangPart(const bool upper, const bool lower, const bool diag)
    {
        const std::size_t N = this->m_elems.size();

        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            std::size_t r;
            std::size_t c;

            // iterator to the final element of the block:
            const typename std::vector<T>::const_iterator final = this->m_elems.begin() + iend;
            // iterator to the first/current element of this->m_elems:
            typename std::vector<T>::iterator it = this->m_elems.begin() + istart;

            for ( std::size_t i = istart;
                    it != final;
                    ++i, ++it )
            {
                T& currElem = *it;

                // row and column of the current index 'i':
                r = i / this->m_cols;
                c = i % this->m_cols;

                // Set the current element to 0, if specified by the input arguments
                if ( !( (r==c && true==diag) ||
                        (r>c && true==lower) ||
                        (r<c && true==upper) ) )
                {
                    currElem = static_cast<T>(0);
                }
            }
        }  // omp parallel
    }



    /*
     * A convenience function that finds the minimum or maximum
     * (depending on 'maxv') value of the given row or column
     * (depending on 'row'). If requested by 'absval', the absolute
     * minmax will be returned.
     *
     * @param rc - row/column number
     * @param row - if TRUE, return minmax value of the row 'rc', otherwise minmax value of the column 'rc'
     * @param maxv - should return row's/column's minimum (FALSE) or maximum (TRUE) value
     * @param absval - should return minmax value of the givenrow/column?
     *
     * @return minimum or maximum value of the given row/column
     *
     * @throw MatrixException if 'rc' is invalid
     */
    T __minmaxRowCol(const std::size_t rc, const bool row, const bool maxv, const bool absval) const
    {
        // Is 'rc' a valid row/column number?
        const std::size_t& ROWCOL = ( true==row ? this->m_rows : this->m_cols );
        if ( rc >= ROWCOL )
        {
            throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
        }

        // Number of elements to check:
        const std::size_t& N = ( true==row ? this->m_cols : this->m_rows );

        /*
         * The first value of the row/column, also
         * the first candidate for the minmax value
         */
        T retVal = ( true==row ? this->at(rc, 0) : this->at(0, rc) );
        if ( true == absval )
        {
            retVal = std::abs(retVal);
        }

        // Coarse grained parallelism
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(N, retVal)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            // The first element of the block, also the first candidate for a minmax value
            T temp = ( true==row ? this->at(rc, istart) : this->at(istart, rc) );
            if ( true == absval )
            {
                temp = std::abs(temp);
            }

            std::size_t idx = istart;
            for (std::size_t cntr = 0;
                 cntr<elems_per_thread && idx<N;
                 ++cntr, ++idx )
            {
                T curr = ( true==row ? this->at(rc, idx): this->at(idx, rc) );
                if ( true == absval )
                {
                     curr = std::abs(curr);
                }

                // Update 'temp' depending on 'maxv'
                temp = ( true==maxv ? std::max(temp, curr) : std::min(temp, curr) );
            }

            // prevent possible race condition when updating retVal
            #pragma omp critical(matrixgeneric_minmaxrowcol)
            {
                retVal = ( true==maxv ? std::max(retVal, temp) : std::min(retVal, temp) );
            }

            (void) iend;
        }  // omp parallel

        return retVal;
    }



    /*
     * A convenience function that obtains a vector of minimum or
     * maximum (depending on 'maxv') values of each row or column
     * (depending on 'row') of the given matrix 'm' and assigns
     * the resulting vector to itself. If requested by 'absval',
     * a vector of rows'/columns' minmax absolute values is assigned
     * to itself.
     *
     * @note 'this' is not required to be of the correct dimensions
     *       as the method takes care of it.
     *
     * @param m - input matrix
     * @param row - if TRUE, return minmax values of each row, otherwise minmax values of each column
     * @param maxv - should return rows'/columns' minimum (FALSE) or maximum (TRUE) values
     * @param absval - should return minmax values of m's absolute values
     *
     * @throw MatrixException if (re)allocation of memory failed
     */
    void __minmaxRowsCols(const MatrixGeneric<T>& m, const bool row, const bool maxv, const bool absval)
    {
        // dimensions of the "output" (i.e 'this'), depending on 'row':
        const std::size_t NROWS = ( true==row ?  m.m_rows : 1 );
        const std::size_t NCOLS = ( false==row ? m.m_cols : 1 );

        // length of the "output" vector
        const std::size_t N = ( true==row ? NROWS : NCOLS );

        // resize itself
        this->__init(NROWS, NCOLS);

        /*
         * Each m's row/column can be processed independently
         * from the others and in parallel.
         */
        #pragma omp parallel for default(none) shared(m)
        for ( std::size_t i=0; i<N; ++i )
        {
            /*
             * Note that 'this' is actually a one dimensional vector, hence
             * each row's/column's minmax can be written directly into the
             * appropriate position of 'm_elems'.
             */

            this->m_elems.at(i) = m.__minmaxRowCol(i, row, maxv, absval);
        }
    }


public:

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
    MatrixGeneric(const std::size_t rows, const std::size_t columns)
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
    MatrixGeneric(const std::size_t n)
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
    MatrixGeneric(const MatrixGeneric<T>& orig)
    {
        this->__copyElems(orig);
    }



    /**
     * @return number of rows
     */
    std::size_t nrRows() const
    {
        return this->m_rows;
    }



    /**
     * @return number of columns
     */
    std::size_t nrColumns() const
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
    T get(const std::size_t row, const std::size_t column) const
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
    T& at(const std::size_t row, const std::size_t column)
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
    const T& at(const std::size_t row, const std::size_t column) const
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
    T& operator()(const std::size_t row, const std::size_t column)
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
    const T& operator()(const std::size_t row, const std::size_t column) const
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
    T& operator()(const std::size_t idx)
    {
        // check if 'idx' is within elems' range
        if ( idx >= this->m_elems.size() )
        {
            throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
        }

        const std::size_t col = idx / this->m_rows;
        const std::size_t row = idx % this->m_rows;

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
    const T& operator()(const std::size_t idx) const
    {
        // check if 'idx' is within elems' range
        if ( idx >= this->m_elems.size() )
        {
            throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
        }

        const std::size_t col = idx / this->m_rows;
        const std::size_t row = idx % this->m_rows;

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
    MatrixGeneric<T>& set(const std::size_t row, const std::size_t column, const T& element)
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
    void display(std::ostream& str = std::cout) const
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

        const std::size_t tabsPerRow = this->m_cols - 1;

        for ( std::size_t r=0; r<(this->m_rows); ++r )
        {
            // display elements of the row r, separated by tabs
            for ( std::size_t c=0; c<(this->m_cols); ++c )
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
     * @return reference to itself
     *
     * @throw MatrixException if memory allocation fails
     */
    MatrixGeneric<T>& operator=(const MatrixGeneric<T>& orig)
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
     * Assignment operator (=) that assigns the scalar value
     * to the only field of a 1x1 matrix.
     *
     * @param scalar - a scalar to be assigned to the only field of the matrix
     *
     * @return reference to itself
     *
     * @throw MatrixException if memory allocation fails
     */
    MatrixGeneric<T>& operator=(const T& scalar)
    {
        this->__init(1, 1);
        this->m_elems.at(0) = scalar;

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
    MatrixGeneric<T>& operator+=(const MatrixGeneric<T>& m)
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
     * @param m - matrix to be subtracted from this one
     *
     * @return reference to this
     *
     * @throw MatrixException if dimensions do not match
     */
    MatrixGeneric<T>& operator-=(const MatrixGeneric<T>& m)
    {
        // Check if dimensions of both matrices match
        if ( this->m_rows != m.m_rows || this->m_cols != m.m_cols )
        {
            throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
        }

        // For a definition of matrix subtraction, see operator-
        math::mtvectadd<T>(this->m_elems, m.m_elems, this->m_elems, false);

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
    MatrixGeneric<T>& operator*=(const MatrixGeneric<T>&m )
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
    MatrixGeneric<T>& operator*=(const T& scalar)
    {
        // Multiply each element by the 'scalar'
        math::mtvectmult<T>(this->m_elems, scalar, this->m_elems);

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
    MatrixGeneric<T>& operator+=(const T& scalar)
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
    MatrixGeneric<T>& operator-=(const T& scalar)
    {
        math::mtvectscalaradd<T>(this->m_elems, scalar, this->m_elems, false, true);

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
    MatrixGeneric<T>& operator/=(const T& scalar)
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
    MatrixGeneric<T>& ewMult(const MatrixGeneric<T>& m)
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
    MatrixGeneric<T>& ewDiv(const MatrixGeneric<T>& m)
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




    /**
     * Matrix transpose operation
     *
     * @return this^T
     *
     * @throw MatrixException if matrix does not contain enough elements
     */
    MatrixGeneric<T> transpose() const
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
    
        const std::size_t& tcols = this->m_cols;

        const std::size_t N = this->m_rows * this->m_cols;

        // Coarse grained parallelism
        const std::vector<T>& els = this->m_elems;

        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(retVal, els, tcols)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            for ( std::size_t idx=istart; idx<iend; ++idx )
            {
                const std::size_t r = idx / tcols;
                const std::size_t c = idx % tcols;

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
    MatrixGeneric<T>& transpose_()
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

            const std::size_t& N = this->m_rows;    // number of rows (and columns)
            const std::size_t Ntr = N * (N-1) / 2;  // number of all elements to be transposed

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
            for ( std::size_t r=0; r<N-1; ++r )
            {
                for ( std::size_t c=r+1; c<N; ++c )
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
     * Conjugation of  all matrix' elements.
     * For complex types, elements' imaginary parts are reversed
     * their signs. For all other types, nothing is done.
     * 
     * @return reference to itself
     */
    MatrixGeneric<T>& conj_()
    {
        math::__matrixprivate::__matconj(*this);
        return *this;
    }



    /**
     * Conjugation of  all matrix' elements.
     * For complex types, elements' imaginary parts are reversed
     * their signs. For all other types, a copy of *this
     * is returned.
     * 
     * @return a new matrix with conjugated elements of *this
     * 
     * @throw MatrixException if allocation of memory fails
     */
    MatrixGeneric<T> conj() const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.conj_();
        return mret;
    }



    /**
     * "Rounds" all small elements (with the absolute value below
     * the default 'eps') to 0.
     * 
     * @return reference to itself
     */
    MatrixGeneric<T>& roundSmallElements_()
    {
        return this->roundSmallElements_( math::NumericUtil::getEPS<T>() );
    }



    /**
     * "Rounds" all small elements (with the absolute value below
     * the given 'eps') to 0.
     * 
     * @param eps - threshold to determine whether each component is "rounded" to 0
     * 
     * @return reference to itself
     */
    MatrixGeneric<T>& roundSmallElements_(const T& eps)
    {
        const std::size_t N = this->m_elems.size();

        // Coarse grained parallelism if OpenMP is enabled
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(eps)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            typename std::vector<T>::iterator it = this->m_elems.begin() + istart;
            for ( std::size_t i = istart;
                  i<iend && it!=this->m_elems.end();
                  ++it, ++i )
            {
                T& currElem = *it;

                currElem = math::NumericUtil::smallValToZero<T>(currElem, eps);
            }
        }  // pragma omp

        return *this;
    }



    /**
     * "Rounds" all small elements (with the absolute value below
     * the default 'eps') to 0.
     * 
     * @return a new matrix with "rounded" small elements
     * 
     * @throw MatrixException if allocation of memory for the new matrix failed
     */
    MatrixGeneric<T> roundSmallElements() const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.roundSmallElements_();
        return mret;
    }



    /**
     * "Rounds" all small elements (with the absolute value below
     * the given 'eps') to 0.
     * 
     * @param eps - threshold to determine whether each component is "rounded" to 0
     * 
     * @return a new matrix with "rounded" small elements
     * 
     * @throw MatrixException if allocation of memory for the new matrix failed
     */
    MatrixGeneric<T> roundSmallElements(const T& eps) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.roundSmallElements_(eps);
        return mret;
    }



    /**
     * Returns the minimum value of the given row. If requested via 'absval',
     * the row's minimum absolute value will be returned.
     *
     * @param row - row number
     * @param absval - should return the row's minimum absolute value? (default: FALSE)
     *
     * @return minimum value of the given row
     *
     * @throw MatrixException if 'row' is invalid
     */
    T minRow(const std::size_t row, const bool absval = false) const
    {
        return this->__minmaxRowCol(row, true, false, absval);
    }



    /**
     * Returns the maximum value of the given row. If requested via 'absval',
     * the row's maximum absolute value will be returned.
     *
     * @param row - row number
     * @param absval - should return the row's maximum absolute value? (default: FALSE)
     *
     * @return maximum value of the given row
     *
     * @throw MatrixException if 'row' is invalid
     */
    T maxRow(const std::size_t row, const bool absval = false) const
    {
        return this->__minmaxRowCol(row, true, true, absval);
    }



    /**
     * Returns the minimum value of the given column. If requested via 'absval',
     * the column's minimum absolute value will be returned.
     *
     * @param column - column number
     * @param absval - should return the column's minimum absolute value? (default: FALSE)
     *
     * @return minimum value of the given column
     *
     * @throw MatrixException if 'column' is invalid
     */
    T minColumn(const std::size_t column, const bool absval = false) const
    {
        return this->__minmaxRowCol(column, false, false, absval);
    }



    /**
     * Returns the maximum value of the given column. If requested via 'absval',
     * the column's maximum absolute value will be returned.
     *
     * @param column - column number
     * @param absval - should return the column's maximum absolute value? (default: FALSE)
     *
     * @return maximum value of the given column
     *
     * @throw MatrixException if 'column' is invalid
     */
    T maxColumn(const std::size_t column, const bool absval = false) const
    {
        return this->__minmaxRowCol(column, false, true, absval);
    }



    /**
     * Returns a nx1 vector with minimum values of each row.
     * The method is equivalent to the Matlab function min(A, [], 2)
     *
     * @param absval - should the function return a vector of minimum absolute values (default: FALSE)
     *
     * @return a column vector with minimum values of each row
     *
     * @throw MatrixException if allocation memory failed
     */
    MatrixGeneric<T> minRows(const bool absval = false) const
    {
        math::MatrixGeneric<T> ret(1, 1);
        ret.__minmaxRowsCols(*this, true, false, absval);
        return ret;
    }



    /**
     * Returns a nx1 vector with maximum values of each row.
     * The method is equivalent to the Matlab function max(A, [], 2)
     *
     * @param absval - should the function return a vector of maximum absolute values (default: FALSE)
     *
     * @return a column vector with maximum values of each row
     *
     * @throw MatrixException if allocation memory failed
     */
    MatrixGeneric<T> maxRows(const bool absval = false) const
    {
        math::MatrixGeneric<T> ret(1, 1);
        ret.__minmaxRowsCols(*this, true, true, absval);
        return ret;
    }



    /**
     * Returns a 1xn vector with minimum values of each column.
     * The method is equivalent to the Matlab function min(A, [], 1)
     *
     * @param absval - should the function return a vector of minimum absolute values (default: FALSE)
     *
     * @return a row vector with minimum values of each column
     *
     * @throw MatrixException if allocation memory failed
     */
    MatrixGeneric<T> minColumns(const bool absval = false) const
    {
        math::MatrixGeneric<T> ret(1, 1);
        ret.__minmaxRowsCols(*this, false, false, absval);
        return ret;
    }



    /**
     * Returns a 1xn vector with maximum values of each column.
     * The method is equivalent to the Matlab function max(A, [], 1)
     *
     * @param absval - should the function return a vector of maximum absolute values (default: FALSE)
     *
     * @return a row vector with maximum values of each column
     *
     * @throw MatrixException if allocation memory failed
     */
    MatrixGeneric<T> maxColumns(const bool absval = false) const
    {
        math::MatrixGeneric<T> ret(1, 1);
        ret.__minmaxRowsCols(*this, false, true, absval);
        return ret;
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
    MatrixGeneric<T>& removeRow_(const std::size_t rowNr)
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
    MatrixGeneric<T>& removeColumn_(const std::size_t colNr)
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
        for ( std::size_t i=1; i<=this->m_rows; ++i )
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
    MatrixGeneric<T>& insertRow_(const std::size_t rowNr, const T& el = static_cast<T>(0))
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
    MatrixGeneric<T>& insertColumn_(const std::size_t colNr, const T& el = static_cast<T>(0))
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
            for ( std::size_t r = 0; r < this->m_rows; ++r )
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
     * Returns a new matrix with the specified row number removed.
     *
     * @param rowNr - the row number to remove
     *
     * @return a new matrix with the specified row removed
     *
     * @throw MatrixException if attempting to remove the nonexistent row or if allocation of memory failed
     */
    MatrixGeneric<T> removeRow(const std::size_t rowNr) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.remoweRow_(rowNr);
        return mret;
    }



    /**
     * Returns a new matrix with the specified column removed.
     *
     * @param colNr - the column number to remove
     *
     * @return a new matrix with the specified column removed
     *
     * @throw MatrixException if attempting to remove the nonexistent column or if allocation of memory failed
     */
    MatrixGeneric<T> removeColumn(const std::size_t colNr) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.remoweColumn_(colNr);
        return mret;
    }



    /**
     * Returns a new matrix with a row inserted in front of the rowNr.th row.
     *
     * @param rowNr - a row will be inserted in front of this row. Valid values between 0 and rows
     * @param el - value to be assigned to all inserted elements (default: 0)
     *
     * @return a new matrix with a row inserted in front of the rowNr.th row
     *
     * @throw MatrixException if invalid rowNr or if allocation of memory failed
     */
    MatrixGeneric<T> insertRow(const std::size_t rowNr, const T& el = static_cast<T>(0)) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.insertRow_(rowNr, el);
        return mret;
    }



    /**
     * Returns a new matrix with a column inserted in front of the colNr.th column.
     *
     * @param colNr - a column will be inserted in front of this column. Valid values between 0 and cols
     * @param el - value to be assigned to all inserted elements (default: 0)
     * 
     * @return a new matrix with a column inserted in front of the colNr.th column
     *
     * @throw MatrixException if invalid colNr or if allocation of memory failed
     */
    MatrixGeneric<T> insertColumn(const std::size_t colNr, const T& el = static_cast<T>(0)) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.insertColumn_(colNr, el);
        return mret;
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
    MatrixGeneric<T>& swapRows_(const std::size_t r1, const std::size_t r2)
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
    MatrixGeneric<T>& swapColumns_(const std::size_t c1, const std::size_t c2)
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

        const std::size_t N = this->m_rows;

        // Coarse grained parallelization
        #pragma omp parallel num_threads(ompIdeal(N)) \
                        if(N>OMP_CHUNKS_PER_THREAD) \
                        default(none)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            for ( std::size_t r=istart; r<iend; ++r )
            {
                std::swap(
                    this->at(r, c1),
                    this->at(r, c2) );
            }
        }  // omp parallel

        return *this;
    }



    /**
     * Returns a new matrix with the specified rows swapped.
     *
     * @param r1 - first row's number
     * @param r2 - second row's number
     *
     * @return a new matrix with the specified rows swapped
     *
     * @throw MatrixException if any input argument is out of range or if allocation of memory failed
     */
    MatrixGeneric<T> swapRows(const std::size_t r1, const std::size_t r2) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.swapRows_(r1, r2);
        return mret;
    }



    /**
     * Returns a new matrix with the specified columns swapped.
     *
     * @param c1 - first column's number
     * @param c2 - second column's number
     *
     * @return a new matrix with the specified columns swapped
     *
     * @throw MatrixException if any input argument is out of range or if allocation of memory failed
     */
    MatrixGeneric<T> swapColumns(const std::size_t c1, const std::size_t c2) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.swapColumns_(c1, c2);
        return mret;
    }



    /**
     * Transforms itself into the upper triangular matrix, optionally including
     * the diagonal.
     * 
     * @param inclDiag - should the diagonal be included as well (default: TRUE)
     * 
     * @return reference to itself
     */
    MatrixGeneric<T>& upperTriangularPart_(const bool inclDiag=true)
    {
        this->__triangPart(true, false, inclDiag);
        return *this;
    }



    /**
     * Transforms itself into the lower triangular matrix, optionally including
     * the diagonal.
     * 
     * @param inclDiag - should the diagonal be included as well (default: TRUE)
     * 
     * @return reference to itself
     */
    MatrixGeneric<T>& lowerTriangularPart_(const bool inclDiag=true)
    {
        this->__triangPart(false, true, inclDiag);
        return *this;
    }



    /**
     * Transforms itself into the diagonal matrix.
     * 
     * @return reference to itself
     */    
    MatrixGeneric<T>& diagPart_()
    {
        this->__triangPart(false, false, true);
        return *this;
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
    MatrixGeneric<T> upperTriangularPart(const bool inclDiag=true) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.upperTriangularPart_(inclDiag);
        return mret;
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
    MatrixGeneric<T> lowerTriangularPart(const bool inclDiag=true) const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.lowerTriangularPart_(inclDiag);
        return mret;
    }



    /**
     * @return the diagonal part of the matrix
     * 
     * @throw MatrixException if allocation of memory fails
     */
    MatrixGeneric<T> diagPart() const
    {
        math::MatrixGeneric<T> mret(*this);
        mret.diagPart_();
        return mret;
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
    MatrixGeneric<T>& setDiag_(const T& scalar)
    {
        // Sanity check
        if ( false == this->isSquare() )
        {
            throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
        }

        // A double for loop will traverse the matrix, its diagonal elements
        // (row == column) will be set to the scalar, others to 0

        const std::size_t& N = this->m_rows;
        const std::size_t N2 = N * N;

        // Coarse grained parallelism:
        std::vector<T>& els = this->m_elems;

        #pragma omp parallel num_threads(ompIdeal(N2)) \
                    if(N2>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(N, els, scalar)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N2);

            typename std::vector<T>::iterator it = els.begin() + istart;
            for ( std::size_t idx = istart;
                  idx<iend && it!=els.end();
                  ++it, ++idx )
            {
                T& currElem = *it;

                const std::size_t r = idx / N;
                const std::size_t c = idx % N;

                currElem = ( r==c ? scalar : static_cast<T>(0) );
            }
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
    MatrixGeneric<T>& setUnit_()
    {
        // Actually this is a diagonal matrix with units (ones)
        // on its diagonal
        this->setDiag_(static_cast<T>(1));

        return *this;
    }



    /**
     * Calculates matrix's determinant.
     * This operation makes sense if T is float, double, Rational, Complex.
     * The result may be wrong if T is any implementation of int !!!!!!
     *
     * @note The method is only supported for square matrices
     * 
     * @param fullp - should the algorithm perform full pivoting (default: TRUE)
     * @param physSwap - should the internal algorithm perform physical swapping of matrix elements (default: FALSE)
     * 
     * @return determinant of the matrix
     *
     * @throw MatrixException if the matrix is not square or if allocation of memory for auxiliary variables fails
     */
    T determinant(const bool fullp=MATRIX_DET_FULL_PIVOT, const bool physSwap=MATRIX_PHYSSWAP_COEF) const;



    /**
     * Matrix inversion.
     * R = A^(-1) if R*A = A*R = I
     * If matrix's determinant equals 0, the matrix is not invertible
     *
     * @note The method is only supported for square matrices
     * 
     * @note The algorithm supports partial and full pivoting. Partial pivoting
     *       should be sufficient in most cases however full pivoting may sometimes be
     *       numerically more stable, on the other hand it introduces additional overhead.
     *
     * @param fullp - should the algorithm perform full pivoting (default: TRUE)
     * @param physSwap - should the internal algorithm perform physical swapping of matrix elements (default: FALSE)
     *
     * @return inverse matrix
     *
     * @throw MatrixException if the matrix is not square or not invertible
     */
    MatrixGeneric<T> inverse(const bool fullp=MATRIX_INVERSE_FULL_PIVOT, const bool physSwap=MATRIX_PHYSSWAP_COEF) const;



    /**
     * @return a logical value indicating whether the matrix is square
     */
    bool isSquare() const
    {
        return ( this->m_rows == this->m_cols );
    }



    /**
     * Calculates rank of a matrix.
     *
     * @param physSwap - should the internal algorithm perform physical swapping of matrix elements (default: FALSE)
     *
     * @return rank of the matrix
     *
     * @throw MatrixException if internal allocation of memory failed
     */
    std::size_t rank(const bool physSwap=MATRIX_PHYSSWAP_COEF) const;



    /**
     * Destructor
     */
    virtual ~MatrixGeneric()
    {
        // Vector's destructors would probably clean up this automatically.
        // Anyway let us clear the vector, just to be aware of allocated resources.
        this->m_elems.clear();

        // Other dynamically allocated memory (via malloc or new) should be freed here.
        // There are no other resources to release.
    }



    // Declaration of friend functions
    friend MatrixGeneric<T> (math::operator+ <>) (const MatrixGeneric<T>& m);
    friend MatrixGeneric<T> (math::operator- <>) (const MatrixGeneric<T>& m);
    friend MatrixGeneric<T> (math::operator+ <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);
    friend MatrixGeneric<T> (math::operator- <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);
    friend MatrixGeneric<T> (math::operator* <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);
    friend MatrixGeneric<T> (math::operator* <>) (const MatrixGeneric<T>& m, const T& sc);
    friend MatrixGeneric<T> (math::operator* <>) (const T& sc, const MatrixGeneric<T>& m);
    friend MatrixGeneric<T> (math::operator+ <>) (const MatrixGeneric<T>& m, const T& sc);
    friend MatrixGeneric<T> (math::operator+ <>) (const T& sc, const MatrixGeneric<T>& m);
    friend MatrixGeneric<T> (math::operator- <>) (const MatrixGeneric<T>& m, const T& sc);
    friend MatrixGeneric<T> (math::operator- <>) (const T& sc, const MatrixGeneric<T>& m);
    friend MatrixGeneric<T> (math::operator/ <>) (const MatrixGeneric<T>& m, const T& sc);
    friend MatrixGeneric<T> (math::matEwMult <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);
    friend MatrixGeneric<T> (math::matEwDiv <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2);
    friend void (math::__matrixprivate::__matconj <>) (MatrixGeneric<T>& m);
    friend void (math::__matrixprivate::__matconj <>) (MatrixGeneric<std::complex<T> >& m);
};  // class MatrixGeneric



/*
 * A separate namespace Pivot with auxilliary functions that perform partial or full pivoting. 
 */
namespace Pivot
{


template <class T>
void rearrangeMatrixRows(MatrixGeneric<T>& x, std::vector<std::size_t>& cols );


namespace __private
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
      )
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

}  // namespace math::Pivot::__private



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
bool getDiagonallyDominantMatrix(
        const MatrixGeneric<T>& A,
        MatrixGeneric<T>* const pB,
        std::vector<std::size_t>& rows,
        std::vector<std::size_t>& cols
      )
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
void rearrangeMatrixRows(
        MatrixGeneric<T>& x,
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
bool solveGaussJordan(
        const MatrixGeneric<T>& A,
        const MatrixGeneric<T>& b,
        MatrixGeneric<T>& x,
        const bool fullp,
        const bool physSwap
      )
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
T getDeterminant(
        const MatrixGeneric<T>& A,
        const bool fullp,
        const bool physSwap
      )
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
std::size_t getRank(
        const MatrixGeneric<T>& A,
        const bool physSwap
      )
{
    math::MatrixGeneric<T> temp(A);
    std::size_t rank;

    math::Pivot::__private::pivot<T>(temp, true, false, NULL, &rank, NULL, physSwap);

    return rank;
}


}  // namespace Pivot



/*
 * The following methods have already been declared inside the class MatrixGeneric.
 * They can only be impleneted after implementation of the namespace Pivot.
 */
 
template <class T>
T MatrixGeneric<T>::determinant(const bool fullp, const bool physSwap) const
{
	// Sanity check
	if ( false == this->isSquare() )
	{
		throw math::MatrixException(math::MatrixException::NONSQUARE_MATRIX);
	}

	return math::Pivot::getDeterminant(*this, fullp, physSwap);
}



template <class T>
MatrixGeneric<T> MatrixGeneric<T>::inverse(const bool fullp, const bool physSwap) const
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
	id.setUnit_();

	// Inverse matrix is a solution (if it exists) of the equation:
	// this * inv = id
	math::MatrixGeneric<T> retVal(id);

	const bool succ = math::Pivot::solveGaussJordan<T>(*this, id, retVal, fullp, physSwap);

	// is *this an uninvertible matrix? (determinant()=0):
	if ( false == succ )
	{
		throw math::MatrixException(math::MatrixException::NON_INVERTIBLE_MATRIX);
	}

	return retVal;
}



template <class T>
std::size_t MatrixGeneric<T>::rank(const bool physSwap) const
{
	return math::Pivot::getRank(*this, physSwap);
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
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m)
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
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m)
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
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2)
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
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2)
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
MatrixGeneric<T> operator*(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2)
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
    for ( std::size_t r=0; r<m1.m_rows; ++r )
    {
        for ( std::size_t c=0; c<m2.m_cols; ++c)
        {
            T sum = static_cast<T>(0);
            for ( std::size_t i=0; i<m1.m_cols; ++i )
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
MatrixGeneric<T> operator*(const MatrixGeneric<T>& m, const T& sc)
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
MatrixGeneric<T> operator*(const T& sc, const MatrixGeneric<T>& m)
{
    math::MatrixGeneric<T> retVal = m * sc;
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
MatrixGeneric<T> operator+(
            const MatrixGeneric<T>& m,
            const T& sc)
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
MatrixGeneric<T> operator+(
            const T& sc,
            const MatrixGeneric<T>& m)
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
MatrixGeneric<T> operator-(
            const MatrixGeneric<T>& m,
            const T& sc)
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
MatrixGeneric<T> operator-(
            const T& sc,
            const MatrixGeneric<T>& m)
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
MatrixGeneric<T> operator/(
            const MatrixGeneric<T>& m,
            const T& sc)
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
MatrixGeneric<T> matEwMult(
            const MatrixGeneric<T>& m1,
            const MatrixGeneric<T>& m2)
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
MatrixGeneric<T> matEwDiv(
            const MatrixGeneric<T>& m1,
            const MatrixGeneric<T>& m2)
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
 * does not modify 'm'
 *
 * @param m - matrix to conjugate
 */
template <class T>
void __matrixprivate::__matconj( MatrixGeneric<T>& m )
{
    // nothing to do if 'T' is not a complex type
    return;
    (void) m;
}



/*
 * Overloading (partial "specialization") of __matconj
 * for complex numbers. It actually conjugates each element of 'm'
 *
 * @param m - matrix to conjugate and write conjugated elements to
 */
template <class T>
void __matrixprivate::__matconj( MatrixGeneric<std::complex<T> >& m )
{
    /*
     * Specialization into another templated class implemented
     * as suggested here:
     * http://www.cplusplus.com/forum/general/68298/
     */
    const std::size_t N = m.m_elems.size();

    //Coarse grained parallelization
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(m)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<std::complex<T> >::iterator it = m.m_elems.begin() + istart;
        for ( std::size_t i=istart;
              i<iend && it!=m.m_elems.end(); ++it, ++i )
        {
            std::complex<T>& currElem = *it;
 
            currElem = std::conj(currElem);
        }
    }  // omp parallel

}



// Matrices with elements of types float, double and long double
// make most sense so these three types are predefined
typedef MatrixGeneric<float>       FMatrix;
typedef MatrixGeneric<double>      Matrix;
typedef MatrixGeneric<long double> LDMatrix;

} // namespace math


#endif  // _MATH_MATRIXGENERIC_HPP_
