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
* @file SqMatrixGeneric.cpp
*
* Implementation of the class SqMatriGeneric.
*
* As the class is templated, this file must not be compiled.
* Instead it must be included after the class declaration in the .h file
*
* @author Jernej Kovacic
*/

// Deliberately there is no #include "SqMatrixGeneric.h" !
#include "NumericUtil.h"
#include "MatrixException.h"
#include "MatrixGeneric.h"

#include <stdexcept>

// Matrix's element of the row r and column c is accessed as "r*N+c"
// (N must be additionaly declared as either this->rows or this->cols).
// This expression would be used frequently inside this file. To improve
// readability, to simplify code reviews and to eliminate possibilities of
// typing errors, this macro has been defined:
#define ELM(r,c)    ( (r) * N + (c) )

/**
 * Constructor.
 * Creates an instance of a square matrix with the specified number of rows and columns.
 *
 * @param dim  - number of rows and columns (default: 1)
 *
 * @throw MatrixException when allocation ofmemory fails or in case of incorrect
 *        input parameters (dim must be at least 1)
 */
template<class T>
SqMatrixGeneric<T>::SqMatrixGeneric(unsigned int dim) throw (MatrixException) :
    MatrixGeneric<T>(dim, dim)
{
    // Square matrices have the same number of rows and columns.
    // This is specified by the initializer list.
    // The class has no extra members, so nothing else to do
}

/**
 * Copy constructor.
 * Creates an instance of a matrix with the same dimensions as orig
 * and copies its elements. Orig must have the same number of rows and columns,
 * otherwise an exception is thrown.
 *
 * @param orig - original matrix to be copied into this one
 *
 * @throw MatrixException if orig has different number of rows and columns
 */
template<class T>
SqMatrixGeneric<T>::SqMatrixGeneric(const MatrixGeneric<T>& orig) throw(MatrixException) :
    MatrixGeneric<T>(orig)
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
 * @throw MatrixException if memory allocation fails or attempting to assign a nonsquare matrix
 */
template<class T>
SqMatrixGeneric<T>& SqMatrixGeneric<T>::operator= (const MatrixGeneric<T>& m) throw (MatrixException)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &m )
    {
        return *this;
    }

    // m must be a "square" matrix, i.e. with equal number of rows and columns
    if ( m.nrRows() != m.nrColumns() )
    {
        throw MatrixException(MatrixException::INVALID_DIMENSION);
    }

    copyElems(m);

    return *this;
}

/**
 * Modifies the matrix into a diagonal matrix.
 * Values of all diagonal elements are set to scalar, the other ones to zero.
 *
 * @param scalar - value of diagonal elements
 *
 * @throw MatrixException
 */
template<class T>
void SqMatrixGeneric<T>::setDiag(const T& scalar) throw(MatrixException)
{
    // A double for loop will traverse the matrix, its diagonal elements
    // (row == column) will be set to the scalar, others to 0
    try
    {
        unsigned int i;
        unsigned int j;
        const unsigned int N = this->rows;

        for ( i=0; i<N; i++ )
        {
            for ( j=0; j<N; j++ )
            {
                this->elems.at(ELM(i, j)) = ( 0 == i-j ? scalar : (T) 0 );
            } // for j
        } // for i
    }  // try
    catch ( std::out_of_range& oor )
    {
        throw MatrixException(MatrixException::OUT_OF_RANGE);
    }
}

/**
 * Modifies the matrix into a unit matrix (a diagonal matrix with ones on the diagonal)
 *
 * @throw MatrixException
 */
template<class T>
void SqMatrixGeneric<T>::setUnit() throw(MatrixException)
{
    // Actually this is a diagonal matrix with units (ones)
    // on its diagonal
    setDiag((T) 1);
}

/**
 * Calculates matrix's determinant.
 * This operation makes sense if T is float, double, Rational, Complex.
 * The result may be wrong if T is any implementation of int !!!!!!
 *
 * @return determinant
 *
 * @throw MatrixException
 */
template<class T>
T SqMatrixGeneric<T>::determinant() const throw(MatrixException)
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
    T retVal = (T) 1;

    try
    {
        // We do not want to modify the matrix, therefore its vector of
        // elements will be copied into a temporary one where any modifications
        // are permitted
        std::vector<T> temp = this->elems;
        unsigned int i;
        unsigned int r;
        unsigned int c;

        const unsigned int N = this->rows;  // number of rows (and columns)

        for ( i=0; i<N-1; i++ )
        {
            // if temp(i,i) equals zero, swap the i^th line with
            // another one (r; r>i) satisfying temp(r,i)!=0
            // Each swap multiplies the determinant by -1

            if ( true == NumericUtil<T>::isZero(temp.at(ELM(i, i))) )
            {

                // Line swap will be necessary.
                // Find the r^th line meeting the criteria above
                // If not found, the determinant will be 0.

                for ( r=i+1; r<N; r++ )
                {
                    if ( false == NumericUtil<T>::isZero(temp.at(ELM(r, i))) )
                    {
                        // Found, no need to search further,
                        // so end the for (r) loop
                        break;  // out of for (r)
                    }
                }  // for r

                // if no appropriate row r was found, the determinant will be 0
                // and the method is finished
                if ( 0 == r - N )
                {
                    return (T) 0;
                }

                // otherwise swap the lines one element by one

                // swap i.th and r.th line by replacing elements one by one

                // BTW, all elements left of (i,i) and (r,i) should already be
                // equal to 0 and it wouldn't be necessary to swap them.
                // But a few extra "operations" shouldn't considerably affect complexity
                for ( c=ELM(i,0); c<ELM(i+1,0); c++ )
                {
                    const T tempElem = temp.at(c);
                    temp.at(c) = temp.at(ELM(r, c));
                    temp.at(ELM(r, c)) = tempElem;
                }

                // finally, if two lines are swaped, det = -det
                retVal = -retVal;
            } // if temp(i,i) equals 0

            // Now temp(i,i) definitely does not equal 0
            // all temp(r,i) will be "set" to 0 where r>i.

            // Note that determinant is not changed if
            // a multiplier of one line (i in this algorithm)
            // is added to another line (r; r>i)

            for ( r=i+1; r<N; r++ )
            {
                // temp(r,i) will be calculated to 0 immediately.
                // However, its initial value is necessary to properly
                // calculate all other elements of the r^th row
                T ri = temp.at(ELM(r, i));

                for ( c=i; c<N; c++ )
                {
                    // temp(r,c) = temp(r,c) - temp(i,c) * temp(r,i) / temp(i,i)
                    temp.at(ELM(r, c)) -= temp.at(ELM(i, c)) * ri / temp.at(ELM(i, i));
                }  // for c
            }  // for r
        }  // for i

        // Now temp is an upper triangular matrix so all its diagonal
        // elements can be multipled
        for ( i=0; i<N; i++)
        {
            retVal *= temp.at(ELM(i, i));
        }
    } // try
    catch ( std::out_of_range& oor )
    {
        throw MatrixException(MatrixException::OUT_OF_RANGE);
    }

    return retVal;
}

/**
 * Matrix inversion.
 * R = A^(-1) if R*A = A*R = I
 * If matrix's determinat equals 0, the matrix is not invertible
 *
 * @return inverse matrix
 * 
 * @throw MatrixException if the matrix is not invertible
 */
// This macro was defined specifically for inversion function. It requires an
// auxilary "matrix" (2*N, N) and its elements are accessed as "r*2*N+c". N must
// be declared and be equal to number of rows of original matrix
#define TMPELM(r,c)    ( (r) * 2*N + (c) )

template<class T>
SqMatrixGeneric<T> SqMatrixGeneric<T>::inverse() const throw(MatrixException)
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
     */

    const unsigned int N = this->rows;  // number of rows (and columns)
    // Result will be a matrix with the same dimensions.
    // Let it be instantiated now and later its elements will be set.
    SqMatrixGeneric<T> retVal(N);

    try
    {
        std::vector<T> temp;
        unsigned int i;
        unsigned int r;
        unsigned int c;
        T el;

        // reserve enough space for the temporary "matrix" N x 2*N
        try
        {
            temp.resize(2*N*N);
        }
        catch (std::bad_alloc &ba)
        {
            throw MatrixException(MatrixException::OUT_OF_MEMORY);
        }

        // Fill the temp "matrix":
        // its left half is the actual matrix, the right half will be an identity matrix

        for ( r=0; r<N; r++ )
        {
            for ( c=0; c<N; c++ )
            {
                temp.at(TMPELM(r, c)) = this->elems.at(ELM(r, c));
                temp.at(TMPELM(r, c + N)) = ( 0 == r-c ? (T) 1 : (T) 0 );
            }  // for c
        }  // for r
        // The temp matrix is filled accordingly

        // Finally try to convert the left half into an identity matrix
        // by appropriate adding multiples of other lines to each line
        for ( i=0; i<N; i++ )
        {
            // first check if the diagonal element equals 0
            if ( true == NumericUtil<T>::isZero(temp.at(TMPELM(i, i))) )
            {
                // if yes, try to find another row r where temp(r,i)!=0
                for ( r=0; r<N; r++ )
                {
                    if ( 0 == r-i )
                    {
                        // it is known in advance, that temp(i,i)==0, so skip it
                        continue;  // for r
                    }

                    if ( false == NumericUtil<T>::isZero(temp.at(TMPELM(r, i))) )
                    {
                        // found, no need to search further
                        break;  // out of for r
                    }
                }  // for r

                if ( 0 == r-N )
                {
                    // No temp(r,i)!=0 was found, the matrix is non-invertible.
                    // Throw an exception
                    throw MatrixException(MatrixException::NON_INVERTIBLE_MATRIX);
                }

                // add the r^th line to the i^th one
                for ( c=0; c<2*N; c++ )
                {
                    temp.at(TMPELM(i, c)) += temp.at(TMPELM(r, c));
                }  // for c
            }  // if temp(i,i)==0

            // Let the diag element be 1. So divide the whole row by temp(i,i)
            //  (columns smaller than i are already 0)
            el = temp.at(TMPELM(i, i));
            for ( c=i; c<2*N; c++ )
            {
                temp.at(TMPELM(i, c)) /= el;
            }  // for c

            // set the i^th column of all other rows (r>i) to 0 by
            // adding the appropriate multiple of the i^th row
            for ( r=i+1; r<N; r++ )
            {
                // Nothing to do if temp(r,i) is already 0.
                if ( true == NumericUtil<T>::isZero(temp.at(TMPELM(r, i))) )
                {
                    continue;  // for r
                }

                // Subtract a multiple of the i^th row. Note that temp(i,i) is already 1.
                el = temp.at(TMPELM(r, i));
                for ( c=i; c<2*N; c++ )
                {
                    temp.at(TMPELM(r, c)) -= el * temp.at(TMPELM(i, c));
                }  // for c
            }  // for r
        }  // for i

        // Now the lower triangle (below diag excl.) is 0, the diagonal consists of 1,
        // The upper triangle (above the diag) must be set to 0 as well.

        for ( r=0; r<N; r++ )
        {
            for ( c=r+1; c<N; c++ )
            {
                // Nothing to do if already 0
                if ( true == NumericUtil<T>::isZero(temp.at(TMPELM(r, c))) )
                {
                    continue;  // for c
                }

                // To set temp(r,c) to 0 it is a good idea to add the c^th row to it.
                // temp(c,i); i<c are already 0 (i.e. will not affect anything left of temp(i,c)
                // and temp(c,c) is already 1.

                el = temp.at(TMPELM(r, c));
                for ( i=c; i<2*N; i++ )
                {
                    temp.at(TMPELM(r, i)) -= el * temp.at(TMPELM(c, i));
                }  // for i
            }  // for c
        }  // for r

        // The right half is now the inverse matrix.
        // Copy it into retVal.
        for ( r=0; r<N; r++ )
        {
            for ( c=0; c<N; c++ )
            {
                retVal.elems.at(ELM(r, c)) = temp.at(TMPELM(r, c + N));
            } // for c
        }  // for r
    } // try
    catch ( std::out_of_range& oor )
    {
        throw MatrixException(MatrixException::OUT_OF_RANGE);
    }

    return retVal;
}
// This macro was defined especially for the inverse function so it should be undef'ed
#undef TMPELM

/*
 * The following functions are inherited from the base class, however they change
 * matrix dimensions so they should not be permitted in a square matrix.
 * (insert or removal of a row/column of a square matrix results in a non-square matrix!)
 * Therefore they will automatically throw an exception if called.
 */
template<class T>
void SqMatrixGeneric<T>::removeRow(unsigned int rowNr) throw (MatrixException)
{
    throw MatrixException(MatrixException::FORBIDDEN);
}

template<class T>
void SqMatrixGeneric<T>::removeColumn(unsigned int colNr) throw (MatrixException)
{
    throw MatrixException(MatrixException::FORBIDDEN);
}

template<class T>
void SqMatrixGeneric<T>::insertRow(unsigned int rowNr) throw (MatrixException)
{
    throw MatrixException(MatrixException::FORBIDDEN);
}

template<class T>
void SqMatrixGeneric<T>::insertColumn(unsigned int colNr) throw (MatrixException)
{
    throw MatrixException(MatrixException::FORBIDDEN);
}

// The macro was defined for implementation in this file only. Undef it now
#undef ELM
