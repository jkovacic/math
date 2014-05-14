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
#include "LinearEquationSolverGeneric.h"
#include "LinearEquationSolverException.h"

#include <stdexcept>

// The function 'pos', defined in the parent class, must be called as
// this->pos. The following macro has been defined to shorten this and
// to make the code a little bit more readable:
#define POS(r, c)   this->pos((r), (c))

// 'Zero' and 'one' constant have already been defined in the class NumericUtil.
// They can only be accessed as math::NumericUtil<T>::ZERO or math::NumericUtil<T>::ONE, respectively
// As this notation is a bit long, these convenience macro are defined:
#define ZERO math::NumericUtil<T>::ZERO
#define ONE  math::NumericUtil<T>::ONE


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

    this->copyElems(m);

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
    try
    {
        const size_t N = this->rows;

        for ( size_t i=0; i<N; ++i )
        {
            for ( size_t j=0; j<N; ++j )
            {
                this->elems.at(POS(i, j)) = ( i==j ? scalar : ZERO );
            } // for j
        } // for i
    }  // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

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
    setDiag(ONE);

    return *this;
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
    T retVal = ONE;

    try
    {
        // We do not want to modify the matrix, therefore its vector of
        // elements will be copied into a temporary one where any modifications
        // are permitted
        std::vector<T> temp = this->elems;
        size_t r;
        
        const size_t N = this->rows;  // number of rows (and columns)

        for ( size_t i=0; i<N-1; ++i )
        {
            // if temp(i,i) equals zero, swap the i^th line with
            // another one (r; r>i) satisfying temp(r,i)!=0
            // Each swap multiplies the determinant by -1

            if ( true == math::NumericUtil<T>::isZero(temp.at(POS(i, i))) )
            {

                // Line swap will be necessary.
                // Find the r^th line meeting the criteria above
                // If not found, the determinant will be 0.

                for ( r=i+1; r<N; ++r )
                {
                    if ( false == NumericUtil<T>::isZero(temp.at(POS(r, i))) )
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
                    return ZERO;
                }

                // otherwise swap the lines one element by one

                // swap i.th and r.th line by replacing elements one by one

                // BTW, all elements left of (i,i) and (r,i) should already be
                // equal to 0 and it wouldn't be necessary to swap them.
                // But a few extra "operations" shouldn't considerably affect complexity
                for ( size_t c=POS(i,0); c<POS(i+1,0); ++c )
                {
                    const T tempElem = temp.at(c);
                    temp.at(c) = temp.at(POS(r, c));
                    temp.at(POS(r, c)) = tempElem;
                }

                // finally, if two lines are swapped, det = -det
                retVal = -retVal;
            } // if temp(i,i) equals 0

            // Now temp(i,i) definitely does not equal 0
            // all temp(r,i) will be "set" to 0 where r>i.

            // Note that determinant is not changed if
            // a multiplier of one line (i in this algorithm)
            // is added to another line (r; r>i)

            for ( size_t r=i+1; r<N; ++r )
            {
                // temp(r,i) will be calculated to 0 immediately.
                // However, its initial value is necessary to properly
                // calculate all other elements of the r^th row
                T ri = temp.at(POS(r, i));

                for ( size_t c=i; c<N; ++c )
                {
                    // temp(r,c) = temp(r,c) - temp(i,c) * temp(r,i) / temp(i,i)
                    temp.at(POS(r, c)) -= temp.at(POS(i, c)) * ri / temp.at(POS(i, i));
                }  // for c
            }  // for r
        }  // for i

        // Now temp is an upper triangular matrix so all its diagonal
        // elements can be multiplied
        for ( size_t i=0; i<N; ++i )
        {
            retVal *= temp.at(POS(i, i));
        }

        // temp not needed anymore, clean it
        temp.clear();
    } // try
    catch ( const std::out_of_range& oor )
    {
        throw math::MatrixException(math::MatrixException::OUT_OF_RANGE);
    }

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
        SqMatrixGeneric<T> retVal = leq.solve();

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
    const size_t N = this->rows;  // number of rows (and columns)
    T temp;

    // Traverse the upper diagonal part of the matrix,
    // no need to reach the final row and
    // no need to transpose elements on the diagonal:
    for ( size_t r=0; r<N-1; ++r )
    {
        for ( size_t c=r+1; c<N; ++c )
        {
            temp = this->elems.at(POS(r, c));
            this->elems.at(POS(r, c)) = this->elems.at(POS(c, r));
            this->elems.at(POS(c, r)) = temp;
        }  // for c
    }  // for r

    return *this;
}

/**
 * "Reimplementation" of operator*= that multiplies a matrix by this one
 * and assigns the product to itself.
 *
 * The product must remain a square matrix, for that reason 'm' must
 * be a square matrix (it can be declared as an "ordinary" matrix, though)
 * with the same dimensions as 'this'
 *
 * @param m
 *
 * @return reference to itself
 *
 * @throw
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
}

template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::removeColumn(size_t colNr) throw (math::MatrixException)
{
    throw math::MatrixException(math::MatrixException::FORBIDDEN);
}

template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::insertRow(size_t rowNr, const T& el) throw (math::MatrixException)
{
    throw math::MatrixException(math::MatrixException::FORBIDDEN);
}

template<class T>
math::MatrixGeneric<T>& math::SqMatrixGeneric<T>::insertColumn(size_t colNr, const T& el) throw (math::MatrixException)
{
    throw math::MatrixException(math::MatrixException::FORBIDDEN);
}

// The macros were defined for implementation in this file only. Undef them now
#undef POS
#undef ZERO
#undef ONE
