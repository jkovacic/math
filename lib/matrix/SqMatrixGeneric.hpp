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
 * @headername{SqMatrixGeneric.h}
 *
 * Declaration of the class SqMatrixGeneric, representing square matrices.
 * The class is derived from MatrixGeneric and has additional functions defined
 * for square matrices only.
 */

#ifndef _MATH_SQMATRIXGENERIC_HPP_
#define	_MATH_SQMATRIXGENERIC_HPP_


#include <cstddef>

#include "matrix/MatrixGeneric.hpp"
#include "exception/MatrixException.hpp"

namespace math
{

// Some compilers may require forward declaration
// of the templated base class:
template <class T> class MatrixGeneric;


/**
 * @brief A class representing square matrices (number of rows is equal to number of columns).
 * 
 * The class is derived from a more general class MatrixGeneric.
 * 
 * Additional functionalities, e.g. determinant, inverse matrix etc., are implemented.
 * Eigen values and vectors are not implemented yet.
 */
template <class T>
class SqMatrixGeneric : public MatrixGeneric<T>
{

public:
    // Constructor
    SqMatrixGeneric(const size_t dim = 1) throw(MatrixException);
    // Copy constructor
    SqMatrixGeneric(const MatrixGeneric<T>& orig) throw(MatrixException);
    // operator= (it must be reimplemented as it is not inherited from the base class)
    SqMatrixGeneric<T>& operator= (const MatrixGeneric<T>& m) throw (MatrixException);

    // Several methods to create diagonal matrices
    SqMatrixGeneric<T>& setDiag(const T& scalar) throw(MatrixException);
    SqMatrixGeneric<T>& setUnit() throw(MatrixException);

    // Determinant of the matrix
    T determinant() const throw(MatrixException);

    // Inverse matrix
    SqMatrixGeneric<T> inverse() const throw(MatrixException);

    // Special (memory efficient) reimplementation of a virtual method
    SqMatrixGeneric<T>& transposed() throw(MatrixException);

    // Inherited as virtual from the parent class
    SqMatrixGeneric<T>& operator*=(const MatrixGeneric<T>& m) throw (MatrixException);

    // The following functions inherited from the parent class
    // are not applicable for square matrices and throw an exception
    MatrixGeneric<T>& removeRow(const size_t rowNr) throw (MatrixException);
    MatrixGeneric<T>& removeColumn(const size_t colNr) throw (MatrixException);
    MatrixGeneric<T>& insertRow(const size_t rowNr, const T& el = static_cast<T>(0)) throw (MatrixException);
    MatrixGeneric<T>& insertColumn(const size_t colNr, const T& el = static_cast<T>(0)) throw (MatrixException);
};

// Matrices with elements of types float, double and long double
// make most sense so these three types are predefined
typedef SqMatrixGeneric<float>       FSqMatrix;
typedef SqMatrixGeneric<double>      SqMatrix;
typedef SqMatrixGeneric<long double> LDSqMatrix;

} // namespace math

// DEFINITION
#include "matrix/SqMatrixGeneric.cpp"

#endif	// _MATH_SQMATRIXGENERIC_HPP_
