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
@file MatrixGeneric.h

Declaration of the class MatrixGeneric, representing generic matrices.

@author Jernej Kovacic
*/

#ifndef _MATH_MATRIXGENERIC_H_
#define	_MATH_MATRIXGENERIC_H_

#include <vector>
#include <iostream>

#include "MatrixException.h"

namespace math
{

// Templates are used to allow several types (T) of matrix elements.
// Types must have implemented basic arithmetic operators (+, -, *, /),
// otherwise build will fail (which is desired). In practice, types as float,
// double, Rational, Complex make most sense. Integer types may be conditionally
// acceptable (e.g. if the class is used to represent mathematical graphs, etc.)
// but arithmetic operations (e.g. inversion of a square matrix) may return
// incorrect results. This is true for unsigned types as well.

// Advance declaration of the class is necessary...
template<class T> class MatrixGeneric;
// to declare the class's friend function:
template<class T>
MatrixGeneric<T> operator* (const T& sc, const MatrixGeneric<T>& m) throw (math::MatrixException);

template <class T>
class MatrixGeneric
{

    // These properties must be accessible in inhereted classes
protected:
    size_t rows;      /// Number of rows
    size_t cols;      /// Number of columns
    // STL Vector has several advantages over arrays, allocated by new[],
    // e.g. elements can be accessed via at() which checks for range and throws
    // an exception when attempting to access sth. outside the allocated range
    // (e.g. in case of a typing error) which may result in a segmentation fault (crash).
    // Operations such as inserting or removing of elements are simplified as well.
    std::vector<T> elems;   /// Elements of the matrix

    // Copy elements from one matrix into another. Used at copy constructors,
    // assignemt operators etc. It s also suitable for use in derived  classes,
    // so it should be 'protected' instead of 'private'
    void copyElems(const MatrixGeneric& orig) throw (MatrixException);

    // position of the element within the internal vector
    inline size_t pos(size_t row, size_t column) const;

public:
    // Constructor
    MatrixGeneric(size_t rows = 1, size_t columns = 1) throw (MatrixException);
    // Copy constructor
    MatrixGeneric(const MatrixGeneric& orig) throw (MatrixException);

    // Number of rows and columns
    size_t nrRows() const;
    size_t nrColumns() const;

    // Get and set the element of the specified row and column
    T get(size_t row, size_t column) const throw (MatrixException);
    T& at(size_t row, size_t column) throw (MatrixException);
    const T& at(size_t row, size_t column) const throw(MatrixException);
    MatrixGeneric<T>& set(size_t row, size_t column, const T& element) throw (MatrixException);

    // Display elements of the matrix
    void display(std::ostream& str = std::cout) const throw (MatrixException);

    // Matrix arithmetics operators
    virtual MatrixGeneric<T>& operator= (const MatrixGeneric<T>& m) throw (MatrixException);
    MatrixGeneric<T> operator+ (const MatrixGeneric<T>& m) const throw (MatrixException);
    MatrixGeneric<T>& operator+= (const MatrixGeneric<T>& m) throw (MatrixException);
    MatrixGeneric<T> operator- (const MatrixGeneric<T>& m) const throw (MatrixException);
    MatrixGeneric<T>& operator-= (const MatrixGeneric<T>& m) throw (MatrixException);
    MatrixGeneric<T> operator* (const MatrixGeneric<T>& m) const throw (MatrixException);
    virtual MatrixGeneric<T>& operator*= (const MatrixGeneric<T>&m ) throw (MatrixException);
    MatrixGeneric<T> operator* (const T& sc) const throw (MatrixException);
    MatrixGeneric<T>& operator*= (const T& scalar);
    MatrixGeneric<T> operator- () const throw (MatrixException);
    // A friend function that multiplies a scalar and a matrix
    friend MatrixGeneric<T> (math::operator* <>) (const T& sc, const MatrixGeneric<T>& m) throw (MatrixException);

    // Transpose the matrix
    MatrixGeneric<T> transpose() const throw (MatrixException);
    virtual MatrixGeneric<T>& transposed() throw (MatrixException);

    // Insert or remove rows/columns.
    virtual MatrixGeneric<T>& removeRow(size_t rowNr) throw (MatrixException);
    virtual MatrixGeneric<T>& removeColumn(size_t colNr) throw (MatrixException);
    virtual MatrixGeneric<T>& insertRow(size_t rowNr) throw (MatrixException);
    virtual MatrixGeneric<T>& insertColumn(size_t colNr) throw (MatrixException);

    // Destructor
    virtual ~MatrixGeneric();

};

// Matrices with elements of types float and double make most sense
// so these two types are predefined
typedef MatrixGeneric<float> FMatrix;
typedef MatrixGeneric<double> Matrix;

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
} // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "MatrixGeneric.cpp"

#endif	/* _MATH_MATRIXGENERIC_H_ */
