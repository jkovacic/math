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
 * @headername{MatrixGeneric.h}
 *
 * Declaration of the class MatrixGeneric, representing generic matrices.
 */

#ifndef _MATH_MATRIXGENERIC_HPP_
#define	_MATH_MATRIXGENERIC_HPP_

#include <vector>
#include <ostream>
#include <cstddef>
#include <complex>

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
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator*(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator*(const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator*(const T& sc, const MatrixGeneric<T>& m) throw (MatrixException);

template <class T>
MatrixGeneric<T> operator+(const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator+(const T& sc, const MatrixGeneric<T>& m) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator-(const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator-(const T& sc, const MatrixGeneric<T>& m) throw(MatrixException);

template <class T>
MatrixGeneric<T> operator/(const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);

template <class T>
MatrixGeneric<T> matEwMult(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);

template <class T>
MatrixGeneric<T> matEwDiv(const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);


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
    size_t m_rows;      /// Number of rows
    size_t m_cols;      /// Number of columns

    /*
     * STL Vector has several advantages over arrays, allocated by new[],
     * e.g. elements can be accessed via at() which checks for range and throws
     * an exception when attempting to access sth. outside the allocated range
     * (e.g. in case of a typing error) which may result in a segmentation fault (crash).
     * Operations such as inserting or removing of elements are simplified as well.
     */
    std::vector<T> m_elems;   /// Elements of the matrix


    // Initializes a new matrix
    void __init(const size_t rows, const size_t cols) throw (MatrixException);

    /*
     * Copy elements from one matrix into another. Used at copy constructors,
     * assignment operators etc. It s also suitable for use in derived  classes,
     * so it should be 'protected' instead of 'private'
     */
    void __copyElems(const MatrixGeneric<T>& orig) throw (MatrixException);

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
    inline size_t __pos(const size_t row, const size_t column) const
    {
        return ( row * this->m_cols + column );
    }

    // a convenience function to extract the upper/lower triangular part
    void __triangPart(const bool upper, const bool lower, const bool diag);

public:
    // Constructors
    MatrixGeneric(const size_t rows, const size_t columns) throw (MatrixException);
    MatrixGeneric(const size_t n) throw (MatrixException);
    // Copy constructor
    MatrixGeneric(const MatrixGeneric<T>& orig) throw (MatrixException);

    // Number of rows and columns
    size_t nrRows() const;
    size_t nrColumns() const;

    // Get and set the element of the specified row and column
    T get(const size_t row, const size_t column) const throw (MatrixException);
    T& at(const size_t row, const size_t column) throw (MatrixException);
    const T& at(const size_t row, const size_t column) const throw(MatrixException);
    T& operator()(const size_t row, const size_t column) throw(MatrixException);
    const T& operator()(const size_t row, const size_t column) const throw(MatrixException);
    T& operator()(const size_t idx) throw(MatrixException);
    const T& operator()(const size_t idx) const throw(MatrixException);
    MatrixGeneric<T>& set(const size_t row, const size_t column, const T& element) throw (MatrixException);  

    // Display elements of the matrix
    void display(std::ostream& str = std::cout) const throw (MatrixException);

    // Matrix arithmetics operators
    MatrixGeneric<T>& operator=(const MatrixGeneric<T>& m) throw(MatrixException);
    MatrixGeneric<T>& operator=(const T& scalar) throw(MatrixException);
    MatrixGeneric<T>& operator+=(const MatrixGeneric<T>& m) throw(MatrixException);
    MatrixGeneric<T>& operator-=(const MatrixGeneric<T>& m) throw(MatrixException);
    MatrixGeneric<T>& operator*=(const MatrixGeneric<T>&m ) throw(MatrixException);
    MatrixGeneric<T>& operator*=(const T& scalar);
    MatrixGeneric<T>& operator+=(const T& scalar);
    MatrixGeneric<T>& operator-=(const T& scalar);
    MatrixGeneric<T>& operator/=(const T& scalar) throw(MatrixException);

    MatrixGeneric<T>& ewMult(const MatrixGeneric<T>& m) throw (MatrixException);
    MatrixGeneric<T>& ewDiv(const MatrixGeneric<T>& m) throw (MatrixException);

    // Transpose the matrix
    MatrixGeneric<T> transpose() const throw (MatrixException);
    MatrixGeneric<T>& transpose_() throw (MatrixException);
    MatrixGeneric<T>& conj_();
    MatrixGeneric<T> conj() const throw (MatrixException);

    // Round very small elements to 0
    MatrixGeneric<T>& roundSmallElements_();
    MatrixGeneric<T>& roundSmallElements_(const T& eps);
    MatrixGeneric<T> roundSmallElements() const throw(MatrixException);
    MatrixGeneric<T> roundSmallElements(const T& eps) const throw(MatrixException);

    // Insert or remove rows/columns.
    MatrixGeneric<T>& removeRow_(const size_t rowNr) throw (MatrixException);
    MatrixGeneric<T>& removeColumn_(const size_t colNr) throw (MatrixException);
    MatrixGeneric<T>& insertRow_(const size_t rowNr, const T& el = static_cast<T>(0)) throw (MatrixException);
    MatrixGeneric<T>& insertColumn_(const size_t colNr, const T& el = static_cast<T>(0)) throw (MatrixException);
    MatrixGeneric<T> removeRow(const size_t rowNr) const throw (MatrixException);
    MatrixGeneric<T> removeColumn(const size_t colNr) const throw (MatrixException);
    MatrixGeneric<T> insertRow(const size_t rowNr, const T& el = static_cast<T>(0)) const throw (MatrixException);
    MatrixGeneric<T> insertColumn(const size_t colNr, const T& el = static_cast<T>(0)) const throw (MatrixException);

    // Swap rows and columns
    MatrixGeneric<T>& swapRows_(const size_t r1, const size_t r2) throw(MatrixException);
    MatrixGeneric<T>& swapColumns_(const size_t c1, const size_t c2) throw(MatrixException);
    MatrixGeneric<T> swapRows(const size_t r1, const size_t r2) const throw(MatrixException);
    MatrixGeneric<T> swapColumns(const size_t c1, const size_t c2) const throw(MatrixException);

    // Triangular parts of this one:
    MatrixGeneric<T>& upperTriangularPart_(const bool inclDiag=true);
    MatrixGeneric<T>& lowerTriangularPart_(const bool inclDiag=true);
    MatrixGeneric<T>& diagPart_();
    MatrixGeneric<T> upperTriangularPart(const bool inclDiag=true) const throw(MatrixException);
    MatrixGeneric<T> lowerTriangularPart(const bool inclDiag=true) const throw(MatrixException);
    MatrixGeneric<T> diagPart() const throw(MatrixException);

    // These methods are only applicable for square matrices:
    MatrixGeneric<T>& setDiag_(const T& scalar) throw(MatrixException);
    MatrixGeneric<T>& setUnit_() throw(MatrixException);
    T determinant(const bool fullp = true) const throw(MatrixException);
    MatrixGeneric<T> inverse(const bool fullp = true) const throw(MatrixException);

    bool isSquare() const;
    size_t rank() const throw(MatrixException);

    // Destructor
    virtual ~MatrixGeneric();


    // Declaration of friend functions
    friend MatrixGeneric<T> (math::operator+ <>) (const MatrixGeneric<T>& m) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator- <>) (const MatrixGeneric<T>& m) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator+ <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator- <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator* <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator* <>) (const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator* <>) (const T& sc, const MatrixGeneric<T>& m) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator+ <>) (const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator+ <>) (const T& sc, const MatrixGeneric<T>& m) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator- <>) (const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator- <>) (const T& sc, const MatrixGeneric<T>& m) throw(MatrixException);
    friend MatrixGeneric<T> (math::operator/ <>) (const MatrixGeneric<T>& m, const T& sc) throw(MatrixException);
    friend MatrixGeneric<T> (math::matEwMult <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);
    friend MatrixGeneric<T> (math::matEwDiv <>) (const MatrixGeneric<T>& m1, const MatrixGeneric<T>& m2) throw(MatrixException);
    friend void (math::__matrixprivate::__matconj <>) (MatrixGeneric<T>& m);
    friend void (math::__matrixprivate::__matconj <>) (MatrixGeneric<std::complex<T> >& m);
};

// Matrices with elements of types float, double and long double
// make most sense so these three types are predefined
typedef MatrixGeneric<float>       FMatrix;
typedef MatrixGeneric<double>      Matrix;
typedef MatrixGeneric<long double> LDMatrix;

} // namespace math

// DEFINITION
#include "matrix/MatrixGeneric.cpp"

#endif	// _MATH_MATRIXGENERIC_HPP_
