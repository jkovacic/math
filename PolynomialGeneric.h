/*
Copyright 2013, Jernej Kovacic

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
 * @file PolynomialGeneric.h
 *
 * Declaration of the class PolynomialGeneric, representing polynomials.
 *
 * @author Jernej Kovacic
 */

#ifndef _MATH_POLYNOMIALGENERIC_H_
#define _MATH_POLYNOMIALGENERIC_H_

#include <vector>
#include <iostream>

#include "NumericUtil.h"
#include "PolynomialException.h"


namespace math
{

// Advance declaration of the class is necessary...
template<class T> class PolynomialGeneric;


// to declare the class's friend functions:
template<class T>
PolynomialGeneric<T> operator+(const PolynomialGeneric<T>& p1, const PolynomialGeneric<T>& p2) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator-(const PolynomialGeneric<T>& p1, const PolynomialGeneric<T>& p2) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator*(const PolynomialGeneric<T>& p1, const PolynomialGeneric<T>& p2) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator+(const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator-(const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator*(const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator/(const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator%(const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator*(const T& sc, const PolynomialGeneric<T>& poly) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator+(const T& sc, const PolynomialGeneric<T>& poly) throw(PolynomialException);

template<class T>
PolynomialGeneric<T> operator-(const T& sc, const PolynomialGeneric<T>& poly) throw(PolynomialException);

// and its friend << operator:
template<class T>
std::ostream& operator<<(std::ostream& output, const PolynomialGeneric<T>& q);


/**
 * A class representing monovariable polynomials and their basic operators.
 * 
 * Additional functionalities, such as polynomial derivation and integration,
 * are implemented.
 * 
 * Root finding algorithm is not implemented yet.
 */
template<class T>
class PolynomialGeneric
{
private:
    /**
     * Coefficients of the polynomial:
     * p(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n
     */
    std::vector<T> coef;

    // A utility function that copies vector's coefficients
    void copyCoefs(const std::vector<T>& cvect) throw (PolynomialException);
    // A utility function that reduces zero-coeeficients from the highest terms
    void reduce();

public:
    // Constructors
    PolynomialGeneric(std::vector<T> cvect) throw (PolynomialException);
    PolynomialGeneric(const T& c0=math::NumericUtil<T>::ZERO) throw (PolynomialException);
    PolynomialGeneric(const T* carray, size_t n) throw (PolynomialException);
    PolynomialGeneric(bool ignored, size_t n = 1) throw (PolynomialException);
    PolynomialGeneric(const PolynomialGeneric<T>& poly) throw (PolynomialException);

    // Getters
    std::vector<T> get() const throw (PolynomialException);
    std::vector<T> getDesc() const throw (PolynomialException);
    T get(size_t pos) const;

    // Setters
    PolynomialGeneric<T>& set(const std::vector<T>& cvect) throw (PolynomialException);
    PolynomialGeneric<T>& setDesc(const std::vector<T>& cvect) throw (PolynomialException);
    PolynomialGeneric<T>& set(size_t pos, const T& c = NumericUtil<T>::ZERO) throw (PolynomialException);

    // Insert and remove coefficients
    PolynomialGeneric<T>& insert(size_t pos, const T& c) throw (PolynomialException);
    PolynomialGeneric<T>& remove(size_t pos);

    // Degree of the polynomial
    size_t degree() const;
    // Value of the polynomial for given x
    T value(const T& x) const;
    // Derivative of the polynomial
    PolynomialGeneric<T> deriv() const throw (PolynomialException);
    // Indefinite integral of the polynomial
    PolynomialGeneric<T> integ(const T& c = NumericUtil<T>::ZERO) const throw (PolynomialException);

    // Display the polynomial
    void display(char arg = 'x', std::ostream& str = std::cout) const;

    // Assignment operator
    PolynomialGeneric<T>& operator=(const PolynomialGeneric<T>& poly) throw (PolynomialException);

    // Operators
    PolynomialGeneric<T>& operator+=(const PolynomialGeneric<T>& poly) throw (PolynomialException);
    PolynomialGeneric<T>& operator-=(const PolynomialGeneric<T>& poly) throw (PolynomialException);
    PolynomialGeneric<T>& operator*=(const PolynomialGeneric<T>& poly) throw (PolynomialException);
    PolynomialGeneric<T> operator-() const throw (PolynomialException);
    PolynomialGeneric<T>& operator+=(const T& sc);
    PolynomialGeneric<T>& operator-=(const T& sc);
    PolynomialGeneric<T>& operator*=(const T& sc);
    PolynomialGeneric<T>& operator/=(const T& sc) throw (PolynomialException);
    PolynomialGeneric<T>& operator%=(const T& sc) throw (PolynomialException);

    // Destructor
    virtual ~PolynomialGeneric();


    // operators in separate functions declared as friend functions
    friend PolynomialGeneric<T> (math::operator+ <>) (const PolynomialGeneric<T>& p1, const PolynomialGeneric<T>& p2) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator- <>) (const PolynomialGeneric<T>& p1, const PolynomialGeneric<T>& p2) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator* <>) (const PolynomialGeneric<T>& p1, const PolynomialGeneric<T>& p2) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator+ <>) (const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator- <>) (const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator* <>) (const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator/ <>) (const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator% <>) (const PolynomialGeneric<T>& poly, const T& sc) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator+ <>) (const T& sc, const PolynomialGeneric<T>& poly) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator- <>) (const T& sc, const PolynomialGeneric<T>& poly) throw(PolynomialException);
    friend PolynomialGeneric<T> (math::operator* <>) (const T& sc, const PolynomialGeneric<T>& poly) throw(PolynomialException);

    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& (math::operator<< <>) (std::ostream& output, const PolynomialGeneric<T>& poly);

};  // class PolynomialGeneric

// Polynomials with elements of types float and double make most sense
// so the following two types are predefined:
typedef PolynomialGeneric<float> FPolynomial;
typedef PolynomialGeneric<double> Polynomial;

}  // namespace math


// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "PolynomialGeneric.cpp"

#endif // _MATH_POLYNOMIALGENERIC_H_
