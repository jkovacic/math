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
 * @file
 * @author Jernej Kovacic
 *
 * An internal header file, it should not be included directly.
 * @headername{PolynomialGeneric.h}
 *
 * Declaration of the class PolynomialGeneric, representing polynomials.
 */

#ifndef _MATH_POLYNOMIALGENERIC_HPP_
#define _MATH_POLYNOMIALGENERIC_HPP_

#include <vector>
#include <iostream>
#include <cstdlib>

#include "util/NumericUtil.hpp"
#include "exception/PolynomialException.hpp"


namespace math
{

// Forward declaration of the class is necessary...
template <typename F> class PolynomialGeneric;


// to declare the class's friend functions:
template <typename F>
PolynomialGeneric<F> operator+(const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator-(const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator*(const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator/(const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator%(const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator+(const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator-(const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator*(const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator/(const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator%(const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator*(const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator+(const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator-(const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator/(const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);

template <typename F>
PolynomialGeneric<F> operator%(const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);

// and its friend << operator:
template <typename F>
std::ostream& operator<<(std::ostream& output, const PolynomialGeneric<F>& q);


/**
 * A class representing monovariable polynomials and their basic operators.
 * 
 * Additional functionalities, such as polynomial derivation and integration,
 * are implemented.
 * 
 * Root finding algorithm is not implemented yet.
 */
template <typename F>
class PolynomialGeneric
{
private:
    /**
     * Coefficients of the polynomial:
     * p(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n
     */
    std::vector<F> coef;

    // A utility function that copies vector's coefficients
    void __copyCoefs(const std::vector<F>& cvect) throw (PolynomialException);

    // A utility function that reduces zero-coefficients from the highest terms
    void __reduce();

    // A utility function for polynomial division
    static void __polyDivision(const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2, PolynomialGeneric<F>* q, PolynomialGeneric<F>* rem) throw (PolynomialException);

    /*
     * @return whether the polynomial is of zero degree and its only coefficient equals 0
     */
    inline bool __isZero() const
    {
        return ( 1==this->coef.size() && true==math::NumericUtil::isZero<F>(this->coef.at(0)) );
    }


public:
    // Constructors
    PolynomialGeneric(const std::vector<F>& cvect) throw (PolynomialException);
    PolynomialGeneric(const F& c0 = static_cast<F>(0)) throw (PolynomialException);
    PolynomialGeneric(const F* carray, size_t n) throw (PolynomialException);
    PolynomialGeneric(bool ignored, size_t n = 1) throw (PolynomialException);
    PolynomialGeneric(const PolynomialGeneric<F>& poly) throw (PolynomialException);

    // Getters
    void get(std::vector<F>& vec) const throw (PolynomialException);
    void getDesc(std::vector<F>& vec) const throw (PolynomialException);
    F get(size_t pos) const;

    // Setters
    PolynomialGeneric<F>& set(const std::vector<F>& cvect) throw (PolynomialException);
    PolynomialGeneric<F>& setDesc(const std::vector<F>& cvect) throw (PolynomialException);
    PolynomialGeneric<F>& set(size_t pos, const F& c = static_cast<F>(0)) throw (PolynomialException);

    // Insert and remove coefficients
    PolynomialGeneric<F>& insert(size_t pos, const F& c) throw (PolynomialException);
    PolynomialGeneric<F>& remove(size_t pos);

    // Degree of the polynomial
    size_t degree() const;

    // Value of the polynomial for given x
    F value(const F& x) const;
    // Derivative of the polynomial

    PolynomialGeneric<F> deriv() const throw (PolynomialException);

    // Indefinite integral of the polynomial
    PolynomialGeneric<F> integ(const F& c = static_cast<F>(0)) const throw (PolynomialException);

    // Display the polynomial
    void display(char arg = 'x', std::ostream& str = std::cout) const;

    // Assignment operator
    PolynomialGeneric<F>& operator=(const PolynomialGeneric<F>& poly) throw (PolynomialException);

    // Operators
    PolynomialGeneric<F>& operator+=(const PolynomialGeneric<F>& poly) throw (PolynomialException);
    PolynomialGeneric<F>& operator-=(const PolynomialGeneric<F>& poly) throw (PolynomialException);
    PolynomialGeneric<F>& operator*=(const PolynomialGeneric<F>& poly) throw (PolynomialException);
    PolynomialGeneric<F>& operator/=(const PolynomialGeneric<F>& poly) throw (PolynomialException);
    PolynomialGeneric<F>& operator%=(const PolynomialGeneric<F>& poly) throw (PolynomialException);
    PolynomialGeneric<F> operator-() const throw (PolynomialException);
    PolynomialGeneric<F>& operator+=(const F& sc);
    PolynomialGeneric<F>& operator-=(const F& sc);
    PolynomialGeneric<F>& operator*=(const F& sc);
    PolynomialGeneric<F>& operator/=(const F& sc) throw (PolynomialException);
    PolynomialGeneric<F>& operator%=(const F& sc) throw (PolynomialException);

    // Destructor
    virtual ~PolynomialGeneric();


    // operators in separate functions declared as friend functions
    friend PolynomialGeneric<F> (math::operator+ <>) (const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator- <>) (const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator* <>) (const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator/ <>) (const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator% <>) (const PolynomialGeneric<F>& p1, const PolynomialGeneric<F>& p2) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator+ <>) (const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator- <>) (const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator* <>) (const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator/ <>) (const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator% <>) (const PolynomialGeneric<F>& poly, const F& sc) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator+ <>) (const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator- <>) (const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator* <>) (const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator/ <>) (const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);
    friend PolynomialGeneric<F> (math::operator% <>) (const F& sc, const PolynomialGeneric<F>& poly) throw(PolynomialException);

    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& (math::operator<< <>) (std::ostream& output, const PolynomialGeneric<F>& poly);

};  // class PolynomialGeneric

// Polynomials with elements of types float, double and long double
// make most sense so the following three types are predefined:
typedef PolynomialGeneric<float>       FPolynomial;
typedef PolynomialGeneric<double>      Polynomial;
typedef PolynomialGeneric<long double> LDPolynomial;


}  // namespace math


// DEFINITION
#include "polynomial/PolynomialGeneric.cpp"

#endif // _MATH_POLYNOMIALGENERIC_HPP_
