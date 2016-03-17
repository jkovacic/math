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
 * @headername{Rational.h}
 *
 * Declaration of the class RationalGeneric, representing mathematical
 * rational numbers, i.e. reduced fractions of two integers.
 */

#ifndef _MATH_RATIONALGENERIC_HPP_
#define _MATH_RATIONALGENERIC_HPP_

#include "util/NumericUtil.hpp"
#include "exception/RationalException.hpp"

#include <ostream>
#include <cstddef>
#include <string>

namespace math
{

// Forward declaration of the class is necessary...
template <typename I> class RationalGeneric;

// to declare the class's friend functions...
template <typename I>
std::ostream& operator<<(std::ostream& output, const RationalGeneric<I>& fraction);

// ... unary operators...
template <typename I>
RationalGeneric<I> operator+(const RationalGeneric<I>& f);

template <typename I>
RationalGeneric<I> operator-(const RationalGeneric<I>& f) throw(RationalException);

// ...and binary operators as friend functions:
template <typename I>
RationalGeneric<I> operator+(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2) throw(RationalException);

template <typename I>
RationalGeneric<I> operator+(const RationalGeneric<I>& f, const I& i) throw(RationalException);

template <typename I>
RationalGeneric<I> operator+(const I& i, const RationalGeneric<I>& f) throw(RationalException);

template <typename I>
RationalGeneric<I> operator-(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2) throw(RationalException);

template <typename I>
RationalGeneric<I> operator-(const RationalGeneric<I>& f, const I& i) throw(RationalException);

template <typename I>
RationalGeneric<I> operator-(const I& i, const RationalGeneric<I>& f) throw(RationalException);

template <typename I>
RationalGeneric<I> operator*(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2) throw(RationalException);

template <typename I>
RationalGeneric<I> operator*(const RationalGeneric<I>& f, const I& i) throw(RationalException);

template <typename I>
RationalGeneric<I> operator*(const I& i, const RationalGeneric<I>& f) throw(RationalException);

template <typename I>
RationalGeneric<I> operator/(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2) throw(RationalException);

template <typename I>
RationalGeneric<I> operator/(const RationalGeneric<I>& f, const I& i) throw(RationalException);

template <typename I>
RationalGeneric<I> operator/(const I& i, const RationalGeneric<I>& f) throw(RationalException);

// Comparison operators:
// - returns true if fractions' values are equal
template <typename I>
bool operator==(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2);

template <typename I>
bool operator==(const RationalGeneric<I>& f, const I& i);

template <typename I>
bool operator==(const I& i, const RationalGeneric<I>& f);

// - returns true if fractions' values are not equal
template <typename I>
bool operator!=(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2);

template <typename I>
bool operator!=(const RationalGeneric<I>& f, const I& i);

template <typename I>
bool operator!=(const I& i, const RationalGeneric<I>& f);

// - returns true if strictly less than fraction
template <typename I>
bool operator<(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2);

template <typename I>
bool operator<(const RationalGeneric<I>& f, const I& i);

template <typename I>
bool operator<(const I& i, const RationalGeneric<I>& f);

// - returns true if less than or equal to fraction
template <typename I>
bool operator<=(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2);

template <typename I>
bool operator<=(const RationalGeneric<I>& f, const I& i);

template <typename I>
bool operator<=(const I& i, const RationalGeneric<I>& f);

// - returns true if strictly greater than fraction
template <typename I>
bool operator>(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2);

template <typename I>
bool operator>(const RationalGeneric<I>& f, const I& i);

template <typename I>
bool operator>(const I& i, const RationalGeneric<I>& f);

// - returns true if greater than or equal to fraction
template <typename I>
bool operator>=(const RationalGeneric<I>& f1, const RationalGeneric<I>& f2);

template <typename I>
bool operator>=(const RationalGeneric<I>& f, const I& i);

template <typename I>
bool operator>=(const I& i, const RationalGeneric<I>& f);

template <typename I>
std::ostream& operator<<(std::ostream& output, const RationalGeneric<I>& fraction);


/**
 * @brief A class representing rational numbers, i.e. fractions of integer numbers.
 * 
 * Rationals are automatically reduced as soon as its nominator or denominator are modified.
 * 
 * I is expected to be an integer type, preferably signed and preferably shorter
 * than '(un)signed long long int'.
 */
template <typename I>
class RationalGeneric
{
    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& (math::operator<< <>) (std::ostream& output, const RationalGeneric<I>& fraction);


private:
    I m_num;                /// fraction's numerator
    I m_denom;              /// fraction's denominator (will be always assigned a positive value, cannot be 0)

public:

    //Constructor, assigns fraction's numerator and denominator
    RationalGeneric(
            const I& numerator = static_cast<I>(0),
            const I& denominator = static_cast<I>(1) ) 
        throw(RationalException);

    // Constructor from a string
    RationalGeneric(
            const std::string& str, 
            const std::size_t repSeqLen = 0 ) 
        throw (RationalException);

    // Copy constructor
    RationalGeneric(const RationalGeneric<I>& orig);
    
    // Destructor
    ~RationalGeneric();


    // Returns simplified fraction's numerator
    I getNumerator() const;

    // Returns simplified fraction's denominator
    I getDenominator() const;

    // Assigns fraction's numerator and denominator and simplifies the fraction
    RationalGeneric<I>& set(
            const I& numerator = static_cast<I>(0), 
            const I& denominator = static_cast<I>(1) ) 
        throw(RationalException);

    // Parses the fraction from its decimal representation
    RationalGeneric<I>& set(
            const std::string& str, 
            const std::size_t repSeqLen = 0 ) 
        throw (RationalException);

    // Outputs the fraction to std::cout, optionally multiplies both members by a factor
    void display(
            const I& factor = static_cast<I>(1), 
            std::ostream& str = std::cout) const;
    
    // Converts the fraction into the specified numeric type
    template <typename F>
    F toNum() const;

    // Returns true if the fraction equals 0
    bool isZero() const;

    // Returns true if the fraction's value is strictly greater than 0
    bool isPositive() const;

    // Returns true if the fraction's value is strictly less than zero
    bool isNegative() const;

    // Inverts the fraction
    RationalGeneric<I>& inv_() throw(RationalException);
    RationalGeneric<I> inv() const throw(RationalException);


    // Operator =, copies values of fraction's numerator and denominator
    RationalGeneric<I>& operator=(const RationalGeneric<I>& fraction);
    RationalGeneric<I>& operator=(const I& sc);

    // Binary operations are implemented as separate friend functions,
    // the following operators remain that actually modify instance of the class:
    RationalGeneric<I>& operator+=(const RationalGeneric<I>& fraction) throw(RationalException);
    RationalGeneric<I>& operator-=(const RationalGeneric<I>& fraction) throw(RationalException);
    RationalGeneric<I>& operator*=(const RationalGeneric<I>& fraction) throw(RationalException);
    RationalGeneric<I>& operator/=(const RationalGeneric<I>& fraction) throw(RationalException);


private:

    // tries to assign a rational from long long int arguments:
    void __setLL( const long long int numerator,
                  const long long int denominator ) 
            throw (RationalException);

    // Reduces the fraction (divides numerator and denominator by their greatest common divisor)
    void __reduce();

    // Friend functions that implement arithmetic operators:

    friend RationalGeneric<I> (math::operator+ <>) (
            const RationalGeneric<I>& f );

    friend RationalGeneric<I> (math::operator- <>) (
            const RationalGeneric<I>& f )
        throw(RationalException);

    friend RationalGeneric<I> (math::operator+ <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator+ <>) (
            const RationalGeneric<I>& f, 
            const I& i ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator+ <>) (
            const I& i, 
            const RationalGeneric<I>& f ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator- <>) (
            const RationalGeneric<I>& f1,
            const RationalGeneric<I>& f2 ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator- <>) (
            const RationalGeneric<I>& f,
            const I& i ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator- <>) (
            const I& i,
            const RationalGeneric<I>& f ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator* <>) (
            const RationalGeneric<I>& f1,
            const RationalGeneric<I>& f2 ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator* <>) (
            const RationalGeneric<I>& f,
            const I& i ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator* <>) (
            const I& i,
            const RationalGeneric<I>& f ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator/ <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator/ <>) (
            const RationalGeneric<I>& f, 
            const I& i ) 
        throw(RationalException);

    friend RationalGeneric<I> (math::operator/ <>) (
            const I& i, 
            const RationalGeneric<I>& f ) 
        throw(RationalException);

    friend bool (math::operator== <>) (
            const RationalGeneric<I>& f1,
            const RationalGeneric<I>& f2 );

    // Friend functions that implement comparison operators:
    friend bool (math::operator== <>) (
            const RationalGeneric<I>& f,
            const I& i );

    friend bool (math::operator== <>) (
            const I& i,
            const RationalGeneric<I>& f);

    friend bool (math::operator!= <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 );

    friend bool (math::operator!= <>) (
            const RationalGeneric<I>& f,
            const I& i );

    friend bool (math::operator!= <>) (
            const I& i,
            const RationalGeneric<I>& f );

    friend bool (math::operator< <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 );

    friend bool (math::operator< <>) (
            const RationalGeneric<I>& f,
            const I& i );

    friend bool (math::operator< <>) (
            const I& i,
            const RationalGeneric<I>& f );

    friend bool (math::operator<= <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 );

    friend bool (math::operator<= <>) (
            const RationalGeneric<I>& f,
            const I& i );

    friend bool (math::operator<= <>) (
            const I& i,
            const RationalGeneric<I>& f);

    friend bool (math::operator> <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 );

    friend bool (math::operator> <>) (
            const RationalGeneric<I>& f,
            const I& i );

    friend bool (math::operator> <>) (
            const I& i,
            const RationalGeneric<I>& f );

    friend bool (math::operator>= <>) (
            const RationalGeneric<I>& f1, 
            const RationalGeneric<I>& f2 );

    friend bool (math::operator>= <>) (
            const RationalGeneric<I>& f,
            const I& i );

    friend bool (math::operator>= <>) (
            const I& i,
            const RationalGeneric<I>& f );

};


typedef RationalGeneric<long int>        Rational;
typedef RationalGeneric<int>             IRational;
typedef RationalGeneric<short int>       SIRational;

} // namespace math


// DEFINITION
#include "rational/RationalGeneric.cpp"

#endif  // _MATH_RATIONALGENERIC_HPP_
