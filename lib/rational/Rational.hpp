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
 * Declaration of the class Rational, representing mathematical rational numbers,
 * i.e. reduced fractions of two integers.
 */

#ifndef _MATH_RATIONAL_HPP_
#define	_MATH_RATIONAL_HPP_

#include "exception/RationalException.hpp"

#include <ostream>
#include <string>

namespace math
{

// Advance declaration of the class is necessary...
class Rational;

// to declare the class's friend function...
std::ostream& operator<<(std::ostream& output, const Rational& fraction);

// ...and binary operators as friend functions:
Rational operator+(const Rational& f1, const Rational& f2) throw(RationalException);
Rational operator-(const Rational& f1, const Rational& f2) throw(RationalException);
Rational operator*(const Rational& f1, const Rational& f2) throw(RationalException);
Rational operator/(const Rational& f1, const Rational& f2) throw(RationalException);

// Comparison operators:
// - returns true if fractions' values are equal
bool operator==(const Rational& f1, const Rational& f2);
// - returns true if fractions' values are not equal
bool operator!=(const Rational& f1, const Rational& f2);
// - returns true if strictly less than fraction
bool operator<(const Rational& f1, const Rational& f2);
// - returns true if less than or equal to fraction
bool operator<=(const Rational& f1, const Rational& f2);
// - returns true if strictly greater than fraction
bool operator>(const Rational& f1, const Rational& f2);
// - returns true if greater than or equal to fraction
bool operator>=(const Rational& f1, const Rational& f2);



/**
 * @brief A class representing rational numbers, i.e. fractions of integer numbers.
 * 
 * Rationals are automatically reduced as soon as its nominator or denominator are modified.
 */
class Rational
{
    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& operator<<(std::ostream& output, const Rational& fraction);


private:
    long int num;                /// fraction's numerator
    unsigned long int denom;     /// fraction's denominator (will be always assigned a positive value, cannot be 0)

public:
    //Constructor, assigns fraction's numerator and denominator
    Rational(long int numerator = 0L, long int denominator = 1L) throw(RationalException);
    // Constructor from a string
    Rational(const std::string& str, unsigned int repSeqLen=0) throw (RationalException);
    // Copy constructor
    Rational(const Rational& orig);
    // Destructor
    ~Rational();

    // Returns simplified fraction's numerator
    long int getNumerator() const;
    // Returns simplified fraction's denominator
    unsigned long int getDenominator() const;
    // Assigns fraction's numerator and denominator and simplifies the fraction
    Rational& set(long int numerator = 0L, long int denominator = 1L) throw(RationalException);
    // Parses the fraction from its decimal representation
    Rational& set(const std::string& str, unsigned int repSeqLen=0) throw (RationalException);

    // Outputs the fraction to std::cout, optionally multiplies both members by a factor
    void display(long int factor = 1L, std::ostream& str = std::cout) const;
    // Converts the fraction into its float value
    float toFloat() const;
    // Converts the fraction into its double value
    double toDouble() const;
    // Converts the fraction into its long double value
    long double toLongDouble() const;

    // Returns true if the fraction equals 0
    bool isZero() const;
    // Returns true if the fraction's value is strictly greater than 0
    bool isPositive() const;
    // Returns true if the fraction's value is strictly less than zero
    bool isNegative() const;

    // Returns inverse value of the fraction
    Rational invert() const throw(RationalException);
    // Inverts (i.e. modifies) the fraction
    Rational& inverse() throw(RationalException);

    // Operator =, copies values of fraction's numerator and denominator
    Rational& operator=(const Rational& fraction);
    // Binary operations are implemented as separate friend functions,
    // the following operators remain that actually modify instance of the class:
    Rational& operator+=(const Rational& fraction) throw(RationalException);
    Rational& operator-=(const Rational& fraction) throw(RationalException);
    Rational& operator*=(const Rational& fraction) throw(RationalException);
    Rational& operator/=(const Rational& fraction) throw(RationalException);
    // Unary operator -
    Rational operator-() const;

private:
    // Reduces the fraction (divides numerator and denominator by their greatest common divisor)
    void reduce();
    // Auxiliary function that calculates unreduced numerator of difference of two fraction
    // Only its sign actually matters, so it returns -1, 0 or 1
    static short int sign(const Rational& f1, const Rational& f2);
    // Absolute value of an integer (just an auxiliary function for others)
    static unsigned long int absolute(long int a);
    static long int auxSum(long int num1, long int denom2, long int num2, long int denom1) throw(RationalException);
    static long int auxProd(long int factor1, long int factor2) throw(RationalException);
    // 10^n
    static unsigned long long int pow10(unsigned int n) throw (RationalException);
    // parses a string into a long long value
    static long long int str2ll(const std::string& str) throw (RationalException);

    // Friend functions that implement operators:
    friend Rational operator+(const Rational& f1, const Rational& f2) throw(RationalException);
    friend Rational operator-(const Rational& f1, const Rational& f2) throw(RationalException);
    friend Rational operator*(const Rational& f1, const Rational& f2) throw(RationalException);
    friend Rational operator/(const Rational& f1, const Rational& f2) throw(RationalException);

    friend bool operator==(const Rational& f1, const Rational& f2);
    friend bool operator!=(const Rational& f1, const Rational& f2);
    friend bool operator<(const Rational& f1, const Rational& f2);
    friend bool operator<=(const Rational& f1, const Rational& f2);
    friend bool operator>(const Rational& f1, const Rational& f2);
    friend bool operator>=(const Rational& f1, const Rational& f2);
};

} // namespace math

#endif	// _MATH_RATIONAL_HPP_