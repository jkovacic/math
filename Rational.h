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
@file Rational.h

Declaration of the class Rational, representing mathematical rational numbers,
i.e. reduced fractions of two integers

@author Jernej Kovacic
*/

#ifndef _RATIONAL_H_
#define	_RATIONAL_H_

#include "RationalException.h"

#include <iostream>

class Rational
{
    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& operator<<(std::ostream& output, const Rational& fraction);

private:
    int num;                /// fraction's numerator
    unsigned int denom;     /// fraction's denominator (will be always assigned a positive value, cannot be 0)

public:
    //Constructor, assigns fraction's numerator and denominator
    Rational(int numerator = 0, int denominator = 1) throw(RationalException);
    // Copy constructor
    Rational(const Rational& orig);
    // Destructor
    ~Rational();

    // Returns simplified fraction's numerator
    int getNumerator() const;
    // Returns simplified fraction's denominator
    unsigned int getDenominator() const;
    // Assigns fraction's numerator and denominator and simplifies the fraction
    Rational& set(int numerator = 0, int denominator = 1) throw(RationalException);

    // Outputs the fraction to std::cout, optionally multiplies both members by a factor
    void display(int factor = 1, std::ostream& str = std::cout) const;
    // Converts the fraction into its float value
    float toFloat() const;
    // Converts the fraction into its double value
    double toDouble() const;

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
    // typical mathematical operators
    Rational operator+(const Rational& fraction) const throw(RationalException);
    Rational& operator+=(const Rational& fraction) throw(RationalException);
    Rational operator-(const Rational& fraction) const throw(RationalException);
    Rational& operator-=(const Rational& fraction) throw(RationalException);
    Rational operator*(const Rational& fraction) const throw(RationalException);
    Rational& operator*=(const Rational& fraction) throw(RationalException);
    Rational operator/(const Rational& fraction) const throw(RationalException);
    Rational& operator/=(const Rational& fraction) throw(RationalException);
    // Unary operator -
    Rational operator-() const;

    // Comparison operators:
    // - returns true if fractions' values are equal
    bool operator==(const Rational& fraction) const;
    // - returns true if fractions' values are not equal
    bool operator!=(const Rational& fraction) const;
    // - returns true if strictly less than fraction
    bool operator<(const Rational& fraction) const;
    // - returns true if less than or equal to fraction
    bool operator<=(const Rational& fraction) const;
    // - returns true if strictly greater than fraction
    bool operator>(const Rational& fraction) const;
    // - returns true if greater than or equal to fraction
    bool operator>=(const Rational& fraction) const;

    // The following two functions are public and static as they
    // may also be useful elsewhere:
    
    //  The greatest common divisor of two integer values
    static unsigned int greatestCommonDivisor(unsigned int first, unsigned int second);
    // The least common multiple of two integer values
    static unsigned int leastCommonMultiple(unsigned int first, unsigned int second);

private:
    // Reduces the fraction (divides numerator and denominator by their greatest common divisor)
    void reduce();
    // Auxiliary function that calculates unreduced numerator of difference of two fraction
    // Only its sign actually matters, so it returns -1, 0 or 1
    int sign(const Rational& fraction) const;
    // Absolute value of an integer (just an auxiliary function for others)
    static unsigned int absolute(int a);
    static int auxSum(int num1, int denom2, int num2, int denom1) throw(RationalException);
    static int auxProd(int factor1, int factor2) throw(RationalException);
};

#endif	/* _RATIONAL_H_ */
