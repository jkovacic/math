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
* @file Rational.cpp
*
* Implementation of the class Rational
*
* @author Jernej Kovacic
*/


#include "Rational.h"

#include <limits.h>

/**
 * Constructor.
 *
 * Creates an instance of a fraction and assigns its numerator and denominator
 * as set by input parameters. The fraction will be automatically reduced to
 * relative prime numbers. Zero denominator is not permitted, any attempt to do
 * so will throw an exception. The numerator and denominator are not constant
 * and can be modified later.
 *
 * Rational(a) will assign a/1
 * Rational()  will assign 0/1 = 0
 *
 * The possibility with one input parameter ensures that int will be
 * automatically correctly 'casted' to a/1 when necessary.
 * 
 * @param numerator - fraction's numerator (default: 0)
 * @param denominator - fraction's denominator (default: 1)
 *
 * @throw RationalException - when attempting to set a zero denominator
 *
 * @see set
 */
Rational::Rational(int numerator, int denominator) throw(RationalException)
{
    // Just call a function that actually sets both members.
    // It will throw a RationalException if the denominator
    // is attempted to be set to 0
    set(numerator, denominator);
}

/**
 * Copy constructor.
 * Creates an instance of a fraction and copies numerator
 * and denominator from orig.
 *
 * @param orig - a valid instance of Rational
 */
Rational::Rational(const Rational& orig) : num(orig.num), denom(orig.denom)
{
    // num and denom are already set, nothing more to do here.
    // It is impossible to pass an invalid input parameter (orig.denom == 0)
    // so there is no need for any check. All fractions are automatically reduced
    // so no need to do this explicitly
}

/**
 * Returns fraction's numerator.
 * Note that setting functions automatically reduce the fraction so this function
 * might return a different value than set by Rational() or set()
 *
 * @return numerator
 *
 * @see set
 */
int Rational::getNumerator() const
{
    return this->num;
}

/**
 * Returns fraction's denominator.
 * Note that setting functions automatically reduce the fraction so this function
 * might return a different value than set by Rational() or set()
 *
 * @return denominator
 *
 * @see set
 */
unsigned int Rational::getDenominator() const
{
    return this->denom;
}

/**
 * Sets the fractions numerator and denominator.
 * It will also check if the denominator is zero (this is not permitted and
 * will throw an exception), it will reduce the fraction (divide the numerator
 * and denominator by their greatest common divisor) and make sure the denominator
 * is never negative.
 *
 * set(a,b) will assign a/b
 * set(a)   will assign a/1 = a
 * set()    will assign 0/1 = 0
 *
 * @param numerator (default: 0)
 * @param denominator (default: 1)
 *
 * @return reference to itself
 *
 * @throw RationalException when attempting to set a zero denominator
 */
Rational& Rational::set(int numerator, int denominator) throw(RationalException)
{
    // Check for zero denominator
    if ( 0 == denominator )
    {
        // this is forbidden so throw an exception
        throw RationalException(RationalException::ZERO_DENOMINATOR);
    }

    // The fraction will be valid, so assign given parameters to fraction's members
    this->num = numerator;

    // The denominator cannot be negative.
    if ( denominator > 0 )
    {
        this->denom = denominator;
    }
    else
    {
        // If this is the case, multiply both values by -1
        this->denom = -denominator;
        this->num = -num;
    }

    // and finally reduce the fraction
    reduce();

    return *this;
} // Rational::set

/**
 * Outputs the fraction to stdout in form 'num/den'
 *
 * @param factor (default: 1): optional, if set, both members will be multiplied by it
 * @param str (default cout): output stream, the fraction will be dislayed
 */
void Rational::display(int factor, std::ostream& str) const
{
    str << num*factor << '/' << denom*factor;
}

/**
 * @return fraction's value converted to float
 */
float Rational::toFloat() const
{
    return ( (float) num)/( (float) denom);
}

/**
 * @return fraction's value converted to double
 */
double Rational::toDouble() const
{
    return ( (double) num)/( (double) denom);
}

/**
 * Checks if the fraction's value equals zero.
 *
 * @return true if equal to zero, false otherwise
 */
bool Rational::isZero() const
{
    bool retVal = false;

    // It is enough to check if the numerator equals 0
    if ( 0 == num )
    {
        retVal = true;;
    }

    return retVal;
} // Rational::isZero

/**
 * Checks if the fraction's value is positive (strictly greater than zero).
 *
 * @return true if positive, false otherwise
 */
bool Rational::isPositive() const
{
    bool retVal = false;

    // the denominator is always positive (see set), so it is enough
    // to check if the numerator is greater than 0
    if ( num > 0 )
    {
        retVal = true;
    }

    return retVal;
}  // Rational::isPositive

/**
 * Checks if the fraction's value is negative (strictly less than zero)
 *
 * @return true if negative, false otherwise
 */
bool Rational::isNegative() const
{
    bool retVal = false;

    // the denominator is always positive (see set), so it is enough
    // to check if the numerator is less than 0
    if ( num < 0 )
    {
        retVal = true;
    }

    return retVal;
} // Rational::isNegative

/**
 * Returns an inverse value (i.e. denom/num) of the fraction.
 * If the fraction's value is zero, an exception will be thrown.
 *
 * @return fraction's inverse fraction
 *
 * @throw RationalException the fraction's value is equal to zero
 */
Rational Rational::invert() const throw(RationalException)
{
    // Check if the numerator is equal to 0
    if ( 0 == num )
    {
        // In this case, inversion cannot be done, throw an exception
        throw RationalException(RationalException::UNINVERTIBLE);
    }

    // inversion will be possible, return a fraction with with exchanged
    // numerator and denominator
    return Rational(denom, num);
} //Rational::invert

/**
 * Modifies the fraction with its inverse value (swaps the numerator and denominator).
 * If the numerator equals zero, this is not possible so an exception will be thrown.
 *
 * @return reference to itself
 * 
 * @throw RationalException if attempting to invert a fraction whose numerator equals zero
 */
Rational& Rational::inverse() throw(RationalException)
{
    // Check if inversion is possible (num!=0)
    if ( 0 == num )
    {
        // inversion is impossible, throw an exception
        throw RationalException(RationalException::UNINVERTIBLE);
    }

    // inversion is possible, swap the fractions num. and denom.
    set(denom, num);

    return *this;
} // Rational::inverse

/**
 * Assignment operator, copies frac's numerator and denominator.
 *
 * @return reference to itself
 */
Rational& Rational::operator=(const Rational& frac)
{
    // The function's behaviour is similar to the copy constructor.
    // With one exception.
    // Nothing to do when attempting to assign it to itself.
    if ( this != &frac )
    {
        // otherwise copy the num and denom
        this->num = frac.num;
        this->denom = frac.denom;
    }

    return *this;
} // Rational::operator=

/**
 * Binary addition operator (+) for addition of two fractions
 *
 * @param frac - a fraction to be added to this
 *
 * @return this + frac
 */
Rational Rational::operator+(const Rational& frac) const throw(RationalException)
{
    // Unreduced sum of two fractions:
    //
    //   a     c     a*d + c*b
    //  --- + --- = -----------
    //   b     d       b * d

    // both fractions are valid which always results in a valid sum, so
    // no check is necesseary. Integer overflow is possible but it is
    // not checked right now.

    // Just construct an unreduced fraction as shown above:
    int numerator;
    int denominator;

    try
    {
        numerator = auxSum(this->num, frac.denom, frac.num, this->denom);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    return Rational(numerator, denominator);
    // the constructor will reduce the result
}  // Rational::operator+

/**
 * Addition operator (+=) that adds frac to itself and assigns the resulting value to itself.
 * a += b  ---> a -> a + b
 *
 * @param frac - a fraction to be added to this
 *
 * @return reference to itself
 */
Rational& Rational::operator+=(const Rational& frac) throw(RationalException)
{
    // See definition of fraction addition in operator+()
    // Result will be assigned to itself so use set

    int numerator;
    int denominator;

    try
    {
        numerator = auxSum(this->num, frac.denom, frac.num, this->denom);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    set(numerator, denominator);
    // set will reduce the result

    return *this;
}  // Rational::operator+=

/**
 * Binary subtraction operator (-) for substraction of two fractions
 *
 * @param frac - a fraction to be substracted from this
 *
 * @return this - frac
 */
Rational Rational::operator-(const Rational& frac) const throw(RationalException)
{
    // Unreduced difference of two fractions:
    //
    //   a     c     a*d - c*b
    //  --- - --- = -----------
    //   b     d       b * d

    // both fractions are valid which always results in a valid difference, so
    // no check is necesseary.

    int numerator;
    int denominator;

    try
    {
        numerator = auxSum(this->num, frac.denom, -frac.num, this->denom);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    return Rational(numerator, denominator);
    // constructor will automatically reduce the result

    // NOTE: it would also be possible to multiply frac by -1 and add it to this
    // (using operator+) but it would require a bit more operations.
}  // Rational::operator-

/**
 * Subtraction operator (-=) that subtracts frac from this and assigns the resulting value to itself.
 * a -= b  ---> a -> a - b
 *
 * @param frac - a fraction to be subtracted from this
 *
 * @return reference to itself
 */
Rational& Rational::operator-=(const Rational& frac) throw(RationalException)
{
    // See definition of fraction subtraction in operator-
    // Result will be assigned to itself so use set().

    int numerator;
    int denominator;

    try
    {
        numerator = auxSum(this->num, frac.denom, -frac.num, this->denom);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }
    
    set(numerator, denominator);
    // set will reduce the result

    return *this;
}  // Rational::operator-=

/**
 * Binary multiplication operator (*) for multiplication of two fractions.
 *
 * @param frac - a fraction, this will be multilied by
 *
 * @return this * frac
 */
Rational Rational::operator*(const Rational& frac) const throw (RationalException)
{
    // Unreduced product of two fractions:
    //
    //   a     c      a * c
    //  --- * --- = ---------
    //   b     d      b * d

    // Multiplication of two valid fractions always results in a valid
    // fraction, so no need for further checks.

    int numerator;
    int denominator;

    try
    {
        numerator = auxProd(this->num, frac.num);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( RationalException ex )
    {
        // The only possible exception is int overflow
        throw ex;
    }

    return Rational(numerator, denominator);
    // constructor will reduce the result
} // Rational::operator*

/**
 * Multiplication operator (*=) that multiplies frac and this and assigns the resulting value to itself.
 * a *= b  ---> a -> a * b
 *
 * @param frac - a fraction, this will be mltiplied by
 *
 * @return reference to itself
 */
Rational& Rational::operator*=(const Rational& frac) throw(RationalException)
{
    // See definition of fraction multiplication in operator*.
    // Result will be assigned to itself so use set().

    int numerator;
    int denominator;

    try
    {
        numerator = auxProd(this->num, frac.num);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    set(numerator, denominator);

    return *this;
} // Rational::operator*=

/**
 * Binary division operator (/) for division of two fractions.
 * Division by zero is not permitted and will throw an exception.
 *
 * @param frac - a fraction, this will be divided by
 *
 * @return this / frac
 *
 * @throw RationalException when attempting to divide by zero
 */
Rational Rational::operator/(const Rational& frac) const throw(RationalException)
{
    // Unreduced quotient of two fractions:
    //
    //   a     c     a     d      a * d
    //  --- / --- = --- * --- = ---------
    //   b     d     b     c      b * c

    // Check if frac's numerator equals 0
    if ( true == frac.isZero() )
    {
        throw RationalException(RationalException::DIVIDE_BY_ZERO);
    }

    int numerator;
    int denominator;

    try
    {
        numerator = auxProd(this->num, frac.denom);
        denominator = auxProd(frac.num, this->denom);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    return Rational(numerator, denominator);
} // Rational::operator/

/**
 * Division operator (/=) that divides this by frac and assigns the resulting value to itself.
 * a /= b  ---> a -> a / b
 * Division by zero is not permitted and will throw an exception.
 *
 * @param frac - a fraction, this will be divided by
 *
 * @return reference to itself
 *
 * @throw RationalException when attempting to divide by zero
 */
Rational& Rational::operator/=(const Rational& frac) throw(RationalException)
{
    // See definition of fraction division in operator/.
    // Check if frac's numerator equals 0
    if ( true == frac.isZero() )
    {
        throw RationalException(RationalException::DIVIDE_BY_ZERO);
    }

    // Result will be assigned to itself so use set()
    int numerator;
    int denominator;

    try
    {
        numerator = auxProd(this->num, frac.denom);
        denominator = auxProd(this->denom, frac.num);
    }
    catch ( RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    set(numerator, denominator);

    return *this;
} // Rational operator/=

/**
 * Unary negation operator (-)
 *
 * @return -this
 */
Rational Rational::operator-() const
{
    // One (doesn't matter which) member of the fraction must be negated.
    // Let it be the numerator
    return Rational(-this->num, this->denom);
}

/**
 * Comparison operator 'equal' (==)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return true if both fractions are equal, false otherwise
 */
bool Rational::operator==(const Rational& frac) const
{
    bool retVal = false;
    if ( 0 == sign(frac) )
    {
        retVal = true;
    }

    return retVal;
}

/**
 * Comparison operator 'not equal' (!=)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return true, if fractions are not equal, false if they are equal
 */
bool Rational::operator!=(const Rational& frac) const
{
    bool retVal = false;
    if ( 0 != sign(frac) )
    {
        retVal = true;
    }

    return retVal;
}

/**
 * Comparison operator 'greater than' (>)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return true if this is strictly greater than frac, false otherwise
 */
bool Rational::operator>(const Rational& frac) const
{
    bool retVal = false;
    if ( 0 < sign(frac) )
    {
        retVal = true;
    }

    return retVal;
}

/**
 * Comparison operator 'greater or equal' (>=)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return true if this is greater than or equal to frac, false otherwise
 */
bool Rational::operator>=(const Rational& frac) const
{
    bool retVal = false;
    if ( 0 <= sign(frac) )
    {
        retVal = true;
    }

    return retVal;
}

/**
 * Comparison operator 'less than' (<)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return true if this is strictly less than frac, false otherwise
 */
bool Rational::operator<(const Rational& frac) const
{
    bool retVal = false;
    if ( 0 > sign(frac) )
    {
        retVal = true;
    }

    return retVal;
}

/**
 * Comparison operator 'less or equal' (<=)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return - true if this is less than or equal to frac, false otherwise
 */
bool Rational::operator<=(const Rational& frac) const
{
    bool retVal = false;
    if ( 0 >= sign(frac) )
    {
        retVal = true;
    }

    return retVal;
}

/**
 * A friend function that outputs the reduced fraction to output stream
 * in form num/denom
 *
 * @param output - stream to write to
 * @param frac - fraction to be displayed
 *
 * @return reference of output stream (same as output)
 */
std::ostream& operator<<(std::ostream& output, const Rational& frac)
{
    // output the num and denom into the stream
    output << frac.num << '/' << frac.denom;
    // and return reference of the stream
    return output;
}

/**
 * An auxiliary function that reduces the fraction:
 * divides the numerator and denominator by their greatest common divisor
 */
void Rational::reduce()
{
    // the GCD is not defined when numerator is equal to 0.
    // In that case just assign denominator to 1
    // (it can be assigned to any integer value except 0)
    if ( 0 == num )
    {
        denom = 1;
        return;
    }

    // The GCD algorithm requires both integers to be positive.
    // The numerator is allowed to be negative, so get its absolute value.
    const unsigned int absNum = absolute(num);  // absolute value of num

    // finally obtain the greatest common divisor...
    const unsigned int gcd = greatestCommonDivisor(absNum, denom);

    // ... and divide both members by it.
    // if both num (handled a few lines above)and denom (not permitted when seting)
    // are different than 0, the GCD is guaranteed to be a non-zero value
    num /= (int) gcd;
    denom /= gcd;
} // Rational::reduce

/**
 * A simple integer implementation of abs
 *
 * @param a
 * @return absolute value of a
 */
unsigned int Rational::absolute(int a)
{
    return ( a>=0 ? a : -a );
}

/**
 * Greatest common divisor of first and second, calculated using the Euclidean algorithm.
 *
 * @param first
 * @param second
 *
 * @return zero if any of the two input parameters equals zero, GCD otherwise
 */
unsigned int Rational::greatestCommonDivisor(unsigned int first, unsigned int second)
{
    // A well known Euclidean algorithm is utilized to find the greatest common divisor.
    // It is known to be efficient, more details about it are available, for instance,
    // at http://en.wikipedia.org/wiki/Euclidean_algorithm


    // If any of both arguments is 0, the algorithm may "end up" in an infinite loop
    // or division by zero can occur. In such a case return 0 immediately.
    if ( 0 == first || 0 == second )
    {
        return 0;
    }

    // Now it is guaranteed to converge towards the GCD (or 1)
    unsigned int a = first;
    unsigned int b = second;
    unsigned int t;

    while ( 0 != b )
    {
        t = b;
        b = a % b;
        a = t;
    }  // while

    return a;
}  // Rational::gretestCommonDivisor

/**
 * Least common multiple of first and second.
 * It only makes sense when both arguments are not zero.
 *
 * @param first
 * @param second
 *
 * @return zero if any of the two input arguments equals zero, the LCM otherwise
 */
unsigned int Rational::leastCommonMultiple(unsigned int first, unsigned int second)
{
    // The following mathematical relation can be proven easily:
    // LCM(a,b) * GCD(a,b) = a * b
    // So, knowing that an efficient algorithm for the GCD exists, the LCM can be derived:
    // LCM(a,b) = (a*b)/GCD(a,b)

    // If any of the two input values equals zero, neither GCD
    // nor LCM doesn't make any sense, return 0 in such cases
    if ( 0 == first || 0 == second )
    {
        return 0;
    }

    // At this point the GCD will be no less than 1, definitely not 0
    const unsigned int gcd = greatestCommonDivisor(first, second);

    // so it's safe to divide
    return (first*second)/gcd;
} // Rational::leastCommonMultiple

/**
 * An auxiliary function that calculates (unreduced) numerator of (this-frac).
 * It is only needed by comparison operators who actually only need to know its sign
 *
 * @param fraction
 * 
 * @return -1 if frac is greater, 0 if they are equal, 1 if this is greater
 */
int Rational::sign(const Rational& frac) const
{
    // One approach to compare two objects is to calculate their difference
    // and check its sign. As denominators are always positive, it is
    // sufficient to calculate just the numerator of the difference.
    // See operator- for definition of fraction difference.

    int retVal = 0;

    // long prevents theoretically possible int overflow
    const long int diff = this->num * frac.denom - frac.num * this->denom;

    if ( diff > 0 )
    {
        retVal = 1;
    }
    else if ( diff < 0 )
    {
        retVal = -1;
    }
    else  // diff == 0
    {
        retVal = 0;
    }

    return retVal;
}

/**
 * Auxiliary function to calculate a sum of two products.
 * It checks for integer overflows and throws an exception in this case.
 *
 * The result is actually an unreduced numerator of a sum of two fractions,
 * which is reflected by "weird" parameter names.
 *
 * @param num1
 * @param denom2
 * @param num2
 * @param denom1
 *
 * @return num1*denom2+num2*denom1
 *
 * @throw RationalException in case of an integer overflow
 */
int Rational::auxSum(int num1, int denom2, int num2, int denom1) throw(RationalException)
{
    const long int sum = num1 * denom2 + num2 * denom1;

    if ( sum > INT_MAX || sum < INT_MIN )
    {
        throw RationalException(RationalException::OVERFLOW);
    }

    return (int) sum;
}

/**
 * Auxiliary function to calcuate a product of two integer numbers.
 * It checks for integer overflows and throws an exception in this case.
 *
 * @param first
 * @param second
 *
 * @return first * second
 *
 * @throw RationalException in case of an integer overflow
 */
int Rational::auxProd(int first, int second) throw(RationalException)
{
    const long int prod = first * second;

    if ( prod > INT_MAX || prod < INT_MIN )
    {
        throw RationalException(RationalException::OVERFLOW);
    }

    return (int) prod;
}

/**
 * Destructor
 */
Rational::~Rational()
{
    // nothing to do, i.e. no dynamically allocated memory to free,
    // no resources to release, etc.
}
