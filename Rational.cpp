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

#include <climits>
#include <limits>
#include <cstdio>
#include <new>


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
math::Rational::Rational(int numerator, int denominator) throw(math::RationalException)
{
    // Just call a function that actually sets both members.
    // It will throw a RationalException if the denominator
    // is attempted to be set to 0
    set(numerator, denominator);
}

/**
 * Constructor from a decimal representation in a string.
 *
 * @see set(std::string, unsigned int) for more details.
 * 
 * @param str - decimal representation of a fraction
 * @param repSeqLen - length of the repeating sequence if applicable (default: 0)
 *
 * @throw RationalException if input argument is invalid
 */
math::Rational::Rational(const std::string& str, unsigned int repSeqLen) throw (math::RationalException)
{
    // Just call a function that actually sets both members.
    // It will throw a RationalException if 'str' is invalid etc.
    set(str, repSeqLen);
}

/**
 * Copy constructor.
 * Creates an instance of a fraction and copies numerator
 * and denominator from orig.
 *
 * @param orig - a valid instance of Rational
 */
math::Rational::Rational(const math::Rational& orig) : num(orig.num), denom(orig.denom)
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
int math::Rational::getNumerator() const
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
unsigned int math::Rational::getDenominator() const
{
    return this->denom;
}

/**
 * Sets the fraction's numerator and denominator.
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
math::Rational& math::Rational::set(int numerator, int denominator) throw(math::RationalException)
{
    // Check for zero denominator
    if ( 0 == denominator )
    {
        // this is forbidden so throw an exception
        throw math::RationalException(math::RationalException::ZERO_DENOMINATOR);
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
 * Parses a string with a decimal number into a fraction.
 * 
 * Each rational number can be represented by a decimal number, either with a
 * finite number of digits or ending with a finite length repeating sequence of digits.
 *
 * @note Only decimal system is supported. The string may start with a sign
 *       ('+' or '-'), it may contain no more than one decimal point (either '.' or ','),
 *       all other characters must be decimal digits (from '0' to '9'). Scientific
 *       and engineering formats are not supported. Any separators of thousands
 *       in the integer part are also not allowed. Additionally, 'repSeqLen' must not
 *       include the decimal point, in other words, the whole repeating sequence
 *       must comprise the fractional part only.
 *   
 *                       __
 * For instance, for 12.345 (45 is the repeating sequence), enter:
 *   set("12.345", 2) 
 * 
 * @param str - string to be parsed
 * @param repSeqLen - length of the repeating sequence if applicable (default: 0)
 *
 * @return reference to itself
 *
 * @throw RationalException if input argument is invalid
 */
math::Rational& math::Rational::set(const std::string& str, unsigned int repSeqLen) throw (math::RationalException)
{
    const unsigned int LEN = str.length();
    unsigned int decPoint = LEN; // position of the decimal point. LEN represents no decimal point.
    bool hasDigits = false; // does the string contain at least one decimal digit

    // Check if 'str' is a proper decimal representation
    for ( unsigned int i=0; i<LEN; i++ )
    {
        const char ch = str.at(i);
        
        if ('+' == ch || '-' == ch)
        {
            // the first (and only the first) character may be a sign
            if ( i>0 )
            {
                throw math::RationalException(math::RationalException::INVALID_INPUT);
            }
        }
        else if ( '.' == ch || ',' == ch )
        {
            // only one decimal point is allowed
            if ( decPoint < LEN )
            {
                // if less than LEN, a decimal point already exists
                throw math::RationalException(math::RationalException::INVALID_INPUT);
            }

            decPoint = i;
        }
        else if ( ch>='0' && ch <='9' )
        {
            // at least one digit found
            hasDigits = true;
        }
        else
        {
            // an invalid character
            throw math::RationalException(math::RationalException::INVALID_INPUT);
        }
    } // for i

    if ( false==hasDigits || (decPoint<LEN && repSeqLen>(LEN-decPoint-1)) )
    {
        throw math::RationalException(math::RationalException::INVALID_INPUT);
    }

    // the string does represent a decimal number, parse it
    try
    {
        // Note: if RationalException is thrown by str2ll or pow10, it will not be
        // caught by this try, instead it will be thrown to the caller. 
        
        std::string buf = str;
        
        /*
         * Temporary terms are declared as longer integers than supported.
         * This still allows the possibility that their difference will
         * probably still be within supported ranges. Of course the differences'
         * values are carefully checked and casted to integers of appropriate length
         * prior to passing them to set(int, unsigned int).
         * 
         * Sensible default values are also assigned to the temporary variables.
         */
        long long int num1 = 0LL;
        long long int num2 = 0LL;
        unsigned int den1 = 1;
        unsigned int den2 = 0;

        if ( 0 == repSeqLen )
        {
            /*
             * No repeating sequence.
             * In that case, simply "multiply" the number by a power of 10
             * (i.e. remove the decimal point from the string) and set the 
             * denominator to that power of 10.
             */
            if ( LEN != decPoint )
            {
                buf.erase(decPoint, 1);
            }

            num1 = str2ll(buf);
            den1 = ( LEN==decPoint ? 1 : pow10(LEN-decPoint-1) );
        }
        else
        {
            /*
             * Repeating sequence.
             * "Multiply" the number (denoted by x) by a power of 10 (10^n) to 
             * eliminate the decimal point and also include one repeating sequence.
             * Then remove the repeating sequence which is equal to multiplying
             * the original number by 10^(n-repSeqLen). When the second product is
             * subtracted from the first one, all repeating sequences except one are
             * eliminated and the fraction can be simply set.
             * 
             *                 ___
             * Example: x = 3.5167
             *                           ___
             *          10^4 * x = 35167.167
             *                      ___
             *          10 * x = 35.167
             * 
             *     (10^4-10)*x = (35167-35)
             * 
             *                   35132
             *              x = -------
             *                   9990   
             */
            buf.erase(decPoint, 1);
            num1 = str2ll(buf);
            buf.erase(LEN-1-repSeqLen, repSeqLen);
            num2 = str2ll(buf);  
            
            den1 = pow10(LEN-1-decPoint);
            den2 = pow10(LEN-1-decPoint-repSeqLen);
        }
        
        const long long int dnum = num1 - num2;
        const unsigned long int dden = den1 - den2;

        if ( dnum<INT_MIN || dnum>INT_MAX || dden>UINT_MAX )
        {
            throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
        }

        set(static_cast<int>(dnum), static_cast<unsigned int>(dden));
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::RationalException(math::RationalException::OUT_OF_MEMORY);
    }

    return *this;
}

/**
 * Outputs the fraction to stdout in form 'num/den'
 *
 * @param factor (default: 1): optional, if set, both members will be multiplied by it
 * @param str (default cout): output stream, the fraction will be displayed
 */
void math::Rational::display(int factor, std::ostream& str) const
{
    str << num*factor << '/' << denom*factor;
}

/**
 * @return fraction's value converted to float
 */
float math::Rational::toFloat() const
{
    return ( static_cast<float>(num)/static_cast<float>(denom) );
}

/**
 * @return fraction's value converted to double
 */
double math::Rational::toDouble() const
{
    return ( static_cast<double>(num)/static_cast<double>(denom) );
}

/**
 * @return fraction's value converted to long double
 */
long double math::Rational::toLongDouble() const
{
    return ( static_cast<long double>(num)/static_cast<long double>(denom) );
}

/**
 * Checks if the fraction's value equals zero.
 *
 * @return true if equal to zero, false otherwise
 */
bool math::Rational::isZero() const
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
bool math::Rational::isPositive() const
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
bool math::Rational::isNegative() const
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
math::Rational math::Rational::invert() const throw(math::RationalException)
{
    // Check if the numerator is equal to 0
    if ( 0 == num )
    {
        // In this case, inversion cannot be done, throw an exception
        throw math::RationalException(math::RationalException::UNINVERTIBLE);
    }

    // inversion will be possible, return a fraction with with exchanged
    // numerator and denominator
    return math::Rational(denom, num);
} //Rational::invert

/**
 * Modifies the fraction with its inverse value (swaps the numerator and denominator).
 * If the numerator equals zero, this is not possible so an exception will be thrown.
 *
 * @return reference to itself
 *
 * @throw RationalException if attempting to invert a fraction whose numerator equals zero
 */
math::Rational& math::Rational::inverse() throw(math::RationalException)
{
    // Check if inversion is possible (num!=0)
    if ( 0 == num )
    {
        // inversion is impossible, throw an exception
        throw math::RationalException(math::RationalException::UNINVERTIBLE);
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
math::Rational& math::Rational::operator=(const math::Rational& frac)
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
math::Rational math::Rational::operator+(const math::Rational& frac) const throw(math::RationalException)
{
    // Unreduced sum of two fractions:
    //
    //   a     c     a*d + c*b
    //  --- + --- = -----------
    //   b     d       b * d

    // both fractions are valid which always results in a valid sum, so
    // no check is necessary. Integer overflow is possible but it is
    // not checked right now.

    // Just construct an unreduced fraction as shown above:
    int numerator;
    int denominator;

    try
    {
        numerator = auxSum(this->num, frac.denom, frac.num, this->denom);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( math::RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    return math::Rational(numerator, denominator);
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
math::Rational& math::Rational::operator+=(const math::Rational& frac) throw(math::RationalException)
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
    catch ( math::RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    set(numerator, denominator);
    // set will reduce the result

    return *this;
}  // Rational::operator+=

/**
 * Binary subtraction operator (-) for subtraction of two fractions
 *
 * @param frac - a fraction to be subtracted from this
 *
 * @return this - frac
 */
math::Rational math::Rational::operator-(const math::Rational& frac) const throw(math::RationalException)
{
    // Unreduced difference of two fractions:
    //
    //   a     c     a*d - c*b
    //  --- - --- = -----------
    //   b     d       b * d

    // both fractions are valid which always results in a valid difference, so
    // no check is necessary.

    int numerator;
    int denominator;

    try
    {
        numerator = auxSum(this->num, frac.denom, -frac.num, this->denom);
        denominator = auxProd(this->denom, frac.denom);
    }
    catch ( math::RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    return math::Rational(numerator, denominator);
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
math::Rational& math::Rational::operator-=(const math::Rational& frac) throw(math::RationalException)
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
    catch ( math::RationalException ex )
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
 * @param frac - a fraction, this will be multiplied by
 *
 * @return this * frac
 */
math::Rational math::Rational::operator*(const math::Rational& frac) const throw (math::RationalException)
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
    catch ( math::RationalException ex )
    {
        // The only possible exception is int overflow
        throw ex;
    }

    return math::Rational(numerator, denominator);
    // constructor will reduce the result
} // Rational::operator*

/**
 * Multiplication operator (*=) that multiplies frac and this and assigns the resulting value to itself.
 * a *= b  ---> a -> a * b
 *
 * @param frac - a fraction, this will be multiplied by
 *
 * @return reference to itself
 */
math::Rational& math::Rational::operator*=(const math::Rational& frac) throw(math::RationalException)
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
    catch ( math::RationalException ex )
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
math::Rational math::Rational::operator/(const math::Rational& frac) const throw(math::RationalException)
{
    // Unreduced quotient of two fractions:
    //
    //   a     c     a     d      a * d
    //  --- / --- = --- * --- = ---------
    //   b     d     b     c      b * c

    // Check if frac's numerator equals 0
    if ( true == frac.isZero() )
    {
        throw math::RationalException(RationalException::DIVIDE_BY_ZERO);
    }

    int numerator;
    int denominator;

    try
    {
        numerator = auxProd(this->num, frac.denom);
        denominator = auxProd(frac.num, this->denom);
    }
    catch ( math::RationalException ex )
    {
        // int overflow is the only possible exception
        throw ex;
    }

    return math::Rational(numerator, denominator);
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
math::Rational& math::Rational::operator/=(const math::Rational& frac) throw(math::RationalException)
{
    // See definition of fraction division in operator/.
    // Check if frac's numerator equals 0
    if ( true == frac.isZero() )
    {
        throw math::RationalException(math::RationalException::DIVIDE_BY_ZERO);
    }

    // Result will be assigned to itself so use set()
    int numerator;
    int denominator;

    try
    {
        numerator = auxProd(this->num, frac.denom);
        denominator = auxProd(this->denom, frac.num);
    }
    catch ( math::RationalException ex )
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
math::Rational math::Rational::operator-() const
{
    // One (doesn't matter which) member of the fraction must be negated.
    // Let it be the numerator
    return math::Rational(-this->num, this->denom);
}

/**
 * Comparison operator 'equal' (==)
 *
 * @param frac - a fraction to be compared to this
 *
 * @return true if both fractions are equal, false otherwise
 */
bool math::Rational::operator==(const math::Rational& frac) const
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
bool math::Rational::operator!=(const math::Rational& frac) const
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
bool math::Rational::operator>(const math::Rational& frac) const
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
bool math::Rational::operator>=(const math::Rational& frac) const
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
bool math::Rational::operator<(const math::Rational& frac) const
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
bool math::Rational::operator<=(const math::Rational& frac) const
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
std::ostream& math::operator<<(std::ostream& output, const math::Rational& frac)
{
    // output the num and denom into the stream
    output << frac.num << '/' << frac.denom;
    // and return reference of the stream
    return output;
}

/*
 * An auxiliary function that reduces the fraction:
 * divides the numerator and denominator by their greatest common divisor
 */
void math::Rational::reduce()
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
    // if both num (handled a few lines above)and denom (not permitted when setting)
    // are different than 0, the GCD is guaranteed to be a non-zero value
    num /= static_cast<int>(gcd);
    denom /= gcd;
} // Rational::reduce

/*
 * A simple integer implementation of abs
 *
 * @param a
 * @return absolute value of a
 */
unsigned int math::Rational::absolute(int a)
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
unsigned int math::Rational::greatestCommonDivisor(unsigned int first, unsigned int second)
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
unsigned int math::Rational::leastCommonMultiple(unsigned int first, unsigned int second)
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
    
    //TODO check ranges
} // Rational::leastCommonMultiple

/*
 * An auxiliary function that calculates (unreduced) numerator of (this-frac).
 * It is only needed by comparison operators who actually only need to know its sign
 *
 * @param fraction
 *
 * @return -1 if frac is greater, 0 if they are equal, 1 if this is greater
 */
int math::Rational::sign(const math::Rational& frac) const
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

/*
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
int math::Rational::auxSum(int num1, int denom2, int num2, int denom1) throw(math::RationalException)
{
    const long int sum = num1 * denom2 + num2 * denom1;

    if ( sum > INT_MAX || sum < INT_MIN )
    {
        throw math::RationalException(math::RationalException::OVERFLOW);
    }

    return static_cast<int>(sum);
}

/*
 * Auxiliary function to calculate a product of two integer numbers.
 * It checks for integer overflows and throws an exception in this case.
 *
 * @param first
 * @param second
 *
 * @return first * second
 *
 * @throw RationalException in case of an integer overflow
 */
int math::Rational::auxProd(int first, int second) throw(math::RationalException)
{
    const long int prod = first * second;

    if ( prod > INT_MAX || prod < INT_MIN )
    {
        throw math::RationalException(math::RationalException::OVERFLOW);
    }

    return static_cast<int>(prod);
}

/*
 * Positive integer power of 10. Additionally, the result's range is checked
 * to prevent an overflow.
 *
 * @param n - exponent (a positive integer number)
 *
 * @return 10^n
 * 
 * @throw RationalException if the result exceeds unsigned long's range
 */
unsigned long int math::Rational::pow10(unsigned int n) throw (math::RationalException)
{
    // simply multiply 10 by itself n times
    unsigned long long int temp = 1;
    for ( unsigned int i=0; i<n; i++ )
    {
        temp *= 10;
        if ( temp>ULONG_MAX )
        {
            throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
        }
    }

    return static_cast<unsigned long int>(temp);
}

/*
 * Parses a string into a long long integer value and also prevents
 * an integer overflow
 * 
 * @param str - string to be parsed
 * 
 * @return long long value of 'str'
 * 
 * @throw RationalException if absolute value of 'str' exceeds long long's range
 */
long long int math::Rational::str2ll(const std::string& str) throw (math::RationalException)
{
    // sanity check already performed by the caller function
    
    unsigned int lmax = std::numeric_limits<long long int>::digits10;
    if ( '+'==str.at(0) || '-'==str.at(0) )
    {
        // if the string starts with a sign, the allowed number
        // of string's characters may be increased by 1.
        lmax++;
    }

    if ( str.length()>lmax )
    {
        throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
    }

    // std::atoll from <cstdlib> would be a more elegant option, however
    // it may not be supported by older compilers not supporting C++11

    long long int retVal = 0LL;
    std::sscanf(str.c_str(), "%lld", &retVal);
    return retVal;
}

/**
 * Destructor
 */
math::Rational::~Rational()
{
    // nothing to do, i.e. no dynamically allocated memory to free,
    // no resources to release, etc.
}
