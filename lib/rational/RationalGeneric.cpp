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
 * Implementation of the class RationalGeneric
 */


// no #include "rational/RationalGeneric.hpp" !!
#include "int_util/IntUtilGeneric.hpp"
#include "int_util/IntFactorizationGeneric.hpp"
#include "util/NumericUtil.hpp"
#include "exception/RationalException.hpp"

#include <climits>
#include <limits>
#include <cstddef>
#include <cstdio>
#include <new>
#include <ostream>

//TODO some major improvements are planned!

// A namespace with "private" functions:
namespace math {  namespace RationalNS {  namespace __private
{

/*
 * Parses a string into an integer value and also prevents
 * an integer overflow
 *
 * @param str - string to be parsed
 *
 * @return long long value of 'str'
 *
 * @throw RationalException if absolute value of 'str' exceeds I's range
 */
long long int __str2ll(const std::string& str) throw (math::RationalException)
{
    // sanity check already performed by the caller function

    size_t lmax = std::numeric_limits<long long int>::digits10;
    if ( '+' == str.at(0) || '-' == str.at(0) )
    {
        // if the string starts with a sign, the allowed number
        // of string's characters may be increased by 1.
        ++lmax;
    }

    if ( str.length() > lmax )
    {
        throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
    }

    // std::atoll from <cstdlib> would be a more elegant option, however
    // it may not be supported by older compilers not supporting C++11

    long long int buf = 0LL;
    std::sscanf(str.c_str(), "%lld", &buf);

    return buf;
}


/*
 * Positive integer power of 10. Additionally, the result's range is checked
 * to prevent an overflow.
 *
 * @param n - exponent (a positive integer number)
 *
 * @return 10^n
 *
 * @throw RationalException if the result exceeds unsigned I's range
 */
 long long int __pow10(const size_t n) throw (math::RationalException)
{

#define POW10_BASE            ( 10ULL )
    // simply multiply 10 by itself n times
    unsigned long long int temp = 1ULL;
    const unsigned long long int MAX_FACTOR = ULLONG_MAX / POW10_BASE;

    for ( size_t i=0; i<n; ++i )
    {
        // prevent a possible integer overflow
        if ( temp > MAX_FACTOR )
        {
            throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
        }

        temp *= POW10_BASE;
    }

    if ( temp > static_cast<unsigned long long int>( LLONG_MAX ) )
    {
        throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
    }
 
    return static_cast<long long int>(temp);
#undef POW10_BASE
}


/*
 * An auxiliary function that returns sign of the difference
 * num1/den1 - num2/den2.
 * 
 * The function expects that both denominators are greater than 0 and
 * that fractions are reduced.
 *
 * @param num1 - first fraction's numerator
 * @param den1 - first fraction's denominator
 * @param num2 - second fraction's numerator
 * @param den2 - second fraction's denominator
 *
 * @return -1 if num1/den1 < num2/den2, 0 if fractions are equal, 1 if the first fraction is greater
 */
template <typename I>
short int __sign(const I& num1, const I& den1, const I& num2, const I& den2) 
{

    /*
     * The only case when the reduced fractions can be equal:
     * their numerators and denominators match or both numerators
     * equal 0:
     */
    if ( (num1==num2 && den1==den2) ||
         (static_cast<I>(0)==num1 && static_cast<I>(0)==num2) )
    {
        return 0;
    }

    /*
     * Handling situations when fractions (numerators) have different signs.
     * Note that the situations when both numerators equal 0 has been handled
     * above.
     */
    if ( num1<=static_cast<I>(0) && false==math::IntUtil::isNegative<I>(num2) )
    {
        return -1;
    }
    
    if ( false==math::IntUtil::isNegative<I>(num1) && num2<=static_cast<I>(0) )
    {
        return 1;
    }

    /*
     * Both numerators have the same sign.
     * In this case it is necessary to obtain the sign of the difference
     * 
     *    num1 * den2  -  num2 * den1
     * 
     * The algorithm takes care that the integer overflow does not occur.
     */

    if ( 
        ( ( false==math::IntUtil::isNegative<I>(num1) && 
            static_cast<long long int>(num1) <= LLONG_MAX/static_cast<long long int>(den2) ) ||
          ( true==math::IntUtil::isNegative<I>(num1)  && 
            static_cast<long long int>(num1) >= LLONG_MIN/static_cast<long long int>(den2) ) ) &&
        ( ( false==math::IntUtil::isNegative<I>(num2) && num2 <= LLONG_MAX/den1 ) ||
          ( true==math::IntUtil::isNegative<I>(num2)  && num2 >= LLONG_MIN/den1 ) ) )
    {
        // Both terms of the difference num1*den2-num2*den1 can be evaluated
        // and compared in long long int's range:

        const long long int p1 = static_cast<long long int>(num1) *
                                 static_cast<long long int>(den2);

        const long long int p2 = static_cast<long long int>(num2) *
                                 static_cast<long long int>(den1);

        return ( p1<p2 ? -1 : 1 );
    }

    /*
     * If integer overflow has occurred in the comparison above, 
     * try the last option: convert both fractions to floats and 
     * compare the values:
     */
    const float diff = static_cast<float>(num1) / static_cast<float>(den2) -
                       static_cast<float>(num2) / static_cast<float>(den2);

    return ( diff<static_cast<float>(0) ? -1 : 1 );
}


/*
 * Auxiliary function to calculate a sum or difference of two products.
 * It checks for integer overflows and throws an exception in this case.
 *
 * The result is actually an unreduced numerator of a sum/difference of two 
 * fractions, which is reflected by "weird" parameter names.
 *
 * @param num1 - first fraction's numerator
 * @param denom2 - second fraction's denominator
 * @param num2 - second fraction's numerator
 * @param denom1 - first fraction's denominator
 *
 * @return num1*denom2+num2*denom1 if add==true, num1*denom2-num2*denom1 if add==false
 *
 * @throw RationalException in case of an integer overflow
 */
template <typename I>
long long int __auxSum(
            const I& num1, 
            const I& denom2, 
            const I& num2, 
            const I& denom1,
            const bool add ) 
        throw (math::RationalException)
{
    // All intermediate results are long long int values:

    // Check for possible long long int overflow at both products:
    if ( ( num1 >= 0LL && 
           num1 > LLONG_MAX/static_cast<long long int>(denom2) ) ||
         ( num1 < 0LL  && 
           num1 < LLONG_MIN/static_cast<long long int>(denom2) ) )
    {
        throw math::RationalException(math::RationalException::INT_OVERFLOW);
    }

    const long long int term1 = 
        static_cast<long long int>(num1) * static_cast<long long int>(denom2);

    if ( ( num2 >= 0LL && 
           num2 > LLONG_MAX/static_cast<long long int>(denom1) ) ||
         ( num2 < 0LL  && 
           num2 < LLONG_MIN/static_cast<long long int>(denom1) ) )
    {
        throw math::RationalException(math::RationalException::INT_OVERFLOW);
    }

    const long long int term2 = 
        static_cast<long long int>(num2) * static_cast<long long int>(denom1);

    /*
     * These rather complicated set of conditions will ensure that
     * the result will not fall out of long long int's range:
     * 
     *   LLONG_MIN <= (term1+term2) <= LLONG_MAX   (for add==true)
     *   LLONG_MIN <= (term1-term2) <= LLONG_MAX   (for add==false)
     */
    if ( ( true == add  && 
         ( ( term2<0LL && term1 < (LLONG_MIN-term2) ) ||
           ( term2>0LL && term1 > (LLONG_MAX-term2) ) ) ) ||
         ( false == add &&
         ( ( term2>0LL && term1 < (LLONG_MIN+term2) ) ||
           ( term2<0LL && term1 > (LLONG_MAX+term2) ) ) ) )
    {
        throw math::RationalException(math::RationalException::INT_OVERFLOW);
    }
    
    const long long int res = ( true==add ? term1 + term2 : term1 - term2 );

    return res;
}


/*
 * Auxiliary function to calculate a product of two integer numbers.
 * It checks for integer overflows and throws an exception in this case.
 *
 * @param a - first factor
 * @param b - second factor
 *
 * @return a * b
 *
 * @throw RationalException in case of an integer overflow
 */
template <typename I>
long long int __auxProd(
            const I& a, 
            const I& b) 
        throw (math::RationalException)
{
    // if any factor equals 0, the product will also be 0:
    if ( static_cast<I>(0)==a || static_cast<I>(0)==b)
    {
        return 0LL;
    }

    /*
     * Special handling of b == -1
     * In this case Imin/b would would result in -Imin which
     * is just out of I's range!
     * In this case just return -a, provided that a does
     * not equal Imin as well.
     */
    if ( true == math::IntUtil::isNegative<I>(b) &&
         static_cast<I>(-1) == b )
    {
        if ( std::numeric_limits<long long int>::min() == a )
        {
            throw math::RationalException(math::RationalException::INT_OVERFLOW);
        }

        return -a;
    }

    /*
     * This rather complicated set of checks will properly ensure that
     * 
     *   MIN_I <= (a * b) <= MAX_I
     */
    const long long int Amin = LLONG_MIN / b;
    const long long int Amax = LLONG_MAX / b;

    if ( ( b > static_cast<I>(0) && 
           ( a < Amin || a > Amax ) ) ||
         ( b < static_cast<I>(0) && 
           ( a > Amin || a < Amax ) ) )
    {
        throw math::RationalException(math::RationalException::INT_OVERFLOW);
    }

    return static_cast<long long int>(a) * static_cast<long long int>(b);
}

}}}  // namespace math::RationalNS::__private




/**
 * Constructor.
 *
 * Creates an instance of a fraction and assigns its numerator and denominator
 * as set by input parameters. The fraction will be automatically reduced to
 * relative prime numbers. Zero denominator is not permitted, any attempt to do
 * so will throw an exception. The numerator and denominator are not constant
 * and can be modified later.
 *
 * RationalGeneric<I>(a) will assign a/1
 * RationalGeneric<I>()  will assign 0/1 = 0
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
template <typename I>
math::RationalGeneric<I>::RationalGeneric(
            const I& numerator, 
            const I& denominator) 
        throw(math::RationalException)
{
    /*
     * Just call a function that actually sets both members.
     * It will throw a RationalException if the denominator
     * is attempted to be set to 0
     */
    this->set(numerator, denominator);
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
template <typename I>
math::RationalGeneric<I>::RationalGeneric(
            const std::string& str, 
            const size_t repSeqLen) 
        throw (math::RationalException)
{
    /*
     * Just call a function that actually sets both members.
     * It will throw a RationalException if 'str' is invalid etc.
     */
    this->set(str, repSeqLen);
}


/**
 * Copy constructor.
 * Creates an instance of a fraction and copies numerator
 * and denominator from orig.
 *
 * @param orig - a valid instance of Rational
 */
template <typename I>
math::RationalGeneric<I>::RationalGeneric(const math::RationalGeneric<I>& orig) : 
                m_num(orig.m_num), m_denom(orig.m_denom)
{
    /*
     * 'm_num' and 'm_denom' are already set, nothing more to do here.
     * It is impossible to pass an invalid input parameter (orig.denom == 0)
     * so there is no need for any check. All fractions are automatically reduced
     * so no need to do this explicitly
     */
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
template <typename I>
I math::RationalGeneric<I>::getNumerator() const
{
    return this->m_num;
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
template <typename I>
I math::RationalGeneric<I>::getDenominator() const
{
    return this->m_denom;
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
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::set(
            const I& numerator, 
            const I& denominator) 
        throw(math::RationalException)
{
    /*
     * Some internal methods may call this function, passing
     * input arguments as references to rational's members.
     * To prevent a confusion in such cases, both input arguments
     * are copied to local variables.
     */
    const I local_num = numerator;
    const I local_den = denominator;

    // Check for zero denominator
    if ( static_cast<I>(0) == local_den )
    {
        // this is forbidden so throw an exception
        throw math::RationalException(math::RationalException::ZERO_DENOMINATOR);
    }

    // denominator will always be positive:
    if ( true == math::IntUtil::isNegative<I>(local_den) )
    {
        this->m_denom = -local_den;
        this->m_num = -local_num;
    }
    else
    {
        this->m_denom = local_den;
        this->m_num = local_num;
    }

    // reduce the fraction
    this->__reduce();

    return *this;
}


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
 *       include the decimal point or any digit of the whole-part.
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
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::set(
            const std::string& str, 
            const size_t repSeqLen) 
        throw (math::RationalException)
{
    const size_t LEN = str.length();
    size_t decPoint = LEN; // position of the decimal point. LEN represents no decimal point.
    bool hasDigits = false; // does the string contain at least one decimal digit

    // Check if 'str' is a proper decimal representation
    for ( size_t i=0; i<LEN; ++i )
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
        /*
         * Note: if RationalException is thrown by str2ll or pow10, it will not be
         * caught by this try, instead it will be thrown to the caller.
         */

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
        long long int den1 = 1LL;
        long long int den2 = 0LL;

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

            num1 = math::RationalNS::__private::__str2ll(buf);
            den1 = ( LEN==decPoint ? 1 : math::RationalNS::__private::__pow10(LEN-decPoint-1) );
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
            num1 = math::RationalNS::__private::__str2ll(buf);
            buf.erase(LEN-1-repSeqLen, repSeqLen);
            num2 = math::RationalNS::__private::__str2ll(buf);

            den1 = math::RationalNS::__private::__pow10(LEN-1-decPoint);
            den2 = math::RationalNS::__private::__pow10(LEN-1-decPoint-repSeqLen);
        }

        /*
         * Make sure, differences will not exceed long long's ranges.
         * Note that 'num1' and 'num2' always have the same signs and that
         * 'num1' is always greater than 'num2' by their absolute values.
         * The same also applies for 'den1' and 'den2'
         * From these facts it is not difficult to conclude, that differences cannot
         * exceed ranges of (unsigned) long long.
         */
        const long long int dnum = num1 - num2;
        const long long int dden = den1 - den2;

        if ( dnum < 0LL && false==std::numeric_limits<I>::is_signed ) 
        {
            throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
        }

        this->__setLL(dnum, dden);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::RationalException(math::RationalException::OUT_OF_MEMORY);
    }

    return *this;
}


/**
 * Tries to assign a rational from long long int input arguments.
 * Both arguments are reduced and if they both are within the I's range,
 * the rational is assigned.
 * 
 * @param numerator - desired rational's numerator
 * @param denominator - desired rational's denominator
 * 
 * @throw RationalException if any argument is invalid or out of I's range
 *        after reduction
 */
template <typename I>
void math::RationalGeneric<I>::__setLL(
        const long long int numerator,
        const long long int denominator )
        throw (math::RationalException)
{
    // Zero denominator is not valid:
    if ( 0LL == denominator )
    {
        throw math::RationalException(math::RationalException::ZERO_DENOMINATOR);
    }

    // Special handling of numerator==0
    if ( 0LL == numerator )
    {
        this->m_num = static_cast<I>(0);
        this->m_denom = static_cast<I>(1);
        return;
    }

    
    // Obtain both arguments' greatest common divisor...
    const long long int gcd = math::IntFactorization::greatestCommonDivisor<long long int>(
                math::IntUtil::absolute<long long int>(numerator),
                math::IntUtil::absolute<long long int>(denominator) );

    // ... divide both arguments by it...   
    const long long int numRed = numerator / gcd;
    const long long int denRed = denominator / gcd;

    // ... check if both reduced arguments are within the I's range...
    if ( numRed < std::numeric_limits<I>::min() || 
         numRed > std::numeric_limits<I>::max() ||
         denRed < std::numeric_limits<I>::min() || 
         denRed > std::numeric_limits<I>::max() )
    {
        throw math::RationalException(math::RationalException::INPUT_OUT_OF_RANGE);
    }

    // ... and set the rational to casted values of both reduced arguments.

    /*
     * Note that both terms are already reduced so they are assigned directly
     * to 'this'. Denominator's sign must be positive.
     */
    if ( denRed < 0LL )
    {
        this->m_num = -static_cast<I>(numRed);
        this->m_denom = -static_cast<I>(denRed);
    }
    else
    {
        this->m_num = static_cast<I>(numRed);
        this->m_denom = static_cast<I>(denRed);
    }
}


/**
 * Outputs the fraction to stdout in form 'num/den'
 *
 * @param factor (default: 1): optional, if set, both members will be multiplied by it
 * @param str (default cout): output stream, the fraction will be displayed
 *
 * @note There is no check whether multiplied members exceed integer range.
 */
template <typename I>
void math::RationalGeneric<I>::display(
            const I& factor, 
            std::ostream& str ) const
{
    str << this->m_num * factor << '/' << this->m_denom * factor;
}


/**
 * @return fraction's value converted to the specified numeric type F
 */
template <typename I> template <typename F>
F math::RationalGeneric<I>::toNum() const
{
    return ( static_cast<F>(this->m_num) /
             static_cast<F>(this->m_denom) );
}


/**
 * Checks if the fraction's value equals zero.
 *
 * @return true if equal to zero, false otherwise
 */
template <typename I>
inline bool math::RationalGeneric<I>::isZero() const
{
    // just check if the numerator equals 0
    return ( static_cast<I>(0) == this->m_num );
}


/**
 * Checks if the fraction's value is positive (strictly greater than zero).
 *
 * @return true if positive, false otherwise
 */
template <typename I>
bool math::RationalGeneric<I>::isPositive() const
{
    /*
     * The denominator is always positive (see set), so it is
     * sufficient to check if the numerator is greater than 0
     */
    return ( this->m_num > static_cast<I>(0) );
}


/**
 * Checks if the fraction's value is negative (strictly less than zero)
 *
 * @return true if negative, false otherwise
 */
template <typename I>
bool math::RationalGeneric<I>::isNegative() const
{
    /*
     * the denominator is always positive (see set), so it is
     * sufficient to check if the numerator is less than 0
     */
    return ( math::IntUtil::isNegative<I>(this->m_num) );
}


/**
 * Returns an inverse value (i.e. denom/num) of the fraction.
 * If the fraction's value is zero, an exception will be thrown.
 *
 * @return fraction's inverse fraction
 *
 * @throw RationalException the fraction's value is equal to zero
 */
template <typename I>
math::RationalGeneric<I> math::RationalGeneric<I>::invert() const throw(math::RationalException)
{
    // Check if the numerator is equal to 0
    if ( static_cast<I>(0) == this->m_num )
    {
        // In this case, inversion cannot be done, throw an exception
        throw math::RationalException(math::RationalException::UNINVERTIBLE);
    }

    /*
     * inversion will be possible, return a fraction with with exchanged
     * numerator and denominator
     */
    return math::RationalGeneric<I>(this->m_denom, this->m_num);
}


/**
 * Modifies the fraction with its inverse value (swaps the numerator and denominator).
 * If the numerator equals zero, this is not possible so an exception will be thrown.
 *
 * @return reference to itself
 *
 * @throw RationalException if attempting to invert a fraction whose numerator equals zero
 */
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::inverse() throw(math::RationalException)
{
    // Check if inversion is possible (m_num!=0)
    if ( static_cast<I>(0) == this->m_num )
    {
        // inversion is impossible, throw an exception
        throw math::RationalException(math::RationalException::UNINVERTIBLE);
    }

    // inversion is possible, swap the fractions 'm_num' and 'm_denom'.
    this->set(this->m_denom, this->m_num);

    return *this;
}


/**
 * Assignment operator, copies frac's numerator and denominator.
 *
 * @return reference to itself
 */
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::operator=(
            const math::RationalGeneric<I>& frac )
{
    /*
     * The function's behaviour is similar to the copy constructor.
     * With one exception.
     * Nothing to do when attempting to assign it to itself.
     */
    if ( this != &frac )
    {
        // otherwise copy the m_num and m_denom
        this->m_num = frac.m_num;
        this->m_denom = frac.m_denom;
    }

    return *this;
}


/**
 * Addition operator (+=) that adds frac to itself and assigns the resulting value to itself.
 * a += b  ---> a -> a + b
 *
 * @param frac - a fraction to be added to this
 *
 * @return reference to itself
 */
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::operator+=(
            const math::RationalGeneric<I>& frac) throw(math::RationalException)
{
    // See definition of fraction addition in operator+()
    // Result will be assigned to itself so use set

    const long long int numerator = 
            math::RationalNS::__private::__auxSum<I>(
                this->m_num, frac.m_denom, frac.m_num, this->m_denom, true);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(
                this->m_denom, frac.m_denom);

    this->__setLL(numerator, denominator);
    // __setLL will reduce the result

    return *this;
}


/**
 * Subtraction operator (-=) that subtracts frac from this and assigns the resulting value to itself.
 * a -= b  ---> a -> a - b
 *
 * @param frac - a fraction to be subtracted from this
 *
 * @return reference to itself
 */
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::operator-=(
            const math::RationalGeneric<I>& frac) throw(math::RationalException)
{
    // See definition of fraction subtraction in operator-
    // Result will be assigned to itself so use set().

    const long long int numerator = 
            math::RationalNS::__private::__auxSum<I>(
                this->m_num, frac.m_denom, frac.m_num, this->m_denom, false);
 
    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(this->m_denom, frac.m_denom);

    this->__setLL(numerator, denominator);
    // __setLL will reduce the result

    return *this;
}


/**
 * Multiplication operator (*=) that multiplies frac and this and assigns the resulting value to itself.
 * a *= b  ---> a -> a * b
 *
 * @param frac - a fraction, this will be multiplied by
 *
 * @return reference to itself
 */
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::operator*=(
            const math::RationalGeneric<I>& frac ) 
        throw(math::RationalException)
{
    // See definition of fraction multiplication in operator*.
    // Result will be assigned to itself so use set().

    const long long numerator = 
            math::RationalNS::__private::__auxProd<I>(this->m_num, frac.m_num);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(this->m_denom, frac.m_denom);

    this->__setLL(numerator, denominator);

    return *this;
}


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
template <typename I>
math::RationalGeneric<I>& math::RationalGeneric<I>::operator/=(
            const math::RationalGeneric<I>& frac ) 
        throw (math::RationalException)
{
    // See definition of fraction division in operator/.
    // Check if frac's numerator equals 0
    if ( true == frac.isZero() )
    {
        throw math::RationalException(math::RationalException::DIVIDE_BY_ZERO);
    }

    // Result will be assigned to itself so use set()

    const long long int numerator = 
            math::RationalNS::__private::__auxProd<I>(this->m_num, frac.m_denom);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(this->m_denom, frac.m_num);

    this->__setLL(numerator, denominator);

    return *this;
}


/*
 * An auxiliary function that reduces the fraction:
 * divides the numerator and denominator by their greatest common divisor
 */
template <typename I>
void math::RationalGeneric<I>::__reduce()
{
    /*
     * the GCD is not defined when numerator is equal to 0.
     * In that case just assign denominator to 1
     * (it can be assigned to any integer value except 0)
     */
    if ( static_cast<I>(0) == this->m_num )
    {
        this->m_denom = static_cast<I>(1);
        return;
    }

    /*
     * The GCD algorithm requires both integers to be positive.
     * The numerator is allowed to be negative, so get its absolute value.
     */
    const I absNum = math::IntUtil::absolute<I>(this->m_num);

    /*
     * Finally obtain the greatest common divisor...
     * Note: since neither 'absNum' nor 'denom' cannot be zero (handled before),
     * IntFactorization::gcd will never throw an exception
     */
    const I gcd = 
        math::IntFactorization::greatestCommonDivisor<I>(absNum, this->m_denom);

    /*
     * ... and divide both members by it.
     * if both m_num (handled a few lines above)and m_denom (not permitted when setting)
     * are different than 0, the GCD is guaranteed to be a non-zero value
     */
    this->m_num /= gcd;
    this->m_denom /= gcd;
}


/**
 * Destructor
 */
template <typename I>
math::RationalGeneric<I>::~RationalGeneric()
{
    // nothing to do, i.e. no dynamically allocated memory to free,
    // no resources to release, etc.
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
template <typename I>
std::ostream& math::operator<<(
            std::ostream& output, 
            const math::RationalGeneric<I>& frac)
{
    // output the num and denom into the stream
    output << frac.m_num << '/' << frac.m_denom;
    // and return reference of the stream
    return output;
}


/**
 * Unary operator '+', returns a copy of the input argument 'f'.
 * 
 * @note Usage of this operator should be avoided
 * 
 * @param f - rational number to be copied
 * 
 * @return copy of 'f'
 */
template <typename I>
math::RationalGeneric<I> math::operator+(const math::RationalGeneric<I>& f)
{
   return f; 
}


/**
 * Unary negation operator (-)
 *
 * @param f - rational number to be negated
 * 
 * @return -f
 * 
 * @throw RationalExcpetion if I is an unsigned type
 */
template <typename I>
math::RationalGeneric<I> math::operator-(const math::RationalGeneric<I>& f)
                         throw (math::RationalException)
{
    // check if I is an unsigned type
    if ( f.m_num != static_cast<I>(0) &&
         false == math::IntUtil::isNegative<I>(-1) )
    {
        throw math::RationalException(math::RationalException::UNSIGNED);
    }

    /*
     * One (doesn't matter which) member of the fraction must be negated.
     * Let it be the numerator
     */
    return math::RationalGeneric<I>(-f.m_num, f.m_denom);
}


/**
 * Binary addition operator (+) for addition of two rational numbers
 *
 * @param f1 - augend
 * @param f2 - addend
 *
 * @return f1 + f2
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator+(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2 ) 
        throw (math::RationalException)
{
    /*
     * Unreduced sum of two fractions:
     *
     *   a     c     a*d + c*b
     *  --- + --- = -----------
     *   b     d       b * d
     */

    // Just construct an unreduced fraction as shown above:

    const long long int numerator = 
            math::RationalNS::__private::__auxSum<I>(
                f1.m_num, f2.m_denom, f2.m_num, f1.m_denom, true);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(f1.m_denom, f2.m_denom);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, denominator);

    return retVal;
}


/**
 * Binary addition operator (+) for addition of a rational
 * and an integer
 *
 * @param f - augend (a rational)
 * @param i - addend (an integer)
 *
 * @return f + i
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator+(
            const math::RationalGeneric<I>& f,
            const I& i )
        throw (RationalException)
{
    const long long int numerator =
            math::RationalNS::__private::__auxSum<I>(
                f.m_num, static_cast<I>(1), i, f.m_denom, true );

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_denom));

    return retVal;
}


/**
 * Binary addition operator (+) for addition of an integer
 * and a rational
 *
 * @param i - augend (an integer)
 * @param f - addend (a rational)
 *
 * @return i + f
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator+(
            const I& i,
            const math::RationalGeneric<I>& f )
        throw (RationalException)
{
    const long long int numerator =
            math::RationalNS::__private::__auxSum<I>(
                i, f.m_denom, f.m_num, static_cast<I>(1), true );

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_denom));

    return retVal;
}


/**
 * Binary subtraction operator (-) for subtraction of two rational numbers
 *
 * @param f1 - minuend
 * @param f2 - subtrahend
 *
 * @return f1 - f2
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator-(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2 ) 
        throw (math::RationalException)
{
    /*
     * Unreduced difference of two fractions:
     *
     *   a     c     a*d - c*b
     *  --- - --- = -----------
     *   b     d       b * d
     */

    const long long int numerator = 
            math::RationalNS::__private::__auxSum<I>(
                f1.m_num, f2.m_denom, f2.m_num, f1.m_denom, false);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(f1.m_denom, f2.m_denom);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, denominator);

    return retVal;

    /*
     * NOTE: it would also be possible to multiply frac by (-1) and add it to this
     * (using operator+) but it would require a bit more operations.
     */
}


/**
 * Binary subtraction operator (-) for subtraction of an integer
 *  from a two rational
 *
 * @param f - minuend (a rational)
 * @param i - subtrahend (an integer)
 *
 * @return f - i
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator-(
            const math::RationalGeneric<I>& f, 
            const I& i ) 
        throw (math::RationalException)
{
    const long long int numerator = 
            math::RationalNS::__private::__auxSum<I>(
                f.m_num, static_cast<I>(1), i, f.m_denom, false);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_denom));

    return retVal;
}


/**
 * Binary subtraction operator (-) for subtraction of a rational
 *  from an integer
 *
 * @param i - minuend (an integer)
 * @param f - subtrahend (a rational)
 *
 * @return i - f
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator-(
            const I& i, 
            const math::RationalGeneric<I>& f ) 
        throw (math::RationalException)
{
    const long long int numerator = 
            math::RationalNS::__private::__auxSum<I>(
                i, f.m_denom, f.m_num, static_cast<I>(1), false);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_denom));

    return retVal;
}


/**
 * Binary multiplication operator (*) for multiplication of two
 * rational numbers.
 *
 * @param f1 - multiplicand
 * @param f2 - multiplier
 *
 * @return f1 * f2
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator*(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2) 
        throw (math::RationalException)
{
    /*
     * Unreduced product of two fractions:
     *
     *   a     c      a * c
     *  --- * --- = ---------
     *   b     d      b * d
     */

    const long long int numerator = 
            math::RationalNS::__private::__auxProd<I>(f1.m_num, f2.m_num);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(f1.m_denom, f2.m_denom);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, denominator);
    
    return retVal;
}


/**
 * Binary multiplication operator (*) for multiplication of a
 * rational by an integer
 *
 * @param f - multiplicand (a rational)
 * @param i - multiplier (an integer)
 *
 * @return f * i
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator*(
            const math::RationalGeneric<I>& f, 
            const I& i) 
        throw (math::RationalException)
{
    const long long int numerator = 
            math::RationalNS::__private::__auxProd<I>(f.m_num, i);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_denom));

    return retVal;
}


/**
 * Binary multiplication operator (*) for multiplication of an
 * integer by a rational
 *
 * @param i - multiplicand (an integer)
 * @param f - multiplier (a rational)
 *
 * @return i * f
 * 
 * @throw RationalException if an integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator*(
            const I& i, 
            const math::RationalGeneric<I>& f ) 
        throw (math::RationalException)
{
    const long long int numerator = 
            math::RationalNS::__private::__auxProd<I>(f.m_num, i);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_denom));

    return retVal;
}


/**
 * Binary division operator (/) for division of two rational numbers.
 * Division by zero is not permitted and will throw an exception.
 *
 * @param f1 - dividend
 * @param f2 - divisor
 *
 * @return f1 / f2
 *
 * @throw RationalException when attempting to divide by zero or an
 *        integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator/(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2 ) 
        throw (math::RationalException)
{
    /*
     * Unreduced quotient of two fractions:
     *
     *   a     c     a     d      a * d
     *  --- / --- = --- * --- = ---------
     *   b     d     b     c      b * c
     */

    // Check if frac's numerator equals 0
    if ( true == f2.isZero() )
    {
        throw math::RationalException(RationalException::DIVIDE_BY_ZERO);
    }

    const long long int numerator = 
            math::RationalNS::__private::__auxProd<I>(f1.m_num, f2.m_denom);

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(f2.m_num, f1.m_denom);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, denominator);

    return retVal;
}


/**
 * Binary division operator (/) for division of a rational
 * by an integer.
 * Division by zero is not permitted and will throw an exception.
 *
 * @param f - dividend (a rational)
 * @param i - divisor (an integer)
 *
 * @return f / i
 *
 * @throw RationalException when attempting to divide by zero or an
 *        integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator/(
            const math::RationalGeneric<I>& f, 
            const I& i ) 
        throw (math::RationalException)
{
    // Check if 'i' equals 0
    if ( static_cast<I>(0) == i )
    {
        throw math::RationalException(RationalException::DIVIDE_BY_ZERO);
    }

    const long long int denominator = 
            math::RationalNS::__private::__auxProd<I>(i, f.m_denom);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(static_cast<long long int>(f.m_num), denominator);

    return retVal;
}


/**
 * Binary division operator (/) for division of an integer
 * by a rational.
 * Division by zero is not permitted and will throw an exception.
 *
 * @param i - dividend (an integer)
 * @param f - divisor (a rational)
 *
 * @return i / f
 *
 * @throw RationalException when attempting to divide by zero or an
 *        integer overflow occurs
 */
template <typename I>
math::RationalGeneric<I> math::operator/(
            const I& i, 
            const math::RationalGeneric<I>& f ) 
        throw (math::RationalException)
{
    // Check if f's numerator equals 0
    if ( true == f.isZero() )
    {
        throw math::RationalException(RationalException::DIVIDE_BY_ZERO);
    }

    const long long int numerator = 
            math::RationalNS::__private::__auxProd<I>(i, f.m_denom);

    math::RationalGeneric<I> retVal;
    retVal.__setLL(numerator, static_cast<long long int>(f.m_num));

    return retVal;
}


/**
 * Comparison operator 'equal' (==)
 *
 * @param f1 - the first rational number to be compared
 * @param f2 - the second rational number to be compared
 *
 * @return true if both rationals are equal, false otherwise
 */
template <typename I>
bool math::operator==(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2)
{
    return ( 0 == math::RationalNS::__private::__sign<I>(
                  f1.m_num, f1.m_denom, f2.m_num, f2.m_denom) );
}


/**
 * Comparison operator 'equal' (==) for comparison of a rational and
 * an integer
 *
 * @param f - rational number to be compared
 * @param i - integer number to be compared to 'f'
 *
 * @return true if both numbers are equal, false otherwise
 */
template <typename I>
bool math::operator==(
            const math::RationalGeneric<I>& f,
            const I& i )
{
    return ( 0 == math::RationalNS::__private::__sign<I>(
                  f.m_num, f.m_denom, i, static_cast<I>(1)) );
}


/**
 * Comparison operator 'equal' (==) for comparison of and integer
 * and a rational
 *
 * @param i - integer number to be compared to 'f'
 * @param f - rational number to be compared
 *
 * @return true if both numbers are equal, false otherwise
 */
template <typename I>
bool math::operator==(
            const I& i,
            const math::RationalGeneric<I>& f )
{
    return ( 0 == math::RationalNS::__private::__sign<I>(
                  i, static_cast<I>(1), f.m_num, f.m_denom) );
}


/**
 * Comparison operator 'not equal' (!=)
 *
 * @param f1 - the first rational number to be compared
 * @param f2 - the second rational number to be compared
 *
 * @return true, if rationals are not equal, false if they are equal
 */
template <typename I>
bool math::operator!=(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2)
{
    return ( 0 != math::RationalNS::__private::__sign<I>(
                  f1.m_num, f1.m_denom, f2.m_num, f2.m_denom) );
}


/**
 * Comparison operator 'not equal' (!=) for comparison of
 * a rational and an integer
 *
 * @param f - rational number to be compared
 * @param i - integer number to be compared to 'f'
 *
 * @return true, if numbers are not equal, false if they are equal
 */
template <typename I>
bool math::operator!=(
            const math::RationalGeneric<I>& f,
            const I& i )
{
    return ( 0 != math::RationalNS::__private::__sign<I>(
                  f.m_num, f.m_denom, i, static_cast<I>(1)) );
}


/**
 * Comparison operator 'not equal' (!=) for comparison of
 * an integer and a rational
 *
 * @param i - integer number to be compared to 'f'
 * @param f - rational number to be compared
 *
 * @return true, if numbers are not equal, false if they are equal
 */
template <typename I>
bool math::operator!=(
            const I& i,
            const math::RationalGeneric<I>& f )
{
    return ( 0 != math::RationalNS::__private::__sign<I>(
                  i, static_cast<I>(1), f.m_num, f.m_denom) );
}


/**
 * Comparison operator 'greater than' (>)
 *
 * @param f1 - the first rational number to be compared
 * @param f2 - the second rational number to be compared
 *
 * @return 'true' if f1>f2, 'false' otherwise
 */
template <typename I>
bool math::operator>(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2)
{
    return ( 0 < math::RationalNS::__private::__sign<I>(
                 f1.m_num, f1.m_denom, f2.m_num, f2.m_denom) );
}


/**
 * Comparison operator 'greater than' (>) for comparison of
 * a rational and an integer
 *
 * @param f - rational number to be compared
 * @param i - integer number to be compared to 'f'
 *
 * @return 'true' if f>i, 'false' otherwise
 */
template <typename I>
bool math::operator>(
            const math::RationalGeneric<I>& f,
            const I& i )
{
    return ( 0 < math::RationalNS::__private::__sign<I>(
                 f.m_num, f.m_denom, i, static_cast<I>(1)) );
}


/**
 * Comparison operator 'greater than' (>) for comparison of
 * an integer and a rational
 *
 * @param i - integer number to be compared to 'f'
 * @param f - rational number to be compared
 *
 * @return 'true' if i>f, 'false' otherwise
 */
template <typename I>
bool math::operator>(
            const I& i,
            const math::RationalGeneric<I>& f )
{
    return ( 0 < math::RationalNS::__private::__sign<I>(
                 i, static_cast<I>(1), f.m_num, f.m_denom) );
}


/**
 * Comparison operator 'greater or equal' (>=)
 *
 * @param f1 - the first rational number to be compared
 * @param f2 - the second rational number to be compared
 *
 * @return 'true' if f1>=f2, 'false' otherwise
 */
template <typename I>
bool math::operator>=(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2)
{
    return ( 0 <= math::RationalNS::__private::__sign<I>(
                  f1.m_num, f1.m_denom, f2.m_num, f2.m_denom) );
}


/**
 * Comparison operator 'greater or equal' (>=) for comparison of
 * a rational and an integer
 *
 * @param f - rational number to be compared
 * @param i - integer number to be compared to 'f'
 *
 * @return 'true' if f>=i, 'false' otherwise
 */
template <typename I>
bool math::operator>=(
            const math::RationalGeneric<I>& f,
            const I& i )
{
    return ( 0 <= math::RationalNS::__private::__sign<I>(
                  f.m_num, f.m_denom, i, static_cast<I>(1)) );
}


/**
 * Comparison operator 'greater or equal' (>=) for comparison of
 * an integer and a rational
 *
 * @param i - integer number to be compared to 'f'
 * @param f - rational number to be compared
 *
 * @return 'true' if i>=f, 'false' otherwise
 */
template <typename I>
bool math::operator>=(
            const I& i,
            const math::RationalGeneric<I>& f)
{
    return ( 0 <= math::RationalNS::__private::__sign<I>(
                  i, static_cast<I>(1), f.m_num, f.m_denom) );
}


/**
 * Comparison operator 'less than' (<)
 *
 * @param f1 - the first rational number to be compared
 * @param f2 - the second rational number to be compared
 *
 * @return 'true' if f1<f2, 'false' otherwise
 */
template <typename I>
bool math::operator<(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2)
{
    return ( 0 > math::RationalNS::__private::__sign<I>(
                 f1.m_num, f1.m_denom, f2.m_num, f2.m_denom) );
}


/**
 * Comparison operator 'less than' (<) for comparison of
 * a rational and an integer
 *
 * @param f - rational number to be compared
 * @param i - integer number to be compared to 'f'
 *
 * @return 'true' if f<i, 'false' otherwise
 */
template <typename I>
bool math::operator<(
            const math::RationalGeneric<I>& f,
            const I& i )
{
    return ( 0 > math::RationalNS::__private::__sign<I>(
                 f.m_num, f.m_denom, i, static_cast<I>(1)) );
}


/**
 * Comparison operator 'less than' (<) for comparison of
 * an integer and a rational
 *
 * @param i - integer number to be compared to 'f'
 * @param f - rational number to be compared
 *
 * @return 'true' if i<f, 'false' otherwise
 */
template <typename I>
bool math::operator<(
            const I& i,
            const math::RationalGeneric<I>& f )
{
    return ( 0 > math::RationalNS::__private::__sign<I>(
                 i, static_cast<I>(1), f.m_num, f.m_denom ) );
}

/**
 * Comparison operator 'less or equal' (<=)
 *
 * @param f1 - the first rational number to be compared
 * @param f2 - the second rational number to be compared
 *
 * @return - 'true' if f1<=f2, 'false' otherwise
 */
template <typename I>
bool math::operator<=(
            const math::RationalGeneric<I>& f1, 
            const math::RationalGeneric<I>& f2)
{
    return ( 0 >= math::RationalNS::__private::__sign<I>(
                  f1.m_num, f1.m_denom, f2.m_num, f2.m_denom) );
}


/**
 * Comparison operator 'less or equal' (<=) for comparison of
 * a rational and an integer
 *
 * @param f - rational number to be compared
 * @param i - integer number to be compared to 'f'
 *
 * @return - 'true' if f<=i, 'false' otherwise
 */
template <typename I>
bool math::operator<=(
            const math::RationalGeneric<I>& f,
            const I& i )
{
    return ( 0 >= math::RationalNS::__private::__sign<I>(
                  f.m_num, f.m_denom, i, static_cast<I>(1)) );
}


/**
 * Comparison operator 'less or equal' (<=) for comparison of
 * an integer and a rational
 *
 * @param i - integer number to be compared to 'f'
 * @param f - rational number to be compared
 *
 * @return - 'true' if i<=f, 'false' otherwise
 */
template <typename I>
bool math::operator<=(
            const I& i,
            const math::RationalGeneric<I>& f )
{
    return ( 0 >= math::RationalNS::__private::__sign<I>(
                  i, static_cast<I>(1), f.m_num, f.m_denom) );
}



/*
 * Specialization of other classes' templated functions for
 * the class RationalGeneric.
 *
 * Note: the specialized functions must be implemented within
 *       classes' corresponding namespaces.
 */

namespace math
{

namespace NumericUtil
{

/*
 * "Specialization" of NumericUtil::isZero()
 */
template <typename I>
bool isZero(const math::RationalGeneric<I>& value, const math::RationalGeneric<I>& eps)
{
    // RationalGeneric already contains its own isZero()...
    return value.isZero();

    (void) eps;
}

}  // namespace NumericUtil

}  // namespace math


namespace std
{

/*
 * "Specialization" of std::abs()
 */
template <typename I>
math::RationalGeneric<I> abs(const math::RationalGeneric<I>& f)
{
    return ( true==f.isNegative() ? -f : f );
}

}
