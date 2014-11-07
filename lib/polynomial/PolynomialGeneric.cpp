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
 * Implementation of the class PolynomialGeneric.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */


// deliberately there is no #include "PolynomialGeneric.h" !
#include "util/NumericUtil.hpp"
#include "../settings/omp_settings.h"

#include <vector>
#include <cstddef>
#include <new>
#include <stdexcept>
#include <limits>


/**
 * Constructor.
 *
 * @param cvect - vector of coefficients in ascending order ([c0, c1, c2, ... cn])
 *
 * @throw PolynomialException if 'cvect' is an empty vector or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>::PolynomialGeneric(std::vector<T> cvect) throw (math::PolynomialException)
{
    // sanity check
    if ( cvect.size() <= 0 )
    {
        // 'cvect' must contain at least one coefficient
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }
    
    try
    {
        this->copyCoefs(cvect);
        // reduce zero-coefficients from the highest order terms
        this->reduce();
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Constructor that assigns a scalar value to a zero-degree polynomial
 * 
 * @param c0 - scalar value of the polynomial (default: 0)
 * 
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>::PolynomialGeneric(const T& c0) throw (math::PolynomialException)
{
    try
    {
        this->coef.resize(1);
        this->coef.at(0) = c0;
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Copy constructor.
 * Creates an exact copy of the original polynomial 'poly'.
 *
 * @param poly - original polynomial to be copied into this one
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>::PolynomialGeneric(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    try
    {
        // Just copy the poly's coefficients
        this->copyCoefs(poly.coef);
        // 'poly' is supposed to be already reduced, so no need to call reduce()
    }
    catch ( const std::bad_alloc &ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Constructor.
 *
 * @note Usage of this constructor may be dangerous if 'n' is larger then the actual size of the array.
 *       Use it carefully and double check both input arguments!
 *
 * @param carray - array of coefficients in ascending order ([c0, c1, c2, ... cn])
 * @param n - number of elements in the array
 *
 * @throw PolynomialException if input arguments are invalid or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>::PolynomialGeneric(const T* carray, size_t n) throw (math::PolynomialException)
{
    // sanity check
    if ( NULL==carray || n<=0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }
    
    // check that carray's size does not exceed maximum allowed vector's size
    if ( n > this->coef.max_size() )
    {
        throw math::PolynomialException(math::PolynomialException::TOO_LARGE);
    }

    try
    {
        // allocate coef:
        this->coef.clear();
        this->coef.resize(n);

        // And copy all elements from the array.
        #pragma omp parallel for if(n>OMP_CHUNKS_PER_THREAD) default(none) shared(carray, n)
        for ( size_t i=0; i<n; ++i )
        {
            this->coef.at(i) = carray[i];
        }

        reduce();
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
    catch ( const std::out_of_range& oor )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_RANGE);
    }
}


/**
 * Constructor.
 * Allocates memory for 'n' coefficients. All coefficients are assigned a zero value except the
 * highest order coefficients that is assigned '1' to prevent immediate reduction.
 *
 * @param ignored - ignored, its sole purpose is to distinguish the constructor from PolynomialGeneric(const T&)
 * @param n - number of all coefficients (default: 1)
 *
 * @throw PolynomialException if 'n' is invalid or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>::PolynomialGeneric(bool ignored, size_t n) throw (math::PolynomialException)
{
    // sanity check:
    if ( n<=0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }
    
    // check that 'n' does not exceed max. allowed vector's size
    if ( n > this->coef.max_size() )
    {
        throw math::PolynomialException(math::PolynomialException::TOO_LARGE);
    }

    try
    {
        this->coef.clear();
        // Populate the whole vector with "zeros":
        this->coef.resize(n, math::NumericUtil<T>::ZERO);

        // Set the highest order coefficients is set to "1" to prevent reductions:
        if ( n>1 )
        {
            // However this is not necessary for 0 - degree polynomials
            this->coef.at(n-1) = math::NumericUtil<T>::ONE;
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }

    (void) ignored;
}


/*
 * A utility function (only used internally) that copies elements of a vector to coef
 *
 * @param cvect - vector to be copied into coef
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
void math::PolynomialGeneric<T>::copyCoefs(const std::vector<T>& cvect) throw (math::PolynomialException)
{
    try
    {
        this->coef.clear();
        // std::vector's assignment operator (=) actually copies all elements from one vector into the other one
        this->coef = cvect;
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/*
 * Zero-values for highest order coefficients are not allowed. For that reason,
 * this utility function (only used internally) is introduced that "reduces" the
 * polynomial, i.e. reduces all highest order zero coefficients until a non-zero
 * coefficient is found:
 *
 * For instance, 1 + x + 3x^2 + 0x^3 + 0x^4 is reduced into:
 *               1 + x + 3x^2, i.e. coef(4) and coef(3) are removed.
 */
template<class T>
void math::PolynomialGeneric<T>::reduce()
{
    const size_t N = this->coef.size() - 1;
    size_t f = 0;
    
    /*
     * Find the first coefficient (i.e. of the highest order) that is not
     * equal to zero. If necessary, delete all coefficients of higher order,
     * excluding the first non-zero coefficient. The first coefficient (coef[0])
     * must never be deleted even if it also equals zero.
     */
    for ( f=N; f>0 && true==math::NumericUtil<T>::isZero(coef.at(f)); --f );
    
    if ( f<N )
    {
        this->coef.erase(this->coef.begin()+f+1, this->coef.end());
    }
}


/**
 * Assignment operator (=)
 *
 * @param poly - a polynomial to be copied into this one
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator=(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &poly )
    {
        return *this;
    }

    // otherwise copy coefficients from 'poly';
    this->copyCoefs(poly.coef);

    // 'poly' is supposed to be already reduced, so no need to call reduce()

    return *this;
}


/**
 * @return vector of coefficients in ascending order ([c0, c1, c2, ..., cn])
 *
 * @throw PolynomialException if allocation of memory for the output vector fails
 */
template<class T>
std::vector<T> math::PolynomialGeneric<T>::get() const throw (math::PolynomialException)
{
    try
    {
        // This will copy contents of coef into the return value:
        return this->coef;
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * @return vector of coefficients in descending order ([cn, ... , c2, c1, c0])
 *
 * @throw PolynomialException if allocation of memory for the output vector fails
 */
template<class T>
std::vector<T> math::PolynomialGeneric<T>::getDesc() const throw (math::PolynomialException)
{
    try
    {
        const size_t N = this->coef.size();
        // Allocate the return vector:
        std::vector<T> retVal(N);

        // copy elements from coef to retVal in reverse order:
        #pragma omp parallel for if(N>OMP_CHUNKS_PER_THREAD) default(none) shared(retVal)
        for ( size_t i=0; i<N; ++i )
        {
            retVal.at(i) = this->coef.at(N-1-i);
        }

        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * @param pos - degree of the desired term
 *
 * @return coefficient of the 'pos'-degree term or zero if pos  is greater than degree of the polynomial
 */
template<class T>
T math::PolynomialGeneric<T>::get(size_t pos) const
{
    // TODO does it make any sense checking if 'pos' exceeds vector:max_size()???
    
    // if 'pos' exceeds the polynomial's degree, return zero
    if ( pos >= this->coef.size() )
    {
        return math::NumericUtil<T>::ZERO;
    }

    // otherwise it is safe to access the desired coefficient:
    return this->coef.at(pos);
}


/**
 * @return degree of the polynomial, i.e. the non-zero coefficient of the highest degree term
 */
template<class T>
size_t math::PolynomialGeneric<T>::degree() const
{
    /*
        If the polynomial: c0 + c1*x + c2*x^2 + ... + cn*x^n
        is reduced then the size of its 'coef' vector
        is  n+1 elements (0 to n), one more than the actual degree (n).
    */
    return this->coef.size()-1;
}


/**
 * Assigns polynomial's coefficients from 'cvect'
 *
 * @param cvect - vector of coefficients in ascending order ([c0, c1, c2, ..., cn])
 *
 * @return reference to itself
 *
 * @throw PolynomialException if 'cvect' is invalid or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::set(const std::vector<T>& cvect) throw (math::PolynomialException)
{
    // sanity check
    if ( cvect.size() <= 0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }

    // If 'cvect' is OK, just copy its values to coef:
    this->copyCoefs(cvect);

    // 'cvect' is an arbitrary vector that does not
    // necessarily represent a reduced polynomial....
    this->reduce();

    return *this;
}


/**
 * Assigns polynomial's coefficients from 'cvect'
 *
 * @param cvect - vector of coefficients in descending order ([cn, ..., c2, c1, c0])
 *
 * @return reference to itself
 *
 * @throw PolynomialException if 'cvect' is invalid or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::setDesc(const std::vector<T>& cvect) throw (math::PolynomialException)
{
    const size_t N = cvect.size();

    // sanity check
    if ( N<=0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }

    try
    {

        this->coef.clear();
        this->coef.resize(N);

        #pragma omp parallel for if(N>OMP_CHUNKS_PER_THREAD) default(none) shared(cvect)
        for ( size_t i=0; i<N; ++i )
        {
            this->coef.at(i) = cvect.at(N-1-i);
        }

        this->reduce();
        return *this;

    }
    catch ( const std::bad_alloc& ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Modifies the coefficient at the position 'pos'. If 'pos' exceeds the polynomial's degree
 * and 'c' does not equal zero, the appropriate number of coefficients will be appended.
 *
 * @param pos - degree of the desired term
 * @param c - value of the coefficient at 'pos' (default: 0)
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory (possible when 'pos' exceeds the polynomial's degree) fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::set(size_t pos, const T& c) throw (math::PolynomialException)
{
    /*
        If 'pos' exceeds the polynomial's degree, the appropriate number of zero-coefficients
        will be inserted between the highest coef and 'pos'.
        If 'c' equals zero, this does not make any sense as reduce would revert the
        polynomial back into its original state.
    */
    if ( pos >= this->coef.size() )
    {
        if ( false==math::NumericUtil<T>::isZero(c) )
        {
            // check that pos does not exceed max. allowed vector's size
            if ( pos > this->coef.max_size() )
            {
                throw math::PolynomialException(math::PolynomialException::OUT_OF_RANGE);
            }
            
            this->insert(pos, c);
        }
    }
    else
    {
        /*
            If 'pos' is less or equal than the polynomial's degree,
            it is safe to access the desired coefficient directly:
        */
        this->coef.at(pos) = c;

        // it is possible that coef(N) is set to zero....
        this->reduce();
    }

    return *this;
}


/**
 * Inserts a coefficient in front of the pos^th coefficient. Degree of the polynomial is increased.
 * If 'pos' exceeds the current polynomial's degree, the appropriate number of coefficients
 * (their values are assigned to zero) will be appended.
 *
 * @param pos - the new coefficient will be inserted in front of this one
 * @param c - value to be assigned to the inserted coefficient (default: 0)
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::insert(size_t pos, const T& c) throw (math::PolynomialException)
{
    const size_t N = this->coef.size();

    try
    {
        //pos must not exceed vector's max. allowed size
        if ( pos > this->coef.max_size() )
        {
            throw math::PolynomialException(math::PolynomialException::OUT_OF_RANGE);    
        }
        
        // If 'pos' exceeds the polynomial's degree...
        if (pos>=N )
        {
            /*
                ... the appropriate number of coefficients will be inserted.

                However, if 'c' equals zero, reduce would revert the polynomial
                into its original state.
            */
            if ( true==math::NumericUtil<T>::isZero(c) )
            {
                // no need to insert a zero as the top coefficient as reduce() would remove it
                return *this;
            }

            // 'c' is not equal to zero, just append the necessary number of zero-coefficients
            this->coef.insert(this->coef.begin()+N, pos-N+1, math::NumericUtil<T>::ZERO);
            // and assign the highest order coefficient the value of 'c':
            this->coef.at(pos) = c;
        }
        else
        {
            // 'pos' does not exceed the polynomial's degree, just insert the 'c' into coef.
            // Prior to that, check that the expanded coef would not exceed max. allowed size:
            if ( N == this->coef.max_size() )
            {
                throw math::PolynomialException(math::PolynomialException::TOO_LARGE);
            }
            
            this->coef.insert(this->coef.begin()+pos, c);
            // as 'c' is not appended at the end of coef, it cannot be set to "reduced" state
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }

    return *this;
}


/**
 * Removes the desired coefficient from the polynomial.
 * If 'pos' exceeds the polynomial's degree, nothing will be done.
 *
 * @param pos - position of the coefficient to be removed
 *
 * @return reference to itself
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::remove(size_t pos)
{
    // Nothing to do if 'pos' exceeds the polynomial's degree
    if ( pos >= this->coef.size() || this->coef.size() <= 1 )
    {
        return *this;
    }

    // now it is safe to erase the desired coefficient:
    this->coef.erase(this->coef.begin() + pos);

    /*
     It is possible that the highest order coefficient is removed and an
     unreduced polynomial remains:
    */
    this->reduce();

    return *this;
}


/**
 * Evaluates the polynomial.
 *
 * @param x - input argument
 *
 * @return value of the polynomial for 'x': c0 + c1*x + c2*x^2 + ... + cn*x^n
 */
template<class T>
T math::PolynomialGeneric<T>::value(const T& x) const
{
    /*
        Horner's method (explained at: http://en.wikipedia.org/wiki/Horner%27s_method)
        is known to be computationally more efficient than calculating powers of x:

        p(x) = cn*x^n + ... + c2*x^2 + c1*x + c0

        can be rewritten into:

        p(x) = ((((cn*x + c[n-1])*x + ... + c3)*x + c2)*x + c1)*x + c0
    */

    size_t i = this->coef.size() - 1;
    T retVal = this->coef.at(i);

    for ( ; i>0; --i )
    {
        retVal = retVal*x + this->coef.at(i-1);
    }

    return retVal;
}


/**
 * Derivation of the polynomial.
 *
 * If p(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n
 *
 * its derivative is another polynomial, expressed as:
 *
 *   dp(x)/dx = c1 + 2*c2*x + ... + n*cn*x^(n-1)
 *
 * @return derivative of 'this'
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::deriv() const throw (math::PolynomialException)
{
    const size_t N = this->coef.size();

    // if 'this' is already a constant (degree==0), its derivative will also be a constant with value 0:
    const size_t DEG = (1==N ? 1 : N-1);
    math::PolynomialGeneric<T> retVal(true, DEG);

    if ( 1==N )
    {
        /*
            Handling of a zero degree polynomial ( p(x) = c )
            Its derivative is 0.
        */
        retVal.coef.at(0) = NumericUtil<T>::ZERO;
        return retVal;
    }

    // For polynomials of higher degree (>0) apply the formula above:
    #pragma omp parallel for default(none) shared(retVal)
    for ( size_t i=0; i<DEG; ++i )
    {
        retVal.coef.at(i) = static_cast<T>(i+1) * this->coef.at(i+1);
    }

    return retVal;
}


/**
 * Indefinite integral of the polynomial.
 *
 * If p(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n
 *
 * its indefinite integral is another polynomial, expressed as:
 *
 *  /                       c1   2    c2    3          cn    n+1
 *  | p(x) dx = c + c0*x + ----*x  + -----*x  + ... + -----*x
 * /                         2         3               n+1
 *
 * Note that 'c' can be an arbitrary value (this is why, the integral is called 'indefinite')
 *
 * @param c - value to be assigned to "indefinite" coefficient at pos=0 (default: 0)
 *
 * @return indefinite integral of 'this'
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::integ(const T& c) const throw (math::PolynomialException)
{
    const size_t N = this->coef.size();
    
    // Note that integration increments degree od the polynomial...
    if ( N == this->coef.max_size() )
    {
        throw math::PolynomialException(math::PolynomialException::TOO_LARGE);
    }
    
    math::PolynomialGeneric<T> retVal(true, N+1);

    // coef(0) is an arbitrary value, given as c:
    retVal.coef.at(0) = c;

    // for other coefficients, apply the formula above:
    #pragma omp parallel for default(none) shared(retVal)
    for ( size_t i=0; i<N; ++i )
    {
        retVal.coef.at(i+1) = this->coef.at(i)/static_cast<T>(i+1);
    }

    return retVal;
}


/*
 * A static utility function that divides twov polynomials and returns
 * their quotient an remainder.
 *
 * @note If 'q' or 'rem' is a NULL pointer, it will not be filled by
 *       quotient or remainder coeficients, respectively.
 * 
 * @param p1 - dividend polynomial
 * @param p2 - divisor polynomial
 * @param q - pointer to a polynomial that will be assigned "p1 / p2"
 * @param rem - pointer to a polynomial that will be assigned "p1 mod p2"
 *
 * @throw PolynomilExcpetion if attempting to divide by a zero polynomial
 */
template<class T>
void math::PolynomialGeneric<T>::polyDivision(
        const math::PolynomialGeneric<T>& p1,
        const math::PolynomialGeneric<T>& p2,
        math::PolynomialGeneric<T>* q,
        math::PolynomialGeneric<T>* rem )
    throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == p2.isZero() )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    // Nothing to do if no output is specified
    if ( NULL==q && NULL==rem )
    {
        return;
    }

    // Degrees of 'p1' and 'p2':
    const size_t Np1 = p1.coef.size() - 1;
    const size_t Np2 = p2.coef.size() - 1;

    // Use specialized operators of 'p2' is a scalar
    if ( Np2==0 && false==math::NumericUtil<T>::isZero(p2.coef.at(0)) )
    {
        if ( NULL != q )
        {
            *q = p1 / p2.get(0);
        }

        if ( NULL != rem )
        {
            rem->coef.resize(1, math::NumericUtil<T>::ZERO);
        }

        return;
    }

    // Handling of situations when p2's degree is higher than p1's
    if ( Np2 > Np1 )
    {
        if ( NULL != q )
        {
            q->coef.resize(1, math::NumericUtil<T>::ZERO );
        }

        if ( NULL != rem )
        {
            *rem = p1;
        }

        return;
    }

    // Quotient's degree:
    const size_t Nq = Np1 - Np2;

    // As 'p1' must remain constant, this vector will store its coefficients
    // during the division procedure:
    std::vector<T> p(p1.coef);

    // The highest degree coefficient of 'p2':
    const T Cdiv = p2.coef.at( Np2 );

    // Preallocate q's vector of coefficients
    if ( NULL != q )
    {
        q->coef.resize(Nq+1, math::NumericUtil<T>::ZERO);
    }

    // This for loop sequentially updates 'p' so it is
    // not suitable for parallelization
    for ( size_t i=0; i<=Nq; ++i )
    {
        // Divide p's and p2's highest order coefficients...
        const T c = p.at(Np1-i) / Cdiv;
        // ... the quotient is also one of q's coefficients...
        if ( NULL != q )
        {
            q->coef.at(Nq-i) = c;
        }

        // If the for loop below were extended by one iteration,
        // it would also calculate this p's coefficients to 0
        p.at(Np1-i) = math::NumericUtil<T>::ZERO;

        /*
         * The for loop is actually equivalent to multiplication
         * of 'p2' by the q's i.th term, subtracting the product
         * from 'p' and assigning the difference to 'p'.
         *
         * Unlike the outer for loop, the inner loop
         * be parallelized.
         */
        #pragma omp paralel for if(Np2>OMP_CHUNKS_PER_THREAD) default(none) shared(p, p2, i, c)
        for ( size_t j=0; j<Np2; ++j )
        {
            p.at(Nq-i+j) -= c * p2.coef.at(j);
        }
    }

    // Finally assign the remainder of 'p' to 'rem'
    if ( NULL != rem )
    {
        rem->coef = p;
        rem->reduce();
    }
}


/**
 * Addition operator (+=) that adds a polynomial to 'this' and assigns the sum to itself.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param poly - polynomial to be added to this one
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator+=(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // For a definition of polynomial addition, see operator+
    try
    {
        const size_t nthis = this->coef.size();
        const size_t npoly = poly.coef.size();

        // If 'poly' is of higher degree,
        // insert the appropriate number of coefficients and set them to 0:
        if ( nthis<npoly )
        {
            this->coef.insert(this->coef.begin()+nthis, npoly-nthis, math::NumericUtil<T>::ZERO);
        }

        // ... and perform addition of same degree terms' coefficients
        #pragma omp parallel for if(npoly>OMP_CHUNKS_PER_THREAD) default(none) shared(poly)
        for ( size_t i=0; i<npoly; ++i )
        {
            this->coef.at(i) += poly.coef.at(i);
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }

    this->reduce();
    return *this;
}


/**
 * Addition operator (+=) that adds a scalar to 'this' and assigns the sum to itself.
 *
 * @param sc - scalar to be added to "this" polynomial
 *
 * @return reference to itself
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator+=(const T& sc)
{
    this->coef.at(0) += sc;
    this->reduce();
    return *this;
}


/**
 * Subtraction operator (-=) that subtracts a polynomial from 'this' and assigns
 * the difference to itself.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param poly - polynomial to be subtracted from this one
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator-=(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // For a definition of polynomial subtraction, see operator-
    try
    {
        const size_t nthis = this->coef.size();
        const size_t npoly = poly.coef.size();

        // If 'poly' is of higher degree, insert appropriate number of coefficients and set them to 0:
        if ( nthis<npoly )
        {
            this->coef.insert(this->coef.begin()+nthis, npoly-nthis, math::NumericUtil<T>::ZERO);
        }

        // ... and perform addition of same degree terms' coefficients
        #pragma omp parallel for if(npoly>OMP_CHUNKS_PER_THREAD) default(none) shared(poly)
        for ( size_t i=0; i<npoly; ++i )
        {
            this->coef.at(i) -= poly.coef.at(i);
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }

    this->reduce();
    return *this;
}


/**
 * Subtraction operator (-=) that subtracts a scalar from 'this' and assigns the difference to itself.
 *
 * @param sc - scalar to be subtracted from "this" polynomial
 *
 * @return reference to itself
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator-=(const T& sc)
{
    this->coef.at(0) -= sc;
    this->reduce();
    return *this;
}


/**
 * Unary negation operator (-)
 *
 * @return -this
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::operator-() const throw (math::PolynomialException)
{
    /*
        Definition of polynomial negation is similar to negation of vectors/matrices:

                                   Np
                                 -----
                                 \            i
                        -p(x) =  |     -pi * x
                                 /
                                 -----
                                  i=0
    */
    try
    {
        math::PolynomialGeneric<T> retVal(*this);
        const size_t N = this->coef.size();

        // Just negate each coefficient:
        #pragma omp parallel for if(N>OMP_CHUNKS_PER_THREAD) default(none) shared(retVal)
        for ( size_t i=0; i<N; ++i )
        {
            retVal.coef.at(i) = -this->coef.at(i);
        }

        // no need to reduce
        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Multiplication operator (*=) that multiplies two polynomials and assigns the product to itself.
 *
 * @note Polynomials can be of different degrees.
 * @note Multiplication of polynomials is commutative.
 *
 * @param poly - polynomial to be multiplied by this one
 *
 * @return reference to itself
 *
 * @throw PolynomialException if the product would be too large or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator*=(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // for a definition of polynomial multiplication, see operator*

    // perform a regular polynomial multiplication and copy the product's coefficients to itself
    math::PolynomialGeneric<T> temp = *this * poly;
    this->copyCoefs(temp.coef);

    // Note that copyCoefs would throw an exception in case of unsuccessful allocation of memory

    // 'reduce' already performed by operator*
    return *this;
}


/**
 * Multiplication operator (*=) that multiplies a polynomial by a scalar
 * and assigns the product to itself.
 *
 * @param sc - scalar
 *
 * @return reference to itself
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator*=(const T& sc)
{
    // Multiply each coefficient by the scalar
    const size_t N = this->coef.size();

    #pragma omp parallel for default(none) shared(sc)
    for ( size_t i=0; i<N; ++i )
    {
        this->coef.at(i) *= sc;
    }

    // applicable when the scalar is 0...
    this->reduce();
    return *this;
}


/**
 * Division operator (/=) that divides two polynomials and assigns
 * the quotient to itself.
 * 
 * @param poly - polynomial to divide this one
 *
 * @return reference to itself
 * 
 * @throw PolynomialException if attempting to divide by a zero polynomial
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator/=(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == poly.isZero() )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    math::PolynomialGeneric<T> quot;

    math::PolynomialGeneric<T>::polyDivision(*this, poly, &quot, NULL);
    *this = quot;

    return *this;
}


/**
 * Modulation operator (%=) that divides two polynomials and assigns
 * the remainder to itself.
 * 
 * @param poly - polynomial to divide this one
 *
 * @return reference to itself
 * 
 * @throw PolynomialException if attempting to divide by a zero polynomial
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator%=(const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == poly.isZero() )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    math::PolynomialGeneric<T> rem;

    math::PolynomialGeneric<T>::polyDivision(*this, poly, NULL, &rem);
    *this= rem;

    return *this;
}


/**
 * Division operator (/=) that divides a polynomial by a scalar
 * and assigns the quotient to itself:
 *    *this /= sc
 *
 * @param sc - scalar; must not be zero
 *
 * @return reference to itself
 *
 * @throw PolynomialException if attempting to divide by zero
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator/=(const T& sc) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == math::NumericUtil<T>::isZero(sc) )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    /*
     * Just divide each polynoial's coefficient by 'sc'.
     * This is equivalent to multiplication by a reciprocal value
     * of 'sc' which is already implemented by operator*.
     */
    *this *= math::NumericUtil<T>::ONE / sc;

    return *this;
}


/**
 * Modulation operator (%=) that divides a polynomial by a scalar
 * and assigns the remainder to itself:
 *    *this %= sc
 * When a polynomial is divided by a scalar, the remainder is always zero.
 *
 * @param sc - scalar; must not be zero
 *
 * @return reference to itself
 *
 * @throw PolynomialException if attempting to divide by zero
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator%=(const T& sc) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == math::NumericUtil<T>::isZero(sc) )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    /*
     * When a polynomial is divided by a scalar, there is
     * no remainder. Hence *this will be set to
     * a zero polynomial.
     */
    this->coef.clear();
    this->coef.resize(1, math::NumericUtil<T>::ZERO);

    return *this;
}


/**
 * Outputs the polynomial to the stream 'str' in form 'c0 +c1*x +c2*x^2 + ...'
 *
 * @param arg - variable to be displayed (default: 'x')
 * @param str - output stream, where the polynomial will be dislayed to (default: cout)
 */
template<class T>
void math::PolynomialGeneric<T>::display(char arg, std::ostream& str) const
{
    /*
        Primarily the method was introduced for brief unit testing purposes
        and not much effort was invested into a visually "nice" output
    */

    // Display coefficients with powers of the variable in ascending order:
    for ( size_t i=0; i < this->coef.size(); ++i )
    {
        /*
            A space will be displayed between terms to better distinguish them.
            The first (lowest order) term is an exception.
        */
        if ( i>0 )
        {
            str << ' ';
            
            /*
              Display signs of all coefficients except of the first one
              (if it is positive).
            */
            str << std::showpos;
        }

        str << this->coef.at(i);
        str << std::noshowpos;

        // Display '*' between a coefficient and a variable (not necessary for i=0)
        if ( i>0 )
        {
            str << '*' << arg;
        }

        // Display power of the variable (where pow>1)
        if ( i>1 )
        {
            str << '^' << i;
        }
    }

}


/**
 * Destructor
 */
template<class T>
math::PolynomialGeneric<T>::~PolynomialGeneric()
{
    // Vector's destructors would probably clean up this automatically.
    // Anyway, let us clear the vector, just to be aware of allocated resources.
    this->coef.clear();

    // Other dynamically allocated memory (via malloc or new) should be freed here.
    // There are no other resources to release.
}




/**
 * Addition operator (+) of two polynomials.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param p1 - augend
 * @param p2 - addend
 *
 * @return p1 + p2
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator+(const math::PolynomialGeneric<T>& p1, const math::PolynomialGeneric<T>& p2) throw (math::PolynomialException)
{

    try
    {
        const size_t Np1 = p1.coef.size();
        const size_t Np2 = p2.coef.size();

        /*
            Addition of polynomials is similar to addition of vectors/matrices:

                                   N
                                 -----
                                 \                   i
                  p(x) + q(x) =  |      (pi + qi) * x
                                 /
                                 -----
                                  i=0

            where N = max(Np, Nq) and pi=0 if i>Np and qi=0 if i>Nq
        */
        const size_t nmax = ( Np1>=Np2 ? Np1 : Np2 );

        math::PolynomialGeneric<T> retVal(true, nmax);

        // Note that the constructor assigns 1 to coef(nmax) which is not desirable at the moment:
        retVal.coef.at(nmax-1) = NumericUtil<T>::ZERO;

        /*
            Add coefficients of the same degree terms. Where 'i' exceeds size of any polynomial,
            consider its i^th coefficient as 0 (already set above)
        */
        #pragma omp parallel for if(nmax>OMP_CHUNKS_PER_THREAD) default(none) shared(retVal, p1, p2)
        for ( size_t i=0; i<nmax; ++i )
        {
            if ( i<Np1 )
            {
                retVal.coef.at(i) = p1.coef.at(i);
            }

            if ( i<Np2 )
            {
                retVal.coef.at(i) += p2.coef.at(i);
            }
        }

        retVal.reduce();
        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Subtraction operator (-) of two polynomials.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param p1 - minuend
 * @param p2 - subtrahend
 *
 * @return p1 - p2
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator-(const math::PolynomialGeneric<T>& p1, const math::PolynomialGeneric<T>& p2) throw (math::PolynomialException)
{
    try
    {
        const size_t Np1 = p1.coef.size();
        const size_t Np2 = p2.coef.size();

        /*
            Subtraction of polynomials is similar to subtraction of vectors/matrices:

                                   N
                                 -----
                                 \                  i
                  p(x) - q(x) =  |     (pi - qi) * x
                                 /
                                 -----
                                  i=0

            where N = max(Np, Nq) and pi=0 if i>Np and qi=0 if i>Nq
        */

        const size_t nmax = ( Np1>=Np2 ? Np1 : Np2 );

        math::PolynomialGeneric<T> retVal(true, nmax);

        // Note that the constructor assigns 1 to coef(nmax) which is not desirable at the moment:
        retVal.coef.at(nmax-1) = math::NumericUtil<T>::ZERO;

        /*
            Subtract coefficients of the same degree terms. Where 'i' exceeds size of any polynomial,
            consider its ith coefficient as 0 (already set above)
        */
        #pragma omp parallel for if(nmax>OMP_CHUNKS_PER_THREAD) default(none) shared(retVal, p1, p2)
        for ( size_t i=0; i<nmax; ++i )
        {
            if ( i<Np1 )
            {
                retVal.coef.at(i) = p1.coef.at(i);
            }

            if ( i<Np2 )
            {
                retVal.coef.at(i) -= p2.coef.at(i);
            }
        }

        retVal.reduce();
        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Multiplication operator (*) of two polynomials.
 *
 * @note Polynomials can be of different degrees.
 * @note Multiplication of polynomials is commutative.
 *
 * @param p1 - multiplicand
 * @param p2 - multiplier
 *
 * @return p1 * p2
 *
 * @throw PolynomialException if the product would be too large or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator*(const math::PolynomialGeneric<T>& p1, const math::PolynomialGeneric<T>& p2) throw (math::PolynomialException)
{
    /*
        As explained at http://www.mathworks.com/help/matlab/ref/conv.html, multiplication of
        polynomials is equivalent to convolution of vectors.

        If Np is size of p(x) and Nq is size of q(x), the product's size will be:
          N = Np + Nq - 1.

        Coefficients of prod(x) = p(x) * q(x), rewritten for the library's order of coefficients,
        can be calculated as follows:

                         N
                       ------
                       \
            prod(k) =  |      p(i) * q(k-i)          (k = 0 .. N-1)
                       /
                       ------
                        i=0

        where prod(i), p(i) and q(i) denote the i^th coefficient of prod, p or q, respectively.

        Where index is out of any polynomial's coefficients' range, consider its coefficient as zero.

        For a unit test, function conv() can be used in Matlab or GNU Octave.
        Note that both programs take polynomial's coefficients in the opposite order as this library!

    */

    try
    {
        const size_t Np1 = p1.coef.size();
        const size_t Np2 = p2.coef.size();

        /*
            At polynomial multiplication it is possible that product's number of coefficients exceeds
            the maximum allowed vector's size. For that reason, this check is performed.
        */
        if ( Np2>(p1.coef.max_size()-Np1 + 1) )
        {
            throw math::PolynomialException(math::PolynomialException::TOO_LARGE);
        }

        // Size of the product polynomial:
        const size_t N = Np1 + Np2 - 1;

        math::PolynomialGeneric<T> retVal(true, N);

        // Each product's coefficient...
        #pragma omp parallel for default(none) shared(retVal, p1, p2)
        for ( size_t i=0; i<N; ++i )
        {
            // ... is a sum of products as specified above
            T temp = math::NumericUtil<T>::ZERO;
            for ( size_t j=0; j<N; ++j )
            {
                // if any index would point out of respective polynomial's range,
                // treat it as a zero (i.e. skip this iteration)
                if ( j>i || j>=Np1 || (i-j)>=Np2 )
                {
                    continue;   // for j
                }

                temp += p1.coef.at(j) * p2.coef.at(i-j);
            }

            retVal.coef.at(i) = temp;
        }

        retVal.reduce();
        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}


/**
 * Division operator (/) of two polynomials.
 *
 * @param p1 - dividend
 * @param p2 - divisor
 *
 * @return p1 / p2
 *
 * @throw PolynomialException if attempting to divide by a zero polynomial or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator/(const math::PolynomialGeneric<T>& p1, const math::PolynomialGeneric<T>& p2) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == p2.isZero() )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    math::PolynomialGeneric<T> retVal;

    math::PolynomialGeneric<T>::polyDivision(p1, p2, &retVal, NULL);
    return retVal;
}


/**
 * Modulation operator (%) of two polynomials.
 *
 * @param p1 - dividend
 * @param p2 - divisor
 *
 * @return p1 mod p2
 *
 * @throw PolynomialException if attempting to divide by a zero polynomial or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator%(const math::PolynomialGeneric<T>& p1, const math::PolynomialGeneric<T>& p2) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == p2.isZero() )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    math::PolynomialGeneric<T> retVal;

    math::PolynomialGeneric<T>::polyDivision(p1, p2, NULL, &retVal);
    return retVal;
}


/**
 * Addition operator (+) of a polynomial and a scalar.
 *
 * @param poly - augend (a polynomial)
 * @param sc - addend (a scalar)
 *
 * @return poly + sc
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator+(const math::PolynomialGeneric<T>& poly, const T& sc) throw (math::PolynomialException)
{
    math::PolynomialGeneric<T> retVal(poly);

    retVal.coef.at(0) += sc;
    retVal.reduce();
    return retVal;
}


/**
 * Subtraction operator (-) of a polynomial and a scalar.
 *
 * @param poly - minuend (a polynomial)
 * @param sc - subtrahend (a scalar)
 *
 * @return poly - sc
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator-(const math::PolynomialGeneric<T>& poly, const T& sc) throw (math::PolynomialException)
{
    math::PolynomialGeneric<T> retVal(poly);

    retVal.coef.at(0) -= sc;
    retVal.reduce();
    return retVal;
}


/**
 * Multiplication operator (*) for multiplication of a polynomial and a scalar.
 *
 * @param poly - multiplicand (a polynomial)
 * @param sc - multiplier (a scalar)
 *
 * @return poly * sc
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator*(const math::PolynomialGeneric<T>& poly, const T& sc) throw (math::PolynomialException)
{
    /*
        Multiplication of a polynomial by a scalar is trivial:
        Each coefficient is multiplied by the scalar.

                          Np
                         -----
                         \               i
            p(x) * sc =  |    sc * pi * x
                         /
                         -----
                          i=0
    */

    const size_t N = poly.coef.size();
    math::PolynomialGeneric<T> retVal(poly);

    #pragma omp parallel for default(none) shared(retVal, sc)
    for ( size_t i=0; i<N; ++i )
    {
        retVal.coef.at(i) *= sc;
    }

    // applicable when sc==o
    retVal.reduce();
    return retVal;
}


/**
 * Division operator (/) for division of a polynomial by a scalar.
 *
 * @param poly - dividend (a polynomial)
 * @param sc - divisor (a scalar), should not be zero
 *
 * @return poly / sc
 *
 * @throw PolynomialException if attempting to divide by zero or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator/(const math::PolynomialGeneric<T>& poly, const T& sc) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == math::NumericUtil<T>::isZero(sc) )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    /*
     * Just divide each polynoial's coefficient by 'sc'.
     * This is equivalent to multiplication by a reciprocal value
     * of 'sc' which is already implemented by operator*.
     */
    return poly * (math::NumericUtil<T>::ONE / sc);
}


/**
 * Modulation operator (%) for division of a polynomial by a scalar.
 *
 * @param poly - dividend (a polynomial)
 * @param sc - divisor (a scalar), should not be zero
 *
 * @return poly mod sc; always zero
 *
 * @throw PolynomialException if attempting to divide by zero or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator%(const math::PolynomialGeneric<T>& poly, const T& sc) throw (math::PolynomialException)
{
    // Division by zero is not permitted
    if ( true == math::NumericUtil<T>::isZero(sc) )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    // There is no remainder when a polynomial is divided by a scalar,
    // return a zero polynomial.
    return math::PolynomialGeneric<T>(math::NumericUtil<T>::ZERO);
    (void) poly;
}


/**
 * Addition operator (+) of a scalar and a polynomial.
 * This operation is commutative and does the same as operator+(scalar).
 *
 * @param sc - augend (a scalar)
 * @param poly - addend (a polynomial)
 *
 * @return sc + poly
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator+(const T& sc, const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    return (poly + sc);
}


/**
 * Subtraction operator (-) of a scalar and a polynomial.
 *
 * @param sc - minuend (a scalar)
 * @param poly - subtrahend (a polynomial)
 *
 * @return sc - poly
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::operator-(const T& sc, const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    return (-poly + sc);
}


/**
  * Multiplication operator (*) of a scalar and a polynomial.
  * This operation is commutative and does the same as operator*(scalar).
  *
  * @param sc - multiplicand (a scalar)
  * @param poly - multiplier (a polynomial)
  *
  * @return sc * poly
  *
  * @throw PolynomialException if allocation of memory fails
  */
template<class T>
math::PolynomialGeneric<T> math::operator*(const T& sc, const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // 'reduce' performed by operator*(poly, sc)
    return (poly * sc);
}


/**
 * Division operator (/) for division of a scalar by a polynomial.
 *
 * @param sc - dividend (a scalar)
 * @param poly - divisor (a polynomial)
 *
 * @return sc / poly
 *
 * @throw PolynomialEsception if attempting to divide by a zero polynomial
 */
template<class T>
math::PolynomialGeneric<T> math::operator/(const T& sc, const math::PolynomialGeneric<T>& poly) throw (PolynomialException)
{
    math::PolynomialGeneric<T> retVal;

    // Division by a zero polynomial is not permitted
    if ( true == poly.isZero() )
    {
        throw math::PolynomialException(math::PolynomialException::DIVIDE_BY_ZERO);
    }

    if ( poly.coef.size() > 1 )
    {
        /*
         * If poly's degree is 1 or higher,
         * return a zero polynomial.
         */
        retVal = math::PolynomialGeneric<T>( math::NumericUtil<T>::ZERO );
    }
    else
    {
        /*
         * Otherwise this is actually division of two scalars.
         * In this case convert the quotient into a 0-degree polynomial.
         */
        retVal = math::PolynomialGeneric<T>( math::PolynomialGeneric<T>(sc / poly.coef.at(0)) );
    }

    return retVal;
}


/**
 * Modulation operator (%) for division of a scalar by a polynomial.
 *
 * @param sc - dividend (a scalar)
 * @param poly - divisor (a polynomial)
 *
 * @return sc mod poly
 *
 * @throw PolynomialEsception if attempting to divide by a zero polynomial
 */
template<class T>
math::PolynomialGeneric<T> math::operator%(const T& sc, const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    math::PolynomialGeneric<T> retVal;

    // Division by a zero polynomial is not permitted
    if ( true == poly.isZero() )
    {
        throw math::PolynomialException( math::PolynomialException::DIVIDE_BY_ZERO );
    }

    if ( poly.coef.size() > 1 )
    {
        /*
         * If poly's degree is 1 or higher, 'sc' is returned,
         * converted into a 0-degree polynomial.
         */
        retVal = math::PolynomialGeneric<T>(sc);
    }
    else
    {
        /*
         * Otherwise this is a division of two scalars.
         * Since the result's degree must always be lower than
         * the divisor's, a zero polynomial is returned.
         */
        retVal = math::PolynomialGeneric<T>( math::NumericUtil<T>::ZERO );
    }

    return retVal;
}


/**
 * A friend function that outputs the polynomial to output stream
 *
 * @param output - stream to write to
 * @param poly - polynomial to be displayed
 *
 * @return reference of output stream (i.e. 'output')
 */
template<class T>
std::ostream& math::operator<<(std::ostream& output, const math::PolynomialGeneric<T>& poly)
{
    // just pass the polynomial to PolynomialGeneric::display()...
    poly.display('x', output);

    // ... and return a reference of the stream
    return output;
}
