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
 * @file PolynomialGeneric.cpp
 *
 * Implementation of the class PolynomialGeneric.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 *
 * @author Jernej Kovacic
 */


// deliberately there is no #include "PolynomialGeneric.h" !
#include "NumericUtil.h"

#include <vector>
#include <cstdlib>
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
        // cvect must contain at least one coefficent
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }

    try
    {
        copyCoefs(cvect);
        // reduce zero-coefficients from the highest order terms
        reduce();
    }
    catch ( std::bad_alloc &ba )
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
math::PolynomialGeneric<T>::PolynomialGeneric(const PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    try
    {
        // Just copy the poly's coefficients
        copyCoefs(poly.coef);
        // 'poly' is supposed to be already reduced, so no need to call reduce()
    }
    catch ( std::bad_alloc &ba )
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
math::PolynomialGeneric<T>::PolynomialGeneric(T* carray, unsigned int n) throw (math::PolynomialException)
{
    // sanity check
    if ( NULL==carray || n<=0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }

    try
    {
        // allocate coef:
        coef.clear();
        coef.resize(n);

        // And copy all elements from the array.
        for ( unsigned int i=0; i<n; i++ )
        {
            coef.at(i) = carray[i];
        }

        reduce();
    }
    catch ( std::bad_alloc &ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
    catch ( std::out_of_range& oor )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_RANGE);
    }
}

/**
 * Constructor.
 * Allocates memory for 'n' coefficients. All coefficients are assigned a zero value except the
 * highest order coefficients that is assigned '1' to prevent immediate reduction.
 *
 * @param n - number of all coefficients (default: 1)
 *
 * @throw PolynomialException if 'n' is invalid or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>::PolynomialGeneric(unsigned int n) throw (math::PolynomialException)
{
    // sanity check:
    if ( n<=0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }

    try
    {
        coef.clear();
        // Populate the whole vector with "zeros":
        coef.resize(n, math::NumericUtil<T>::ZERO);

        // Set the highest order coefficients is set to "1" to prevent reductions:
        if ( n>1 )
        {
            // However this is not necessary for 0 - degree polynomials
            coef.at(n-1) = math::NumericUtil<T>::ONE;
        }
    }
    catch ( std::bad_alloc &ba )
    {
        // Memory allocation failed
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
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
        coef.clear();
        // std::vector's assignment operator (=) actually copies all elements from one vector into the other one
        coef = cvect;
    }
    catch ( std::bad_alloc &ba )
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
    const unsigned int N = coef.size();

    /*
        If all coefficients equal zero, coef(0) will be preserved.
        This is allowed.
    */
    for ( unsigned int i=N-1; i>0 && true==NumericUtil<T>::isZero(coef.at(i)); i-- )
    {
        coef.erase(coef.begin()+i);
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
    copyCoefs(poly.coef);

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
    catch ( std::bad_alloc &ba )
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
        const unsigned int N = this->coef.size();
        // Allocate the return vector:
        std::vector<T> retVal(N);

        // copy elements from coef to retVal in reverse order:
        for ( unsigned int i=0; i<N; i++ )
        {
            retVal.at(i) = this->coef.at(N-1-i);
        }

        return retVal;
    }
    catch ( std::bad_alloc& ba )
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
T math::PolynomialGeneric<T>::get(unsigned int pos) const
{
    // if 'pos' exceeds the polynomial's degree, return zero
    if ( pos>=coef.size() )
    {
        return math::NumericUtil<T>::ZERO;
    }

    // otherwise it is safe to access the desired coefficient:
    return coef.at(pos);
}

/**
 * @return degree of the polynomial, i.e. the non-zero coefficient of the highest degree term
 */
template<class T>
unsigned int math::PolynomialGeneric<T>::degree() const
{
    /*
        If the polynomial: c0 + c1*x + c2*x^2 + ... + cn*x^n
        is reduced then the size of its 'coef' vector
        is  n+1 elements (0 to n), one more than the actual degree (n).
    */
    return coef.size()-1;
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
    copyCoefs(cvect);

    // 'cvect' is an arbitrary vector that does not
    // necessarily represent a reduced polynomial....
    reduce();

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
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::setDesc(const std::vector<T> cvect) throw (math::PolynomialException)
{
    const unsigned int N = cvect.size();

    // sanity check
    if ( N<=0 )
    {
        throw math::PolynomialException(math::PolynomialException::INVALID_ARGUMENT);
    }

    try
    {

        coef.clear();
        coef.resize(N);

        for ( unsigned int i=0; i<N; i++ )
        {
            coef.at(i) = cvect.at(N-1-i);
        }

        reduce();
        return *this;

    }
    catch ( std::bad_alloc& ba )
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
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::set(unsigned int pos, const T& c) throw (math::PolynomialException)
{
    /*
        If 'pos' exceeds the polynomial's degree, the appropriate number of zero-coefficients
        will be inserted between the highest coef and 'pos'.
        If 'c' equals zero, this does not make any sense as reduce would revert the
        polynomial back into its original state.
    */
    if ( pos>=coef.size() )
    {
        if ( false==math::NumericUtil<T>::isZero(c) )
        {
            insert(pos, c);
        }
    }
    else
    {
        /*
            If 'pos' is less or equal than the polynomial's degree,
            it is safe to access the desired coefficient directly:
        */
        coef.at(pos) = c;

        // it is possible that coef(N) is set to zero....
        reduce();
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
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::insert(unsigned int pos, const T& c) throw (math::PolynomialException)
{
    const unsigned int N = this->coef.size();

    try
    {
        // If 'pos' excees the polynomial's degre...
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
            coef.insert(coef.begin()+N, pos-N+1, math::NumericUtil<T>::ZERO);
            // and assign the highest order coefficient the value of 'c':
            coef.at(pos) = c;
        }
        else
        {
            //    'pos' does not exceed the polynomial's degree, just insert the 'c' into coef:
            coef.insert(coef.begin()+pos, c);
            // as 'c' is not appended at the end of coef, it cannot be set to "reduced" state
        }
    }
    catch ( std::bad_alloc& ba )
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
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::remove(unsigned int pos)
{
    // Nothing to do if 'pos' exceeds the polynomial's degree
    if ( pos<0 || pos>=coef.size() || coef.size()<=1 )
    {
        return *this;
    }

    // now it is safe to erase the desired coefficient:
    coef.erase(coef.begin() + pos);

    /*
     It is possible that the highest order coefficient is removed and an
     unreduced polynomial remains:
    */
    reduce();

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

    unsigned int i = coef.size() - 1;
    T retVal = coef.at(i);

    for ( ; i>0; i-- )
    {
        retVal = retVal*x + coef.at(i-1);
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
    const unsigned int N = coef.size();

    // if 'this' is already a constant (degree==0), its derivative will also be a constant with value 0:
    const unsigned int DEG = (1==N ? 1 : N-1);
    math::PolynomialGeneric<T> retVal(DEG);

    if ( 1==N )
    {
        /*
            Handling of a zero degree polynimial ( p(x) = c )
            Its derivative is 0.
        */
        retVal.coef.at(0) = NumericUtil<T>::ZERO;
        return retVal;
    }

    // For polynomials of higher degree (>0) apply the formula above:
    for ( unsigned int i=0; i<DEG; i++ )
    {
        retVal.coef.at(i) = static_cast<T>(i+1)*coef.at(i+1);
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
    const unsigned int N = coef.size();
    math::PolynomialGeneric<T> retVal(N+1);

    // coef(0) is an arbitrary value, given as c:
    retVal.coef.at(0) = c;

    // for other coefficients, apply the formula above:
    for ( unsigned int i=0; i<N; i++ )
    {
        retVal.coef.at(i+1) = coef.at(i)/static_cast<T>(i+1);
    }

    return retVal;
}

/**
 * Addition operator (+) of two polynomials.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param poly - polynomial to be added to this one
 *
 * @return *this + poly
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::operator+(const math::PolynomialGeneric<T>& poly) const throw (math::PolynomialException)
{

    try
    {
        const unsigned int nthis = this->coef.size();
        const unsigned int npoly = poly.coef.size();

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
        const unsigned int nmax = ( nthis>=npoly ? nthis : npoly );

        math::PolynomialGeneric<T> retVal(nmax);

        // Note that the constructor assigns 1 to coef(nmax) which is not desirable at the moment:
        retVal.coef.at(nmax-1) = NumericUtil<T>::ZERO;

        /*
            Add coefficients of the same degree terms. Where 'i' exceeds size of any polynomial,
            consider its i^th coefficient as 0 (already set above)
        */
        for ( unsigned int i=0; i<nmax; i++ )
        {
            if ( i<nthis )
            {
                retVal.coef.at(i) = coef.at(i);
            }

            if ( i<npoly )
            {
                retVal.coef.at(i) += poly.coef.at(i);
            }
        }

        retVal.reduce();
        return retVal;
    }
    catch ( std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}

/**
 * Addition operator (+=) that adds a polynomial to 'this' and assigns the sum to itself.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param poly - polynomaial to be added to this one
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator+=(const math::PolynomialGeneric<T> poly) throw (math::PolynomialException)
{
    // For a definiton of polynomial addition, see operator+
    try
    {
        const unsigned int nthis = this->coef.size();
        const unsigned int npoly = poly.coef.size();

        // If 'poly' is of higher degree,
        // insert the appropriate number of coefficients and set them to 0:
        if ( nthis<npoly )
        {
            coef.insert(coef.begin()+nthis, npoly-nthis, math::NumericUtil<T>::ZERO);
        }

        // ... and perform addition of same degree terms' coefficients
        for ( unsigned int i=0; i<npoly; i++ )
        {
            coef.at(i) += poly.coef.at(i);
        }
    }
    catch ( std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }

    reduce();
    return *this;
}

/**
 * Subtraction operator (-) of two polynomials.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param poly - polynomial to be subtracted from this one
 *
 * @return *this - poly
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::operator-(const math::PolynomialGeneric<T>& poly) const throw (math::PolynomialException)
{
    try
    {
        const unsigned int nthis = this->coef.size();
        const unsigned int npoly = poly.coef.size();

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

        const unsigned int nmax = ( nthis>=npoly ? nthis : npoly );

        math::PolynomialGeneric<T> retVal(nmax);

        // Note that the constructor assigns 1 to coef(nmax) which is not desirable at the moment:
        retVal.coef.at(nmax-1) = math::NumericUtil<T>::ZERO;

        /*
            Subtract coefficients of the same degree terms. Where 'i' exceeds size of any polynomial,
            consider its ith coefficient as 0 (already set above)
        */
        for ( unsigned int i=0; i<nmax; i++ )
        {
            if ( i<nthis )
            {
                retVal.coef.at(i) = coef.at(i);
            }

            if ( i<npoly )
            {
                retVal.coef.at(i) -= poly.coef.at(i);
            }
        }

        retVal.reduce();
        return retVal;
    }
    catch ( std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}

/**
 * Subtraction operator (-=) that subtracts a polynomial from 'this' and assigns
 * the difference to itself.
 *
 * @note Polynomials can be of different degrees.
 *
 * @param poly - polynomaial to be subtracted from this one
 *
 * @return reference to itself
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T>& math::PolynomialGeneric<T>::operator-=(const math::PolynomialGeneric<T> poly) throw (math::PolynomialException)
{
    // For a definiton of polynomial subtraction, see operator-
    try
    {
        const unsigned int nthis = this->coef.size();
        const unsigned int npoly = poly.coef.size();

        // If poly is of higher degree, insert appropriate number of coefficients and set them to 0:
        if ( nthis<npoly )
        {
            coef.insert(coef.begin()+nthis, npoly-nthis, math::NumericUtil<T>::ZERO);
        }

        // ... and perform addition of same degree terms' coefficients
        for ( unsigned int i=0; i<npoly; i++ )
        {
            coef.at(i) -= poly.coef.at(i);
        }
    }
    catch ( std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }

    reduce();
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

        // Just negate each coefficient:
        for ( unsigned int i=0; i<coef.size(); i++ )
        {
            retVal.coef.at(i) = -coef.at(i);
        }

        // no need to reduce
        return retVal;
    }
    catch ( std::bad_alloc& ba )
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
 * @param poly - polynomial to be multiplied by this one
 *
 * @return *this * poly
 *
 * @throw PolynomialException if the product would be too large or if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::operator*(const math::PolynomialGeneric<T>& poly) const throw (math::PolynomialException)
{
    /*
        As explained at http://www.mathworks.com/help/matlab/ref/conv.html, multiplication of
        polynomials is equivalent to convolution of vectors.

        If Np is size of p(x) and Nq is size of q(x), the product's size will be:
          N = Np + Nq - 1.

        Coeffiecients of prod(x) = p(x) * q(x), rewritten for the library's order of coefficients,
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
        const unsigned int nthis = this->coef.size();
        const unsigned int npoly = poly.coef.size();

        // Size of the product polynomial:
        const unsigned long int N = nthis + npoly - 1;

        /*
            Unsigned int is used internally for manipulation of polynomials. This should be more
            than enough for most real life applications. At polynomial multiplication it is possible
            that product's number of coefficients exceeds the upper limit of unsigned int. Behaviour
            of most methods in this class would be wrong in this case. For that reason, this check of
            unsigned int's limit is performed.
        */
        if ( N > std::numeric_limits<unsigned int>::max() )
        {
            throw math::PolynomialException(math::PolynomialException::TOO_LARGE);
        }

        math::PolynomialGeneric<T> retVal(N);

        // Each product's coefficient...
        for ( unsigned int i=0; i<N; i++ )
        {
            // ... is a sum of products as specified above
            T temp = math::NumericUtil<T>::ZERO;
            for ( unsigned int j=0; j<N; j++ )
            {
                // if any index would point out of respective polynomial's range,
                // treat it as a zero (i.e. skip this iteration)
                if ( j>i || j>=nthis || (i-j)>=npoly )
                {
                    continue;   // for j
                }

                temp += coef.at(j) * poly.coef.at(i-j);
            }

            retVal.coef.at(i) = temp;
        }

        retVal.reduce();
        return retVal;
    }
    catch ( std::bad_alloc& ba )
    {
        throw math::PolynomialException(math::PolynomialException::OUT_OF_MEMORY);
    }
}

/**
 * Multiplication operator (*=) that multiplies two polynomials and assigns the product to itself.
 *
 * @note Polynomials can be of different degrees.
 * @note Multiplication of polyniomials is commutative.
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

    // perform a regular polynomial multiplication and copy the products's coefficients to itself
    math::PolynomialGeneric<T> temp = *this * poly;
    copyCoefs(temp.coef);

    // Note that copyCoefs would throw an exception in case of unsuccessful allocation of memory

    // 'reduce' already performed by operator*
    return *this;
}

/**
 * Multiplication operator (*) for multiplication of a polynomial and a scalar.
 *
 * @param sc - scalar to be multiplied by the polynomial
 *
 * @return *this * sc
 *
 * @throw PolynomialException if allocation of memory fails
 */
template<class T>
math::PolynomialGeneric<T> math::PolynomialGeneric<T>::operator*(const T& sc) const throw (math::PolynomialException)
{
    /*
        Multiplication of a polynomial by a scalar is trivial:
        Each cofficient is multiplied by the scalar.

                         Np
                        -----
                        \               i
            p(x) *sc =  |    sc * pi * x
                        /
                        -----
                         i=0
    */

    const unsigned int N = this->coef.size();
    math::PolynomialGeneric<T> retVal(*this);

    for ( unsigned int i=0; i<N; i++ )
    {
        retVal.coef.at(i) *= sc;
    }

    // applicable when sc==o
    retVal.reduce();
    return retVal;
}

/**
  * Multiplication operator (*) of a scalar and a polynomial.
  * This operation is commutative and does the same as operator*(scalar).
  * Since the first operand is not a polynomial, it must be implemented as
  * a friend function.
  *
  * @param sc - scalar
  * @param poly - polynomial
  *
  * @return sc * poly
  *
  * @throw PolynomialException if allocation of memory fails
  */
template<class T>
math::PolynomialGeneric<T> math::operator* (const T& sc, const math::PolynomialGeneric<T>& poly) throw (math::PolynomialException)
{
    // 'reduce' performed by operator*(poly, sc)
    return (poly * sc);
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
    for (unsigned int i=0; i<coef.size(); i++ )
    {
        /*
            A space will be displayed between terms to better distinguish them.
            The first (lowest order) term is an exception.
        */
        if ( i>0 )
        {
            str << ' ';
        }

        /*
            If the coefficient is negative, '-' will be displayed automatically.
            This is not true for '+' that must must be displayed explicitly.
        */
        if ( i>0 && coef.at(i)>=0 )
        {
            str << '+';
        }

        str << coef.at(i);

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

/**
 * Destructor
 */
template<class T>
math::PolynomialGeneric<T>::~PolynomialGeneric()
{
    // Vector's destructors would probably clean up this automatically.
    // Anyway, let us clear the vector, just to be aware of allocated resources.
    coef.clear();

    // Other dynamically allocated memory (via malloc or new) should be freed here.
    // There are no other resources to release.
}
