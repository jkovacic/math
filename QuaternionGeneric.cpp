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
 * @file QuaternionGeneric.cpp
 *
 * Implementation of the class QuaternionGeneric.
 *
 * Please note that quaternions make sense if T represents a real-number type
 * (e.g. float or double), in case of any other type, results of many methods
 * do not make much sense, some methods even throw an exception imediately.
 * If you wish to template it with any "officially" unsupported type, you should
 * specialize those methods (i.e. introduce your own implementation of sqrt for the
 * desired type).
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 *
 * @author Jernej Kovacic
 */


#include <cmath>

// Deliberately there is no #include "QuaternionGeneric.h" !
#include "NumericUtil.h"


/**
 * Constructor.
 * Creates an instance of a quaternion.
 *
 * @param o - scalar component of the quaternion (default: 0)
 * @param i - component i of the quaternion (default: 0)
 * @param j - component j of the quaternion (default: 0)
 * @param k - component j of the quaternion (default: 0)
 */
template<class T>
math::QuaternionGeneric<T>::QuaternionGeneric(const T& o, const T& i, const T& j, const T& k) :
    quat_o(o), quat_i(i), quat_j(j), quat_k(k)
{
    // Nothing else to do
}

/**
 * Copy constructor.
 * Creates an exact copy of the original quaternion q.
 *
 * @param q - original quaternion to be copied into this one
 */
template<class T>
math::QuaternionGeneric<T>::QuaternionGeneric(const math::QuaternionGeneric<T>& q) :
    quat_o(q.quat_o), quat_i(q.quat_i), quat_j(q.quat_j), quat_k(q.quat_k)
{
    // Nothing else to do
}

/**
 * Assignment operator (=)
 *
 * @param q - a quaternion to be copied into this one
 *
 * @return reference to this
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::operator=(const math::QuaternionGeneric<T>& q)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &q )
    {
        return *this;
    }

    // otherwise copy the components:
    this->quat_o = q.quat_o;
    this->quat_i = q.quat_i;
    this->quat_j = q.quat_k;
    this->quat_k = q.quat_k;

    // ... and return the reference to itself:
    return *this;
}

/**
 * @return scalar componet ('1') of the quaternion
 */
template<class T>
T math::QuaternionGeneric<T>::getOne() const
{
    return this->quat_o;
}

/**
 * @return component 'i' of the quaternion
 */
template<class T>
T math::QuaternionGeneric<T>::getI() const
{
    return this->quat_i;
}

/**
 * @return component 'j' of the quaternion
 */
template<class T>
T math::QuaternionGeneric<T>::getJ() const
{
    return this->quat_j;
}

/**
 * @return component 'k' of the quaternion
 */
template<class T>
T math::QuaternionGeneric<T>::getK() const
{
    return this->quat_k;
}


/**
 * Assigns values of all quaternion's components.
 *
 * @param o - scalar component of the quaternion (default: 0)
 * @param i - component i of the quaternion (default: 0)
 * @param j - component j of the quaternion (default: 0)
 * @param k - component j of the quaternion (default: 0)
 *
 * @return reference to itself
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::set(const T& o, const T& i, const T& j, const T& k)
{
    this->quat_o = o;
    this->quat_i = i;
    this->quat_j = j;
    this->quat_k = k;

    // return the reference to itself:
    return *this;
}

/**
 * Assigns the scalar component of the quaternion.
 *
 * @param o - value of the scalar component ('1') (default: 0)
 *
 * @return reference to itself
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::setOne(const T& o)
{
    this->quat_o = o;

    // return the reference to itself:
    return *this;
}

/**
 * Assigns the component 'i' of the quaternion.
 *
 * @param i - value of the component 'i' (default: 0)
 *
 * @return reference to itself
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::setI(const T& i)
{
    this->quat_i = i;

    // return the reference to itself:
    return *this;
}

/**
 * Assigns the component 'j' of the quaternion.
 *
 * @param j - value of the component 'j' (default: 0)
 *
 * @return reference to itself
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::setJ(const T& j)
{
    this->quat_j = j;

    // return the reference to itself:
    return *this;
}

/**
 * Assigns the component 'k' of the quaternion.
 *
 * @param k - value of the component 'k' (default: 0)
 *
 * @return reference to itself
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::setK(const T& k)
{
    this->quat_k = k;

    // return the reference to itself:
    return *this;
}

/**
 * Outputs the quaternion to the stream 'str' in form '(a+bi+cj+dk)'
 *
 * @param str (default cout): output stream, where the quaternion will be dislayed
 */
template<class T>
void math::QuaternionGeneric<T>::display(std::ostream& str) const
{
    /*
        Primarily the method was introduced for brief unit testing purposes
        and not much effort was invested into a visually "nice" output
    */

    // start with an opening bracket:
    str << '(';

    // output he first component ('1')
    str << quat_o;

    /*
        For the other components, display '+' if the component is positive (ot zero).
        If the component's value is negative, '-' will be displayed automatically.
        After the numeric value append the component's name ('i', 'j' or 'k').
    */

    // Start with 'i':
    if ( quat_i >= math::NumericUtil<T>::ZERO )
    {
        str << '+';
    }
    str << quat_i << 'i';

    // do the same for 'j':
    if ( quat_j >= math::NumericUtil<T>::ZERO )
    {
        str << '+';
    }
    str << quat_j << 'j';

    // and finally the same for 'k':
    if ( quat_k >= math::NumericUtil<T>::ZERO )
    {
        str << '+';
    }
    str << quat_k << 'k';

    // finish with a closing bracket
    str << ')';
}

/**
 * A friend function that outputs the quaternion to an output stream
 *
 * @param output - stream to write to
 * @param q - quaternion to be displayed
 *
 * @return reference of output stream (i.e. 'output')
 */
template<class T>
std::ostream& math::operator<<(std::ostream& output, const math::QuaternionGeneric<T>& q)
{
    // just pass the quaternion to QuaternionGeneric::display()...
    q.display(output);

    // ... and return reference of the stream
    return output;
}

/**
 * Addition operator (+) of two quaternions.
 *
 * @param q - quaternion to be added to this one
 *
 * @return *this + q
 */
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::operator+(const math::QuaternionGeneric<T>& q) const
{
    /*
        Addition of quaternions is trivial:
        (a1 + b1*i + c1*j + d1*k) + (a2 + b2*i + c2*j + d2*k) =
        = ( (a1+a2) + (b1+b2)*i + (c1+c2)*j + (d1+d2)*k )
    */

    return math::QuaternionGeneric<T>( quat_o + q.quat_o,
                                       quat_i + q.quat_i,
                                       quat_j + q.quat_j,
                                       quat_k + q.quat_k );
}

/**
 * Addition operator (+=) that adds a quaternion to this and assigns the sum to itself.
 *
 * @param q - quaternion to be added to this one
 *
 * @return reference to this
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::operator+=(const math::QuaternionGeneric<T>& q)
{
    // For a definition of quaternion addition, see operator+

    // assign the updated values directly to relevant components:
    this->quat_o += q.quat_o;
    this->quat_i += q.quat_i;
    this->quat_j += q.quat_j;
    this->quat_k += q.quat_k;

    // return the reference to itself
    return *this;
}

/**
 * Subtraction operator (-) of two quaternions.
 *
 * @param q - quaternion to be subtracted to this one
 *
 * @return *this - q
 */
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::operator-(const math::QuaternionGeneric<T>& q) const
{
    /*
        Subtraction of quaternions is trivial:
        (a1 + b1*i + c1*j + d1*k) - (a2 + b2*i + c2*j + d2*k) =
        = ( (a1-a2) + (b1-b2)*i + (c1-c2)*j + (d1-d2)*k )
    */

    return math::QuaternionGeneric<T>( quat_o - q.quat_o,
                                       quat_i - q.quat_i,
                                       quat_j - q.quat_j,
                                       quat_k - q.quat_k );
}

/**
 * Subtraction operator (-=) that subtracts a quaternion from this and assigns the difference to itself.
 *
 * @param q - quaternion to be subtracted from this one
 *
 * @return reference to this
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::operator-=(const math::QuaternionGeneric<T>& q)
{
    // For a definition of quaternion subtraction, see operator-

    // assign the updated values directly to relevant components:
    this->quat_o -= q.quat_o;
    this->quat_i -= q.quat_i;
    this->quat_j -= q.quat_j;
    this->quat_k -= q.quat_k;

    // return the reference to itself
    return *this;
}

/**
 * Multiplication operator (*) of two quaternions.
 *
 * Note that quaternion multiplication is not comutative (q*p != p*q).
 *
 * @param q- quaternion to be multiplied by this
 *
 * @return this * q
*/
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::operator*(const math::QuaternionGeneric<T>& q) const
{
    /*
        From the following definitions:
          i*i = j*j = k*k = -1,
          i*j = k, j*i = -k, j*k = i, k*j = -i, k*i = j and i*k = -j,
        the following formula can be quickly derived:

        (a1 + b1*i + c1*j + d1*k) * (a2 + b2*i + c2*j + d2*k) =
        =  (a1*a2 - b1*b2 - c1*c2 - d1*d2)     +
        +  (a1*b2 + b1*a2 + c1*d2 - d1*c2) * i +
        +  (a1*c2 - b1*d2 + c1*a2 + d1*b2) * j +
        +  (a1*d2 + b1*c2 - c1*b2 + d1*a2) * k

        Note: The following script for GNU Octave or Matlab can be used
        for a quick unit test of the function:
        http://mind.cog.jhu.edu/courses/680/octave/Installers/Octave/Octave.OSX10.6/Applications/MATLAB_R2009b.app/toolbox/aero/aero/quatmultiply.m
    */

    return math::QuaternionGeneric<T>(
            quat_o*q.quat_o - quat_i*q.quat_i - quat_j*q.quat_j - quat_k*q.quat_k,
            quat_o*q.quat_i + quat_i*q.quat_o + quat_j*q.quat_k - quat_k*q.quat_j,
            quat_o*q.quat_j - quat_i*q.quat_k + quat_j*q.quat_o + quat_k*q.quat_i,
            quat_o*q.quat_k + quat_i*q.quat_j - quat_j*q.quat_i + quat_k*q.quat_o );
}

/**
 * Multiplication operator (*=) that multiplies a quaternion to this and assigns the product to itself.
 *
 * Note that quaternion multiplication is not comutative (q*p != p*q).
 *
 * @param q - quaternion to be multiplied by this one
 *
 * @return reference to this (=this * q)
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::operator*=(const math::QuaternionGeneric<T>& q)
{
    // For a definition of quaternion multiplication, see operator*

    /*
        Any method would require storing componets' values in 4 variables. Just performing
        a general quaternion multiplication (by operator*) and instantiation of a temporary
        variable consumes the same amount of memory. Additionally it improves maintainability
        and reduces the risk of typing errors.
    */

    const math::QuaternionGeneric<T> p( *this * q );

    //copy the temporary variable's values into this:
    set(p.quat_o, p.quat_i, p.quat_j, p.quat_k);

    // return the reference to itself
    return *this;
}

/**
 * Multiplication operator (*) for multiplication of a quaternion and a scalar.
 * It multiplies each quaternion's component by the scalar.
 *
 * @param scalar - scalar to be multiplied by the quaternion
 *
 * @return this * scalar
 */
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::operator*(const T& scalar) const
{
    /*
        From the definition of quaternion multiplication (see operator*),
        one can quickly derive the following simplified formula:
          (a+ b*i + c*j + d*k) * s = (a*s + (b*s)*i + (c*s)*j + (d*s)*k))
    */

    return math::QuaternionGeneric<T>( quat_o * scalar,
                                       quat_i * scalar,
                                       quat_j * scalar,
                                       quat_k * scalar );
}

/**
 * Multiplication operator (*) of a scalar and a quaternion.
 * In general ( if T represents a real number) this operation is comutative
 * and does the same as operator*(scalar).
 * Since the first operand is not a quaternion, it must be implemented as
 * a friend function.
 *
 * @param scalar
 * @param q - quaternion
 *
 * @return scalar * q
 */
template<class T>
math::QuaternionGeneric<T> math::operator*(const T& scalar, const math::QuaternionGeneric<T>& q)
{
    /*
        In general multiplication of a scalar and quaternion is comutative.
        If this is not a case, implement a specialization.
    */
    return (q * scalar);
}

/**
 * Unary negation operator (-)
 *
 * @return -this
*/
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::operator-() const
{
    /*
        Negation of a quaternion is trivial: just negate each component.
    */

    return math::QuaternionGeneric<T>( -quat_o,
                                       -quat_i,
                                       -quat_j,
                                       -quat_k );
}

/**
 * Conjugation of the quaternion.
 *
 * @return conjugation of this
 */
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::conj() const
{
    /*
        From the definition of quaternion conjugation:
        q* = 0.5 * (q + i*q*i + j*q*j + k*q*k)
        one can quickly derive the following simplified formula:
        (a + b*i + c*j + d*k)* = (a - b*i - c*j -d*k)

        In other words, the "scalar" part remains unmodified, the other three
        "vector" components are negated.
    */

    return math::QuaternionGeneric<T>(quat_o, -quat_i, -quat_j, -quat_k);
}

/**
 * Conjugate the quaternion and write all changes into it
 * (disregard the original quaternion).
 *
 * @return reference to this
 */
template<class T>
math::QuaternionGeneric<T>& math::QuaternionGeneric<T>::conjugated()
{
    // For a definition of quaternion conjugation, see conj().

    // quat_o remains unmodified

    // the other members are assigned their opposite values:
    quat_i = -quat_i;
    quat_j = -quat_j;
    quat_k = -quat_k;

    // return the reference to itself
    return *this;
}

/*
 * A utilty function that calculates a sum of components' squares.
 * The function is caled  by other public methods, such as norm, reciprocal, unit, etc.
 *
 * As the function is hort, it is declared as inline to slightly reduce overhead
 *
 * @return sum of all components' squares
 */
template<class T>
inline T math::QuaternionGeneric<T>::sqsum() const
{
    return (quat_o*quat_o + quat_i*quat_i + quat_j*quat_j + quat_k*quat_k);
}

/**
 * Norm of the quaternion.
 *
 * @return ||this||
 *
 * @throw QuaternionException if operation is not supported for type T
 */
template<class T>
T math::QuaternionGeneric<T>::norm() const throw (math::QuaternionException)
{
    /*
        Norm of a quaternion is defined as:

        ||q|| = sqrt( q*conj(q) ) = sqrt( conj(q)*q )

        From the definition of quatenion multiplication,
        one can quickly derive the following simplified formula:

        ||(a+b*i+c*j+d*k)|| = sqrt(a*a + b*b + c*c + d*d)
    */

    /*
        This operation is only supported for T=float or T=double.
        Specialized implementation is provided below for these two types.

        For any other type, the operation is (probably) not supported
        as sqrt may not be defined for the type. In such a case throw an
        exception immediately.
    */

    throw math::QuaternionException(math::QuaternionException::UNSUPPORTED_TYPE);

    // will never execute, but some compilers may produce a warning if nothing is returned
    return math::NumericUtil<T>::ZERO;
}

/*
    Specialization of norm() for float, double and long double.
    All three specializations are very similar and only differ in types of the
    output value and internal variables.
    For easier maintainability, the specialization will be implemented
    only once using a parameterized #define.
*/

#define _MATH_QUATERNIONGENERIC_SPECIALIZED_NORM(FD) \
template<> \
FD math::QuaternionGeneric<FD>::norm() const throw (math::QuaternionException) \
{ \
    return std::sqrt(sqsum()); \
}
// end of #define

// the actual specialization for float:
_MATH_QUATERNIONGENERIC_SPECIALIZED_NORM(float)

// for double:
_MATH_QUATERNIONGENERIC_SPECIALIZED_NORM(double)

// and for long double:
_MATH_QUATERNIONGENERIC_SPECIALIZED_NORM(long double)

// definition of _MATH_QUATERNIONGENERIC_SPECIALIZED_NORM not needed anymore, #undef it
#undef _MATH_QUATERNIONGENERIC_SPECIALIZED_NORM


/**
 * Transforms 'this' into a unit quaternion (its norm is equal to 1).
 *
 * @return U(this)
 *
 * @throw QuaternionException if operation is not supported for type T or if 'this' is a zero-quaternion
 */
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::unit() const throw (math::QuaternionException)
{
    /*
        A unit quaternion is the quaternion, divided by its norm:

        U(q) = q / ||q||

        Its norm is equal to 1.

        A unit quaternion for the zero quaternion (norm=0) cannot be calculated.
    */

    /*
        This operation is only supported for T=float or T=double.
        Specialized implementation is provided below for these two types.

        For any other type, the operation is (probably) not supported
        as norm may not be defined for the type. In such a case throw an
        exception immediately.
    */
    throw math::QuaternionException(math::QuaternionException::UNSUPPORTED_TYPE);

    // will never execute, but some compilers may produce a warning if nothing is returned
    return math::QuaternionGeneric<T>();
}

/*
    Specialization of unit() for float, double and long double.
    All three specializations are very similar and only differ in types of the
    output value and internal variables.
    For easier maintainability, the specialization will be implemented
    only once using a parameterized #define.
*/

#define _MATH_QUATERNIONGENERIC_SPECIALIZED_UNIT(FD) \
template<> \
math::QuaternionGeneric<FD> math::QuaternionGeneric<FD>::unit() const throw (math::QuaternionException) \
{ \
    const FD norm = this->norm(); \
    if ( true == math::NumericUtil<FD>::isZero(norm) ) \
    { \
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO); \
    } \
    return ( *this * (math::NumericUtil<FD>::ONE/norm) ); \
}
// end of #define

// specialization for float:
_MATH_QUATERNIONGENERIC_SPECIALIZED_UNIT(float)

// for double:
_MATH_QUATERNIONGENERIC_SPECIALIZED_UNIT(double)

// and for long double:
_MATH_QUATERNIONGENERIC_SPECIALIZED_UNIT(long double)

// definition of _MATH_QUATERNIONGENERIC_SPECIALIZED_UNIT not needed anymore, #undef it
#undef _MATH_QUATERNIONGENERIC_SPECIALIZED_UNIT


/**
 * Reciprocal quaternion (q^(-1)), satisfying the condition:
 * q * q^(-1)  =  q^(-1) * q  = 1
 *
 * Note: even though reciprocal should only be defined for types T that do support norm
 * (float and type), norm's square can be calculated for many other types as well. For that
 * reason, the operation is implemented for any type. However, T must support division (operator/),
 * otherwise the class will not compile.
 *
 * @return this^(-1)
 *
 * @throw QuaternionException if 'this' is a zero-quaternion
 */
template<class T>
math::QuaternionGeneric<T> math::QuaternionGeneric<T>::reciprocal() const throw (math::QuaternionException)
{
    /*
        q^(-1) is a reciprocal quaternion of q if the following condition is satisfied;

            q * q^(-1) = q^(-1) * q = 1

        Reciprocal of q is defined as:

            q^(-1) = q* / ||q||^2

        The following formula can be derived from it:

                                      a - b*i - c*j - d*k
            (a+b*i+c*j+d*k)^(-1) = -------------------------
                                     a^2 + b^2 + c^2 + d^2
    */

    const T nsq = sqsum();

    // avoid possible division by zero
    if ( true == math::NumericUtil<T>::isZero(nsq) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    return math::QuaternionGeneric<T>(quat_o/nsq, -quat_i/nsq, -quat_j/nsq, -quat_k/nsq);
}
