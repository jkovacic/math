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
 * Implementation of the class QuaternionGeneric.
 *
 * Please note that quaternions make sense if F represents a real-number type
 * (e.g. float, double or long double), in case of any other type, results of
 * many methods do not make much sense, some methods even throw an exception
 * immediately. If you wish to template it with any "officially" unsupported
 * type, you should specialize those methods (i.e. introduce your own
 * implementation of sqrt for the desired type).
 */


#include <cmath>

// no #include "QuaternionGeneric.hpp" !!!
#include "util/NumericUtil.hpp"
#include "../settings/omp_settings.h"


/*
 * Notes about parallelization:
 *
 * In most functions, each quaternion's element can be calculated
 * independently from the others, hence it would be possible to
 * apply OpenMP's directive "sections" to parallelize most algorithms into
 * 4 threads, one for each element. However, this is implemented only in
 * operator*() which requires a little bit more operations for each element.
 * It could be done in most other functions in a very similar manner, however
 * it doesn't make much sense as each element's operations are very basic.
 *
 * Furthermore, parallelization of any algorithm in this file can be
 * enabled/disabled by setting the macro OMP_QUAT_PARALLELIZE in omp_settings.h
 * If its value is set to 0, any parallelization within this file is completely
 * disabled, otherwise it is enabled. By default this macro is set
 * to disable parallelization.
 */


/**
 * Constructor.
 * Creates an instance of a quaternion.
 *
 * @param o - scalar component of the quaternion (default: 0)
 * @param i - component i of the quaternion (default: 0)
 * @param j - component j of the quaternion (default: 0)
 * @param k - component j of the quaternion (default: 0)
 */
template <typename F>
math::QuaternionGeneric<F>::QuaternionGeneric(const F& o, const F& i, const F& j, const F& k) :
    m_o(o), m_i(i), m_j(j), m_k(k)
{
    // Nothing else to do
}


/**
 * Copy constructor.
 * Creates an exact copy of the original quaternion q.
 *
 * @param q - original quaternion to be copied into this one
 */
template <typename F>
math::QuaternionGeneric<F>::QuaternionGeneric(const math::QuaternionGeneric<F>& q) :
    m_o(q.m_o), m_i(q.m_i), m_j(q.m_j), m_k(q.m_k)
{
    // Nothing else to do
}


/**
 * Assignment operator (=)
 *
 * @param q - a quaternion to be copied into this one
 *
 * @return reference to itself
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator=(const math::QuaternionGeneric<F>& q)
{
    // Nothing to do, if attempting to assign itself
    if ( this == &q )
    {
        return *this;
    }

    // otherwise copy the components:
    this->m_o = q.m_o;
    this->m_i = q.m_i;
    this->m_j = q.m_k;
    this->m_k = q.m_k;

    // ... and return the reference to itself:
    return *this;
}


/**
 * Assignment operator (=) that copies the scalar to
 * the "real" component and assigns the other 3 "complex"
 * components to 0.
 *
 * @param sc - scalar to assign to the quaternion
 *
 * @return reference to itself
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator=(const F& sc)
{
    this->m_o = sc;
    this->m_i = static_cast<F>(0);
    this->m_j = static_cast<F>(0);
    this->m_k = static_cast<F>(0);

    return *this;
}


/**
 * @return scalar component ('1') of the quaternion
 */
template <typename F>
F math::QuaternionGeneric<F>::getOne() const
{
    return this->m_o;
}


/**
 * @return component 'i' of the quaternion
 */
template <typename F>
F math::QuaternionGeneric<F>::getI() const
{
    return this->m_i;
}


/**
 * @return component 'j' of the quaternion
 */
template <typename F>
F math::QuaternionGeneric<F>::getJ() const
{
    return this->m_j;
}


/**
 * @return component 'k' of the quaternion
 */
template <typename F>
F math::QuaternionGeneric<F>::getK() const
{
    return this->m_k;
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
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::set(const F& o, const F& i, const F& j, const F& k)
{
    this->m_o = o;
    this->m_i = i;
    this->m_j = j;
    this->m_k = k;

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
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::setOne(const F& o)
{
    this->m_o = o;

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
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::setI(const F& i)
{
    this->m_i = i;

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
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::setJ(const F& j)
{
    this->m_j = j;

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
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::setK(const F& k)
{
    this->m_k = k;

    // return the reference to itself:
    return *this;
}


/**
 * Outputs the quaternion to the stream 'str' in form '(a+bi+cj+dk)'
 *
 * @param str (default cout): output stream, where the quaternion will be displayed
 */
template <typename F>
void math::QuaternionGeneric<F>::display(std::ostream& str) const
{
    /*
     * Primarily the method was introduced for brief unit testing purposes
     * and not much effort was invested into a visually "nice" output
     */

    // start with an opening bracket:
    str << '(';

    // output he first component ('1')
    str << this->m_o;

    // For the other components, display each component's sign
    
    str << std::showpos;

    str << this->m_i << 'i';
    str << this->m_j << 'j';
    str << this->m_k << 'k';

    str << std::noshowpos;
    // finish with a closing bracket
    str << ')';
}


/**
 * Addition operator (+=) that adds a quaternion to this and assigns the sum to itself.
 *
 * @param q - quaternion to be added to this one
 *
 * @return reference to this
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator+=(const math::QuaternionGeneric<F>& q)
{
    // For a definition of quaternion addition, see operator+

    // assign the updated values directly to relevant components:
    this->m_o += q.m_o;
    this->m_i += q.m_i;
    this->m_j += q.m_j;
    this->m_k += q.m_k;

    // return the reference to itself
    return *this;
}


/**
 * Addition operator (+=) that adds a scalar to this and assigns the sum to itself.
 *
 * @param scalar - scalar to be added to this one
 *
 * @return reference to this 
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator+=(const F& scalar)
{
    this->m_o += scalar;
    
    return *this;
}


/**
 * Subtraction operator (-=) that subtracts a quaternion from this and assigns the difference to itself.
 *
 * @param q - quaternion to be subtracted from this one
 *
 * @return reference to this
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator-=(const math::QuaternionGeneric<F>& q)
{
    // For a definition of quaternion subtraction, see operator-

    // assign the updated values directly to relevant components:
    this->m_o -= q.m_o;
    this->m_i -= q.m_i;
    this->m_j -= q.m_j;
    this->m_k -= q.m_k;

    // return the reference to itself
    return *this;
}


/**
 * Subtraction operator (-=) that subtracts a scalar from this and assigns 
 * the difference to itself.
 *
 * @param scalar - scalar to be subtracted from this quaternion
 *
 * @return reference to this 
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator-=(const F& scalar)
{
    this->m_o -= scalar;
    
    return *this;
}


/**
 * Multiplication operator (*=) that multiplies a quaternion to this and assigns the product to itself.
 *
 * Note that quaternion multiplication is not commutative (q*p != p*q).
 *
 * @param q - quaternion to be multiplied by this one
 *
 * @return reference to this (=this * q)
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator*=(const math::QuaternionGeneric<F>& q)
{
    // For a definition of quaternion multiplication, see operator*

    /*
     * Any method would require storing components' values in 4 variables. Just performing
     * a general quaternion multiplication (by operator*) and instantiation of a temporary
     * variable consumes the same amount of memory. Additionally it improves maintainability
     * and reduces the risk of typing errors.
     */

    const math::QuaternionGeneric<F> p( *this * q );

    //copy the temporary variable's values into this:
    set(p.m_o, p.m_i, p.m_j, p.m_k);

    // return the reference to itself
    return *this;
}


/**
 * Multiplication operator (*=) that multiplies a quaternion by a scalar
 * and assigns the product to itself.
 * 
 * @param sc - scalar
 *
 * @return reference to itself
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator*=(const F& sc)
{
    // Multiply each component by the scalar
    this->m_o *= sc;
    this->m_i *= sc;
    this->m_j *= sc;
    this->m_k *= sc;

    return *this;
}


/**
 * Division operator (/=) that divides the quaternion by a scalar
 * and assigns the quotient to itself.
 * 
 * @param scalar
 * 
 * @return reference to itself
 *  
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::operator/=(
        const F& scalar
      ) throw(math::QuaternionException)
{
    if ( true == math::NumericUtil::isZero<F>(scalar) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    const F rec = static_cast<F>(1) / scalar;

    this->m_o *= rec;
    this->m_i *= rec;
    this->m_j *= rec;
    this->m_k *= rec;

    return *this;
}


/**
 * Right division "operator" that performs right division of 'this'
 * (this * rec(q)) and assigns the quotient to itself.
 *
 * @param q - quaternion to divide 'this'
 *
 * @return reference to itself
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::rdiv_(const math::QuaternionGeneric<F>& q) throw(math::QuaternionException)
{
    if ( true == math::NumericUtil::isZero<F>(q.__sqsum()) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    *this *= q.reciprocal();

    return *this;
}


/**
 * Right division "operator" that divides 'this' by a scalar
 * and assigns the quotient to itself.
 *
 * @note This operator is actually equivalent to left division by a scalar
 *
 * @param sc - scalar to divide "'his'
 *
 * @return reference to itself
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::rdiv_(const F& sc) throw(math::QuaternionException)
{
    return this->operator/=(sc);
}


/**
 * Left division "operator" that performs left division of 'this'
 * (rec(q) * this) and assigns the quotient to itself.
 *
 * @param q - quaternion to divide 'this'
 *
 * @return reference to itself
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::ldiv_(const math::QuaternionGeneric<F>& q) throw(math::QuaternionException)
{
    math::QuaternionGeneric<F> lquot = math::ldiv(*this, q);

    this->m_o = lquot.m_o;
    this->m_i = lquot.m_i;
    this->m_j = lquot.m_j;
    this->m_k = lquot.m_k;

    return *this;
}


/**
 * Left division "operator" that divides 'this' by a scalar
 * and assigns the quotient to itself.
 *
 * @note This operator is actually equivalent to right division by a scalar
 *
 * @param sc - scalar to divide 'this'
 *
 * @return reference to itself
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::ldiv_(const F& sc) throw(math::QuaternionException)
{
    return this->operator/=(sc);
}


/**
 * Conjugate the quaternion and write all changes into it
 * (disregard the original quaternion).
 *
 * @return reference to itself
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::conj_()
{
    /*
     * From the definition of quaternion conjugation:
     *   q* = 0.5 * (q + i*q*i + j*q*j + k*q*k)
     * one can quickly derive the following simplified formula:
     *   (a + b*i + c*j + d*k)* = (a - b*i - c*j -d*k)
     *
     * In other words, the "scalar" part remains unmodified, the other three
     * "vector" components are negated.
     */

    // 'm_o' remains unmodified

    // the other members are assigned their opposite values:
    this->m_i = -this->m_i;
    this->m_j = -this->m_j;
    this->m_k = -this->m_k;

    // return the reference to itself
    return *this;
}


/**
 * Conjugation of the quaternion.
 *
 * @return conjugation of this
 */
template <typename F>
math::QuaternionGeneric<F> math::QuaternionGeneric<F>::conj() const
{
    math::QuaternionGeneric<F> qret(*this);
    qret.conj_();
    return qret;
}


/**
 * Norm of the quaternion.
 *
 * @return ||this||
 *
 * @throw QuaternionException if operation is not supported for type T
 */
template <typename F>
F math::QuaternionGeneric<F>::norm() const throw (math::QuaternionException)
{
    /*
     * Norm of a quaternion is defined as:
     *
     *   ||q|| = sqrt( q*conj(q) ) = sqrt( conj(q)*q )
     *
     * From the definition of quaternion multiplication,
     * one can quickly derive the following simplified formula:
     *
     *   ||(a+b*i+c*j+d*k)|| = sqrt(a*a + b*b + c*c + d*d)
     */

    return std::sqrt( this->__sqsum() );;
}


/**
 * Transforms itself into a unit quaternion (its norm is equal to 1).
 *
 * @return reference to itself
 *
 * @throw QuaternionException if 'this' is a zero-quaternion
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::unit_() throw (math::QuaternionException)
{
    /*
     * A unit quaternion is the quaternion, divided by its norm:
     *
     *   U(q) = q / ||q||
     *
     * Its norm is equal to 1.
     *
     * A unit quaternion for the zero quaternion (norm=0) cannot be calculated.
     */

    const F norm = this->norm();

    if ( true == math::NumericUtil::isZero<F>(norm) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    const F normrec = static_cast<F>(1) / norm;

    this->m_o *= normrec;
    this->m_i *= normrec;
    this->m_j *= normrec;
    this->m_k *= normrec;

    return *this;
}


/**
 * Transforms 'this' into a unit quaternion (its norm is equal to 1).
 *
 * @return U(this)
 *
 * @throw QuaternionException if 'this' is a zero-quaternion
 */
template <typename F>
math::QuaternionGeneric<F> math::QuaternionGeneric<F>::unit() const throw (math::QuaternionException)
{
    math::QuaternionGeneric<F> qret(*this);
    qret.unit_();
    return qret;
}


/**
 * Transforms itself to its reciprocal quaternion (q^(-1)),
 * satisfying the condition:
 * q * q^(-1)  =  q^(-1) * q  = 1
 *
 * @return reference to itself
 *
 * @throw QuaternionException if 'this' is a zero-quaternion
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::reciprocal_() throw (math::QuaternionException)
{
    /*
     * q^(-1) is a reciprocal quaternion of q if the following condition is satisfied;
     *
     *   q * q^(-1) = q^(-1) * q = 1
     *
     * Reciprocal of q is defined as:
     *
     *   q^(-1) = q* / ||q||^2
     *
     * The following formula can be derived from it:
     *
     *                            a - b*i - c*j - d*k
     *   (a+b*i+c*j+d*k)^(-1) = -------------------------
     *                            a^2 + b^2 + c^2 + d^2
     */

    const F nsq = this->__sqsum();

    // avoid possible division by zero
    if ( true == math::NumericUtil::isZero<F>(nsq) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    const F nsqRec = static_cast<F>(1) / nsq;

    this->m_o *= nsqRec,
    this->m_i *= -nsqRec,
    this->m_j *= -nsqRec,
    this->m_k *= -nsqRec;

    return *this;
}


/**
 * Reciprocal quaternion (q^(-1)), satisfying the condition:
 * q * q^(-1)  =  q^(-1) * q  = 1
 *
 * @return this^(-1)
 *
 * @throw QuaternionException if 'this' is a zero-quaternion
 */
template <typename F>
math::QuaternionGeneric<F> math::QuaternionGeneric<F>::reciprocal() const throw (math::QuaternionException)
{
    math::QuaternionGeneric<F> qret(*this);
    qret.reciprocal_();
    return qret;
}


/**
 * "Rounds" all small quaternion's components (with the absolute value below
 * the given 'eps') to 0.
 * 
 * @param eps - threshold to determine whether each component is "rounded" to 0
 * 
 * @return reference to itself
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::roundSmallElements_(const F& eps)
{
    this->m_o = math::NumericUtil::smallValToZero<F>(this->m_o, eps);
    this->m_i = math::NumericUtil::smallValToZero<F>(this->m_i, eps);
    this->m_j = math::NumericUtil::smallValToZero<F>(this->m_j, eps);
    this->m_k = math::NumericUtil::smallValToZero<F>(this->m_k, eps);

    return *this;
}


/**
 * "Rounds" all small quaternion's components (with the absolute value below
 * the default 'eps') to 0.
 * 
 * @return reference to itself
 */
template <typename F>
math::QuaternionGeneric<F>& math::QuaternionGeneric<F>::roundSmallElements_()
{
    return this->roundSmallElements_(math::NumericUtil::getEPS<F>());
}


/**
 * "Rounds" all small quaternion's components (with the absolute value below
 * the given 'eps') to 0.
 * 
 * @param eps - threshold to determine whether each component is "rounded" to 0
 * 
 * @return "rounded" quaternion
 */
template <typename F>
math::QuaternionGeneric<F> math::QuaternionGeneric<F>::roundSmallElements(const F& eps) const
{
    math::QuaternionGeneric<F> qret(*this);
    qret.roundSmallElements_(eps);
    return qret;
}


/**
 * "Rounds" all small quaternion's components (with the absolute value below
 * the default 'eps') to 0.
 * 
 * @return "rounded" quaternion
 */
template <typename F>
math::QuaternionGeneric<F> math::QuaternionGeneric<F>::roundSmallElements() const
{
    math::QuaternionGeneric<F> qret(*this);
    qret.roundSmallElements_();
    return qret;
}



/**
 * Addition operator (+) of two quaternions.
 *
 * @param q1 - augend
 * @param q2 - addend
 *
 * @return q1 + q2
 */
template <typename F>
math::QuaternionGeneric<F> math::operator+(const math::QuaternionGeneric<F>& q1, const math::QuaternionGeneric<F>& q2)
{
    /*
     * Addition of quaternions is trivial:
     *   (a1 + b1*i + c1*j + d1*k) + (a2 + b2*i + c2*j + d2*k) =
     *    = ( (a1+a2) + (b1+b2)*i + (c1+c2)*j + (d1+d2)*k )
     */

    return math::QuaternionGeneric<F>(
            q1.m_o + q2.m_o,
            q1.m_i + q2.m_i,
            q1.m_j + q2.m_j,
            q1.m_k + q2.m_k );
}


/**
 * Subtraction operator (-) of two quaternions.
 *
 * @param q1 - minuend
 * @param q2 - subtrahend
 *
 * @return q1 - q2
 */
template <typename F>
math::QuaternionGeneric<F> math::operator-(const math::QuaternionGeneric<F>& q1, const math::QuaternionGeneric<F>& q2)
{
    /*
     * Subtraction of quaternions is trivial:
     *   (a1 + b1*i + c1*j + d1*k) - (a2 + b2*i + c2*j + d2*k) =
     *    = ( (a1-a2) + (b1-b2)*i + (c1-c2)*j + (d1-d2)*k )
     */

    return math::QuaternionGeneric<F>(
    		q1.m_o - q2.m_o,
    		q1.m_i - q2.m_i,
    		q1.m_j - q2.m_j,
    		q1.m_k - q2.m_k );
}


/**
 * Multiplication operator (*) of two quaternions.
 *
 * Note that quaternion multiplication is not commutative (q*p != p*q).
 *
 * @param q1 - multiplicand
 * @param q2 - multiplier
 *
 * @return q1 * q2
 */
template <typename F>
math::QuaternionGeneric<F> math::operator*(const math::QuaternionGeneric<F>& q1, const math::QuaternionGeneric<F>& q2)
{
    /*
     * From the following definitions:
     *   i*i = j*j = k*k = -1,
     *   i*j = k, j*i = -k, j*k = i, k*j = -i, k*i = j and i*k = -j,
     * the following formula can be quickly derived:
     *
     *   (a1 + b1*i + c1*j + d1*k) * (a2 + b2*i + c2*j + d2*k) =
     *    =  (a1*a2 - b1*b2 - c1*c2 - d1*d2)     +
     *    +  (a1*b2 + b1*a2 + c1*d2 - d1*c2) * i +
     *    +  (a1*c2 - b1*d2 + c1*a2 + d1*b2) * j +
     *    +  (a1*d2 + b1*c2 - c1*b2 + d1*a2) * k
     *
     * Note: The following script for GNU Octave or Matlab can be used
     *       for a quick unit test of the function:
     * http://mind.cog.jhu.edu/courses/680/octave/Installers/Octave/Octave.OSX10.6/Applications/MATLAB_R2009b.app/toolbox/aero/aero/quatmultiply.m
     */

	math::QuaternionGeneric<F> retVal;

    // Calculation of the product can be parallelized into 4 mutually independent sections:
    #pragma omp parallel sections if (OMP_QUAT_PARALLELIZE!=0)
    {
        #pragma omp section
        {
            retVal.m_o = q1.m_o * q2.m_o -
                         q1.m_i * q2.m_i -
                         q1.m_j * q2.m_j -
                         q1.m_k * q2.m_k;
        }

        #pragma omp section
        {
           retVal.m_i = q1.m_o * q2.m_i +
                        q1.m_i * q2.m_o +
                        q1.m_j * q2.m_k -
                        q1.m_k * q2.m_j;
        }

        #pragma omp section
        {
            retVal.m_j = q1.m_o * q2.m_j -
                         q1.m_i * q2.m_k +
                         q1.m_j * q2.m_o +
                         q1.m_k * q2.m_i;
        }

        #pragma omp section
        {
            retVal.m_k = q1.m_o * q2.m_k +
                         q1.m_i * q2.m_j -
                         q1.m_j * q2.m_i +
                         q1.m_k * q2.m_o;
        }
    }  // omp parallel sections

    return retVal;
}


/**
 * Addition operator (+) of a quaternion and a scalar
 *
 * @param q - augend (a quaternion)
 * @param sc - addend (a scalar)
 *
 * @return q + sc
 */
template <typename F>
math::QuaternionGeneric<F> math::operator+(const math::QuaternionGeneric<F>& q, const F& sc)
{
    return math::QuaternionGeneric<F>(
            q.m_o + sc,
            q.m_i,
            q.m_j,
            q.m_k );
}


/**
 * Subtraction operator (-) of a quaternion and a scalar
 *
 * @param q - minuend (a quaternion)
 * @param sc - subtrahend (a scalar)
 *
 * @return q - sc
 */
template <typename F>
math::QuaternionGeneric<F> math::operator-(const math::QuaternionGeneric<F>& q, const F& sc)
{
    return QuaternionGeneric<F>(
    		q.m_o - sc,
    		q.m_i,
    		q.m_j,
    		q.m_k );
}


/**
 * Multiplication operator (*) for multiplication of a quaternion and a scalar.
 * It multiplies each quaternion's component by the scalar.
 *
 * @param q - multiplicand (a quaternion)
 * @param sc - multiplier (a scalar)
 *
 * @return q * sc
 */
template <typename F>
math::QuaternionGeneric<F> math::operator*(const math::QuaternionGeneric<F>& q, const F& sc)
{
    /*
     * From the definition of quaternion multiplication (see operator*),
     * one can quickly derive the following simplified formula:
     *   (a+ b*i + c*j + d*k) * s = (a*s + (b*s)*i + (c*s)*j + (d*s)*k))
     */

    return math::QuaternionGeneric<F>(
            q.m_o * sc,
            q.m_i * sc,
            q.m_j * sc,
            q.m_k * sc );
}


/**
 * Addition operator (+) of a scalar and a quaternion.
 * In general ( if T represents a real number) this operation is commutative
 * and does the same as operator+(scalar).
 *
 * @param scalar - augend (a scalar)
 * @param q - addend (a quaternion)
 *
 * @return scalar + q
 */
template <typename F>
math::QuaternionGeneric<F> math::operator+(const F& scalar, const math::QuaternionGeneric<F>& q)
{
    // Addition is commutative
    return (q + scalar);
}


/**
 * Subtraction operator (-) of a scalar and a quaternion.
 *
 * @param scalar - minuend (a scalar)
 * @param q - subtrahend (a quaternion)
 *
 * @return scalar - q
 */
template <typename F>
math::QuaternionGeneric<F> math::operator-(const F& scalar, const math::QuaternionGeneric<F>& q)
{
    // Subtraction is not commutative!
    return math::QuaternionGeneric<F>
            ( scalar - q.m_o, -q.m_i, -q.m_j, -q.m_k );
}


/**
 * Multiplication operator (*) of a scalar and a quaternion.
 * In general ( if T represents a real number) this operation is commutative
 * and does the same as operator*(scalar).
 *
 * @param scalar - multiplicand (a scalar)
 * @param q - multiplier (a quaternion)
 *
 * @return scalar * q
 */
template <typename F>
math::QuaternionGeneric<F> math::operator*(const F& scalar, const math::QuaternionGeneric<F>& q)
{
    /*
     * In general, multiplication of a scalar and quaternion is commutative.
     * If this is not a case, implement a specialization.
     */
    return (q * scalar);
}


/**
 * Division operator (/) that divides a quaternion by a scalar.
 * 
 * @param q - dividend (quaternion)
 * @param sc - divisor (scalar)
 * 
 * @return q / sc
 * 
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::operator/(
        const math::QuaternionGeneric<F>&q,
        const F& sc
      ) throw(math::QuaternionException)
{
    math::QuaternionGeneric<F> retVal(q);
    retVal /= sc;
    return retVal;
}


/**
 * Division operator (/) that divides a scalar by a quaternion
 * 
 * @param sc - dividend (scalar)
 * @param q - divisor (quaternion)
 * 
 * @return sc / q = sc * rec(q)
 * 
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::operator/(
        const F& sc,
        const math::QuaternionGeneric<F>& q
      ) throw(math::QuaternionException)
{
    if ( true == math::NumericUtil::isZero<F>(q.__sqsum()) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    return q.reciprocal() * sc;
}


/**
 * Right division "operator" (q1 * rec(q2))
 *
 * @param q1 - dividend quaternion
 * @param q2 - divisor quaternion
 *
 * @return q1 * rec(q2)
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::rdiv(
        const math::QuaternionGeneric<F>& q1,
        const math::QuaternionGeneric<F>& q2
      ) throw(math::QuaternionException)
{
    math::QuaternionGeneric<F> retVal(q1);
    retVal.rdiv_(q2);

    return retVal;
}


/**
 * Right division "operator" (q * rec(sc))
 *
 * @param q - dividend (quaternion)
 * @param sc - divisor (scalar)
 *
 * @return q * rec(sc)
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::rdiv(
        const math::QuaternionGeneric<F>& q,
        const F& sc
      ) throw(math::QuaternionException)
{
    return q / sc;
}


/**
 * Right division "operator" (sc * rec(q))
 *
 * @param sc - dividend (scalar)
 * @param q - divisor (quaternion)
 *
 * @return sc * rec(q)
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::rdiv(
        const F& sc,
        const math::QuaternionGeneric<F>& q
      ) throw(math::QuaternionException)
{
    return math::rdiv(math::QuaternionGeneric<F>(sc), q);
}


/**
 * Left division "operator" (rec(q2) * q1)
 *
 * @param q1 - dividend quaternion
 * @param q2 - divisor quaternion
 *
 * @return rec(q2) * q1
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::ldiv(
        const math::QuaternionGeneric<F>& q1,
        const math::QuaternionGeneric<F>& q2
      ) throw(math::QuaternionException)
{
    if ( true == math::NumericUtil::isZero<F>(q2.__sqsum()) )
    {
        throw math::QuaternionException(math::QuaternionException::DIVIDE_BY_ZERO);
    }

    return q2.reciprocal() * q1;
}


/**
 * Left division "operator" (rec(sc) * q)
 *
 * @param q - dividend (quaternion)
 * @param sc - divisor (scalar)
 *
 * @return rec(sc) * q
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::ldiv(
        const math::QuaternionGeneric<F>& q,
        const F& sc
      ) throw(math::QuaternionException)
{
    return q / sc;
}


/**
 * Left division "operator" (rec(q) * sc)
 *
 * @param sc - dividend (quaternion)
 * @param q - divisor (scalar)
 *
 * @return rec(q) * sc
 *
 * @throw QuaternionException if attempting to divide by zero
 */
template <typename F>
math::QuaternionGeneric<F> math::ldiv(
        const F& sc,
        const math::QuaternionGeneric<F>& q
      ) throw(math::QuaternionException)
{
    return math::ldiv(math::QuaternionGeneric<F>(sc), q);
}


/**
 * Unary operator '+', returns a copy of the input argument 'q'.
 * 
 * @note Usage of this operator should be avoided
 * 
 * @param q - quaternion to be copied
 * 
 * @return copy of 'q'
 */
template <typename F>
math::QuaternionGeneric<F> math::operator+(const math::QuaternionGeneric<F>& q)
{
    return q;
}


/**
 * Unary negation operator (-)
 *
 * @param q - quaternion to be negated
 * 
 * @return -q
*/
template <typename F>
math::QuaternionGeneric<F> math::operator-(const math::QuaternionGeneric<F>& q)
{
    /*
     * Negation of a quaternion is trivial: just negate each component.
     */

    return math::QuaternionGeneric<F>(
            -q.m_o,
            -q.m_i,
            -q.m_j,
            -q.m_k );
}


/**
 * A friend function that outputs the quaternion to an output stream
 *
 * @param output - stream to write to
 * @param q - quaternion to be displayed
 *
 * @return reference of output stream (i.e. 'output')
 */
template <typename F>
std::ostream& math::operator<<(std::ostream& output, const math::QuaternionGeneric<F>& q)
{
    // just pass the quaternion to QuaternionGeneric::display()...
    q.display(output);

    // ... and return reference of the stream
    return output;
}



/*
 * Specialization of other classes'/namespaces' templated functions for
 * the class QuternionGeneric.
 *
 * Note: the specialized functions must be implemented within
 *       corresponding namespaces.
 */

namespace std
{

/*
 * "Specialization" of std::abs()
 */
template <typename F>
F abs(const math::QuaternionGeneric<F>& q)
{
    // Absolute value of a quaternion is actually equal to its norm
    return q.norm();
}

}
