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
 * An internal header file, it should not be included directly.
 * @headername{Quaternion.h}
 *
 * Declaration of the class QuaternionGeneric, representing quaternions.
 */

#ifndef _MATH_QUATERNIONGENERIC_HPP_
#define _MATH_QUATERNIONGENERIC_HPP_


#include "exception/QuaternionException.hpp"

namespace math
{

// Forward declaration of the class is necessary...
template <typename F> class QuaternionGeneric;

// to declare its friend functions:
template <typename F>
QuaternionGeneric<F> operator+(const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);

template <typename F>
QuaternionGeneric<F> operator-(const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);

template <typename F>
QuaternionGeneric<F> operator*(const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);

template <typename F>
QuaternionGeneric<F> operator+(const QuaternionGeneric<F>& q, const F& sc);

template <typename F>
QuaternionGeneric<F> operator-(const QuaternionGeneric<F>& q, const F& sc);

template <typename F>
QuaternionGeneric<F> operator*(const QuaternionGeneric<F>& q, const F& sc);

template <typename F>
QuaternionGeneric<F> operator+(const F& scalar, const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> operator-(const F& scalar, const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> operator*(const F& scalar, const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> operator+(const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> operator-(const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> operator/(const QuaternionGeneric<F>&q, const F& sc);

template <typename F>
QuaternionGeneric<F> operator/(const F& sc, const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> rdiv(const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);

template <typename F>
QuaternionGeneric<F> rdiv(const QuaternionGeneric<F>& q, const F& sc);

template <typename F>
QuaternionGeneric<F> rdiv(const F& sc, const QuaternionGeneric<F>& q);

template <typename F>
QuaternionGeneric<F> ldiv(const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);

template <typename F>
QuaternionGeneric<F> ldiv(const QuaternionGeneric<F>& q, const F& sc);

template <typename F>
QuaternionGeneric<F> ldiv(const F& sc, const QuaternionGeneric<F>& q);

// and its friend << operator:
template <typename F>
std::ostream& operator<<(std::ostream& output, const QuaternionGeneric<F>& q);


/**
 * A class representing quaternions and their basic operations.
 * 
 * Quaternions are a convenient tool to calculate spatial rotations and 
 * as such very handy at computer graphics, robotics, etc.
 * 
 * For more details about quaternions, see:
 * http://en.wikipedia.org/wiki/Quaternion
 */
template <typename F>
class QuaternionGeneric
{
private:
    F m_o;     /// element '1'
    F m_i;     /// element 'i'
    F m_j;     /// element 'j'
    F m_k;     /// element 'k'


    /*
     * A utility function that calculates a sum of components' squares.
     * The function is called  by other public methods, such as norm, reciprocal, unit, etc.
     *
     * As the function is short, it is declared as inline to slightly reduce overhead
     *
     * @return sum of all components' squares
     */
    inline F __sqsum() const
    {
        return (
                this->m_o * this->m_o +
                this->m_i * this->m_i +
                this->m_j * this->m_j +
                this->m_k * this->m_k );
    }


public:
    // Constructor
    QuaternionGeneric(
            const F& o = static_cast<F>(0),
            const F& i = static_cast<F>(0),
            const F& j = static_cast<F>(0),
            const F& k = static_cast<F>(0) );

    // Copy constructor
    QuaternionGeneric(const QuaternionGeneric& q);

    // Getters:
    F getOne() const;
    F getI() const;
    F getJ() const;
    F getK() const;

    // Setters:
    QuaternionGeneric<F>& set(
            const F& o = static_cast<F>(0),
            const F& i = static_cast<F>(0),
            const F& j = static_cast<F>(0),
            const F& k = static_cast<F>(0) );

    QuaternionGeneric<F>& setOne(const F& o = static_cast<F>(0));
    QuaternionGeneric<F>& setI(const F& i = static_cast<F>(0));
    QuaternionGeneric<F>& setJ(const F& j = static_cast<F>(0));
    QuaternionGeneric<F>& setK(const F& k = static_cast<F>(0) );

    // Display the quaternion:
    void display(std::ostream& str = std::cout) const;

    // Quaternion arithmetics operators:
    QuaternionGeneric<F>& operator=(const QuaternionGeneric<F>& q);
    QuaternionGeneric<F>& operator=(const F& sc);
    QuaternionGeneric<F>& operator+=(const QuaternionGeneric<F>& q);
    QuaternionGeneric<F>& operator-=(const QuaternionGeneric<F>& q);
    QuaternionGeneric<F>& operator*=(const QuaternionGeneric<F>& q);
    QuaternionGeneric<F>& operator*=(const F& scalar);
    QuaternionGeneric<F>& operator+=(const F& scalar);
    QuaternionGeneric<F>& operator-=(const F& scalar);
    QuaternionGeneric<F>& operator/=(const F& scalar);
    QuaternionGeneric<F>& rdiv_(const QuaternionGeneric<F>& q);
    QuaternionGeneric<F>& rdiv_(const F& sc);
    QuaternionGeneric<F>& ldiv_(const QuaternionGeneric<F>& q);
    QuaternionGeneric<F>& ldiv_(const F& sc);

    // Conjugate the quaternion:
    QuaternionGeneric<F>& conj_();
    QuaternionGeneric<F> conj() const;

    // Norm of the quaternion (||q||):
    F norm() const;

    // Unit quaternion (q/||q||):
    QuaternionGeneric<F>& unit_();
    QuaternionGeneric<F> unit() const;

    // Reciprocal quaternion (q^(-1)):
    QuaternionGeneric<F>& reciprocal_();
    QuaternionGeneric<F> reciprocal() const;

    // Round very small values of elements to 0
    QuaternionGeneric<F>& roundSmallElements_(const F& eps);
    QuaternionGeneric<F>& roundSmallElements_();
    QuaternionGeneric<F> roundSmallElements(const F& eps) const;
    QuaternionGeneric<F> roundSmallElements() const;
    
    // operators in separate functions declared as friend functions
    friend QuaternionGeneric<F> (math::operator+ <>) (const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);
    friend QuaternionGeneric<F> (math::operator- <>) (const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);
    friend QuaternionGeneric<F> (math::operator* <>) (const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);
    friend QuaternionGeneric<F> (math::operator+ <>) (const QuaternionGeneric<F>& q, const F& sc);
    friend QuaternionGeneric<F> (math::operator- <>) (const QuaternionGeneric<F>& q, const F& sc);
    friend QuaternionGeneric<F> (math::operator* <>) (const QuaternionGeneric<F>& q, const F& sc);
    friend QuaternionGeneric<F> (math::operator* <>) (const F& scalar, const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::operator+ <>) (const F& scalar, const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::operator- <>) (const F& scalar, const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::operator+ <>) (const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::operator- <>) (const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::operator/ <>) (const QuaternionGeneric<F>& q, const F& sc);
    friend QuaternionGeneric<F> (math::operator/ <>) (const F& sc, const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::rdiv <>) (const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);
    friend QuaternionGeneric<F> (math::rdiv <>) (const QuaternionGeneric<F>& q, const F& sc);
    friend QuaternionGeneric<F> (math::rdiv <>) (const F& sc, const QuaternionGeneric<F>& q);
    friend QuaternionGeneric<F> (math::ldiv <>) (const QuaternionGeneric<F>& q1, const QuaternionGeneric<F>& q2);
    friend QuaternionGeneric<F> (math::ldiv <>) (const QuaternionGeneric<F>& q, const F& sc);
    friend QuaternionGeneric<F> (math::ldiv <>) (const F& sc, const QuaternionGeneric<F>& q);

    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& (math::operator<< <>) (std::ostream& output, const QuaternionGeneric<F>& q);

};   // class QuaternionGeneric



// Quaternions with elements of types float, double and long double
// make most sense so the following three types are predefined:
typedef QuaternionGeneric<float>       FQuaternion;
typedef QuaternionGeneric<double>      Quaternion;
typedef QuaternionGeneric<long double> LDQuaternion;

}  // namespace math

// DEFINITION
#include "quaternion/QuaternionGeneric.cpp"

#endif  // _MATH_QUATERNIONGENERIC_HPP_
