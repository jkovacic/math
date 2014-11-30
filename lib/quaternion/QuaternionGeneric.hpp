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
 * @headername{QuaternionGeneric.h}
 *
 * Declaration of the class QuaternionGeneric, representing quaternions.
 */

#ifndef _MATH_QUATERNIONGENERIC_HPP_
#define _MATH_QUATERNIONGENERIC_HPP_


#include "util/NumericUtil.hpp"
#include "exception/QuaternionException.hpp"

namespace math
{

// Advance declaration of the class is necessary...
template<class T> class QuaternionGeneric;

// to declare its friend functions:
template<class T>
QuaternionGeneric<T> operator+(const QuaternionGeneric<T>& q1, const QuaternionGeneric<T>& q2);

template<class T>
QuaternionGeneric<T> operator-(const QuaternionGeneric<T>& q1, const QuaternionGeneric<T>& q2);

template<class T>
QuaternionGeneric<T> operator*(const QuaternionGeneric<T>& q1, const QuaternionGeneric<T>& q2);

template<class T>
QuaternionGeneric<T> operator+(const QuaternionGeneric<T>& q, const T& sc);

template<class T>
QuaternionGeneric<T> operator-(const QuaternionGeneric<T>& q, const T& sc);

template<class T>
QuaternionGeneric<T> operator*(const QuaternionGeneric<T>& q, const T& sc);

template<class T>
QuaternionGeneric<T> operator+(const T& scalar, const QuaternionGeneric<T>& q);

template<class T>
QuaternionGeneric<T> operator-(const T& scalar, const QuaternionGeneric<T>& q);

template<class T>
QuaternionGeneric<T> operator*(const T& scalar, const QuaternionGeneric<T>& q);

// and its friend << operator:
template<class T>
std::ostream& operator<<(std::ostream& output, const QuaternionGeneric<T>& q);


/**
 * A class representing quaternions and their basic operations.
 * 
 * Quaternions are a convenient tool to calculate spatial rotations and 
 * as such very handy at computer graphics, robotics, etc.
 * 
 * For more details about quaternions, see:
 * http://en.wikipedia.org/wiki/Quaternion
 */
template <class T>
class QuaternionGeneric
{
private:
    T quat_o;    /// element '1'
    T quat_i;    /// element 'i'
    T quat_j;    /// element 'j'
    T quat_k;    /// element 'k'


    /*
     * A utility function that calculates a sum of components' squares.
     * The function is called  by other public methods, such as norm, reciprocal, unit, etc.
     *
     * As the function is short, it is declared as inline to slightly reduce overhead
     *
     * @return sum of all components' squares
     */
    inline T __sqsum() const
    {
        return (
                this->quat_o * this->quat_o +
                this->quat_i * this->quat_i +
                this->quat_j * this->quat_j +
                this->quat_k * this->quat_k );
    }


public:
    // Constructor
    QuaternionGeneric(const T& o = static_cast<T>(0),
                      const T& i = static_cast<T>(0),
                      const T& j = static_cast<T>(0),
                      const T& k = static_cast<T>(0) );

    // Copy constructor
    QuaternionGeneric(const QuaternionGeneric& q);

    // Getters:
    T getOne() const;
    T getI() const;
    T getJ() const;
    T getK() const;

    // Setters:
    QuaternionGeneric<T>& set(const T& o = static_cast<T>(0),
             const T& i = static_cast<T>(0),
             const T& j = static_cast<T>(0),
             const T& k = static_cast<T>(0) );

    QuaternionGeneric<T>& setOne(const T& o = static_cast<T>(0));
    QuaternionGeneric<T>& setI(const T& i = static_cast<T>(0));
    QuaternionGeneric<T>& setJ(const T& j = static_cast<T>(0));
    QuaternionGeneric<T>& setK(const T& k = static_cast<T>(0) );

    // Display the quaternion:
    void display(std::ostream& str = std::cout) const;

    // Quaternion arithmetics operators:
    QuaternionGeneric<T>& operator=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T>& operator+=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T>& operator-=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T>& operator*=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T>& operator*=(const T& scalar);
    QuaternionGeneric<T>& operator+=(const T& scalar);
    QuaternionGeneric<T>& operator-=(const T& scalar);
    QuaternionGeneric<T> operator-() const;

    // Conjugate the quaternion:
    QuaternionGeneric<T> conj() const;
    QuaternionGeneric<T>& conjugated();

    // Norm of the quaternion (||q||):
    T norm() const throw (QuaternionException);

    // Unit quaternion (q/||q||):
    QuaternionGeneric<T> unit() const throw (QuaternionException);

    // Reciprocal quaternion (q^(-1)):
    QuaternionGeneric<T> reciprocal() const throw (QuaternionException);


    // operators in separate functions declared as friend functions
    friend QuaternionGeneric<T> (math::operator+ <>) (const QuaternionGeneric<T>& q1, const QuaternionGeneric<T>& q2);
    friend QuaternionGeneric<T> (math::operator- <>) (const QuaternionGeneric<T>& q1, const QuaternionGeneric<T>& q2);
    friend QuaternionGeneric<T> (math::operator* <>) (const QuaternionGeneric<T>& q1, const QuaternionGeneric<T>& q2);
    friend QuaternionGeneric<T> (math::operator+ <>) (const QuaternionGeneric<T>& q, const T& sc);
    friend QuaternionGeneric<T> (math::operator- <>) (const QuaternionGeneric<T>& q, const T& sc);
    friend QuaternionGeneric<T> (math::operator* <>) (const QuaternionGeneric<T>& q, const T& sc);
    friend QuaternionGeneric<T> (math::operator* <>) (const T& scalar, const QuaternionGeneric<T>& q);
    friend QuaternionGeneric<T> (math::operator+ <>) (const T& scalar, const QuaternionGeneric<T>& q);
    friend QuaternionGeneric<T> (math::operator- <>) (const T& scalar, const QuaternionGeneric<T>& q);

    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& (math::operator<< <>) (std::ostream& output, const QuaternionGeneric<T>& q);

};   // class QuaternionGeneric


// Declaration of specialized methods inside the name space declaration
// is essential if implemented elsewhere:
template<> float QuaternionGeneric<float>::norm() const throw (QuaternionException);
template<> double QuaternionGeneric<double>::norm() const throw (QuaternionException);
template<> long double QuaternionGeneric<long double>::norm() const throw (QuaternionException);
template<> QuaternionGeneric<float> QuaternionGeneric<float>::unit() const throw (QuaternionException);
template<> QuaternionGeneric<double> QuaternionGeneric<double>::unit() const throw (QuaternionException);
template<> QuaternionGeneric<long double> QuaternionGeneric<long double>::unit() const throw (QuaternionException);

// Quaternions with elements of types float, double and long double
// make most sense so the following three types are predefined:
typedef QuaternionGeneric<float>       FQuaternion;
typedef QuaternionGeneric<double>      Quaternion;
typedef QuaternionGeneric<long double> LDQuaternion;


// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "quaternion/QuaternionGeneric.cpp"

#endif // _MATH_QUATERNIONGENERIC_HPP_
