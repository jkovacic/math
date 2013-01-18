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
@file QuaternionGeneric.h

Declaration of the class QuaternionGeneric, representing quaternions.

@author Jernej Kovacic
*/

#ifndef _MATH_QUATERNIONGENERIC_H_
#define _MATH_QUATERNIONGENERIC_H_


#include "NumericUtil.h"
#include "QuaternionException.h"

namespace math
{

// Advance declaration of the class is necessary...
template<class T> class QuaternionGeneric;
// to declare the class's friend function:
template<class T>
QuaternionGeneric<T> operator* (const T& scalar, const QuaternionGeneric<T>& q);
// and its friend << operator:
template<class T>
std::ostream& operator<<(std::ostream& output, const QuaternionGeneric<T>& q);


template <class T>
class QuaternionGeneric
{
private:
    T quat_o;    /// element '1'
    T quat_i;    /// element 'i'
    T quat_j;    /// element 'j'
    T quat_k;    /// element 'k'

    // A utility function to calculate a sum of all elements' squares,
    // only called by some methods of the class
    inline T sqsum() const;

public:
    // Constructor
    QuaternionGeneric(const T& o=math::NumericUtil<T>::ZERO,
                      const T& i=math::NumericUtil<T>::ZERO,
                      const T& j=math::NumericUtil<T>::ZERO,
                      const T& k=math::NumericUtil<T>::ZERO);

    // Copy constructor
    QuaternionGeneric(const QuaternionGeneric& q);

    // Getters:
    T getOne() const;
    T getI() const;
    T getJ() const;
    T getK() const;

    // Setters:
    QuaternionGeneric<T>& set(const T& o=math::NumericUtil<T>::ZERO,
             const T& i=math::NumericUtil<T>::ZERO,
             const T& j=math::NumericUtil<T>::ZERO,
             const T& k=math::NumericUtil<T>::ZERO );

    QuaternionGeneric<T>& setOne(const T& o=math::NumericUtil<T>::ZERO);
    QuaternionGeneric<T>& setI(const T& i=math::NumericUtil<T>::ZERO);
    QuaternionGeneric<T>& setJ(const T& j=math::NumericUtil<T>::ZERO);
    QuaternionGeneric<T>& setK(const T& k=math::NumericUtil<T>::ZERO);

    // Display the quaternion:
    void display(std::ostream& str = std::cout) const;

    // Quaternion aruthmetics operators:
    QuaternionGeneric<T>& operator=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T> operator+(const QuaternionGeneric<T>& q) const;
    QuaternionGeneric<T>& operator+=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T> operator-(const QuaternionGeneric<T>& q) const;
    QuaternionGeneric<T>& operator-=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T> operator*(const QuaternionGeneric<T>& q) const;
    QuaternionGeneric<T>& operator*=(const QuaternionGeneric<T>& q);
    QuaternionGeneric<T> operator*(const T& scalar) const;
    QuaternionGeneric<T> operator-() const;
    // A friend function that multiplies a scalar and a quaternion
    friend QuaternionGeneric<T> (math::operator* <>) (const T& scalar, const QuaternionGeneric<T>& q);

    // a friend function to overload the operator << (used by std::cout and std::cerr)
    friend std::ostream& (math::operator<< <>) (std::ostream& output, const QuaternionGeneric<T>& q);

    // Conjugate the quaternion:
    QuaternionGeneric<T> conj() const;
    QuaternionGeneric<T>& conjugated();

    // Norm of the quaternion (||q||):
    T norm() const throw (QuaternionException);

    // Unit quaternion (q/||q||):
    QuaternionGeneric<T> unit() const throw (QuaternionException);

    // Reciprocal quaternion (q^(-1)):
    QuaternionGeneric<T> reciprocal() const throw (QuaternionException);
};   // class QuaternionGeneric


// Declaration of specialized methods inside the name space declaration
// is essential if implemented elsewhere:
template<> float QuaternionGeneric<float>::norm() const throw (QuaternionException);
template<> double QuaternionGeneric<double>::norm() const throw (QuaternionException);
template<> QuaternionGeneric<float> QuaternionGeneric<float>::unit() const throw (QuaternionException);
template<> QuaternionGeneric<double> QuaternionGeneric<double>::unit() const throw (QuaternionException);

// Quaternions with elements of types float and double make most sense
// so the following two types are predefined:
typedef QuaternionGeneric<float> FQuaternion;
typedef QuaternionGeneric<double> Quaternion;

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "QuaternionGeneric.cpp"

#endif // _MATH_QUATERNIONGENERIC_H_
