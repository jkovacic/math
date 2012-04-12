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
@file NumericUtil.cpp

Implementation  of the class NumericUtil, a collection of some useful
numerical utilities. This is a templated class and must not be compiled.
Instead it must be included after the class declaration in the .h file

@author Jernej Kovacic
*/

// Delibarately there is no #include "NumericUtil.h"
#include "Rational.h"

// Note that the optimal EPS depends on application requirements
/*
 * Definition of EPS for type float
 */
template<>
float NumericUtil<float>::EPS = 1e-9f;

/*
 * Double is a more accurate type...
 */
template<>
double NumericUtil<double>::EPS = 1e-16;

/*
 * In int and other types, EPS doesn't make sense, so set it to 0
 */
template<class T>
T NumericUtil<T>::EPS = (T) 0;

/*
 * This "general" implementation is valid for float and double.
 * In case of these two types, the equality operator (==) is useless.
 * In numerical mathematics, two numbers are considered "equal", when
 * absolute value of their difference does not exceed a reasonably set EPS.
 */
template<class T>
bool NumericUtil<T>::isZero(T value)
{
    bool retVal = false;
    // quick definition of an absolute value
    T absValue = ( value>=0 ? value : -value );

    retVal = (absValue < EPS ? true : false );

    return retVal;
}

/*
 * And a "specialized" implementation for integers.
 * In this case no comparing to EPS is necessary.
 */
template<>
bool NumericUtil<int>::isZero(int value)
{
    // Now the operator == does make sense....
    bool retVal = ( 0==value ? true : false );

    return retVal;
}

/*
 * Similar for Rational
 */
template<>
bool NumericUtil<Rational>::isZero(Rational value)
{
    // Rational already contains its own isZero()...
    return value.isZero();
}
