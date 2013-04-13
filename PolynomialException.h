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
@file PolynomialException.h

Declaration of the class PolynomialException

@author Jernej Kovacic
*/

#ifndef _MATH_POLYNOMIALEXCEPTION_H_
#define _MATH_POLYNOMIALEXCEPTION_H_

#include <iostream>

#include "IMathException.h"

namespace math
{

/**
 * @brief An exception, typically  thrown by functions of PolynomialGeneric
 */
struct PolynomialException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        OUT_OF_MEMORY,              /// Could not allocate enough memory
        OUT_OF_RANGE,               /// Attempt to access an element out of range
        INVALID_ARGUMENT,           /// Invalid argument passed to constructor
        TOO_LARGE                   /// Too many polynomial's terms
    };

    err_codes error;     /// Type of an error
    // Constructor
    PolynomialException(err_codes err);
    // Output a short description of the error
    void display(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif // _MATH_POLYNOMIALEXCEPTION_H_
