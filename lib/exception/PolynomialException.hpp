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
 * @headername{PolynomialException.h}
 *
 * Declaration of the class PolynomialException.
 */

#ifndef _MATH_POLYNOMIALEXCEPTION_HPP_
#define _MATH_POLYNOMIALEXCEPTION_HPP_

#include <ostream>

#include "IMathException.hpp"

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
        TOO_LARGE,                  /// Too many polynomial's terms
        DIVIDE_BY_ZERO              /// Attempt of division by zero
    };

    err_codes error;     /// Type of an error
    // Constructor
    PolynomialException(err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif // _MATH_POLYNOMIALEXCEPTION_H_
