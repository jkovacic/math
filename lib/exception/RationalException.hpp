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
 * @file
 * @author Jernej Kovacic
 *
 * An internal header file, it should not be included directly.
 * @headername{RationalException.h}
 *
 * Declaration of the class RationalException.
 */

#ifndef _MATH_RATIONALEXCEPTION_HPP_
#define	_MATH_RATIONALEXCEPTION_HPP_

#include <ostream>

#include "exception/IMathException.hpp"

namespace math
{

/**
 * @brief An exception, typically  thrown by functions of Rational.
 */
struct RationalException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        ZERO_DENOMINATOR,   /// attempted to set a zero denominator
        UNINVERTIBLE,       /// the fraction cannot be inverted
        DIVIDE_BY_ZERO,     /// attempt of division by zero
        INT_OVERFLOW,       /// integer overflow
        INVALID_INPUT,      /// invalid input argument
        INPUT_OUT_OF_RANGE, /// string input out of integer's boundaries
        OUT_OF_MEMORY,      /// (re)allocation of string's memory failed
        UNSIGNED,           /// operation impossible for unsigned type
    };

    err_codes error;    /// type of an error
    // Constructor
    RationalException(const RationalException::err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif	// _MATH_RATIONALEXCEPTION_HPP_
