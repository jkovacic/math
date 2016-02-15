/*
Copyright 2016, Jernej Kovacic

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
 * @headername{IntFunctionException.h}
 *
 * Declaration of the class IntFunctionException
 */

#ifndef _MATH_INTFUNCTIONEXCEPTION_HPP_
#define _MATH_INTFUNCTIONEXCEPTION_HPP_


#include "exception/IMathException.hpp"

#include <ostream>

namespace math
{

/**
 * @brief An exception, typically thrown by functions of IntFunction
 */
struct IntFunctionException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        NEGATIVE_ARG,       /// negative input argument is not allowed
        NONPOSITIVE_ARG,    /// only strictly positive arguments are allowed
        INVALID_TYPE,       /// invalid type of an argument
        INTERNAL_ERROR      /// internal error, probably unusual implementation of the type 'I'
    };

    err_codes error;    /// type of an error
    // Constructor
    IntFunctionException(const IntFunctionException::err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif	// _MATH_INTFUNCTIONEXCEPTION_HPP_
