/*
Copyright 2014, Jernej Kovacic

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
 * @headername{RootFindException.h}
 *
 * Declaration of the class RootFindException.
 */


#ifndef _MATH_ROOT_FIND_EXCEPTION_HPP_
#define _MATH_ROOT_FIND_EXCEPTION_HPP_

#include "exception/IMathException.hpp"

#include <ostream>

namespace math
{

/**
 * @brief An exception, typically thrown by root finding algorithms
 */
struct RootFindException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes
    {
        UNDEFINED,            /// found a point where the function is not defined
        INVALID_ARGS,         /// invalid input arguments
        TOO_MANY_ROOTS,       /// too many roots in the search interval
        NO_CONVERGENCE,       /// algorithm could not converge
        ZERO_SLOPE,           /// encountered a point with a zero slope
    };

    err_codes error;    /// type of an error
    // Constructor
    RootFindException(const RootFindException::err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif	// _MATH_ROOT_FIND_EXCEPTION_HPP_
