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
 * @headername{CombinatoricsException.h}
 *
 * Declaration of the class CombinatoricsException
 */

#ifndef _MATH_COMBINATORICSEXCEPTION_HPP_
#define _MATH_COMBINATORICSEXCEPTION_HPP_

#include "exception/IMathException.hpp"

#include <ostream>

namespace math
{

/**
 * @brief An exception, typically thrown by functions of combinatorics
 *        related classes: IntCombinatorics, CombinationGeneric, PermutationGeneric.
 */
struct CombinatoricsException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes
    {
        OUT_OF_RANGE,       /// result is out of integer range
        INVALID_INPUT,      /// operation does not make any sense with this input argument
        OUT_OF_MEMORY       /// allocation of memory failed
    };

    err_codes error;    /// type of an error
    // Constructor
    CombinatoricsException(const CombinatoricsException::err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif  // _MATH_COMBINATORICSEXCEPTION_HPP_
