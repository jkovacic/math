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
 * @headername{CalculusException.h}
 *
 * Declaration of the class CalculusException.
 */


#ifndef _MATH_CALCULUSEXCEPTION_HPP_
#define _MATH_CALCULUSEXCEPTION_HPP_


#include <ostream>

#include "exception/IMathException.hpp"


namespace math
{

/**
 * @brief An exception, typically  thrown by functions of IntegGeneric,
 * DerivGeneric and their related classes.
 */
struct CalculusException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        UNDEFINED,              /// function undefined at given 'x'
        NOT_ENOUGH_STEPS,       /// not enough integration steps
        INVALID_STEP,           /// step size negative or too small
        UNSUPPORTED_ALGORITHM,  /// unsupported algorithm
    };

    err_codes error;     /// Type of an error
    // Constructor
    CalculusException(err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

}

#endif  // _MATH_CALCULUSEXCEPTION_HPP_
