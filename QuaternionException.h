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
 * @file QuaternionException.h
 *
 * Declaration of the class QuaternionException
 *
 * @author Jernej Kovacic
 */

#ifndef _MATH_QUATERNIONEXCEPTION_H_
#define _MATH_QUATERNIONEXCEPTION_H_

#include <ostream>
#include "IMathException.h"

namespace math
{

/**
 * An exception, typically thrown by functions of QuaternionGeneric
 */
struct QuaternionException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        UNSUPPORTED_TYPE,           /// Operation not supported for this type
        DIVIDE_BY_ZERO              /// Attempt of division by zero
    };

    err_codes error;     /// Type of an error
    // Constructor
    QuaternionException(err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif // _MATH_QUATERNIONEXCEPTION_H_
