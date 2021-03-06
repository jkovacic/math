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
 * @headername{MatrixException.h}
 *
 * Declaration of the class MatrixException.
 */


#ifndef _MATH_MATRIXEXCEPTION_HPP_
#define _MATH_MATRIXEXCEPTION_HPP_

#include <ostream>

#include "exception/IMathException.hpp"

namespace math
{

/**
 * @brief An exception, typically  thrown by functions of MatrixGeneric and
 * its derived classes.
 */
struct MatrixException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        FORBIDDEN,                  /// Operation is forbidden for this type of a matrix
        OUT_OF_MEMORY,              /// Could not allocate enough memory
        TOO_LARGE,                  /// Matrix is too large
        INVALID_DIMENSION,          /// Dimension is not valid for the operation
        OUT_OF_RANGE,               /// Attempt to access an element out of defined range
        NONSQUARE_MATRIX,           /// Operation only allowed on square matrices 
        NOT_ENOUGH_ELEMENTS,        /// Not enough elements
        NON_INVERTIBLE_MATRIX,      /// Matrix cannot be inverted
    };

    err_codes error;     /// Type of an error
    // Constructor
    MatrixException(const MatrixException::err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif  // _MATH_MATRIXEXCEPTION_HPP_
