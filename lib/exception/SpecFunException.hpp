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
 * @headername{SpecFunException.h}
 *
 * Declaration of the class SpecFunException.
 */


#ifndef _MATH_SPECFUNEXCEPTION_HPP_
#define	_MATH_SPECFUNEXCEPTION_HPP_

#include "exception/IMathException.hpp"

#include <ostream>

namespace math
{

/**
 * @brief An exception, typically thrown by functions that
 *        evaluate special functions
 */
struct SpecFunException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes
    {
        UNDEFINED,           /// function not defined for the given input argument
        NO_CONVERGENCE,      /// evaluation algorithm did not converge
    };

    err_codes error;    /// type of an error
    // Constructor
    SpecFunException(err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif	// _MATH_SPECFUNEXCEPTION_HPP_
