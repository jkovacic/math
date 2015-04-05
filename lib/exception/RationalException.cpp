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
 * Implementation of the class RationalException
 */


#include <ostream>

#include "exception/RationalException.hpp"


/**
 * Constructor
 *
 * @param error code
 */
math::RationalException::RationalException(const math::RationalException::err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::RationalException::what(std::ostream& str) const
{
    switch (this->error)
    {
        case math::RationalException::ZERO_DENOMINATOR :
            str << "Zero denominator is forbidden.";
            break;
        case math::RationalException::UNINVERTIBLE :
            str << "Uninvertible fraction.";
            break;
        case math::RationalException::DIVIDE_BY_ZERO :
            str << "Attempt of division by zero.";
            break;
        case math::RationalException::INT_OVERFLOW :
            str << "Operation caused an integer overflow.";
            break;
        case math::RationalException::INVALID_INPUT :
            str << "Invalid input argument.";
            break;
        case math::RationalException::INPUT_OUT_OF_RANGE :
            str << "Input argument out of range.";
            break;
        case math::RationalException::OUT_OF_MEMORY :
            str << "Allocation of memory failed.";
            break;
        case math::RationalException::UNSIGNED :
            str << "Operation is not possible for unsigned types.";
            break;
        default:
            // Should not occur but handle it anyway.
            // Maybe a code was inserted into err_codes and
            // this function hasn't been updated?
            str << "Strange, unspecified error.";
    }  // switch

    // Output a newline character at the end.
    str << std::endl;
}
