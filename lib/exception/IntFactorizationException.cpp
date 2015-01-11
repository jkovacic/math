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
 * Implementation of the class IntFactorizationException
 */


#include <ostream>

#include "exception/IntFactorizationException.hpp"


/**
 * Constructor
 *
 * @param error code
 */
math::IntFactorizationException::IntFactorizationException(err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::IntFactorizationException::what(std::ostream& str) const
{
    switch (error)
    {
        case OUT_OF_RANGE:
            str << "Result of the operation is out of the integer range.";
            break;
        case INVALID_INPUT:
            str << "Operation not defined for this input argument.";
            break;
        case NEGATIVE_ARG :
            str << "Negative input arguments are not allowed.";
            break;
        case OUT_OF_MEMORY:
            str << "Allocation of memory failed.";
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
