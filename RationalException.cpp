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
@file RationalException.cpp

Implementation of the class RationalException

@author Jernej Kovacic
*/


#include "RationalException.h"

using math::RationalException;

/**
 * Constructor
 *
 * @param error code
 */
RationalException::RationalException(err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void RationalException::display(std::ostream& str) const
{
    switch (error)
    {
        case ZERO_DENOMINATOR:
            str << "Zero denominator is forbidden.";
            break;
        case UNINVERTIBLE:
            str << "Uninvertible fraction.";
            break;
        case DIVIDE_BY_ZERO:
            str << "Attempt of division by zero.";
            break;
        case OVERFLOW:
            str << "Operation caused an integer overflow.";
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
