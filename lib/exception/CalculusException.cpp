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
 * Implementation of the class CalculusException
 */


#include <ostream>

#include "exception/CalculusException.hpp"


/**
 * Constructor.
 *
 * @param error code
 */
math::CalculusException::CalculusException(err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::CalculusException::what(std::ostream& str) const
{
    switch (error)
    {
        case UNDEFINED :
            str << "Function undefined at the given input";
            break;
        case NOT_ENOUGH_STEPS :
            str << "Desired number of integration steps is too small";
            break;
        case INVALID_STEP :
            str << "Step size is negative or too small";
            break;
        case UNSUPPORTED_ALGORITHM :
            str << "Unsupported algorithm";
            break;
        case INVALID_BREAKPOINT :
            str << "Invalid sign of a breakpoint";
            break;
        default:
            // Should not occur but handle it anyway.
            // Maybe a code was inserted into err_codes and
            // this function hasn't been updated?
            str << "Strange, unspecified error";
    }  // switch

    str << std::endl;
}
