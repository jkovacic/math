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
 * Implementation of the class RootFindException
 */


#include <ostream>

#include "exception/RootFindException.hpp"


/**
 * Constructor
 *
 * @param error code
 */
math::RootFindException::RootFindException(const math::RootFindException::err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::RootFindException::what(std::ostream& str) const
{
    switch (error)
    {
        case math::RootFindException::UNDEFINED :
            str << "Encountered a point where the function is not defined";
            break;
        case math::RootFindException::INVALID_ARGS :
            str << "Invalid input arguments";
            break;
        case math::RootFindException::TOO_MANY_ROOTS :
            str << "There are too many roots in the specified search interval";
            break;
        case math::RootFindException::NO_CONVERGENCE :
            str << "Algorithm did not converge";
            break;
        case math::RootFindException::ZERO_SLOPE :
            str << "Encountered a point with the slope equal to zero";
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
