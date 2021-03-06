/*
Copyright 2015, Jernej Kovacic

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
 * Implementation of the class SelectionException
 */


#include <ostream>

#include "exception/SelectionException.hpp"


/**
 * Constructor
 *
 * @param error code
 */
math::SelectionException::SelectionException(const math::SelectionException::err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::SelectionException::what(std::ostream& str) const
{
    switch (this->error)
    {
        case math::SelectionException::OUT_OF_MEMORY :
            str << "Could not allocate enough memory.";
            break;
        case math::SelectionException::SAMPLE_EMPTY :
            str << "The sample is empty.";
            break;
        case math::SelectionException::ARG_OUT_OF_RANGE :
            str << "Input argument is out of range.";
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
