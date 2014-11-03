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
 * @file LinearEquationSolverException.cpp
 *
 * Implementation of the class LinearEquationSolverException
 *
 * @author Jernej Kovacic
 */


#include <ostream>

#include "LinearEquationSolverException.h"


/**
 * Constructor.
 *
 * @param error code
 */
math::LinearEquationSolverException::LinearEquationSolverException(err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::LinearEquationSolverException::what(std::ostream& str) const
{
    switch (error)
    {
        case OUT_OF_MEMORY:
            str << "Could not allocate enough memory";
            break;
        case INVALID_DIMENSION:
            str << "Invalid dimensions of input matrices";
            break;
        case NO_UNIQUE_SOLUTION:
            str << "Unique solution of the system of linear equations does not exist";
            break;
        default:
            // Should not occur but handle it anyway.
            // Maybe a code was inserted into err_codes and
            // this function hasn't been updated?
            str << "Strange, unspecified error";
    }  // switch

    str << std::endl;
}
