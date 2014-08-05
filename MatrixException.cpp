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
 * @file MatrixException.cpp
 *
 * Implementation of the class MatrixException
 *
 * @author Jernej Kovacic
 */

#include <ostream>

#include "MatrixException.h"


/**
 * Constructor.
 *
 * @param error code
 */
math::MatrixException::MatrixException(err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::MatrixException::display(std::ostream& str) const
{
    switch (error)
    {
        case FORBIDDEN:
            str << "Forbidden operation";
            break;
        case OUT_OF_MEMORY:
            str << "Could not allocate enough memory";
            break;
        case TOO_LARGE:
            str << "Too many rows or columns";
            break;
        case INVALID_DIMENSION:
            str << "Invalid dimension of matrix";
            break;
        case OUT_OF_RANGE:
            str << "Attempted to access elements out of valid range";
            break;
        case NOT_ENOUGH_ELEMENTS:
            str << "Not enough elements";
            break;
        case NON_INVERTIBLE_MATRIX:
            str << "Matrix is not invertible";
            break;
        default:
            // Should not occur but handle it anyway.
            // Maybe a code was inserted into err_codes and
            // this function hasn't been updated?
            str << "Strange, unspecified error";
    }  // switch

    str << std::endl;
}
