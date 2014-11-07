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
 * Implementation of the class StatisticsException
 */


#include <ostream>

#include "StatisticsException.hpp"


/**
 * Constructor
 *
 * @param error code
 */
math::StatisticsException::StatisticsException(err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::StatisticsException::what(std::ostream& str) const
{
    switch (error)
    {
        case SAMPLE_NOT_PROCESSED_YET :
            str << "The sample has not been processed yet.";
            break;
        case SAMPLE_EMPTY :
            str << "The sample is empty.";
            break;
        case SAMPLE_TOO_SMALL :
            str << "The sample is too small to perform the operation.";
            break;
        case DF_SUBTRAHEND_TOO_LARGE :
            str << "DF subtrahend exceeds the sample's size.";
            break;
        case UNSUPPORTED_TYPE :
            str << "Operation is not supported for the specified type T.";
            break;
        case OUT_OF_MEMORY :
            str << "Could not allocate enough memory.";
            break;
        case INVALID_PROBABILTY :
            str << "Probabilty is not within the valid range.";
            break;
        case UNSUPPORTED_QUANTILE_METHOD :
            str << "The method to estimate the quantile is not supported yet.";
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