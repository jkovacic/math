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

#include "exception/StatisticsException.hpp"


/**
 * Constructor
 *
 * @param error code
 */
math::StatisticsException::StatisticsException(const math::StatisticsException::err_codes err)
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
    switch (this->error)
    {
        case math::StatisticsException::SAMPLE_EMPTY :
            str << "The sample is empty.";
            break;
        case math::StatisticsException::SAMPLE_TOO_SMALL :
            str << "The sample is too small to perform the operation.";
            break;
        case math::StatisticsException::DF_SUBTRAHEND_TOO_LARGE :
            str << "DF subtrahend exceeds the sample's size.";
            break;
        case math::StatisticsException::UNEQUAL_SAMPLE_SIZES :
            str << "Sample sizes are not equal.";
            break;
        case math::StatisticsException::OUT_OF_MEMORY :
            str << "Could not allocate enough memory.";
            break;
        case math::StatisticsException::INVALID_PROBABILTY :
            str << "Probability is not within the valid range.";
            break;
        case math::StatisticsException::UNSUPPORTED_QUANTILE_METHOD :
            str << "The method to estimate the quantile is not supported yet.";
            break;
        case math::StatisticsException::INVALID_STDEV :
            str << "Standard deviation is negative or zero.";
            break;
        case math::StatisticsException::INVALID_DF :
            str << "Zero or negative number of degrees of freedom.";
            break;
        case math::StatisticsException::INVALID_ARG :
            str << "Invalid value of the input argument.";
            break;
        case math::StatisticsException::UNDEFINED :
            str << "Operation not defined for the given combination of arguments.";
            break;
        case math::StatisticsException::OPERATION_FAILED :
            str << "Unable to perform an operation.";
            break;
        case math::StatisticsException::INTEGER_OUT_OF_RANGE :
            str << "The solution is out of integer range.";
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
