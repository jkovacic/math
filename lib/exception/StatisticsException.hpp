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
 * An internal header file, it should not be included directly.
 * @headername{StatisticsException.h}
 *
 * Declaration of the class StatisticsException.
 */


#ifndef _MATH_STATISTICSEXCEPTION_HPP_
#define	_MATH_STATISTICSEXCEPTION_HPP_

#include "exception/IMathException.hpp"

#include <ostream>

namespace math
{

/**
 * @brief An exception, typically thrown by functions of statistics
 *        related classes: SampleStatGeneric, etc.
 */
struct StatisticsException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes
    {
        SAMPLE_EMPTY,                   /// empty sample
        SAMPLE_TOO_SMALL,               /// sample too small
        DF_SUBTRAHEND_TOO_LARGE,        /// DF subtrahend too large
        UNSUPPORTED_TYPE,               /// operation not supported for this type
        UNEQUAL_SAMPLE_SIZES,           /// samples' sizes are not equal
        OUT_OF_MEMORY,                  /// allocation of memory failed
        INVALID_PROBABILTY,             /// probability is invalid
        UNSUPPORTED_QUANTILE_METHOD,    /// unsupported method to estimate the quantile
        INVALID_STDEV,                  /// negative or zero standard deviation
        INVALID_DF,                     /// negative or zero degrees of freedom
        OPERATION_FAILED,               /// unable to perform an operation
    };

    err_codes error;    /// type of an error
    // Constructor
    StatisticsException(err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif	// _MATH_STATISTICSEXCEPTION_HPP_
