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
 * Implementation of the class CurveFittingException
*/


#include <ostream>

#include "exception/CurveFittingException.hpp"

/**
 * Constructor.
 *
 * @param error code
 */
math::CurveFittingException::CurveFittingException(const math::CurveFittingException::err_codes err)
{
    this->error = err;
}

/**
 * Outputs a short error description to stdout
 *
 * @param str - stream, the error description will be written in (default: cerr)
 */
void math::CurveFittingException::what(std::ostream& str) const
{
    switch (error)
    {
        case math::CurveFittingException::OUT_OF_MEMORY :
            str << "Could not allocate enough memory";
            break;

        case math::CurveFittingException::ADD_POINT_NOT_ALLOWED :
            str << "Prohibited to enter more points";
            break;

        case math::CurveFittingException::NO_POINTS :
            str << "No points have been entered so far";
            break;

        case math::CurveFittingException::DUPLICATE_POINTS :
            str << "Duplicate points (with the same value of abscise) have been entered";
            break;

        case math::CurveFittingException::CURVE_NOT_GENERATED :
            str << "Curve has not been generated yet";
            break;

        case math::CurveFittingException::CURVE_ALREADY_GENERATED :
            str << "Curve is already generated";
            break;

        case math::CurveFittingException::OUT_OF_BOUNDS :
            str << "Input is out of the definition range";
            break;

        case math::CurveFittingException::CURVE_GENERATION_FAILED :
            str << "Could not generate a curve";
            break;

        case math::CurveFittingException::INVALID_ARGUMENT :
            str << "Invalid argument passed";
            break;

        default:
            // Should not occur but handle it anyway.
            // Maybe a code was inserted into err_codes and
            // this function hasn't been updated?
            str << "Strange, unspecified error";
    }  // switch

    str << std::endl;
}
