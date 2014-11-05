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
 * Declaration of the class CurveFittingException
 */

#ifndef _MATH_CURVEFITTINGEXCEPTION_H_
#define _MATH_CURVEFITTINGEXCEPTION_H_

#include <iostream>


namespace math
{

/**
 * @brief An exception, typically thrown by functions of various curve 
 * fitting/interpolation classes, implemented by inherited classes of
 * CurveFittingGenericAb.
 */
struct CurveFittingException
{
    /// Enum with possible error codes
    enum err_codes {
        OUT_OF_MEMORY,              /// Could not allocate enough memory
        ADD_POINT_NOT_ALLOWED,      /// Adding of a point is not allowed
        NO_POINTS,                  /// No points added yet
        DUPLICATE_POINTS,           /// Duplicate points (with the same abscise value)
        CURVE_NOT_GENERATED,        /// Curve not generated yet
        CURVE_ALREADY_GENERATED,    /// Curve is already generated
        OUT_OF_BOUNDS,              /// Input value out of bounds
        CURVE_GENERATION_FAILED,    /// Curve generation was not successful
        INVALID_ARGUMENT            /// Invalid argument passed to constructor
    };

    err_codes error;     /// Type of an error
    // Constructor
    CurveFittingException(err_codes err);
    // Output a short description of the error
    void what(std::ostream& str = std::cerr) const;
};

} // namespace math

#endif // _MATH_CURVEFITTINGEXCEPTION_H_

