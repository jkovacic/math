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
 * An internal header file, it should not be included directly.
 *
 * Declaration of auxiliary functions that handle integer variables.
 */


#ifndef _MATH_INTUTILGENERIC_HPP_
#define	_MATH_INTUTILGENERIC_HPP_


namespace math
{

/*
 * @brief A namespace with utilities to handle integer variables
 */
namespace IntUtil
{

    template <typename I>
    bool isNegative(const I& n);

    template <typename I>
    I absolute(const I& n);

}  // namespace IntUtil

}  // namespace math

// DEFINITION
#include "int_util/IntUtilGeneric.cpp"

#endif	 // _MATH_INTUTILGENERIC_HPP_
