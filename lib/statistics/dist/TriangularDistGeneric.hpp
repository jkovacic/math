/*
Copyright 2016, Jernej Kovacic

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
 * @headername{TriangularDist.h}
 *
 * Declaration of functions within the namespace TriangularDist
 * that perform various triangular distribution related operations,
 * such as calculation of upper and lower tail probabilities
 * and quantiles, probability distribution function, etc.
 */


#ifndef _MATH_TRIANGULARDISTGENERIC_HPP_
#define _MATH_TRIANGULARDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform various triangular
 * distribution related operations, such as calculation of upper and
 * lower tail probabilities and quantiles, probability distribution
 * function, etc.
 */
namespace TriangularDist
{

    template <typename F>
    F pdf(
          const F& x,
          const F& a,
          const F& b,
          const F& c
        );


    template <typename F>
    F probInt(
          const F& from,
          const F& to,
          const F& a,
          const F& b,
          const F& c
        );


    template <typename F>
    F prob(
          const F& x,
          const F& a,
          const F& b,
          const F& c,
          const bool lowerTail = true
        );


    template <typename F>
    F quant(
          const F& p,
          const F& a,
          const F& b,
          const F& c,
          const bool lowerTail = true
        );


}  // namespace TriangularDist

}  // namespace math

// DEFINITION
#include "statistics/dist/TriangularDistGeneric.cpp"

#endif  // _MATH_TRIANGULARDISTGENERIC_HPP_
