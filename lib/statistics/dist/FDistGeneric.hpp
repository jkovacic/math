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
 * @headername{FDistGeneric.h}
 *
 * Declaration of functions within the namespace FDist
 * that perform various F-distribution (a.k.a.Fisher - Snedecor distribution)
 * related operations, such as calculation of upper and lower tail
 * probabilities and quantiles, probability distribution function, etc.
 */


#ifndef _MATH_FDISTGENERIC_HPP_
#define _MATH_FDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various F-distribution (a.k.a.Fisher-Snedecor distribution)
 * related operations, such as calculation of upper and lower
 * tail probabilities and quantiles, probability distribution
 * function, etc.
 */
namespace FDist
{

    template <class T>
    T pdf(
          const T& x,
          const T& d1,
          const T& d2
        ) throw (StatisticsException);


    template <class T>
    T probInt(
          const T& a,
          const T& b,
          const T& d1,
          const T& d2
        ) throw (StatisticsException);


    template <class T>
    T prob(
          const T& x,
          const T& d1,
          const T& d2,
          bool lowerTail = true
       ) throw (StatisticsException);


    template <class T>
    T quant(
          const T& p,
          const T& d1,
          const T& d2,
          bool lowerTail = true
        ) throw (math::StatisticsException);

}  // namespace FDist

}  // namespace math


// DEFINITION
#include "statistics/dist/FDistGeneric.cpp"

#endif  // _MATH_FDISTGENERIC_HPP_
