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
 * @headername{NormalDistGeneric.h}
 *
 * Declaration of functions within the namespace NormalDist
 * that perform various normal distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, calculation of z statistics, probability
 * distribution function, etc.
 */


#ifndef _MATH_NORMALDISTGENERIC_HPP_
#define _MATH_NORMALDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various normal distribution related operations, such as
 * calculation of upper and lower tail probabilities and
 * quantiles, calculation of z statistics, probability
 * distribution function, etc.
 */
namespace NormalDist
{

    template <class T>
    T getZ(
          const T& x,
          const T& mu = static_cast<T>(0),
          const T& sigma = static_cast<T>(1)
        ) throw (StatisticsException);


    template <class T>
    T getX(
          const T& z,
          const T& mu = static_cast<T>(0),
          const T& sigma = static_cast<T>(1)
        ) throw (StatisticsException);


    template <class T>
    T pdf(
          const T& x,
          const T& mu = static_cast<T>(0),
          const T& sigma = static_cast<T>(1)
        ) throw (StatisticsException);


    template <class T>
    T probInt(
          const T& a,
          const T& b,
          const T& mu = static_cast<T>(0),
          const T& sigma = static_cast<T>(1)
        ) throw (StatisticsException);


    template <class T>
    T prob(
          const T& x,
          const T& mu = static_cast<T>(0),
          const T& sigma = static_cast<T>(1),
          bool lowerTail = true
        ) throw (StatisticsException);


    template <class T>
    T quant(
          const T& p,
          const T& mu = static_cast<T>(0),
          const T& sigma = static_cast<T>(1),
          bool lowerTail = true
        ) throw (StatisticsException);


}  // namespace NormalDist

}  // namespace math

// DEFINITION
#include "statistics/dist/NormalDistGeneric.cpp"

#endif  // _MATH_NORMALDISTGENERIC_HPP_
