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
 * @headername{StudentDistGeneric.h}
 *
 * Declaration of functions within the namespace StudentDist
 * that perform various Student's (Gosset's) distribution related
 * operations, such as calculation of upper and lower tail
 * probabilities and quantiles, calculation of t statistics,
 * probability distribution function, etc.
 */


#ifndef _MATH_STUDENTDISTGENERIC_HPP_
#define _MATH_STUDENTDISTGENERIC_HPP_


#include <cstddef>

#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various Student's distribution related operations, such as
 * calculation of upper and lower tail probabilities and
 * quantiles, calculation of t-statistics, probability
 * distribution function, etc.
 */
namespace StudentDist
{

    template <class T>
    T getT(
        const T& x,
        size_t n = 1,
        const T& mu = static_cast<T>(0),
        const T& s = static_cast<T>(1)
      ) throw (StatisticsException);


    template <class T>
    T getX(
        const T& t,
        size_t n = 1,
        const T& mu = static_cast<T>(0),
        const T& s = static_cast<T>(1)
      ) throw (StatisticsException);



    template <class T>
    T pdf(
        const T& x,
        const T& df,
        const T& mu = static_cast<T>(0),
        const T& sigma = static_cast<T>(1)
      ) throw (StatisticsException);


    template <class T>
    T probInt(
        const T& a,
        const T& b,
        const T& df,
        const T& mu = static_cast<T>(0),
        const T& sigma = static_cast<T>(1)
      ) throw (StatisticsException);



    template <class T>
    T prob(
        const T& x,
        const T& df,
        bool lowerTail = true,
        const T& mu = static_cast<T>(0),
        const T& sigma = static_cast<T>(1)
      ) throw (StatisticsException);



    template <class T>
    T quant(
        const T& p,
        const T& df,
        bool lowerTail = true,
        const T& mu = static_cast<T>(0),
        const T& sigma = static_cast<T>(1)
      ) throw (StatisticsException);

}  // namespace StudentDist

}  // namespace math


// DEFINITION
#include "statistics/StudentDistGeneric.cpp"

#endif  // _MATH_STUDENTDISTGENERIC_HPP_
