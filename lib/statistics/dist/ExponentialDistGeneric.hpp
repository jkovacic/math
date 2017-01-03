/*
Copyright 2017, Jernej Kovacic

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
 * @headername{ExponentialDist.h}
 *
 * Declaration of functions within the namespace ExponentialDist
 * that perform various exponential distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability distribution function, etc.
 */

#ifndef _MATH_EXPONENTIALDISTGENERIC_HPP_
#define _MATH_EXPONENTIALDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various exponential distribution related operations, such as
 * calculation of upper and lower tail probabilities and
 * quantiles, probability distribution function, etc.
 */
namespace ExponentialDist
{

    template <typename F>
    F pdf(
        const F& x,
        const F& lambda
      ) throw(StatisticsException);


    template <typename F>
    F probInt(
        const F& a,
        const F& b,
        const F& lambda
      ) throw (StatisticsException);


    template <typename F>
    F prob(
        const F& x,
        const F& lambda,
        const bool lowerTail = true
      ) throw(StatisticsException);


    template <typename F>
    F quant(
        const F& prob,
        const F& lambda,
        const bool lowerTail = true
      ) throw(StatisticsException);

}  // namespace ExponentialDist

}  // namespace math

// DEFINITION
#include "statistics/dist/ExponentialDistGeneric.cpp"

#endif  // _MATH_EXPONENTIALDISTGENERIC_HPP_
