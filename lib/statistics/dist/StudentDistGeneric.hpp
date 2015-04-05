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

    template <typename F, typename I>
    F getT(
        const F& x,
        const I& n = static_cast<I>(1),
        const F& mu = static_cast<F>(0),
        const F& s = static_cast<F>(1)
      ) throw (StatisticsException);


    template <typename F, typename I>
    F getX(
        const F& t,
        const I& n = static_cast<I>(1),
        const F& mu = static_cast<F>(0),
        const F& s = static_cast<F>(1)
      ) throw (StatisticsException);



    template <typename F>
    F pdf(
        const F& x,
        const F& df,
        const F& mu = static_cast<F>(0),
        const F& sigma = static_cast<F>(1)
      ) throw (StatisticsException);


    template <typename F>
    F probInt(
        const F& a,
        const F& b,
        const F& df,
        const F& mu = static_cast<F>(0),
        const F& sigma = static_cast<F>(1)
      ) throw (StatisticsException);



    template <typename F>
    F prob(
        const F& x,
        const F& df,
        const bool lowerTail = true,
        const F& mu = static_cast<F>(0),
        const F& sigma = static_cast<F>(1)
      ) throw (StatisticsException);



    template <typename F>
    F quant(
        const F& p,
        const F& df,
        const bool lowerTail = true,
        const F& mu = static_cast<F>(0),
        const F& sigma = static_cast<F>(1)
      ) throw (StatisticsException);

}  // namespace StudentDist

}  // namespace math


// DEFINITION
#include "statistics/dist/StudentDistGeneric.cpp"

#endif  // _MATH_STUDENTDISTGENERIC_HPP_
