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
 * @headername{PoissonDistGeneric.h}
 *
 * Declaration of functions within the namespace PoissonDist
 * that perform various Poisson distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability mass function, etc.
 */

#ifndef _MATH_POISSONDISTGENERIC_HPP_
#define _MATH_POISSONDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various Poisson distribution related operations, such as
 * calculation of upper and lower tail probabilities and
 * quantiles, probability mass function, etc.
 */
namespace PoissonDist
{

    template <typename F, typename I>
    F pmf(
        const I& k,
        const F& lambda
      ) throw(StatisticsException);


    template <typename F, typename I>
    F probInt(
        const I& a,
        const I& b,
        const F& lambda,
        bool incLower = true,
        bool incUpper = true
      ) throw (StatisticsException);


    template <typename F, typename I>
    F prob(
        const I& k,
        const F& lambda,
        bool incl = true,
        bool lowerTail = true
      ) throw(StatisticsException);


    template <typename F, typename I>
    I quant(
        const F& prob,
        const F& lambda,
        bool smallest = true,
        bool lowerTail = true
      ) throw(StatisticsException);

}  // namespace PoissonDist

}  // namespace math

// DEFINITION
#include "statistics/dist/PoissonDistGeneric.cpp"

#endif	// _MATH_POISSONDISTGENERIC_HPP_
