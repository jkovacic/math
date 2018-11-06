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
 * @headername{BinomDist.h}
 *
 * Declaration of functions within the namespace BinomDist
 * that perform various binomial distribution related operations,
 * such as calculation of upper and lower tail probabilities and
 * quantiles, probability mass function, etc.
 */


#ifndef _MATH_BINOMDISTGENERIC_HPP_
#define _MATH_BINOMDISTGENERIC_HPP_


#include "../settings/probdist_settings.h"
#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various binomial distribution related operations, such as
 * calculation of upper and lower tail probabilities and
 * quantiles, probability mass function, etc.
 */
namespace BinomDist
{

    template <typename F, typename I>
    F mean(
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2)
      );


    template <typename F, typename I>
    F var(
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2)
      );


    template <typename F, typename I>
    F stdev(
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2)
      );


    template <typename F, typename I>
    bool normalApprox(
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2),
        const F& th = static_cast<F>(PROBDIST_BINOM_DIST_NORMAL_APPROX_NUM) / static_cast<F>(PROBDIST_BINOM_DIST_NORMAL_APPROX_DEN)
      );


    template <typename F, typename I>
    F pmf(
        const I& k,
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2)
      );


    template <typename F, typename I>
    F probInt(
        const I& a,
        const I& b,
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2),
        const bool incLower = true,
        const bool incUpper = true
      );


    template <typename F, typename I>
    F prob(
        const I& k,
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2),
        const bool incl = true,
        const bool lowerTail = true
      );


    template <typename F, typename I>
    I quant(
        const F& prob,
        const I& n,
        const F& p = static_cast<F>(1) / static_cast<F>(2),
        const bool smallest = true,
        const bool lowerTail = true
      );

}  // namespace BinomDist

}  // namespace math

// DEFINITION
#include "statistics/dist/BinomDistGeneric.cpp"

#endif  // _MATH_BINOMDISTGENERIC_HPP_
