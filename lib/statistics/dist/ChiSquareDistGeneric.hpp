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
 * @headername{ChiSquareDistGeneric.h}
 *
 * Declaration of functions within the namespace ChiSquareDist
 * that perform various chi-squared distribution related
 * operations, such as calculation of upper and lower tail
 * probabilities and quantiles, probability distribution function, etc.
 */


#ifndef _MATH_CHISQUAREDISTGENERIC_HPP_
#define _MATH_CHISQUAREDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform
 * various chi-squared distribution related operations, such as
 * calculation of upper and lower tail probabilities and
 * quantiles, probability distribution function, etc.
 */
namespace ChiSquareDist
{

    template <typename F>
    F pdf(
          const F& x,
          const F& df
        ) throw (StatisticsException);


    template <typename F>
    F probInt(
          const F& a,
          const F& b,
          const F& df
        ) throw (StatisticsException);


    template <typename F>
    F prob(
          const F& x,
          const F& df,
          bool lowerTail = true
       ) throw (StatisticsException);


    template <typename F>
    F quant(
          const F& p,
          const F& df,
          bool lowerTail = true
        ) throw (math::StatisticsException);

}  // namespace ChiSquareDist

}  // namespace math


// DEFINITION
#include "statistics/dist/ChiSquareDistGeneric.cpp"

#endif  // _MATH_CHISQUAREDISTGENERIC_HPP_
