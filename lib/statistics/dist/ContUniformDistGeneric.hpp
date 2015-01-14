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
 * @headername{ContUniformDistGeneric.h}
 *
 * Declaration of functions within the namespace ContUniformDist
 * that perform various continuous unform (a.k.a. rectangular) distribution
 * related operations, such as calculation of upper and lower tail probabilities
 * and quantiles, probability distribution function, etc.
 */


#ifndef _MATH_CONTUNIFORMDISTGENERIC_HPP_
#define _MATH_CONTUNIFORMDISTGENERIC_HPP_


#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace with functions that perform various continuous uniform 
 * (rectangular) distribution related operations, such as calculation of upper
 * and lower tail probabilities and quantiles, probability distribution function, 
 * etc.
 */
namespace ContUniformDist
{

    template <typename F>
    F pdf(
          const F& x,
          const F& min = static_cast<F>(0),
          const F& max = static_cast<F>(1)
        ) throw (StatisticsException);


    template <typename F>
    F probInt(
          const F& a,
          const F& b,
          const F& min = static_cast<F>(0),
          const F& max = static_cast<F>(1)
        ) throw (StatisticsException);


    template <typename F>
    F prob(
          const F& x,
          const F& min = static_cast<F>(0),
          const F& max = static_cast<F>(1),
          bool lowerTail = true
        ) throw (StatisticsException);


    template <typename F>
    F quant(
          const F& p,
          const F& min = static_cast<F>(0),
          const F& max = static_cast<F>(1),
          bool lowerTail = true
        ) throw (StatisticsException);


}  // namespace ContUniformDist

}  // namespace math

// DEFINITION
#include "statistics/dist/ContUniformDistGeneric.cpp"

#endif  // _MATH_CONTUNIFORMDISTGENERIC_HPP_
