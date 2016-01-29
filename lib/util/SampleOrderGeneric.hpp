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
 * @headername{SampleOrderGeneric.h}
 *
 * Declaration of functions within the namespace SampleOrder
 * that return indices of the given vector's elements in a stably
 * sorted vector.
 */

#ifndef _MATH_SAMPLEORDERGENERIC_HPP_
#define _MATH_SAMPLEORDERGENERIC_HPP_

#include <cstddef>
#include <vector>

#include "exception/SampleOrderException.hpp"


namespace math
{


/**
 * @brief A namespace with functions that search indices of
 *        stably sorted vectors of elements.
 */
namespace SampleOrder
{

    template <typename F>
    std::vector<std::size_t>& order(
            const std::vector<F>& x,
            std::vector<std::size_t>& dest,
            const bool asc=true
          ) throw(SampleOrderException);

    template <typename F>
    F min(const std::vector<F>& x) throw(SampleOrderException);

    template <typename F>
    F max(const std::vector<F>& x) throw(SampleOrderException);

    template <typename F>
    std::size_t whichMin(const std::vector<F>& x) throw(SampleOrderException);

    template <typename F>
    std::size_t whichMax(const std::vector<F>& x) throw(SampleOrderException);

}  // namespace SampleOrder
}  // namespace math


// DEFINITION
#include "SampleOrderGeneric.cpp"

#endif  // _MATH_SAMPLEORDERGENERIC_HPP_
