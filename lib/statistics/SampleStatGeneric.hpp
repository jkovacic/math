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
 * @headername{SampleStatGeneric.h}
 *
 * Declaration of functions in the namespace SampleStat that calculate sample's
 * sum, mean, variance, standard deviation, covariance, Pearson's r
 * (correlation), r squared, etc.
 */

#ifndef _MATH_SAMPLESTATGENERIC_HPP_
#define _MATH_SAMPLESTATGENERIC_HPP_

#include <cstddef>
#include <vector>

#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief A namespace containing functions that calculate sample's basic statistics, 
 * such as sum, mean, variance, standard deviation, covariance, Pearson'r
 * (correlation), etc.
 */
namespace SampleStat
{

    template <typename F>
    F min(const std::vector<F>& x) throw(StatisticsException);

    template <typename F>
    F max(const std::vector<F>& x) throw(StatisticsException);

    template <typename F>
    F sum(const std::vector<F>& x);

    template <typename F>
    F mean(const std::vector<F>& x) throw(StatisticsException);

    template <typename F>
    F var(const std::vector<F>& x, const bool sample=true) throw(StatisticsException);

    template <typename F>
    F var(const std::vector<F>& x, const size_t df_sub) throw(StatisticsException);

    template <typename F>
    F stdev(const std::vector<F>& x, const bool sample=true) throw(StatisticsException);

    template <typename F>
    F stdev(const std::vector<F>& x, const size_t df_sub) throw(StatisticsException);

    template <typename F>
    F cov(const std::vector<F>& x1, const std::vector<F>& x2, const size_t df_sub) throw(StatisticsException);

    template <typename F>
    F cov(const std::vector<F>& x1, const std::vector<F>& x2, const bool sample=true) throw(StatisticsException);

    template <typename F>
    F cor(const std::vector<F>& x1, const std::vector<F>& x2) throw(StatisticsException);

    template <typename F>
    F r2(const std::vector<F>& x1, const std::vector<F>& x2) throw(StatisticsException);

    template <typename F, typename I>
    F moment(const std::vector<F>& x, const I& n, const F& about=static_cast<F>(0)) throw(StatisticsException);

    template <typename F, typename I>
    F centralMoment(const std::vector<F>& x, const I& n) throw(StatisticsException);

    template <typename F>
    F skewness(const std::vector<F>& x) throw(StatisticsException);

    template <typename F>
    F kurtosis(const std::vector<F>& x) throw(StatisticsException);

}  // namespace SampleStat

}  // namespace math


// DEFINITION
#include "statistics/SampleStatGeneric.cpp"

#endif  // _MATH_SAMPLESTATGENERIC_HPP_
