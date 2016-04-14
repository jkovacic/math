/*
Copyright 2016, Jernej Kovacic

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
 * @headername{SampleQuantileSelection.h}
 *
 * Declaration of the class SampleQuantileSelectionGeneric that
 * estimates quantiles of a sample, based on implemented functionality
 * for selection of the n.th smallest element.
 */


#ifndef _MATH_SAMPLEQUANTILESELECTIONGENERIC_HPP_
#define _MATH_SAMPLEQUANTILESELECTIONGENERIC_HPP_

#include <cstddef>
#include <vector>
#include <set>

#include "statistics/SampleQuantileGenericAb.hpp"

#include "../settings/stat_settings.h"

#include "exception/StatisticsException.hpp"


namespace math
{


/**
 * @brief A class that estimates sample's quantile for any probability within
 * a valid range. Several estimation methods are supported.
 *
 * The class does not sort samples (which might be a costly operation).
 * Instead it relies on the selection functionality (functions that efficiently
 * select the n.th smallest element).
 *
 * It is advisable to instantiate this class to perform selection on the original
 * sample vector if no elements are inserted to or removed from the original vector
 * (rearrangement of elements is permitted) during the effective usage of this
 * class. When this is not possible, the class can be instantiated to hold an internal
 * copy of the original sample vector. Please note that this option requires additional
 * memory.
 *
 * This class is convenient when the sample is large and/or
 * a handful of sample quantiles need to be obtained.
 */
template <typename F>
class SampleQuantileSelectionGeneric : public math::SampleQuantileGenericAb<F>
{

private:
    // vector to store internal copy of the sample when applicable
    std::vector<F> m_stor;
    // reference to the actual vector of observations, referred by methods of this class
    const std::vector<F>& m_v;

private:
    virtual F _select(const std::size_t n) const throw(StatisticsException);

    virtual void _select2(
            const std::size_t n1,
            const std::size_t n2,
            F& val1,
            F& val2
          ) const throw(StatisticsException);

public:
    // Constructor
    SampleQuantileSelectionGeneric(const std::vector<F>& sample, const bool copy=false) throw(StatisticsException);

    virtual F ecdf(const F& t) const;

    virtual F min() const;

    virtual F max() const;

    virtual F elem(
           const std::size_t n,
           const bool largest = true,
           const bool zerobase = STAT_DEFAULT_ZERO_BASE
         ) const throw(StatisticsException);

    virtual void outliers(
           std::set<F>& outl,
           const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
         ) const throw (StatisticsException);

    // Destructor
    virtual ~SampleQuantileSelectionGeneric();
};


// Samples with elements of types float, double and long double
// make most sense, therefore the following types are predefined

typedef SampleQuantileSelectionGeneric<float>              FSampleQuantileSelection;
typedef SampleQuantileSelectionGeneric<double>             SampleQuantileSelection;
typedef SampleQuantileSelectionGeneric<long double>        LDSampleQuantileSelection;

}  // namespace math


// DEFINITION
#include "statistics/SampleQuantileSelectionGeneric.cpp"

#endif  // _MATH_SAMPLEQUANTILESELECTIONGENERIC_HPP_
