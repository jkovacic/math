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
 * @headername{SampleQuantileSortedArray.h}
 *
 * Declaration of the class SampleQuantileSortedArrayGeneric that
 * estimates quantiles of a sample, based on a sorted copy of the
 * sample vector.
 */


#ifndef _MATH_SAMPLEQUANTILESORTEDARRAYGENERIC_HPP_
#define _MATH_SAMPLEQUANTILESORTEDARRAYGENERIC_HPP_

#include <cstddef>
#include <vector>
#include <set>

#include "statistics/SampleQuantileGenericAb.hpp"

#include "exception/StatisticsException.hpp"


namespace math
{


/**
 * @brief A class that estimates sample's quantile for any probability within
 * a valid range. Several estimation methods are supported.
 *
 * The class creates its own copy of the sample vector and sorts the copy in
 * ascending order. This operation requires additional space, however, after
 * the instantiation, the original sample vector is not needed anymore and can
 * be manipulated freely. The initial sorting might be a costly operation,
 * however, this "price" is only "paid" at the instantiation of the class and
 * all other operations are performed quickly.
 *
 * This class is convenient when the sample is not too large and/or
 * many sample quantiles need to be obtained.
 */
template <typename F>
class SampleQuantileSortedArrayGeneric : public math::SampleQuantileGenericAb<F>
{

private:
    // internal storage to hold a sorted copy of the sample:
    std::vector<F> m_v;

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
    SampleQuantileSortedArrayGeneric(const std::vector<F>& sample) throw(StatisticsException);

    virtual F ecdf(const F& t) const;

    virtual F min() const;

    virtual F max() const;

    virtual F elem(
           const std::size_t n,
           const bool largest = true,
           const bool zerobase = STAT_DEFAULT_ZERO_BASE
         ) const throw(StatisticsException);

    virtual bool isOutlier(
           const F& val,
           const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
        )  const;

    virtual void outliers(
           std::set<F>& outl,
           const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
         ) const throw (StatisticsException);

    // Destructor
    virtual ~SampleQuantileSortedArrayGeneric();
};


// Samples with elements of types float, double and long double
// make most sense, therefore the following types are predefined

typedef SampleQuantileSortedArrayGeneric<float>          FSampleQuantileSortedArray;
typedef SampleQuantileSortedArrayGeneric<double>         SampleQuantileSortedArray;
typedef SampleQuantileSortedArrayGeneric<long double>    LDSampleQuantileSortedArray;

}  // namespace math


// DEFINITION
#include "statistics/SampleQuantileSortedArrayGeneric.cpp"

#endif  // _MATH_SAMPLEQUANTILESORTEDARRAYGENERIC_HPP_
