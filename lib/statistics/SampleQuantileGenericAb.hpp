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
 *
 * Declaration of the class SampleQuantileGenericAb, a superclass
 * of all classes that estimate quantiles of a sample.
 */


#ifndef _MATH_SAMPLEQUANTILEGENERICAB_HPP_
#define _MATH_SAMPLEQUANTILEGENERICAB_HPP_

#include <cstddef>
#include <vector>
#include <set>

#include "../settings/stat_settings.h"
#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief Supported methods to estimate a sample quantile.
 */
struct EQntlType
{
    enum type
    {
        R1,                       /// R quantile, type 1
        R2,                       /// R quantile, type 2
        R3,                       /// R quantile, type 3
        R4,                       /// R quantile, type 4
        R5,                       /// R quantile, type 5
        R6,                       /// R quantile, type 6
        R7,                       /// R quantile, type 7
        R8,                       /// R quantile, type 8
        R9,                       /// R quantile, type 9
        SAS1,                     /// SAS QNTL, method 1
        SAS2,                     /// SAS QNTL, method 2
        SAS3,                     /// SAS QNTL, method 3
        SAS4,                     /// SAS QNTL, method 4
        SAS5,                     /// SAS QNTL, method 5
        EXCEL,                    /// MS Excel, function PERCENTILE
        SCIPY_0_0,                /// SciPy scipy.stats.mstats.mquantiles, alphap=0, betap=0
        SCIPY_0_1,                /// SciPy scipy.stats.mstats.mquantiles, alphap=0, betap=1
        SCIPY_05_05,              /// SciPy scipy.stats.mstats.mquantiles, alphap=0.5, betap=0.5
        SCIPY_1_1,                /// SciPy scipy.stats.mstats.mquantiles, alphap=1, betap=1
        SCIPY_13_13,              /// SciPy scipy.stats.mstats.mquantiles, alphap=1/3, betap=1/3
        SCIPY_38_38,              /// SciPy scipy.stats.mstats.mquantiles, alphap=3/8, betap=3/8
        SCIPY_N05_N05,            /// SciPy scipy.stats.mstats.mquantiles, alphap=-1/2, betap=-1/2
    };
};


/**
 * @brief A base abstarct class for all classes that estimate sample's quantile
 * for any probability within a valid range. It supports several estimation
 * methods, however, it relies on several pure virtual methods that are implemented
 * in derived classes.
 *
 * As this class is abstract, it cannot be instantiated directly.
 */
template <typename F>
class SampleQuantileGenericAb
{

protected:
    const std::size_t m_N;

protected:
   // Constructor (called by derived classes' constructors)
   SampleQuantileGenericAb(const std::size_t N);

private:
    // internally used linear interpolation to estimate noninteger ordinals
    F __linIntrp(const F& h) const;

private:
    // pure virtual methods that select (depending on implementation) n.th smallest element(s)
    virtual F _select(const std::size_t n) const throw(StatisticsException) = 0;

    virtual void _select2(
            const std::size_t n1,
            const std::size_t n2,
            F& val1,
            F& val2
          ) const throw(StatisticsException) = 0;

protected:
    void _outlierBounds(
            F& lower,
            F& upper,
            const F& iqrs,
            const EQntlType::type method
          ) const;

public:

    std::size_t sampleSize() const;

    template <typename I>
    F quantile(
           const I& num,
           const I& den, 
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG 
         ) const throw (StatisticsException);

    F qntl(
           const F& p,
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG 
         ) const throw (StatisticsException);

    F median(const bool approx = false) const;

    F quartile(
           const unsigned short int q,
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG 
         ) const throw (StatisticsException);

    F iqr(const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG) const;

    bool isOutlier(
               const F& val,
               const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
               const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
            )  const;

    // The following methods are pure virtual, i.e. implemented in derived classes
    virtual F ecdf(const F& t) const = 0;

    virtual F min() const = 0;

    virtual F max() const = 0;

    virtual F elem(
           const std::size_t n, 
           const bool largest = true,
           const bool zerobase = STAT_DEFAULT_ZERO_BASE
         ) const throw(StatisticsException) = 0;

    virtual void outliers(
           std::set<F>& outl,
           const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
         ) const throw (StatisticsException) = 0;

    // Destructor
    virtual ~SampleQuantileGenericAb() = 0;

};  // class SampleQuantileGenericAb

}  // namespace math


// DEFINITION
#include "statistics/SampleQuantileGenericAb.cpp"

#endif  // _MATH_SAMPLEQUANTILEGENERICAB_HPP_
