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
 * @headername{SampleQuantileGeneric.h}
 *
 * Declaration of the class SampleQuantileGeneric that estimates
 * quantiles of a sample.
 */

#ifndef _MATH_SAMPLEQUANTILEGENERIC_HPP_
#define _MATH_SAMPLEQUANTILEGENERIC_HPP_

#include <cstddef>
#include <vector>

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
 * @brief A class that estimates sample's quantile for any probabilty within
 * a valid range. Several estimation methods are supported.
 *
 * The class creates its own copy of the sample vector.
 */
template <class T>
class SampleQuantileGeneric
{

private:
    std::vector<T> m_v;
    size_t m_N;

public:
    // Constructor
    SampleQuantileGeneric(const std::vector<T>& sample) throw(StatisticsException);

    // Methods to obtain quantiles of the sample:
    size_t sampleSize();
    T quantile(size_t num, size_t den, EQntlType::type method=EQntlType::R7) throw(StatisticsException);
    T qntl(double p, EQntlType::type method=EQntlType::R7) throw(StatisticsException);
    T median();
    T iqr(EQntlType::type method=EQntlType::R7);

    // Destructor
    virtual ~SampleQuantileGeneric();

private:
    // internally used functions:

    /*
     * Converts a (positive) real number to an integer (size_t).
     * Note that a rounded float/double may be actually a bit less
     * than expected, so addition of 0.5 ensures that the number is
     * properly converted (truncated) to an integer.
     */
    static inline size_t dbl2int(double n) { return static_cast<size_t>(n+0.5); }

    T linIntrp(double h);

};  // class SampleQuantileGeneric


// Samples with elements of types float, double and long double
// make most sense, therefore the following types are predefined

typedef SampleQuantileGeneric<float>        FSampleQuantile;
typedef SampleQuantileGeneric<double>       SampleQuantile;
typedef SampleQuantileGeneric<long double>  LDSampleQuantile;

}  // namespace math


// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "SampleQuantileGeneric.cpp"

#endif   // _MATH_SAMPLEQUANTILEGENERIC_HPP_