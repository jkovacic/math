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
 * Declaration of the class SampleStatGeneric that calculates sample's
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
 * @brief A class that calculates sample's basic statistics, such as
 * sum, mean, variance, standard deviation, covariance, Pearson'r
 * (correlation), etc.
 *
 * The whole class is static, so no instantiation is necessary.
 */
template <class T>
class SampleStatGeneric
{

public:
    static T min(const std::vector<T>& x) throw(StatisticsException);
    static T max(const std::vector<T>& x) throw(StatisticsException);
    static T sum(const std::vector<T>& x);
    static T mean(const std::vector<T>& x) throw(StatisticsException);
    static T var(const std::vector<T>& x, bool sample=true) throw(StatisticsException);
    static T var(const std::vector<T>& x, size_t df_sub) throw(StatisticsException);
    static T stdev(const std::vector<T>& x, bool sample=true) throw(StatisticsException);
    static T stdev(const std::vector<T>& x, size_t df_sub) throw(StatisticsException);
    static T cov(const std::vector<T>& x1, const std::vector<T>& x2, size_t df_sub) throw(StatisticsException);
    static T cov(const std::vector<T>& x1, const std::vector<T>& x2, bool sample=true) throw(StatisticsException);
    static T cor(const std::vector<T>& x1, const std::vector<T>& x2) throw(StatisticsException);
    static T r2(const std::vector<T>& x1, const std::vector<T>& x2) throw(StatisticsException);

private:
    static T __getShift(const std::vector<T>& x, size_t Nmax=5);
    static T __minmax(const std::vector<T>&x, bool min) throw(StatisticsException);

};  // class SampleStatGeneric


// Declaration of specialized methods inside the name space declaration
// is essential if implemented elsewhere:
template<> float SampleStatGeneric<float>::stdev(const std::vector<float>& x, size_t df_sub) throw (StatisticsException);
template<> double SampleStatGeneric<double>::stdev(const std::vector<double>& x, size_t df_sub) throw (StatisticsException);
template<> long double SampleStatGeneric<long double>::stdev(const std::vector<long double>& x, size_t df_sub) throw (StatisticsException);


// Samples with elements of types float, double and long double
// make most sense, therefore the following types are predefined
typedef SampleStatGeneric<float>        FSampleStat;
typedef SampleStatGeneric<double>       SampleStat;
typedef SampleStatGeneric<long double>  LDSampleStat;

}  // namespace math


// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "statistics/SampleStatGeneric.cpp"

#endif  // _MATH_SAMPLESTATGENERIC_HPP_
