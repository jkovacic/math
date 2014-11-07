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
 * mean, variance, standard deviation, sum, etc.
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
 * mean, variance, standard deviation, sum, etc.
 *
 * The class is stateful. A vector of samples can only be passed to the class's
 * constructor. Before any getter methods are called, process() must be executed
 * on an instance of the class. This method calculates and internally keeps all
 * values that are necessary by getter methods. The method processed() may be used
 * to verify whether the results are available to getter methods.
 *
 * Please note that an instance of this class does not create its own copy
 * of the sample vector. Instead it only keeps a reference to a vector. The
 * vector of samples should not be modified until process() is finished!
 * Hence the class cannot be considered as thread safe, however this behavior
 * is acceptable for majority of applications. After process() is finished,
 * the instance of this class will not need the vector of samples anymore.
 *
 * The class itself does not modify the vector of samples at all.
 */
template <class T>
class SampleStatGeneric
{

private:
    const std::vector<T>* m_pData;
    size_t m_N;
    T m_sum;
    T m_sumSqDev;
    T m_mean;
    bool m_processed;

public:
    // Constructor
    SampleStatGeneric(const std::vector<T>& sample) throw(StatisticsException);

    // Process
    void process();

    // Methods to obtain results of the processed sample
    bool processed();
    size_t sampleSize();
    T sum() throw(StatisticsException);
    T mean() throw(StatisticsException);
    T var(bool sample=true) throw(StatisticsException);
    T var(size_t df_sub) throw(StatisticsException);
    T stdev(bool samplen=true) throw(StatisticsException);
    T stdev(size_t df_sub) throw(StatisticsException);
};  // class SampleStatGeneric


// Declaration of specialized methods inside the name space declaration
// is essential if implemented elsewhere:
template<> float SampleStatGeneric<float>::stdev(size_t df_sub) throw (StatisticsException);
template<> double SampleStatGeneric<double>::stdev(size_t df_sub) throw (StatisticsException);
template<> long double SampleStatGeneric<long double>::stdev(size_t df_sub) throw (StatisticsException);


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
