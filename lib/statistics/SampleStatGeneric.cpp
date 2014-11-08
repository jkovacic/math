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
 * Implementation of the class SampleStatGeneric that calculates sample's
 * mean, variance, standard deviation, sum, etc.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */


// deliberately there is no #include "SampleStatGeneric.hpp" !
#include <cstddef>
#include <cmath>
#include <vector>

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"

#include "exception/StatisticsException.hpp"
#include "util/NumericUtil.hpp"


/**
 * @param x - vector of sample elements
 *
 * @return sum of all sample values
 */
template <class T>
T math::SampleStatGeneric<T>::sum(const std::vector<T>& x)
{
    const size_t N = x.size();

    if ( 0 == N )
    {
        return math::NumericUtil<T>::ZERO;
    }

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    // Ideal number of threads
    // (if each one processes approx. OMP_CHUNKS_PER_THREAD items):
    const size_t ideal = N / OMP_CHUNKS_PER_THREAD +
                         ( 0 == N % OMP_CHUNKS_PER_THREAD ? 0 : 1 );

    /*
     * In the first step, each thread calculates the sum of its block
     */
    T sum = math::NumericUtil<T>::ZERO;
    #pragma omp parallel num_threads(ideal) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(x) \
                    reduction(+ : sum)
    {
        // Depending on the number of available threads,
        // determine the ideal nr.of samples per thread,
        // and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (N + nthreads -1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;
        const size_t iend = istart + samples_per_thread;

        // Calculate the sum of the assigned block...
        T partsum = math::NumericUtil<T>::ZERO;
        for ( size_t i=istart; i<iend && i<N; ++i )
        {
            partsum += x.at(i);
        }

        // ... and add it to the total sum in a thread safe manner.
        sum += partsum;
    }

    return sum;

    // In serial mode this variable is never used.
    (void) ideal;
}


/**
 * Arithmetical mean (or average) of the sample.
 *
 * @param x - vector of sample elements
 *
 * @return mean value of the sample
 *
 * @throw StatisticsException if 'x' is empty
 */
template <class T>
T math::SampleStatGeneric<T>::mean(const std::vector<T>& x) throw(math::StatisticsException)
{
    /*
     * Arithmetical mean is calculated as:
     *
     *                   N
     *                 -----
     *             1   \
     * mean(X) = -----  >  X[i]
     *             N   /
     *                 -----
     *                  i=1
     *
     */

    const size_t N = x.size();

    if ( 0 == N )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    return sum(x) / static_cast<T>(N);
}


/**
 * Variance of the sample.
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * the sample size.
 *
 * @param x - vector of sample elements
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return variance of the sample, depending on the given 'df_sub'
 *
 * @throw StatisticsException if 'x' is empty or if 'df_sub' exceeds sample's size
 */
template <class T>
T math::SampleStatGeneric<T>::var(const std::vector<T>& x, size_t df_sub) throw(math::StatisticsException)
{
    /*
     * The best known algorithm to calculate a variance is:
     *
     *                              N
     *                            -----
     *                  1         \                 2
     *   var(X) = ------------     >   (X[i] - mean)
     *             N - df_sub     /
     *                            -----
     *                             i=1
     *
     * It is called a two step algorithm as it must calculate the sample mean
     * in the first step, followed by the algorithm above to calculate the variance.
     *
     * When the mean is not required, the expression above can be replaced by the
     * one step algorithm:
     *
     *
     *               N                       /   N           \  2
     *             -----                     | -----         |
     *             \              2      1   | \             |
     *              >     (X[i]-K)   -  ---  |  >   (X[i]-K) |
     *             /                     N   | /             |
     *             -----                     | -----         |
     *              i=1                      \  i=1          /
     *   var(X) = --------------------------------------------------
     *                               N - df_sub
     *
     * Where K may be an arbitrary value. Typically it is recommended to
     * not equal 0 to avoid the catastrophic cancellation (both terms may be
     * of very similar). Typically it can be assigned any element's value,
     * ideally close to the sample's mean.
     *
     * For more details, see:
     * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
     */

    const size_t N = x.size();

    if ( 0 == N )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    if ( df_sub >= N )
    {
        throw math::StatisticsException(math::StatisticsException::DF_SUBTRAHEND_TOO_LARGE);
    }

    // Let K be equal to the first element:
    const T K = x.at(0);

    T sum  = math::NumericUtil<T>::ZERO;
    T sum2 = math::NumericUtil<T>::ZERO;

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    // Ideal number of threads
    // (if each one processes approx. OMP_CHUNKS_PER_THREAD items):
    const size_t ideal = N / OMP_CHUNKS_PER_THREAD +
                         ( 0 == N % OMP_CHUNKS_PER_THREAD ? 0 : 1 );

    #pragma omp parallel num_threads(ideal) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x) \
                reduction(+ : sum, sum2)
    {
    	// Depending on the number of available threads,
        // determine the ideal nr.of samples per thread,
        // and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (N + nthreads -1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;
        const size_t iend = istart + samples_per_thread;

        // Calculate both sums of the assigned block...
        T partsum  = math::NumericUtil<T>::ZERO;
        T partsum2 = math::NumericUtil<T>::ZERO;
        for ( size_t i=istart; i<iend && i<N; ++i )
        {
            const T diff = x.at(i) - K;
            partsum  += diff;
            partsum2 += diff * diff;
        }

        // ... and add them to the total sums in a thread safe manner.
        sum  += partsum;
        sum2 += partsum2;
    }

    return (sum2 - (sum*sum)/static_cast<T>(N)) / static_cast<T>(N - df_sub);

    // In serial mode this variable is never used.
    (void) ideal;
}


/**
 * Variance of the sample.
 * Calculates either variance of a sample (sum of square deviations from
 * the mean is divided by N-1) or of a population (sum of squares divided by N).
 *
 * @param x - vector of sample elements
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return variance of the sample depending on 'sample'
 *
 * @throw StatisticsException if the sample is empty or too small
 */
template <class T>
T math::SampleStatGeneric<T>::var(const std::vector<T>& x, bool sample) throw(math::StatisticsException)
{
    return var( x, static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Standard deviation of the sample.
 * Calculates either standard deviation of a sample (sum of square deviations from
 * the mean is divided by N-1) or of a population (sum of squares divided by N).
 *
 * @param x - vector of sample elements
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return standard deviation of the sample depending on 'sample'
 *
 * @throw StatisticsException if the sample is empty or too small
 */
template <class T>
T math::SampleStatGeneric<T>::stdev(const std::vector<T>& x, bool sample) throw(math::StatisticsException)
{
    return stdev( x, static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Standard deviation of the sample.
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * the sample size.
 *
 * @param x - vector of sample elements
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return standard deviation of the sample, depending on the given 'df_sub'
 *
 * @throw StatisticsException if 'x' is empty or 'df_sub' exceeds sample's size
 */
template <class T>
T math::SampleStatGeneric<T>::stdev(const std::vector<T>& x, size_t df_sub) throw(math::StatisticsException)
{
    /*
     * Standard deviation is calculated as square root
     * of the variance.
     */

    /*
     *  This operation is only supported for T=float, T=double or T=long double.
     *  Specialized implementation is provided below for these three types.
     *
     *  For any other type, the operation is (probably) not supported
     *  as sqrt may not be defined for the type. In such a case throw an
     *  exception immediately.
     */

    throw math::StatisticsException(math::StatisticsException::UNSUPPORTED_TYPE);

    // will never execute, but some compilers may produce a warning if nothing is returned
    return math::NumericUtil<T>::ZERO;
}



/*
 * Specialization of stdev() for float, double and long double.
 * All three specializations are very similar and only differ in types of the
 * returned  value.
 * For easier maintainability, the specialization will be implemented
 * only once using a parameterized #define.
 */

#define _MATH_SAMPLESTATGENERIC_SPECIALIZED_STDEV(FD) \
template<> \
FD math::SampleStatGeneric<FD>::stdev(const std::vector<FD>& x, size_t df_sub) throw (math::StatisticsException) \
{ \
    return std::sqrt( var(x, df_sub) ); \
}
// end of #define

// the actual specialization for float:
_MATH_SAMPLESTATGENERIC_SPECIALIZED_STDEV(float)

// for double:
_MATH_SAMPLESTATGENERIC_SPECIALIZED_STDEV(double)

// and for long double:
_MATH_SAMPLESTATGENERIC_SPECIALIZED_STDEV(long double)

// definition of _MATH_QUATERNIONGENERIC_SPECIALIZED_NORM not needed anymore, #undef it
#undef _MATH_SAMPLESTATGENERIC_SPECIALIZED_STDEV
