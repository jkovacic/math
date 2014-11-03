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
 * @file SampleStatGeneric.cpp
 *
 * Implementation of the class SampleStaGeneric that calculates sample's
 * mean, variance, standard deviation, sum, etc.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 *
 * @author Jernej Kovacic
 */


// deliberately there is no #include "SampleStatGeneric.h" !
#include <cstddef>
#include <cmath>
#include <vector>

#include "omp_header.h"
#include "omp_settings.h"

#include "StatisticsException.h"
#include "NumericUtil.h"


/**
 * Constructor. Initializes all internal variables but does not
 * calculate anything yet.
 *
 * @param sample - a vector of samples
 *
 * @throw StatisticsException if 'sample' is an empty vector
 */
template <class T>
math::SampleStatGeneric<T>::SampleStatGeneric(const std::vector<T>& sample) throw(math::StatisticsException)
{
    this->m_processed = false;
    this->m_N = 0;
    this->m_sum = math::NumericUtil<T>::ZERO;
    this->m_mean = math::NumericUtil<T>::ZERO;
    this->m_sumSqDev = math::NumericUtil<T>::ZERO;
    this->m_pData = NULL;

    if ( 0 == sample.size() )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    this->m_N = sample.size();
    this->m_pData = &sample;
}


/**
 * Calculates values of all necessary internal variables which become
 * available to be returned by getter methods.
 *
 * Nothing is done if process() has been called before on the same object.
 *
 * After this method is finished, the original vector of samples may be freely
 * modified as other methods do not need it anymore.
 */
template <class T>
void math::SampleStatGeneric<T>::process()
{
    // Check whether the sample has already been processed.
    if ( true == this->m_processed )
    {
        return;
    }

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    // Ideal number of threads
    // (if each one processes approx. OMP_CHUNKS_PER_THREAD items):
    const size_t ideal = this->m_N / OMP_CHUNKS_PER_THREAD +
                         ( 0 == this->m_N % OMP_CHUNKS_PER_THREAD ? 0 : 1 );

    /*
     * In the first step, each thread calculates the sum of its block
     */
    T tempSum = NumericUtil<T>::ZERO;
    #pragma omp parallel num_threads(ideal) \
                    if(this->m_N>OMP_CHUNKS_PER_THREAD) \
                    default(none) \
                    reduction(+ : tempSum)
    {
        // Depending on the number of available threads,
    	// determine the ideal nr.of samples per thread,
    	// and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (this->m_N + nthreads -1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;
        const size_t iend = istart + samples_per_thread;

        // Calculate the sum of the assigned block...
        T partsum = math::NumericUtil<T>::ZERO;
        for ( size_t i=istart; i<iend && i<this->m_N; ++i )
        {
            partsum += this->m_pData->at(i);
        }

        // ... and add it to the total sum in a thread safe manner.
        tempSum += partsum;
    }

    // When the sum is known, it is trivial to obtain the sample's mean:
    this->m_sum = tempSum;
    this->m_mean = this->m_sum / this->m_N;

    // Fork threads once more, this time calculate the sum of
    // squared deviations from the mean.
    // Apply exactly the same approach as above.
    tempSum = NumericUtil<T>::ZERO;
    #pragma omp parallel num_threads(ideal) \
                 if(this->m_N>OMP_CHUNKS_PER_THREAD) \
                 default(none) \
                 reduction(+ : tempSum)
    {
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (this->m_N + nthreads -1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;
        const size_t iend = istart + samples_per_thread;

        T partsum = math::NumericUtil<T>::ZERO;
        for ( size_t i=istart; i<iend && i<this->m_N; ++i )
        {
            T diff = this->m_pData->at(i) - this->m_mean;
            partsum += diff * diff;
        }

        tempSum += partsum;

    }

    // The sum of squared deviations will be used by other functions to
    // calculate the sample's variance and standard deviation.
    this->m_sumSqDev = tempSum;

    // Flag the sample as processed
    this->m_processed = true;

    // In serial mode this variable is never used.
    (void) ideal;
}


/**
 * @return a logical value indicating whether the samples have already been processed
 */
template <class T>
bool math::SampleStatGeneric<T>::processed() const
{
    return this->m_processed;
}


/**
 * @return number of samples or 0 if process() has not been called yet
 */
template <class T>
size_t math::SampleStatGeneric<T>::sampleSize() const
{
    return this->m_N;
}


/**
 * @return sum of all sample values
 *
 * @throw StatisticsException if process() has not been called yet
 */
template <class T>
T math::SampleStatGeneric<T>::sum() const throw(math::StatisticsException)
{
    if ( false == this->m_processed )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_NOT_PROCESSED_YET);
    }

    return this->m_sum;
}


/**
 * Arithmetical mean (or average) of the sample.
 *
 * @return mean value of the sample
 *
 * @throw StatisticsException if process() has not been called yet
 */
template <class T>
T math::SampleStatGeneric<T>::mean() const throw(math::StatisticsException)
{
    /*
     * Arithmetical mean is calculated as:
     *
     *                   N
     *                 -----
     *             1   \
     * mean(X) = -----  >  X(i)
     *             N   /
     *                 -----
     *                  i=1
     *
     */

    if ( false == this->m_processed )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_NOT_PROCESSED_YET);
    }

    return this->m_mean;
}


/**
 * Variance of the sample.
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * the sample size.
 *
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return variance of the sample, depending on the given 'df_sub'
 *
 * @throw StatisticsException if process() has not been called yet or if 'df_sub' exceeds sample's size
 */
template <class T>
T math::SampleStatGeneric<T>::var(size_t df_sub) const throw(math::StatisticsException)
{
    /*
     * Variance of the sample is calculated as:
     *
     *                              N
     *                            -----
     *                  1         \                 2
     *   var(X) = ------------     >   (X(i) - mean)
     *             N - df_sub     /
     *                            -----
     *                             i=1
     */

    if ( false == this->m_processed )
	{
        throw math::StatisticsException(math::StatisticsException::SAMPLE_NOT_PROCESSED_YET);
    }

    if ( df_sub >= this->m_N )
    {
        throw math::StatisticsException(math::StatisticsException::DF_SUBTRAHEND_TOO_LARGE);
    }

    return this->m_sumSqDev / static_cast<T>(this->m_N - df_sub);
}


/**
 * Variance of the sample.
 * Calculates either variance of a sample (sum of square deviations from
 * the mean is divided by N-1) or of a population (sum of squares divided by N).
 *
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return variance of the sample depending on 'sample'
 *
 * @throw StatisticsException if process() has not been called yet or if the sample is too small
 */
template <class T>
T math::SampleStatGeneric<T>::var(bool sample) const throw(math::StatisticsException)
{
    /*
     * Variance of the sample, calculated as:
     *
     *                         N
     *                       -----
     *                1      \                 2
     *   var(X) = --------    >   (X(i) - mean)
     *             N - BC    /
     *                       -----
     *                        i=1
     *
     * where BC (Bessel's correction) equals 1 if 'sample' is 'true', otherwise it equals 0.
     */

    if ( false == this->m_processed )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_NOT_PROCESSED_YET);
    }

    if ( true==sample && 1==this->m_N )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_TOO_SMALL);
    }

    return this->var( static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Standard deviation of the sample.
 * Calculates either standard deviation of a sample (sum of square deviations from
 * the mean is divided by N-1) or of a population (sum of squares divided by N).
 *
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return standard deviation of the sample depending on 'sample'
 *
 * @throw StatisticsException if process() has not been called yet or if the sample is too small
 */
template <class T>
T math::SampleStatGeneric<T>::stdev(bool sample) const throw(math::StatisticsException)
{
    /*
     * Standard deviation is calculated as:
     *
     *                            ---------------------------------
     *                           /            N                   |
     *                          /           -----
     *              ---        /    1       \                 2
     *   stdev(X) =    \      /  --------    >   (X(i) - mean)
     *                  \    /    N - BC    /
     *                   \  /               -----
     *                    \/                 i=1
     *
     * where BC (Bessel's correction) equals 1 if 'sample' is 'true', otherwise it equals 0.
     */

    if ( true==sample && 1==this->m_N )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_TOO_SMALL);
    }

    return this->stdev( static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Standard deviation of the sample.
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * the sample size.
 *
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return standard deviation of the sample, depending on the given 'df_sub'
 *
 * @throw StatisticsException if process() has not been called yet or if 'df_sub' exceeds sample's size
 */
template <class T>
T math::SampleStatGeneric<T>::stdev(size_t df_sub) const throw(math::StatisticsException)
{
    /*
     * (Generalized) standard deviation is calculated as:
     *
     *                            -------------------------------------
     *                           /               N                    |
     *                          /               -----
     *              ---        /     1          \                 2
     *   stdev(X) =    \      / ------------     >   (X(i) - mean)
     *                  \    /   N - df_sub     /
     *                   \  /                   -----
     *                    \/                     i=1
     *
     */

    /*
     *  This operation is only supported for T=float or T=double.
     *  Specialized implementation is provided below for these two types.
     *
     *  For any other type, the operation is (probably) not supported
     *  as sqrt may not be defined for the type. In such a case throw an
     *  exception immediately.
     */

    throw math::StatisticsException(math::StatisticsException::UNSUPPORTED_TYPE);

    // will never execute, but some compilers may produce a warning if nothing is returned
    return math::NumericUtil<T>::ONE;
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
FD math::SampleStatGeneric<FD>::stdev(size_t df_sub) const throw (math::StatisticsException) \
{ \
    return std::sqrt( this->var(df_sub) ); \
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
