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
#include <algorithm>

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"

#include "exception/StatisticsException.hpp"
#include "util/NumericUtil.hpp"


/*
 * Finds an appropriate shift value (necessary to obtain variance or covariance)
 * among the first elements of 'x'.
 * 
 * @param x - vector of samples 
 * @param Nmax - the highest number of elements to check (default: 5)
 * 
 * @return element with the highest absolute value among the first 'Nmax' elements of 'x'
 */
template <class T>
T math::SampleStatGeneric<T>::__getShift(const std::vector<T>& x, size_t Nmax)
{
    T retVal = x.at(0);
    T absRetVal = ( retVal<math::NumericUtil<T>::ZERO ? -retVal : retVal );
    const size_t N = std::min<size_t>(x.size(), Nmax);

    size_t cntr = 1;
    for ( typename std::vector<T>::const_iterator it = x.begin()+cntr;
          cntr<N; ++it, ++cntr)
    {
        const T el = *it;
        const T absx = ( el<math::NumericUtil<T>::ZERO ? -el : el );

        if ( absx > absRetVal )
        {
            retVal = el;
            absRetVal = absx;
        }
    }

    return retVal;
}


/*
 * Finds sample's either minimum or maximum value, depending
 * on 'min'.
 *
 * @param x - vector of sample elements
 * @param min - a logical value indicating whether minimum or maximum value should be returned
 *
 * @return
 *
 * @throw StatisticsException if 'x' is empty
 */
template <class T>
T math::SampleStatGeneric<T>::__minmax(const std::vector<T>& x, bool min) throw(math::StatisticsException)
{
    const size_t N = x.size();

    // sanity check
    if ( 0 == N )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    // the first element is the first candidate for the extreme value...
    T retVal = x.at(0);

    // Coarse grained parallelism:
    const size_t ideal = N / OMP_CHUNKS_PER_THREAD +
                 ( 0 == N % OMP_CHUNKS_PER_THREAD ? 0 : 1 );

    #pragma omp parallel num_threads(ideal) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x, retVal, min)
    {
    	const size_t thnr = omp_get_thread_num();
        const size_t nthreads  = omp_get_num_threads();
        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thnr;

        typename std::vector<T>::const_iterator it = x.begin() + istart;
        // the first value of the block is the first candidate for the local extreme
        T temp = *it;

        for ( size_t cntr=0; cntr<elems_per_thread && it!=x.end(); ++it, ++cntr)
        {
            // update 'temp' depending on 'min'
            temp = ( true==min ? std::min(temp, *it) : std::max(temp, *it) );
        }

        // prevent possible race condition when updating retVal
        #pragma omp critical(samplestatgeneric_minmax)
        {
            retVal = ( true==min ? std::min(retVal, temp) : std::max(retVal, temp) );
        }
    }  // omp parallel

    return retVal;

    // this variable is never used in serial mode
    (void) ideal;
}


/**
 * @param x - vector of sample elements
 *
 * @return minimum value of the sample
 *
 * @throw StatisticsException if 'x' is empty
 */
template <class T>
T math::SampleStatGeneric<T>::min(const std::vector<T>& x) throw(math::StatisticsException)
{
    return __minmax(x, true);
}


/**
 * @param x - vector of sample elements
 *
 * @return maximum value of the sample
 *
 * @throw StatisticsException if 'x' is empty
 */
template <class T>
T math::SampleStatGeneric<T>::max(const std::vector<T>& x) throw(math::StatisticsException)
{
    return __minmax(x, false);
}


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
     * Each thread calculates the sum of its block
     */
    T sum = math::NumericUtil<T>::ZERO;
    #pragma omp parallel num_threads(ideal) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(x) \
                    reduction(+ : sum)
    {
        // Depending on the number of available threads,
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;

        // Calculate the sum of the assigned block...
        T partsum = math::NumericUtil<T>::ZERO;
        size_t cntr = 0;
        for ( typename std::vector<T>::const_iterator it = x.begin() + istart;
              cntr<samples_per_thread && it!=x.end(); ++it, ++cntr )
        {
            partsum += *it;
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
     * not equal 0 to avoid the catastrophic cancellation (both terms can be
     * of very similar values). Typically it can be assigned any element's
     * value, ideally close to the sample's mean.
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
    const T K = __getShift(x);

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
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;

        // Calculate both sums of the assigned block...
        T partsum  = math::NumericUtil<T>::ZERO;
        T partsum2 = math::NumericUtil<T>::ZERO;
        size_t cntr = 0;
        for ( typename std::vector<T>::const_iterator it = x.begin() + istart; 
              cntr<samples_per_thread && it!=x.end(); ++it, ++cntr )
        {
            const T diff = *it - K;
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



/**
 * Covariance of two equally sized samples.
 *
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * a single sample's size.
 *
 * @param x1 - first vector of sample elements
 * @param x2 - second vector of sample elements
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return covariance of both samples, depending on the given 'df_sub'
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or 'df_sub' exceeds single sample's size
 */
template <class T>
T math::SampleStatGeneric<T>::cov(const std::vector<T>& x1, const std::vector<T>& x2, size_t df_sub) throw(math::StatisticsException)
{
    /*
     * Covariance of two equally sized samples (X1 and X2) can be
     * calculated as:
     *
     *                              N
     *                            -----
     *                     1      \          --           --
     * cov(X1, X2) = ------------  >  (X1[i]-X1) * (X2[i]-X2)
     *                N - df_sub  /
     *                            -----
     *                             i=1
     *
     * It is called a two step algorithm as it must calculate the sample means
     * in the first step, followed by the algorithm above to calculate the variance.
     *
     * When the means are not required, the expression above can be replaced by the
     * one step algorithm:
     *
     *                  N                               N             N
     *                -----                           -----         -----
     *                \                            1  \             \
     *                 >  (X1[i]-K1)*(X2[i]-K2) - ---  > (X1[i]-K1)  > (X2[i]-K2)
     *                /                            N  /             /
     *                -----                           -----         -----
     *                 i=1                             i=1           i=1
     * cov(X1, X2) = --------------------------------------------------------------
     *                                        N - df_sub
     *
     * Where K1 and K2 may be arbitrary values. Typically it is recommended to
     * not equal 0 to avoid the catastrophic cancellation (both terms may be
     * be of very similar values). Typically they can be assigned any sample
     * element's value, ideally close to the samples' means.
     *
     * For more details, see:
     * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
     */

    const size_t N1 = x1.size();
    const size_t N2 = x2.size();

    if ( 0==N1 || 0==N2 )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    if ( N1 != N2 )
    {
        throw math::StatisticsException(math::StatisticsException::UNEQUAL_SAMPLE_SIZES);
    }

    if ( df_sub >= N1 )
    {
        throw math::StatisticsException(math::StatisticsException::DF_SUBTRAHEND_TOO_LARGE);
    }

    // K's are equal to the first elements of both samples
    const T K1 = __getShift(x1);
    const T K2 = __getShift(x2);

    T sum  = math::NumericUtil<T>::ZERO;
    T sum1 = math::NumericUtil<T>::ZERO;
    T sum2 = math::NumericUtil<T>::ZERO;

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    // Ideal number of threads
    // (if each one processes approx. OMP_CHUNKS_PER_THREAD items):
    const size_t ideal = N1 / OMP_CHUNKS_PER_THREAD +
                         ( 0 == N1 % OMP_CHUNKS_PER_THREAD ? 0 : 1 );

    #pragma omp parallel num_threads(ideal) \
                if(N1>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x1, x2) \
                reduction(+ : sum, sum1, sum2)
    {
    	// Depending on the number of available threads,
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (N1 + nthreads - 1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;

        // Calculate both sums of the assigned block...
        T partsum  = math::NumericUtil<T>::ZERO;
        T partsum1 = math::NumericUtil<T>::ZERO;
        T partsum2 = math::NumericUtil<T>::ZERO;
        size_t cntr = 0;
        typename std::vector<T>::const_iterator it1;
        typename std::vector<T>::const_iterator it2;
        for ( it1 = x1.begin() + istart, it2 = x2.begin() + istart;
              cntr<samples_per_thread && it1!=x1.end(); ++it1, ++it2, ++cntr )
        {
            const T d1 = *it1 - K1;
            const T d2 = *it2 - K2;
            partsum  += d1 * d2;
            partsum1 += d1;
            partsum2 += d2;
        }

        // ... and add them to the total sums in a thread safe manner.
        sum  += partsum;
        sum1 += partsum1;
        sum2 += partsum2;
    }

    return (sum - (sum1*sum2)/static_cast<T>(N1)) / static_cast<T>(N1 - df_sub);

    // In serial mode this variable is never used.
    (void) ideal;
}


/**
 * Covariance of two equally sized samples.
 *
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * a single sample's size.
 *
 * @param x1 - first vector of sample elements
 * @param x2 - second vector of sample elements
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return covariance of both samples, depending on 'sample'
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or if they are too small
 */
template <class T>
T math::SampleStatGeneric<T>::cov(const std::vector<T>& x1, const std::vector<T>& x2, bool sample) throw(math::StatisticsException)
{
    return cov( x1, x2, static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Correlation a.k.a. "Pearson's r" of two equally sized samples.
 *
 * It equals covariance, divided by the product of
 * both samples' standard deviations.
 *
 * @param x1 - first vector of sample elements
 * @param x2 - second vector of sample elements
 *
 * @return correlation of both samples, always between -1 and 1
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or if they are too small
 */
template <class T>
T math::SampleStatGeneric<T>::cor(const std::vector<T>& x1, const std::vector<T>& x2) throw(math::StatisticsException)
{
    /*
     * Correlation can be calculated as:
     *
     *                             cov(X1, X2)
     *     r = cor(X1, X2) = -----------------------
     *                        stdev(X1) * stdev(X2)
     *
     * It turns out that if the same divisor (N-df_sub) is applied at all operations,
     * it cancels out, hence the correlation is independent on 'df_sub'.
     */

    // Apply the population covariance and standard deviations to
    // allow smaller samples (at least one element each):
    return cov(x1, x2, false) / ( stdev(x1, false) * stdev(x2, false) );
}


/**
 * Coefficient of determination a.k.a. "r squared".
 *
 *     r^2 = cor(X1, X2)^2
 *
 * @param x1 - first vector of sample elements
 * @param x2 - second vector of sample elements
 *
 * @return samples' coefficient of determination
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or if they are too small
 */
template <class T>
T math::SampleStatGeneric<T>::r2(const std::vector<T>& x1, const std::vector<T>& x2) throw(math::StatisticsException)
{
    /*
     * R squared can be calculated by squaring the samples' correlation:
     *
     *                                            2
     *          2              2       cov(X1, X2)
     * r(X1, X2)  = cor(X1, X2)  = -------------------
     *                              var(X1) * var(X2)
     *
     * It can easily be proven that divisors (N - df_sub) cancel out
     * if the same divisor is applied at each operation, hence
     * r2 does not depend on 'df_sub'.
     */

    // The second formula above supports more types T.
    // Additionally apply population covariance and variances
    // to allow smaller samples (at least one element each).

    const T cv = cov(x1, x2, false);

    return (cv * cv) / ( var(x1, false) * var(x2, false) );
}
