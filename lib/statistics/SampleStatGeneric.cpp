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
 * Implementation of functions in the namespace SampleStat that calculate sample's
 * sum, mean, variance, standard deviation, covariance, Pearson's r
 * (correlation), r squared, etc.
 */


// no #include "SampleStatGeneric.hpp" !!!
#include <cstddef>
#include <cmath>
#include <vector>
#include <algorithm>

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"
#include "int_util/IntUtilGeneric.hpp"
#include "int_util/IntExponentiatorGeneric.hpp"
#include "util/NumericUtil.hpp"

#include "exception/StatisticsException.hpp"


// Implementation of "private" functions
namespace math
{

namespace SampleStat
{

namespace __private
{

/*
 * Finds an appropriate shift value (necessary to obtain variance or covariance)
 * among the first elements of 'x'.
 * 
 * @param x - vector of observations 
 * @param Nmax - the highest number of observations to check (default: 5)
 * 
 * @return observation with the highest absolute value among the first 'Nmax' elements of 'x'
 */
template <typename F>
F __getShift(const std::vector<F>& x, const size_t Nmax = 5)
{
    F retVal = x.at(0);
    F absRetVal = ( retVal<static_cast<F>(0) ? -retVal : retVal );
    const size_t N = std::min<size_t>(x.size(), Nmax);
    size_t cntr = 1;
    typename std::vector<F>::const_iterator it = x.begin() + cntr;
    
    for ( ; cntr<N; ++it, ++cntr )
    {
        const F el = *it;
        const F absx = ( el<static_cast<F>(0) ? -el : el );

        if ( absx > absRetVal )
        {
            retVal = el;
            absRetVal = absx;
        }
    }

    return retVal;
}


/*
 * Finds sample's either minimum or maximum observation, depending
 * on 'min'.
 *
 * @param x - vector of observations
 * @param min - a logical value indicating whether minimum or maximum value should be returned
 *
 * @return
 *
 * @throw StatisticsException if 'x' is empty
 */
template <typename F>
F __minmax(const std::vector<F>& x, const bool min) throw(math::StatisticsException)
{
    const size_t N = x.size();

    // sanity check
    if ( 0 == N )
    {
        throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    // the first observation is the first candidate for the extreme value...
    F retVal = x.at(0);

    // Coarse grained parallelism:
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x, retVal, min)
    {
    	const size_t thnr = omp_get_thread_num();
        const size_t nthreads  = omp_get_num_threads();
        const size_t elems_per_thread = (N + nthreads - 1) / nthreads;
        const size_t istart = elems_per_thread * thnr;

        typename std::vector<F>::const_iterator it = x.begin() + istart;
        // the first value of the block is the first candidate for the local extreme
        F temp = *it;

        for ( size_t cntr=0; 
              cntr<elems_per_thread && it!=x.end(); 
              ++it, ++cntr)
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
}

}  // namespace __private

}  // namespace SampleStat

}  // namespace math


/**
 * @param x - vector of observations
 *
 * @return minimum observation of the sample
 *
 * @throw StatisticsException if 'x' is empty
 */
template <typename F>
F math::SampleStat::min(const std::vector<F>& x) throw(math::StatisticsException)
{
    return math::SampleStat::__private::__minmax<F>(x, true);
}


/**
 * @param x - vector of observations
 *
 * @return maximum observation of the sample
 *
 * @throw StatisticsException if 'x' is empty
 */
template <typename F>
F math::SampleStat::max(const std::vector<F>& x) throw(math::StatisticsException)
{
    return math::SampleStat::__private::__minmax<F>(x, false);
}


/**
 * @param x - vector of observations
 *
 * @return sum of all observations
 */
template <typename F>
F math::SampleStat::sum(const std::vector<F>& x)
{
    const size_t N = x.size();

    if ( 0 == N )
    {
        return static_cast<F>(0);
    }

    /*
     * Note: the function 'moment' (for n=1 and about=0) could
     * be used to obtain the sum. However, in a typical application
     * 'sum' and 'mean' are expected to be called much more frequently
     * than 'moment', hence an optimized "specialization" is implemented
     * in this function.
     */

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    /*
     * Each thread calculates the sum of its block
     */
    F sum = static_cast<F>(0);
    #pragma omp parallel num_threads(ompIdeal(N)) \
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
        F partsum = static_cast<F>(0);
        typename std::vector<F>::const_iterator it = x.begin() + istart;
        for ( size_t cntr = 0;
              cntr<samples_per_thread && it!=x.end(); 
              ++it, ++cntr )
        {
            partsum += *it;
        }

        // ... and add it to the total sum in a thread safe manner.
        sum += partsum;
    }  // omp parallel

    return sum;
}


/**
 * Arithmetical mean (or average) of the sample.
 *
 * @param x - vector of observations
 *
 * @return mean value of the sample
 *
 * @throw StatisticsException if 'x' is empty
 */
template <typename F>
F math::SampleStat::mean(const std::vector<F>& x) throw(math::StatisticsException)
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

    return math::SampleStat::sum<F>(x) / static_cast<F>(N);
}


/**
 * Variance of the sample.
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * the sample size.
 *
 * @param x - vector of observations
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return variance of the sample, depending on the given 'df_sub'
 *
 * @throw StatisticsException if 'x' is empty or if 'df_sub' exceeds sample's size
 */
template <typename F>
F math::SampleStat::var(const std::vector<F>& x, const size_t df_sub) throw(math::StatisticsException)
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
     * of very similar values). Typically it can be assigned any observation's
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

    // Let K be equal to the first observation:
    const F K = math::SampleStat::__private::__getShift<F>(x);

    F sum  = static_cast<F>(0);
    F sum2 = static_cast<F>(0);

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    #pragma omp parallel num_threads(ompIdeal(N)) \
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
        F partsum  = static_cast<F>(0);
        F partsum2 = static_cast<F>(0);
        typename std::vector<F>::const_iterator it = x.begin() + istart;
        for ( size_t cntr = 0;
              cntr<samples_per_thread && it!=x.end(); 
              ++it, ++cntr )
        {
            const F diff = *it - K;
            partsum  += diff;
            partsum2 += diff * diff;
        }

        // ... and add them to the total sums in a thread safe manner.
        sum  += partsum;
        sum2 += partsum2;
    }

    return (sum2 - (sum*sum)/static_cast<F>(N)) / static_cast<F>(N - df_sub);
}


/**
 * Variance of the sample.
 * Calculates either variance of a sample (sum of square deviations from
 * the mean is divided by N-1) or of a population (sum of squares divided by N).
 *
 * @param x - vector of observations
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return variance of the sample depending on 'sample'
 *
 * @throw StatisticsException if the sample is empty or too small
 */
template <typename F>
F math::SampleStat::var(const std::vector<F>& x, const bool sample) throw(math::StatisticsException)
{
    return math::SampleStat::var<F>( x, static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Standard deviation of the sample.
 * Calculates either standard deviation of a sample (sum of square deviations from
 * the mean is divided by N-1) or of a population (sum of squares divided by N).
 *
 * @param x - vector of observations
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return standard deviation of the sample depending on 'sample'
 *
 * @throw StatisticsException if the sample is empty or too small
 */
template <typename F>
F math::SampleStat::stdev(const std::vector<F>& x, const bool sample) throw(math::StatisticsException)
{
    return math::SampleStat::stdev<F>( x, static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Standard deviation of the sample.
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * the sample size.
 *
 * @param x - vector of observations
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return standard deviation of the sample, depending on the given 'df_sub'
 *
 * @throw StatisticsException if 'x' is empty or 'df_sub' exceeds sample's size
 */
template <typename F>
F math::SampleStat::stdev(const std::vector<F>& x, const size_t df_sub) throw(math::StatisticsException)
{
    /*
     * Standard deviation is calculated as square root
     * of the variance.
     */

    return std::sqrt( math::SampleStat::var<F>(x, df_sub) );
}


/**
 * Covariance of two equally sized samples.
 *
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * a single sample's size.
 *
 * @param x1 - first vector of observations
 * @param x2 - second vector of observations
 * @param df_sub - generalized Bessel's correction value (typically 1 or 0)
 *
 * @return covariance of both samples, depending on the given 'df_sub'
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or 'df_sub' exceeds single sample's size
 */
template <typename F>
F math::SampleStat::cov(const std::vector<F>& x1, const std::vector<F>& x2, const size_t df_sub) throw(math::StatisticsException)
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
     * be of very similar values). Typically they can be assigned any observation's 
     * value, ideally close to the samples' means.
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
    const F K1 = math::SampleStat::__private::__getShift<F>(x1);
    const F K2 = math::SampleStat::__private::__getShift<F>(x2);

    F sum  = static_cast<F>(0);
    F sum1 = static_cast<F>(0);
    F sum2 = static_cast<F>(0);

    /*
     * Coarse grained parallelism will be applied, i.e. each thread will be
     * assigned an (approximately) equally sized contiguous block of data
     * to be processed.
     */

    #pragma omp parallel num_threads(ompIdeal(N1)) \
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
        F partsum  = static_cast<F>(0);
        F partsum1 = static_cast<F>(0);
        F partsum2 = static_cast<F>(0);
        typename std::vector<F>::const_iterator it1 = x1.begin() + istart;
        typename std::vector<F>::const_iterator it2 = x2.begin() + istart;
        for ( size_t cntr = 0;
              cntr<samples_per_thread && it1!=x1.end(); ++it1, ++it2, ++cntr )
        {
            const F d1 = *it1 - K1;
            const F d2 = *it2 - K2;
            partsum  += d1 * d2;
            partsum1 += d1;
            partsum2 += d2;
        }

        // ... and add them to the total sums in a thread safe manner.
        sum  += partsum;
        sum1 += partsum1;
        sum2 += partsum2;
    }

    return (sum - (sum1*sum2)/static_cast<F>(N1)) / static_cast<F>(N1 - df_sub);
}


/**
 * Covariance of two equally sized samples.
 *
 * The method allows the divisor (sample size) to be reduced by an arbitrary
 * positive integer number 'df_sub' as long as it is strictly smaller than
 * a single sample's size.
 *
 * @param x1 - first vector of observations
 * @param x2 - second vector of observations
 * @param sample - if 'true', sum of squared deviations is divided by (N-1), otherwise by N
 *
 * @return covariance of both samples, depending on 'sample'
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or if they are too small
 */
template <typename F>
F math::SampleStat::cov(const std::vector<F>& x1, const std::vector<F>& x2, const bool sample) throw(math::StatisticsException)
{
    return math::SampleStat::cov<F>( x1, x2, static_cast<size_t>( (false==sample ? 0 : 1) ) );
}


/**
 * Correlation a.k.a. "Pearson's r" of two equally sized samples.
 *
 * It equals covariance, divided by the product of
 * both samples' standard deviations.
 *
 * @param x1 - first vector of observations
 * @param x2 - second vector of observations
 *
 * @return correlation of both samples, always between -1 and 1
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or if they are too small
 */
template <typename F>
F math::SampleStat::cor(const std::vector<F>& x1, const std::vector<F>& x2) throw(math::StatisticsException)
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
    // allow smaller samples (at least one observation each):
    return math::SampleStat::cov<F>(x1, x2, false) /
           ( math::SampleStat::stdev<F>(x1, false) * math::SampleStat::stdev<F>(x2, false) );
}


/**
 * Coefficient of determination a.k.a. "r squared".
 *
 *     r^2 = cor(X1, X2)^2
 *
 * @param x1 - first vector of observations
 * @param x2 - second vector of observations
 *
 * @return samples' coefficient of determination
 *
 * @throw StatisticsException if any vector is empty, if they are not of equal sizes or if they are too small
 */
template <typename F>
F math::SampleStat::r2(const std::vector<F>& x1, const std::vector<F>& x2) throw(math::StatisticsException)
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

    const F cv = math::SampleStat::cov<F>(x1, x2, false);

    return (cv * cv) / 
           ( math::SampleStat::var<F>(x1, false) * math::SampleStat::var<F>(x2, false) );
}


/**
 * Returns the sample's n-th moment about the given value (if provided).
 *
 * @param x - vector of observations
 * @param n - order of the moment
 * @param about - center of the moment (default: 0)
 *
 * @return sample's n-th central moment about the value 'about'
 *
 * @throw StatisticsException if 'n' is negative or 'x' is an empty vector
 */
template <typename F, typename I>
F math::SampleStat::moment(const std::vector<F>& x, const I& n, const F& about) throw(math::StatisticsException)
{
    /*
     * A moment about an arbitrary value is defined as:
     *
     *              N-1
     *             -----
     *          1  \     /              \n
     *   m_n = ---  >    | x[i] - about |
     *          N  /     \              /
     *             -----
     *              i=0
     *
     * More details about the mathematical moment:
     * https://en.wikipedia.org/wiki/Moment_%28mathematics%29
     */

    const size_t NN = x.size();

    // 'n' must be non-negative
    if ( true == math::IntUtil::isNegative<I>(n) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    // 'x' must not be empty
    if ( NN <= 0 )
    {
    	throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    /*
     * Handle a special situation:
     */

    /*
     * If n==0, the expression above turns into a sum of ones (a^0=1)
     */
    if ( 0 == n )
    {
        return static_cast<F>(NN);
    }


    /*
     * If n>=1, just apply the expression above
     */

    F sum = static_cast<F>(0);

    // Coarse grained parallelism
    #pragma omp parallel num_threads(ompIdeal(NN)) \
                if(NN>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x, n, about) \
                reduction(+ : sum)
    {
        // Depending on the number of available threads,
        // determine the ideal nr. of samples per thread,
        // and start and sample of a block that each thread will process.
        const size_t thrnr = omp_get_thread_num();
        const size_t nthreads = omp_get_num_threads();

        const size_t samples_per_thread = (NN + nthreads - 1) / nthreads;
        const size_t istart = samples_per_thread * thrnr;

        F partsum = static_cast<F>(0);

        typename std::vector<F>::const_iterator it = x.begin() + istart;
        for ( size_t cntr = 0;
              cntr<samples_per_thread && it!=x.end(); ++it,  ++cntr )
        {
            /*
             * TODO when does it make sense to calculate powers
             * "manually" and when using the provided function?
             */
            partsum += math::IntExponentiator::power((*it - about), n);
        }

        sum += partsum;
    }  // pragma omp parallel

    return sum / static_cast<F>(NN);
}


/**
 * Returns the sample's n-th central moment about the mean.
 *
 * @param x - vector of observations
 * @param n - order of the moment
 *
 * @return sample's n-th central moment about the mean
 *
 * @throw StatisticsException if 'n' is negative or 'x' is an empty vector
 */
template <typename F, typename I>
F math::SampleStat::centralMoment(const std::vector<F>& x, const I& n) throw(math::StatisticsException)
{
    /*
     * Central moment about the mean is defined as:
     *
     *              N-1
     *             -----
     *          1  \     /        _ \n
     *   m_n = ---  >    | x[i] - x |
     *          N  /     \          /
     *             -----
     *              i=0
     *
     * More details about the central moment:
     * https://en.wikipedia.org/wiki/Central_moment
     * http://mathworld.wolfram.com/CentralMoment.html
     */

    const size_t NN = x.size();

    // 'n' must be non-negative
    if ( true == math::IntUtil::isNegative<I>(n) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    // 'x' must not be empty
    if ( NN <= 0 )
    {
    	throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
    }

    /*
     * Handle some special situations:
     */

    /*
     * If n==0, the expression above turns into a sum of ones
     * (a^0=1), their average equals 1:
     */
    if ( 0 == n )
    {
        return static_cast<F>(1);
    }

    /*
     * If n==1, the sum in the expression above will equal 0
     * (this is actually a definition of a mean value).
     */
    if ( 1 == n )
    {
        return static_cast<F>(0);
    }

    /*
     * If n==2, the expression above is actually a definition of
     * variance, without Bessel's correction:
     */
    if ( 2 == n )
    {
        return math::SampleStat::var<F>(x, false);
    }

    /*
     * If n>=3, just apply the expression above
     */

    const F xbar = math::SampleStat::mean<F>(x);

    return math::SampleStat::moment<F>(x, n, xbar);
}


/**
 * Sample skewness
 *
 * @param x - vector of observations
 *
 * @return skewness of observations in 'x'
 *
 * @throw StatisticsExcpetion if 'x' is empty or all observations are equal
 */
template <typename F>
F math::SampleStat::skewness(const std::vector<F>& x) throw(math::StatisticsException)
{
    /*
     * Sample skewness is defined as:
     *
     *           m3
     *    b1 = ------
     *          sd^3
     *
     * where 'm3' denotes the sample's 3rd central moment and
     * and 'sd' denotes the sample'sstandard deviation, with Bessel's correction.
     *
     * More details:
     * https://en.wikipedia.org/wiki/Skewness
     */

    const F m3 = math::SampleStat::centralMoment<F>(x, 3);
    const F sd = math::SampleStat::stdev<F>(x, true);

    const F sd3 = sd * sd * sd;

    /*
     * Prevent very unlikely division by zero
     * (when all observations in 'x' are equal)
     */
    if ( true == math::NumericUtil::isZero<F>(sd3) )
    {
        // By convention, skewness equals 0 in such cases
        return static_cast<F>(0);
    }

    return m3 / sd3;
}


/**
 * Sample excess kurtosis
 *
 * @param x - vector of observations
 *
 * @return kurtosis of observations in 'x'
 *
 * @throw StatisticsExcpetion if 'x' is empty or all observations are equal
 */
template <typename F>
F math::SampleStat::kurtosis(const std::vector<F>& x) throw(math::StatisticsException)
{
    /*
     * Sample's excess kurtosis is defined as:
     *
     *          m4
     *   g2 = ------ - 3
     *         m2^2
     *
     * where 'm4' denotes the 4th central moment and
     * 'm2' denotes the 2nd central moment a.k.a.
     * the sample's variance w/o Bessel's correction.
     *
     * More details:
     * https://en.wikipedia.org/wiki/Kurtosis
     */

    const F m4 = math::SampleStat::centralMoment<F>(x, 4);
    const F m2 = math::SampleStat::var<F>(x, false);
    const F m2_2 = m2 * m2;

    /*
     * Prevent very unlikely division by zero
     * (when all observations in 'x' are equal)
     */
    if ( true == math::NumericUtil::isZero<F>(m2_2) )
    {
        // By convention, kurtosis equals 0 in such cases
        return static_cast<F>(0);
    }

    return m4/m2_2 - static_cast<F>(3);
}
