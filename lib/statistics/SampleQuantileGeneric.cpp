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
 * Implementation of the class SampleQuantileGeneric that estimates
 * quantiles of a sample.
 */


// no #include "SampleQuantileGeneric.hpp" !

#include <vector>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <set>

#include "util/mtcopy.hpp"
#include "exception/StatisticsException.hpp"
#include "int_util/IntUtilGeneric.hpp"
#include "util/NumericUtil.hpp"


/**
 * Constructor.
 * Creates its own copy of the vector of observations.
 *
 * @param sample - a vector of observations
 *
 * @throw StatisticsException if 'sample' is empty or allocation of memory for its copy failed
 */
template <typename F>
math::SampleQuantileGeneric<F>::SampleQuantileGeneric(const std::vector<F>& sample) throw(math::StatisticsException)
{
    try
    {
        if ( 0==sample.size() )
        {
            throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
        }

        // copy the 'sample' into an internal vector:
        math::mtcopy(sample, this->m_v);
        this->m_N = this->m_v.size();

        // and sort it in ascending order
        std::sort(this->m_v.begin(), this->m_v.end());
    }
    catch ( const std::bad_alloc& ex )
    {
        throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
    }
}


/**
 * @return size of the sample
 */
template <typename F>
size_t math::SampleQuantileGeneric<F>::sampleSize() const
{
    return this->m_N;
}


/**
 * Estimates a quantile for a probability expressed as a ratio
 * of two integers.
 *
 * Note that 'den' must be greater or equal as 2 and that 'num'
 * must be strictly greater than 0 and strictly less than 'den':
 *   den >= 2
 *   1 <= num <= den-1
 *
 * @param num - probability's numerator
 * @param den - probabilty's denominator
 * @param method - one of the supported methods to estimate the quantile (default: R7)
 *
 * @return quantile for the probabilty 'num'/'den'
 *
 * @throw StatisticsException if any input argument is invalid
 */
template <typename F> template <typename I>
F math::SampleQuantileGeneric<F>::quantile(
                    const I& num, 
                    const I& den, 
                    const math::EQntlType::type method ) 
                const throw(math::StatisticsException)
{
    if ( true == math::IntUtil::isNegative<I>(num) ||
         true == math::IntUtil::isNegative<I>(den) ||
         static_cast<I>(0) == num || 
         den < static_cast<I>(2) || num >= den )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    return this->qntl( static_cast<F>(num) / static_cast<F>(den), method);
}


/*
 * Internal (and private) method that estimates the quantile from the given value 'h'.
 * The function was introduced because several methods actually use the same formula
 * to estimate a quantile from 'h' (calculated differently, depending on a method).
 *
 * @param h - h value, calculated depending on a method
 *
 * @return estimated quantile as a function of 'h'
 */
template <typename F>
F math::SampleQuantileGeneric<F>::linIntrp(const F& h) const
{
    // "rename" the vector as it is referred in statistical publications:
    const std::vector<F>& x = this->m_v;
    const F hf = std::floor(h);

    /*
     * if floor(h) is referred as 'hf', the quantile is estimated as:
     *
     *   qp = x[hf] + (h - hf) * (x[hf+1] - x[hf])
     *
     * Note that the expression assumes one-based arrays.
     * In C++ all indices must be additionally decremented by 1
     * prior accessing elements of x.
     *
     * Hence:
     *
     *   qp = x[hf-1] + (h - hf) * (x[hf] - x[hf-1])
     */

    return ( x.at(math::NumericUtil::intRound<F, size_t>(hf) - 1) +
            (h - hf) * (x.at(math::NumericUtil::intRound<F, size_t>(hf)) -
             x.at(math::NumericUtil::intRound<F, size_t>(hf) - 1)) );
}


/**
 * Estimates a quantile for the given probability 'p'.
 *
 * See https://en.wikipedia.org/wiki/Quantile for more
 * details about supported methods.
 *
 * Note that 'p' must be greater than 0 and less than 1:
 *   0 <= p <= 1
 *
 * @param p - probabilty
 * @param method - one of the supported methods to estimate the quantile (default: R7)
 *
 * @return quantile for the probabilty 'p'
 *
 * @throw StatisticsException if any input argument is invalid
 */
template <typename F>
F math::SampleQuantileGeneric<F>::qntl(const F& p, const math::EQntlType::type method) const throw(math::StatisticsException)
{
    F retVal;

    // Use references to "rename" the variables into names,
    // typically used in statistical publications.
    const size_t& N = this->m_N;
    const std::vector<F>& x = this->m_v;

    // N casted to T, used often in nonint expressions
    const F NT = static_cast<F>(N);

    if ( p<static_cast<F>(0) || p>static_cast<F>(1) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    // Frequently used numerical constants:
    const F ONE = static_cast<F>(1);
    const F TWO = static_cast<F>(2);
    const F HALF = ONE / TWO;
    const F QUARTER = ONE / static_cast<F>(4);

    /*
     * Implementation of all methods described in:
     * https://en.wikipedia.org/wiki/Quantile
     * Note that the Wikipedia article assumes one-based arrays.
     * In a language with zero-based arrays (e.g. C++), the index must
     * be decremented by 1 before an array element is accessed.
     */

    switch (method)
    {
    case math::EQntlType::R1 :
    case math::EQntlType::SAS3 :
    {
        // Inverse of empirical distribution function
        if ( true == math::NumericUtil::isZero<F>(p) )
        {
            retVal = x.at(0);
        }
        else
        {
            const F h = NT * p + HALF;
            retVal = x.at(math::NumericUtil::intCeil<F, size_t>(h-HALF) - 1);
        }

        break;
    }

    case math::EQntlType::R2 :
    case math::EQntlType::SAS5 :
    {
        // Inverse of empirical distribution function
        // with averaging at discontinuities
        if ( true == math::NumericUtil::isZero<F>(p) )
        {
            retVal = x.at(0);
        }
        else if ( true == math::NumericUtil::isZero<F>(ONE-p) )
        {
            retVal = x.at(N - 1);
        }
        else
        {
            const F h = NT * p + HALF;
            retVal = ( x.at(math::NumericUtil::intCeil<F, size_t>(h-HALF) - 1) +
                       x.at(math::NumericUtil::intFloor<F, size_t>(h+HALF) - 1) ) * HALF;
        }

        break;
    }

    case math::EQntlType::R3 :
    case math::EQntlType::SAS2 :
    {
        // Observation closest to N*p
        if ( p < HALF/NT )
        {
            retVal = x.at(0);
        }
        else
        {
            const F h = NT * p;
            const F frac = h - std::floor(h);
            size_t idx = math::NumericUtil::intFloor<F, size_t>(h);

            /*
             * The index is h, rounded to the nearest integer.
             * In case of a tie, round to the even integer,
             * assuming one-based indexing.
             */

            if ( true == math::NumericUtil::isZero<F>(frac-HALF) )
            {
                // in this case round to the nearest even (divisible by 2) integer idx
                if ( 0 != idx%2 )
                {
                    ++idx;
                }
            }
            else if ( frac > HALF )
            {
                ++idx;
            }
            // no need to to do anything with 'idx' if frac<0.5

            retVal = x.at(idx-1);
        }

        break;
    }

    case math::EQntlType::R4 :
    case math::EQntlType::SAS1 :
    case math::EQntlType::SCIPY_0_1 :
    {
        // Linear interpolation of the empirical distribution function
        if ( p < static_cast<F>(1)/NT )
        {
            retVal = x.at(0);
        }
        else if ( true == math::NumericUtil::isZero<F>(ONE-p) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp(NT * p);
        }

        break;
    }

    case math::EQntlType::R5 :
    case math::EQntlType::SCIPY_05_05 :
    {
        // Piecewise linear function where the knots are the values midway
        // through the steps of the empirical distribution function
        if ( p < HALF/NT )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT-HALF)/NT )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp(NT * p + HALF);
        }

        break;
    }

    case math::EQntlType::R6 :
    case math::EQntlType::SAS4 :
    case math::EQntlType::SCIPY_0_0 :
    {
        // Linear interpolation of the expectations for the order statistics
        // for the uniform distribution on [0,1]
        if ( p < ONE/(NT+ONE) )
        {
            retVal = x.at(0);
        }
        else if ( p >= NT/(NT+ONE) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + ONE) * p);
        }

        break;
    }

    case math::EQntlType::R7 :
    case math::EQntlType::SCIPY_1_1 :
    case math::EQntlType::EXCEL :
    {
        // Linear interpolation of the modes for the order statistics for
        // the uniform distribution on [0,1]
        if ( true == math::NumericUtil::isZero<F>(ONE-p) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT - ONE) * p + ONE);
        }

        break;
    }

    case math::EQntlType::R8 :
    case math::EQntlType::SCIPY_13_13 :
    {
        const F third = ONE / static_cast<F>(3);

        // Linear interpolation of the approximate medians for order statistics
        if ( p < (TWO * third) / (NT+third) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT-third) / (NT+third) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + third) * p + third);
        }

        break;
    }

    case math::EQntlType::R9 :
    case math::EQntlType::SCIPY_38_38 :
    {
        const F F_3_8 = static_cast<F>(3) / static_cast<F>(8);
        const F F_5_8 = static_cast<F>(5) / static_cast<F>(8);

        // The resulting quantile estimates are approximately unbiased for the
        // expected order statistics if x is normally distributed
        if ( p < F_5_8/(NT+QUARTER) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT-F_3_8)/(NT+QUARTER) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + QUARTER) * p + F_3_8);
        }

        break;
    }

    case math::EQntlType::SCIPY_N05_N05 :
    {
        // If h were rounded, this would give the order statistic with the least
        // expected square deviation relative to p
        if ( p < (static_cast<F>(3) / TWO) / (NT + TWO) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT + HALF) / (NT + TWO) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + TWO) * p - HALF);
        }

        break;
    }

    default:
    	throw math::StatisticsException(math::StatisticsException::UNSUPPORTED_QUANTILE_METHOD);
    }  // switch

    return retVal;
}


/**
 * Median of the sample.
 *
 * If the number of observations is odd, the middle element
 * is returned:
 *
 *   median = sorted_vector[(N-1)/2]
 *
 * If the number of observations is even, the mean of the middle
 * two elements is returned:
 *
 *             sorted_vector[N/2 - 1] + sorted_vector[N/2]
 *   median = ---------------------------------------------
 *                                  2
 *
 *
 * @return median of the sample
 */
template <typename F>
F math::SampleQuantileGeneric<F>::median() const
{
    F retVal;
    const size_t& N = this->m_N;
    const std::vector<F>& x = this->m_v;

    if ( 0 == N % 2 )
    {
        // even number of observations
    	size_t h = N / 2;
        retVal = (x.at(h-1) + x.at(h)) / static_cast<F>(2);
    }
    else
    {
        // odd number of observations
        retVal = x.at((N-1) / 2);
    }

    return retVal;
}


/**
 * Estimation of a quartile of the sample.
 *
 * @param q - desired quartile (between 0 and 4 )
 * @param method - one of the supported methods to estimate the quartile (default: R7)
 *
 * @return estimation of the desired quartile
 *
 * @throw StatisticsException  if any input argument is invalid
 */
template <typename F>
F math::SampleQuantileGeneric<F>::quartile(
            const unsigned short int q, 
            const math::EQntlType::type method) 
        const throw(math::StatisticsException)
{
    F retVal;

    if ( 0 == q )
    {
        retVal = this->m_v.at(0);
    }
    else if ( 4 == q )
    {
        retVal = this->m_v.at(this->m_N - 1);
    }
    else if ( q>=1 && q<=3 )
    {
        retVal = this->quantile<short int>(q, 4, method);
    }
    else
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    return retVal;
}


/**
 * @param method - one of the supported methods to estimate both quartiles (default: R7)
 *
 * @return sample's interquartile range (difference between the 3rd and the 1st quartile)
 */
template <typename F>
F math::SampleQuantileGeneric<F>::iqr(const math::EQntlType::type method) const
{
    return ( this->quartile(3, method) - this->quartile(1, method) );
}


/**
 * Empirical cumulative distribution function.
 * Returns a proportion of all sample's elements that
 * are less than or equal to 't'.
 *  
 * @param t - input value
 * 
 * @return empirical cumulative distribution function for 't' 
 */
template <typename F>
F math::SampleQuantileGeneric<F>::ecdf(const F& t) const
{
    const size_t WMAX = 5;

    const size_t& N = this->m_N;
    const std::vector<F>& v = this->m_v;

    // Handling two "corner cases"
    if ( t < v.at(0) )
    {
        return static_cast<F>(0);
    }

    if ( t > v.at(N-1) )
    {
        return static_cast<F>(1);
    }
    
    size_t k = v.at(N / 2);
    size_t kl, ku;

    // Determine the initial search interval, either the lower or the upper half
    if ( v.at(k) < t )
    {
        kl = k;
        ku = N - 1;
    }
    else
    {
        kl = 0;
        ku = k;
    }

    /*
     * Narrow the search interval using the bisection method.
     * 
     * Note: as the final index will be further adjusted (and this
     * operation is very fast), it is not necessary to narrow down the
     * search interval to the width of 1.
     */
    while ( (ku-kl) > WMAX )
    {
        // A more robust version (w.r.t. int. range) of (ku+kl)/2
        k = kl + (ku - kl) / 2;

        if ( v.at(k) < t  )
        {
            kl = k;
        }
        else
        {
            ku = k;
        }
    }

    /*
     * Final adjustment of the index 'kl' that properly handles the
     * situation when several sample's values are equal to 't'.
     * As 'kl' is guaranteed to be less than than 't' and already
     * reasonably close to the final value, it is sensible to iteratively
     * increment it by 1.
     */
    for ( ; kl<N && v.at(kl)<=t; ++kl );

    /*
     * In zero based indexed arrays, the adjusted 'kl' will
     * denote the number of elements that are less than or equal to 't'
     */

    return static_cast<F>(kl) / static_cast<F>(N);
}


/**
 * @return minimum observation of the sample
 */
template <typename F>
F math::SampleQuantileGeneric<F>::min() const
{
    // Note that this->m_v is already sorted in ascending order
    return this->m_v.at(0);
}


/**
 * @return maximum observation of the sample
 */
template <typename F>
F math::SampleQuantileGeneric<F>::max() const
{
    // Note that this->m_v is already sorted in ascending order
    return this->m_v.at(this->m_N-1);
}


/**
 * Returns sample's n.th largest or smallest element
 * 
 * @note Index origin depends on 'zerobase'
 * 
 * @param n - desired index
 * @param largest - if TRUE, the n.th largest element is returned, otherwise the n.th smallest (default: TRUE)
 * @param zerobase - does 'n' denote a zero-based index (default: TRUE)
 * 
 * @return sample's n.th largest or smallest element, depending on 'largest'
 * 
 * @throw StatisticsException if 'n' is invalid
 */
template <typename F>
F math::SampleQuantileGeneric<F>::elem(
             const size_t n, 
             const bool largest,
             const bool zerobase
           ) const throw(math::StatisticsException)
{
    // sanity check
    if ( (true==zerobase && n>=this->m_N) ||
         (false==zerobase && (n<=0 || n>this->m_N) ) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    const size_t nn = ( true==zerobase ? n : n-1 );


    // this->m_v is already sorted in ascending order

    const size_t idx = ( false==largest ? nn : (this->m_N - 1) - nn );

    return this->m_v.at(idx);
}


/**
 * Test if an observation is an outlier regarding the sample.
 * 
 * An observation is an outlier if it lies either below the first quartile - iqrs * IQR
 * or above the third quartile + iqrs * IQR.
 * 
 * The exact value of 'iqrs' can be passed as an argument and is typically
 * equal to 1.5.
 *
 * Note that quartiles are estimated from the original sample that was passed to
 * the constructor. In other words, if 'val' is not a member of the sample,
 * quartiles will not be reestimated to include this sample.
 * 
 * @param val - value to be checked
 * @param iqrs - number of interquartile ranges below/above Q1 and Q3 (default 1.5)
 * @param method - method to estimate quartiles (default R7)
 * 
 * @return a logical value indicating whether 'val' is an outlier or not 
 */
template <typename F>
bool math::SampleQuantileGeneric<F>::isOutlier(
                const F& val,
                const F& iqrs,
                const math::EQntlType::type method) const
{
    const F q1 = this->quartile(1, method);
    const F q3 = this->quartile(3, method);
    const F diff = q3 - q1;   // IQR

    return ( ( val < (q1-iqrs * diff) ) ||
             ( val > (q3+iqrs * diff) ) );
}


/**
 * Fills all sample's outliers into the set 'outl'.
 * 
 * An observation is an outlier if it lies either below the first quartile - iqrs * IQR
 * or above the third quartile + iqrs * IQR.
 * 
 * The exact value of 'iqrs' can be passed as an argument and is typically
 * equal to 1.5.
 * 
 * @param outl - a reference to a set to be filled with outliers
 * @param iqrs - number of interquartile ranges below/above Q1 and Q3 (default 1.5)
 * @param method - method to estimate quartiles (default R7)
 * 
 * @throw StatisticsException if allocation of memory fails
 */
template <typename F>
void math::SampleQuantileGeneric<F>::outliers(
                std::set<F>& outl,
                const F& iqrs,
                const math::EQntlType::type method) const
            throw (math::StatisticsException)
{
    try
    {
        const F q1 = this->quartile(1, method);
        const F q3 = this->quartile(3, method);
        const F diff = q3 - q1;   // IQR
        const F lowerBound = q1 - iqrs * diff;
        const F upperBound = q3 + iqrs * diff;

        /*
         * Vector this->m_v is already sorted in ascending order. That said,
         * it is sufficient to start iterating from the vector's start as long
         * as the elements are below the lower bound. Then another iteration is
         * performed from the vector's end as long as the elements are above
         * the upper bound. 
         */

        typename std::vector<F>::const_iterator it;
        for ( it=this->m_v.begin(); 
              it!=this->m_v.end() && *it<lowerBound;
              ++it )
        {
            outl.insert(*it);
        }

        typename std::vector<F>::const_reverse_iterator rit;
        for ( rit = this->m_v.rbegin();
              rit!=this->m_v.rend() && *rit>upperBound;
              ++rit )
        {
            outl.insert(*rit);
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
    }
}


/**
 * Destructor
 */
template <typename F>
math::SampleQuantileGeneric<F>::~SampleQuantileGeneric()
{
    // Vector's destructors would probably clean up this automatically.
    // Anyway let us clear the vector, just to be aware of allocated resources.
    this->m_v.clear();

    // No other resources to release.
}
