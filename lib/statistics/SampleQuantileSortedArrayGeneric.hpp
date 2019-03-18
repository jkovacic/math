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
 * Declaration and implementation  of the class SampleQuantileSortedArrayGeneric
 * that estimates quantiles of a sample, based on a sorted copy of the
 * sample vector.
 */


#ifndef _MATH_SAMPLEQUANTILESORTEDARRAYGENERIC_HPP_
#define _MATH_SAMPLEQUANTILESORTEDARRAYGENERIC_HPP_


#include <cstddef>
#include <vector>
#include <set>
#include <new>
#include <algorithm>

#include "statistics/SampleQuantileGenericAb.hpp"
#include "util/mtcopy.hpp"

#include "../settings/stat_settings.h"

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

    /*
     * @param n - ordinal
     *
     * @return n.th smallest element of the sample
     *
     * @throw StatisticsExcpetion if 'n' is larger or equal to the number of elements
     */
    virtual F _select(const std::size_t n) const
    {
        if ( n >= this->m_N )
        {
            throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
        }

        return this->m_v.at(n);
    }



    /*
     * Finds the n1.th and n2.th smallest value of the sample and assigns them
     * to 'val1' and 'val2', respectively.
     * 
     * @param n1 - first ordinal
     * @param n2 - second ordinal
     * @param val1 - reference to the variable to assign the n1.th smallest value
     * @param val2 - reference to the variable to assign the n2.th smallest value
     * 
     * @throw StatisticsExcpetion if either 'n1' or 'n2' is larger or equal to the number of elements
     */
    virtual void _select2(
            const std::size_t n1,
            const std::size_t n2,
            F& val1,
            F& val2
          ) const
    {
        if ( n1>=this->m_N || n2>=this->m_N )
        {
            throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
        }

        val1 = this->m_v.at(n1);
        val2 = this->m_v.at(n2);
    }


public:

    /**
     * Constructor.
     * Creates its own copy of the vector of observations and sorts
     * it in ascending order.
     *
     * @param sample - a vector of observations
     *
     * @throw StatisticsException if 'sample' is empty or allocation of memory for its copy failed
     */
    SampleQuantileSortedArrayGeneric(const std::vector<F>& sample)
      : SampleQuantileGenericAb<F>( sample.size() )
    {
        try
        {
            if ( 0 == this->m_N )
            {
                throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
            }

            // copy the 'sample' into an internal vector:
            math::mtcopy(sample, this->m_v);

            // and sort it in ascending order
            std::sort(this->m_v.begin(), this->m_v.end());
        }
        catch ( const std::bad_alloc& ex )
        {
            throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
        }
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
    virtual F ecdf(const F& t) const
    {
        const std::size_t WMAX = 5;

        const std::size_t& N = this->m_N;
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

        std::size_t k = v.at(N / 2);
        std::size_t kl, ku;

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
    virtual F min() const
    {
        // Note that this->m_v is already sorted in ascending order
        return this->m_v.at(0);
    }



    /**
     * @return maximum observation of the sample
     */
    virtual F max() const
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
    virtual F elem(
           const std::size_t n,
           const bool largest = true,
           const bool zerobase = STAT_DEFAULT_ZERO_BASE
         ) const
    {
        // sanity check
        if ( (true==zerobase && n>=this->m_N) ||
             (false==zerobase && (n<=0 || n>this->m_N) ) )
        {
            throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
        }

        const std::size_t nn = ( true==zerobase ? n : n-1 );


        // this->m_v is already sorted in ascending order

        const std::size_t idx = ( false==largest ? nn : (this->m_N - 1) - nn );

        return this->m_v.at(idx);
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
    virtual void outliers(
           std::set<F>& outl,
           const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
         ) const
   {
        try
        {
            F lowerBound,  upperBound;
            this->_outlierBounds(lowerBound, upperBound, iqrs, method);

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
                const F& vcur = *it;
                outl.insert(vcur);
            }

            typename std::vector<F>::const_reverse_iterator rit;
            for ( rit = this->m_v.rbegin();
                  rit!=this->m_v.rend() && *rit>upperBound;
                  ++rit )
            {
                const F& vcur = *rit;
                outl.insert(vcur);
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
    virtual ~SampleQuantileSortedArrayGeneric()
    {
        // just clean the vector 'm_v'
        this->m_v.clear();
    }
    
};


// Samples with elements of types float, double and long double
// make most sense, therefore the following types are predefined

typedef SampleQuantileSortedArrayGeneric<float>          FSampleQuantileSortedArray;
typedef SampleQuantileSortedArrayGeneric<double>         SampleQuantileSortedArray;
typedef SampleQuantileSortedArrayGeneric<long double>    LDSampleQuantileSortedArray;

}  // namespace math


#endif  // _MATH_SAMPLEQUANTILESORTEDARRAYGENERIC_HPP_
