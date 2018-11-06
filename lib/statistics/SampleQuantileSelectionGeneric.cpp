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
 * Implementation of the class SampleQuantileSelectionGeneric that
 * estimates quantiles of a sample, based on implemented functionality
 * for selection of the n.th smallest element.
 */


// no #include "SampleQuantileSelectionGeneric.hpp" !
#include <cstddef>
#include <vector>
#include <set>
#include <new>
#include <algorithm>

#include "util/mtcopy.hpp"
#include "util/SelectionGeneric.hpp"
#include "statistics/SampleStatGeneric.hpp"

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/StatisticsException.hpp"
#include "exception/SelectionException.hpp"


/**
 * Constructor.
 * Creates its own copy of the vector of observations and sorts
 * it in ascending order.
 *
 * @note Please see a detailed description of this class (documentation
 * comments in SampleQuantileSelectionGeneric.hpp) for more details when
 * instantiation of this class with an internal copy of 'sample' ('copy'
 * set to TRUE) makes sense.
 *
 * @param sample - a vector of observations
 * @param copy - create and use an internal copy of 'sample' (default: FALSE)
 *
 * @throw StatisticsException if 'sample' is empty or allocation of memory for its copy failed
 */
template <typename F>
math::SampleQuantileSelectionGeneric<F>::SampleQuantileSelectionGeneric(
            const std::vector<F>& sample,
            const bool copy
          )
  : math::SampleQuantileGenericAb<F>( sample.size() ),
    m_v( true==copy ? m_stor : sample )
{
    /*
     * When a copy of the sample vector is requested, 'm_v' will refer to
     * the vector 'm_stor' that stores an internal copy. When a copy is
     * not requested, 'm_v' will simply refer to the original sample vector.
     */

    try
    {
        if ( 0 == this->m_N )
        {
            throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
        }

        if ( true == copy )
        {
            // copy the 'sample' into an internal vector:
            math::mtcopy(sample, this->m_stor);
        }
        else
        {
            // 'm_stor' is not needed in this case, it will be cleared
            this->m_stor.clear();
        }
    }
    catch ( const std::bad_alloc& ex )
    {
        throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
    }
}


/*
 * @param n - ordinal
 *
 * @return n.th smallest element of the sample
 *
 * @throw StatisticsExcpetion if 'n' is larger or equal to the number of elements or if selection algorithm's internal allocation of memory fails
 */
template <typename F>
F math::SampleQuantileSelectionGeneric<F>::_select(const std::size_t n) const
{
    if ( n >= this->m_N )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    try
    {
        // An efficient selection algorithm is implemented by this function
        return math::Selection::select(this->m_v, n);
    }
    catch ( const math::SelectionException& sox )
    {
        throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
    }
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
 * @throw StatisticsExcpetion if either 'n1' or 'n2' is larger or equal to the number of elements or if selection algorithm's internal allocation of memory fails
 */
template <typename F>
void math::SampleQuantileSelectionGeneric<F>::_select2(
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

    try
    {
        // An efficient selection algorithm is implemented by this function
        math::Selection::select2(this->m_v, n1, n2, val1, val2);
    }
    catch ( const math::SelectionException& sox )
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
template <typename F>
F math::SampleQuantileSelectionGeneric<F>::ecdf(const F& t) const
{
    // This functionality has been efficiently implemented by this function:
    return math::SampleStat::ecdf(this->m_v, t, true);
}


/**
 * @return minimum observation of the sample
 */
template <typename F>
F math::SampleQuantileSelectionGeneric<F>::min() const
{
    // This functionality has been efficiently implemented by this function:
    return math::Selection::min(this->m_v);
}


/**
 * @return maximum observation of the sample
 */
template <typename F>
F math::SampleQuantileSelectionGeneric<F>::max() const
{
    // This functionality has been efficiently implemented by this function:
    return math::Selection::max(this->m_v);
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
F math::SampleQuantileSelectionGeneric<F>::elem(
             const std::size_t n,
             const bool largest,
             const bool zerobase
           ) const
{
    // sanity check
    if ( (true==zerobase && n>=this->m_N) ||
         (false==zerobase && (n<=0 || n>this->m_N) ) )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_ARG);
    }

    // the actual ordinal
    const std::size_t nn = ( true==zerobase ? n : n-1 );

    // An efficient selection algorithm is implemented by this function
    return math::Selection::select(this->m_v, nn, !largest);
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
void math::SampleQuantileSelectionGeneric<F>::outliers(
                std::set<F>& outl,
                const F& iqrs,
                const math::EQntlType::type method
              ) const
{
    try
    {
        const std::size_t& N = this->m_N;

        // Both boundaries of the range to determine outliers
        F lowerBound,  upperBound;
        this->_outlierBounds(lowerBound, upperBound, iqrs, method);

        /*
         * Coarse grained parallelization.
         * Each thread will inspect all observations within "its" block.
         * If the observation is an outlier, it will be inserted into
         * the output set.
         */
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(N, outl, lowerBound, upperBound)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            typename std::vector<F>::const_iterator it = this->m_v.begin() + istart;
            for ( std::size_t cntr = 0;
                  cntr<elems_per_thread && it!=this->m_v.end();
                  ++it, ++cntr )
            {
                const F& xcur = *it;

                if ( xcur<lowerBound || xcur>upperBound )
                {
                    /*
                     * Note that modification (including insertion) of
                     * STL data structures is not thread safe, hence this
                     * should only be done within a critical section.
                     * Also note that the total number of outliers should
                     * be relatively low so the "serialization" of insertions
                     * should not affect performance significantly.
                     */

                    #pragma omp critical(samplequantileselectiongeneric_outliers)
                    {
                        outl.insert(xcur);
                    }
                }
            }  // for

            (void) iend;
        }  // omp parallel
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
math::SampleQuantileSelectionGeneric<F>::~SampleQuantileSelectionGeneric()
{
    // just clear the vector 'm_stor'
    this->m_stor.clear();
}
