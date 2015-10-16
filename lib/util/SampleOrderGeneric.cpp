/*
Copyright 2015, Jernej Kovacic

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
 * Implementation of functions within the namespace SampleOrder
 * that return indices of the given vector's elements in a stably
 * sorted vector.
 */


// no #include "SampleOrderGeneric.hpp" !
#include <cstddef>
#include <vector>
#include <algorithm>
#include <new>

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/SampleOrderException.hpp"


namespace math {  namespace SampleOrder {  namespace __private
{


/*
 * An internal functor class that implements comparison of two elements, specified by
 * their indices. It implements the operator () that performs either "less than"
 * or "greater than" comparison, depending on the parameter 'asc'.
 */
template <typename F>
class IndexCmp
{
private:
    const std::vector<F>* m_pvec;   // pointer to a vector of sample elements
    bool m_asc;                     // should order in ascending order?

public:

    /*
     * Constructor that accepts the pointer to a vector of elements and
     * type of sorting (ascending or descending)
     *
     * @note: as the class is used internally only, it is reasonable to assume
     *        that'pvec' will never be NULL.
     *
     * @param pvec - pointer to a vector of sample elements (should not be NULL)
     * @param asc - logical value indicating whether elements' indices should be sorted in ascending (true) or descending (false) order
     */
    IndexCmp(const std::vector<F>* pvec, const bool asc=true) :
        m_pvec(pvec), m_asc(asc)
    {
        // nothing else to do
    }


    /*
     * Operator () that performs comparison of elements at positions 'a' and 'b'.
     *
     * @param a - position of the first element
     * @param b - position of the second element
     *
     * @return logical value, indicating either v[a]<v[b] (if 'asc'==true) or v[a]>v[b]
     */
    bool operator()(const size_t a, const size_t b)
    {
        // Note that sorting algorithms only require that F defines the operator '<'

        return ( true == this->m_asc ?
                 ( this->m_pvec->at(a) < this->m_pvec->at(b) ) :
                 ( -this->m_pvec->at(a) < -this->m_pvec->at(b) )  );
    }

};  // class IndexCmp


/*
 * Finds sample's either minimum or maximum observation, depending
 * on 'min'.
 *
 * @param x - vector of observations
 * @param min - a logical value indicating whether minimum or maximum value should be returned
 *
 * @return sample's minimum/maximum observation
 *
 * @throw SampleOrdersException if 'x' is empty
 */
template <typename F>
F __minmax(const std::vector<F>& x, const bool min) throw(math::SampleOrderException)
{
    const size_t N = x.size();

    // sanity check
    if ( 0 == N )
    {
        throw math::SampleOrderException(math::SampleOrderException::SAMPLE_EMPTY);
    }

    // the first observation is the first candidate for the extreme value...
    F retVal = x.at(0);

    // Coarse grained parallelism:
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x, retVal)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

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
        #pragma omp critical(sampleordergeneric_minmax)
        {
            retVal = ( true==min ? std::min(retVal, temp) : std::max(retVal, temp) );
        }

        (void) iend;
    }  // omp parallel

    return retVal;
}


/*
 * Finds index of sample's minimum or maximum observation, depending
 * on 'min'. In case of ties, the smallest index is returned.
 *
 * @param x - vector of observations
 * @param min - a logical value indicating whether minimum or maximum value should be returned
 *
 * @return the smallest index of sample's minimum/maximum observation
 *
 * @throw SampleOrdersException if 'x' is empty
 */
template <typename F>
F __whichMinMax(const std::vector<F>& x, const bool min) throw(math::SampleOrderException)
{
    const size_t N = x.size();

    // sanity check
    if ( 0 == N )
    {
        throw math::SampleOrderException(math::SampleOrderException::SAMPLE_EMPTY);
    }

    // the first observation is the first candidate for the extreme value...
    F extVal = x.at(0);
    F extIdx = 0;

    // Coarse grained parallelism:
    #pragma omp parallel num_threads(ompIdeal(N)) \
                if(N>OMP_CHUNKS_PER_THREAD) \
                default(none) shared(x, extVal, extIdx)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        typename std::vector<F>::const_iterator it = x.begin() + istart;
        size_t i = istart;

        // the first value of the block is the first candidate for the local extreme
        F temp = *it;
        size_t idx = istart;

        for ( size_t cntr=0;
              cntr<elems_per_thread && it!=x.end();
              ++it, ++cntr, ++i )
        {
            // compare the current element with 'temp' and
            // update 'temp' and 'idx' if necessary
            if ( ( true==min && *it<temp ) ||
                 ( false==min && *it>temp ) )
            {
                temp = *it;
                idx = i;
            }
        }

        // prevent possible race condition when updating 'extVal' and 'extIdx'
        #pragma omp critical(sampleordergeneric_whichminmax)
        {
            /*
             * Compare each thread's extreme to the current global extreme.
             *
             * Note: this rather complicated combination of conditions ensures
             * that the smallest index will be returned in case of ties
             */
            if ( ( true==min && (temp<extVal || (temp==extVal && idx<extIdx) ) ) ||
                ( false==min && (temp>extVal || (temp==extVal && idx<extIdx) ) ) )
            {
                extVal = temp;
                extIdx = idx;
            }
        }

        (void) iend;
    }  // omp parallel

    return extIdx;
}

}}}  // namespace math::SampleOrder::__private



/**
 * Returns indices of each element in the stably sorted vector.
 * In case of ties, the order of indices remains unmodified.
 *
 * @note The input vector 'x' is not modified.
 *
 * @note If 'x' is empty, an empty vector is returned.
 *
 * @note The function is equivalent to the function 'order' in R.
 *
 * @param x - vector of elements
 * @param dest - reference to a vector to write indices of elements of 'x' in the sorted vector
 * @param asc - logical value indicating whether to sort in ascending order or not (default: TRUE)
 *
 * @return reference to 'dest'
 *
 * @throw SampleOrderException if allocation of the return vector fails
 */
template <typename F>
std::vector<size_t>& math::SampleOrder::order(
        const std::vector<F>& x,
        std::vector<size_t>& dest,
        const bool asc
      ) throw(math::SampleOrderException)
{
    try
    {
        const size_t N = x.size();

        dest.resize(N);

        if ( 0 == N )
        {
            return dest;
        }


        /*
         * Assign 'dest' with initial positions.
         */
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(dest)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            typename std::vector<size_t>::iterator it = dest.begin() + istart;
            size_t i = istart;
            for ( size_t cntr=0;
                  cntr<elems_per_thread && it!=dest.end();
                  ++cntr, ++i, ++it )
            {
                *it = i;
            }

            (void) iend;
        }  // omp parallel

        math::SampleOrder::__private::IndexCmp<F> idxCmp(&x, asc);

        /*
         * Performs stable sort of indices, assuring that order
         * of indices remains unmodified in case of ties.
         */
        std::stable_sort(dest.begin(), dest.end(), idxCmp);

    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::SampleOrderException(math::SampleOrderException::OUT_OF_MEMORY);
    }

    return dest;
}


/**
 * @param x - vector of observations
 *
 * @return minimum observation of the sample
 *
 * @throw SampleOrderException if 'x' is empty
 */
template <typename F>
F math::SampleOrder::min(const std::vector<F>& x) throw(math::SampleOrderException)
{
    return math::SampleOrder::__private::__minmax<F>(x, true);
}


/**
 * @param x - vector of observations
 *
 * @return maximum observation of the sample
 *
 * @throw SampleOrderException if 'x' is empty
 */
template <typename F>
F math::SampleOrder::max(const std::vector<F>& x) throw(math::SampleOrderException)
{
    return math::SampleOrder::__private::__minmax<F>(x, false);
}


/**
 * @param x - vector of observations
 *
 * @return the smallest zero based index of the minimum observation of the sample
 *
 * @throw SampleOrderException if 'x' is empty
 */
template <typename F>
size_t math::SampleOrder::whichMin(const std::vector<F>& x) throw(math::SampleOrderException)
{
    return math::SampleOrder::__private::__whichMinMax(x, true);
}


/**
 * @param x - vector of observations
 *
 * @return the smallest zero based index of the maximum observation of the sample
 *
 * @throw SampleOrderException if 'x' is empty
 */
template <typename F>
size_t math::SampleOrder::whichMax(const std::vector<F>& x) throw(SampleOrderException)
{
	return math::SampleOrder::__private::__whichMinMax(x, false);
}
