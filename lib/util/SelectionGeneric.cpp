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
 * Implementation of functions within the namespace Selection
 * that select the i.th largest/smallest element of a vector.
 */


// no #include "SelectionGeneric.hpp" !
#include <cstddef>
#include <vector>
#include <algorithm>
#include <new>

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/SelectionException.hpp"


namespace math {  namespace Selection {  namespace __private
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
    bool operator()(const std::size_t a, const std::size_t b)
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
 * @throw SelectionException if 'x' is empty
 */
template <typename F>
F __minmax(const std::vector<F>& x, const bool min) throw(math::SelectionException)
{
    const std::size_t N = x.size();

    // sanity check
    if ( 0 == N )
    {
        throw math::SelectionException(math::SelectionException::SAMPLE_EMPTY);
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

        for ( std::size_t cntr=0;
              cntr<elems_per_thread && it!=x.end();
              ++it, ++cntr)
        {
            const F& xcur = *it;

            // update 'temp' depending on 'min'
            temp = ( true==min ? std::min(temp, xcur) : std::max(temp, xcur) );
        }

        // prevent possible race condition when updating retVal
        #pragma omp critical(selectiongeneric_minmax)
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
 * @throw SelectionException if 'x' is empty
 */
template <typename F>
F __whichMinMax(const std::vector<F>& x, const bool min) throw(math::SelectionException)
{
    const std::size_t N = x.size();

    // sanity check
    if ( 0 == N )
    {
        throw math::SelectionException(math::SelectionException::SAMPLE_EMPTY);
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
        std::size_t i = istart;

        // the first value of the block is the first candidate for the local extreme
        F temp = *it;
        std::size_t idx = istart;

        for ( std::size_t cntr=0;
              cntr<elems_per_thread && it!=x.end();
              ++it, ++cntr, ++i )
        {
            const F& xcur = *it;

            // compare the current element with 'temp' and
            // update 'temp' and 'idx' if necessary
            if ( ( true==min && xcur<temp ) ||
                 ( false==min && xcur>temp ) )
            {
                temp = xcur;
                idx = i;
            }
        }

        // prevent possible race condition when updating 'extVal' and 'extIdx'
        #pragma omp critical(selectiongeneric_whichminmax)
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


/*
 * Allocates the vector 'dest' to 'N' elements and fills in with integers
 * between 0 and N-1 in ascending order.
 *
 * @note If 'N' is zero, an empty vector is returned.
 *
 * @param dest - reference to a vector to be filled with integers
 * @param N - desired number of elements in 'dest'.
 *
 * @return reference to 'dest'
 *
 * @throw SelectionException if allocation of the return vector fails
 */
std::vector<std::size_t>& fillIndices(
        std::vector<std::size_t>& dest,
        const std::size_t N )
        throw(math::SelectionException)
{
    try
    {
        // (re)allocates 'dest' to 'N' elements
        dest.resize(N);

        // Nothing else to do if N==0
        if ( 0 == N )
        {
            return dest;
        }

        /*
         * Assign 'dest' with initial positions 0..(N-1)
         */
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(dest)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            typename std::vector<std::size_t>::iterator it = dest.begin() + istart;
            std::size_t i = istart;
            for ( std::size_t cntr=0;
                  cntr<elems_per_thread && it!=dest.end();
                  ++cntr, ++i, ++it )
            {
                std::size_t& destCurr = *it;
                destCurr = i;
            }

            (void) iend;
        }  // omp parallel
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::SelectionException(math::SelectionException::OUT_OF_MEMORY);
    }

    return dest;
}

}}}  // namespace math::Selection::__private



/**
 * Returns the index table of 'x', i.e. a vector of indices of each
 * element in the stably sorted vector. In case of ties, the order of indices
 * remains unmodified.
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
 * @throw SelectionException if allocation of the return vector fails
 */
template <typename F>
std::vector<std::size_t>& math::Selection::order(
        const std::vector<F>& x,
        std::vector<std::size_t>& dest,
        const bool asc
      ) throw(math::SelectionException)
{
    const std::size_t N = x.size();

    // allocates 'dest' and fills it with initial positions 0..(N-1)
    math::Selection::__private::fillIndices(dest, N);

    if ( 0 == N )
    {
        return dest;
    }

    math::Selection::__private::IndexCmp<F> idxCmp(&x, asc);

    /*
     * Performs stable sort of indices, assuring that order
     * of indices remains unmodified in case of ties.
     */
    std::stable_sort(dest.begin(), dest.end(), idxCmp);

    return dest;
}


/**
 * Returns the rank table of 'x' The table indicates each element's
 * position if 'x' were sorted in ascending or descending order,
 * depending on 'asc'. In case of ties, all equal elements are assigned
 * the same (i.e the smallest possible) rank, as it is typical in sports
 * ranking.
 *
 * @note The input vector 'x' is not modified.
 *
 * @note If 'x' is empty, an empty vector is returned.
 *
 * @note The function is equivalent to the function 'rank(x, ties.method="min")' in R.
 *
 * @param x - vector of elements
 * @param dest - reference to a vector to write ranks of elements of 'x' to
 * @param asc - logical value indicating whether ranks indicate  ascending order or not (default: TRUE)
 *
 * @return reference to 'dest'
 *
 * @throw SelectionException if allocation of the return vector fails
 */
template <typename F>
std::vector<std::size_t>& math::Selection::rank(
            const std::vector<F>& x,
            std::vector<std::size_t>& dest,
            const bool asc = true
          ) throw(math::SelectionException)
{
    try
    {
        const std::size_t N = x.size();
        std::vector<std::size_t> indexTable;

        dest.resize(N);

        if ( 0 == N )
        {
            return dest;
        }

        // The most convenient algorithm requires the
        // index table of 'x'.
        math::Selection::order(x, indexTable, asc);

        // indexTable should not be modified anymore:
        const std::vector<std::size_t>& idx = indexTable;

        /*
         * The algorithm is actually quite simple. It relies on
         * the fact, that idx[i] denotes the position of the i.th
         * largest/smallest element of 'x'. Hence the most basic
         * algorithm would be as follows:
         *
         * for i in 0 (N-1):
         *   rank[ idx[i] ] = i
         *
         * However, this algorithm does not handle ties, therefore it
         * requires a few modifications as explained below. The modified
         * algorithm, however, must only run sequentially and it is
         * not possible to parallelize it.
         */

        /*
         * But in any case it is perfectly safe to assign the value 0
         * to (one of) the smallest/largest element(s)...
         */
        dest.at(idx.at(0)) = 0;

        /*
         * For all the remaining elements, the algorithm above is slightly
         * modified. x[idx[i]] is compared to the previously ranked element
         * (x[idx[i-1]]). If they are equal, just copy the previously ranked
         * element's rank (rank[idx[i-1]]). Otherwise assign the lement x[idx[i]]
         * the value 'i'.
         */
        for (std::size_t i = 1; i<N; ++i )
        {
            dest.at(idx.at(i)) =
                ( x.at(idx.at(i-1)) == x.at(idx.at(i)) ?
                  dest.at(idx.at(i-1)) :  i );
        }
    }
    catch ( const std::bad_alloc& ex )
    {
        throw math::SelectionException(math::SelectionException::OUT_OF_MEMORY);
    }

    return dest;
}


/**
 * @param x - vector of observations
 *
 * @return minimum observation of the sample
 *
 * @throw SelectionException if 'x' is empty
 */
template <typename F>
F math::Selection::min(const std::vector<F>& x) throw(math::SelectionException)
{
    return math::Selection::__private::__minmax<F>(x, true);
}


/**
 * @param x - vector of observations
 *
 * @return maximum observation of the sample
 *
 * @throw SelectionException if 'x' is empty
 */
template <typename F>
F math::Selection::max(const std::vector<F>& x) throw(math::SelectionException)
{
    return math::Selection::__private::__minmax<F>(x, false);
}


/**
 * @param x - vector of observations
 *
 * @return the smallest zero based index of the minimum observation of the sample
 *
 * @throw SelectionException if 'x' is empty
 */
template <typename F>
std::size_t math::Selection::whichMin(const std::vector<F>& x) throw(math::SelectionException)
{
    return math::Selection::__private::__whichMinMax(x, true);
}


/**
 * @param x - vector of observations
 *
 * @return the smallest zero based index of the maximum observation of the sample
 *
 * @throw SelectionException if 'x' is empty
 */
template <typename F>
std::size_t math::Selection::whichMax(const std::vector<F>& x) throw(SelectionException)
{
    return math::Selection::__private::__whichMinMax(x, false);
}
