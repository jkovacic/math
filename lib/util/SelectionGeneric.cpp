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

#include "util/FillVectors.h"

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
    const std::vector<F>& m_vec;    // reference to a vector of sample elements
    const bool m_asc;               // should order in ascending order?

public:

    /*
     * Constructor that accepts the reference to a vector of elements and
     * type of sorting (ascending or descending)
     *
     * @param vec - vector of sample elements
     * @param asc - logical value indicating whether elements' indices should be sorted in ascending (true) or descending (false) order
     */
    IndexCmp(const std::vector<F>& vec, const bool asc=true) :
        m_vec(vec), m_asc(asc)
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
                 ( this->m_vec.at(a) < this->m_vec.at(b) ) :
                 ( -this->m_vec.at(a) < -this->m_vec.at(b) )  );
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


/**
 * Allocates the vector 'dest' and fills in with pointers to all
 * elements of 'x'.
 * 
 * @param x - vector of elements
 * @param dest - reference to a vector to be filled with pointers to elements of 'x'
 * 
 * @return reference to 'dest'
 * 
 * @throw SelectionException if allocation of return vector fails
 */
template <typename F>
std::vector<const F*>& __fillPointers(
        const std::vector<F>& x,
        std::vector<const F*>& dest
      ) throw (math::SelectionException)
{
    try
    {
        const std::size_t N = x.size();
        dest.resize(N, NULL);

        // Nothing else to do if N==0
        if ( 0 == N )
        {
            return dest;
        }

        /*
         * Assign 'dest' with pointers to x[0] .. x[N-1]
         */
        #pragma omp parallel num_threads(ompIdeal(N)) \
                    if(N>OMP_CHUNKS_PER_THREAD) \
                    default(none) shared(x, dest)
        {
            OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

            typename std::vector<F>::const_iterator srcit = x.begin() + istart;
            typename std::vector<const F*>::iterator it = dest.begin() + istart;
            for ( std::size_t cntr=0;
                  cntr<elems_per_thread && it!=dest.end();
                  ++cntr, ++srcit, ++it )
            {
                *it = &(*srcit);
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


/*
 * Median of 3.
 * 
 * Finds index of the median element of *ptrs[idx1], *ptrs[idx2] and *ptrs[idx3].
 * 
 * @param ptrs - vector of pointers to elements
 * @param idx1 - index of the first element
 * @param idx2 - index of the second element
 * @param idx3 - index of the third element
 * 
 * @return index of median(*ptrs[idx1], *ptrs[idx2], *ptrs[idx3])
 */
template <typename F>
std::size_t __medianOf3Idx(
                const std::vector<F*>& ptrs,
                const std::size_t idx1,
                const std::size_t idx2,
                const std::size_t idx3 )
{
    #define MATH_SELECTIONGENERIC_MEDIANOF3_BUFSIZE        ( 3 )

    std::size_t buf[ MATH_SELECTIONGENERIC_MEDIANOF3_BUFSIZE ] =
                   { idx1, idx2, idx3 };

    /*
     * The algorithm is generalized to compute the median of 'N' for small
     * values of 'N' and 'N' ideally being an odd number (not divisible by 2).
     * The values are sorted using the insertion sort algorithm which is known
     * to be very convenient for small values of 'N'. For N=3, the algorithm
     * requires no more than 3 comparisons and exchanges.
     */
    for ( std::size_t i=1; i<MATH_SELECTIONGENERIC_MEDIANOF3_BUFSIZE; ++i )
    {
        for ( std::size_t j = i;
              j>0 && ( *ptrs.at(buf[j]) < *ptrs.at(buf[j-1]) );
              --j )
        {
            std::swap(buf[j], buf[j-1]);
        }  // for j
    }  // for i

    return buf[MATH_SELECTIONGENERIC_MEDIANOF3_BUFSIZE / 2];

    #undef MATH_SELECTIONGENERIC_MEDIANOF3_BUFSIZE
}


/*
 * "Partitions" the given subarray 'px' between 'left' and 'right, rearranges the
 * array of pointers 'px' int he given range and returns an index 'k' s.t. all
 * elements *px[i<k] are guaranteed to be less than *px[k] and all elements *px[i>k]
 * are guaranteed to be greater than *px[k].
 *
 * @note This function is only called internally hence it is assumed that all
 *       arguments are valid, i.e., 'px' is preallocated, left<right, etc.
 *
 * @note The function only rerranges elements of 'px' in the range left <= i <= right
 *
 * @note The same function is also a part of the quicksort algorithm
 *
 * @param px - vector of pointers to be rearranged
 * @param left - left margin of the partition interval
 * @param right - right margin of the partition interval
 *
 * @return 'k' that satisfies *px[i<k] <= *px[k] and *px[k] <= *px[i>k]
 */
template <typename F>
std::size_t __partition(
        std::vector<const F*>& px,
        const std::size_t left,
        const std::size_t right )
{
    /*
     * The algorithm is described in detail at:
     *
     * Robert Sedgewick, Kevin Wayne
     * Algorithms (4th Edition)
     * Addison-Wesley Professional, 2011
     *
     * http://algs4.cs.princeton.edu
     */


    /*
     * The pivot is chosen by the "median of 3" algorithm. It finds the median
     * value of the leftmost, rightmost and the middle element.
     */

    const std::size_t med3 = math::Selection::__private::__medianOf3Idx(
                    px,
                    left,
                    (left + right) / 2,
                    right );

    /*
     * If the selected pivot is not the leftmost element,
     * they must be exchanged. This way the pivot remains
     * the leftmost element of the range.
     */
    if ( med3 != left )
    {
        std::swap(px.at(med3), px.at(left));
    }

    std::size_t i = left + 1;
    std::size_t j = right;
    const F pval = *(px.at(left));

    for ( ; ; )
    {
        // Find both elements to swap
        for ( ; i<right && *(px.at(i)) < pval; ++i );
        for ( ; j>left && pval < *(px.at(j)); --j );

        // If the pointers have crossed, terminate the "infinite" loop
        if ( i >= j )
        {
            break;  // out of for(;;)
        }

        // ...and swap the elements of 'px'
        std::swap( px.at(i), px.at(j) );

    }

    // Move the "pivot" to its appropriate position within 'px'.
    if ( j != left )
    {
        std::swap( px.at(left), px.at(j) );
    }

    return j;
}


/*
 * Finds the P.th smallest element of 'pvec' in the range between 'from' and 'to'.
 *
 * The function is only called internally, hence several assumptions can be made:
 * - from < P < to
 * - *pvec[i] < *pvec[from]  for each i < 'from'
 * - *pvec[j] > * pvec[to]   for each j > 'to'
 * - 'pvec' is only rearranged by previous calls of this function and the range is
 *   narrowed down at each iteration
 *
 * @param pvec - vector of pointers to array's elements
 * @param P - the ordinal
 * @param from - lower bound of the range
 * @param to - upper bound of the range
 */
template <typename F>
void __selectRange(
            std::vector<const F*>& pvec,
            const std::size_t P,
            const std::size_t from,
            const std::size_t to )
{
    // nothing to do if P 'P' is not in the range from..to
    if ( from>to || P<from || P>to )
    {
        return;
    }


    /*
     * The algorithm is based on the well known quicksort algorithm.
     * It will partition the current subarray of 'pvec' and narrow it
     * down either to its left or right subarray, depending whether the
     * partition's return value is less or greater than 'P'
     * The procedure will repeat until the returned value equals 'P'.
     */

    std::size_t lo;
    std::size_t hi;

    for ( lo = from, hi = to;
          hi > lo; )
    {
        const std::size_t pivot_pos = math::Selection::__private::__partition(pvec, lo, hi);

        if ( pivot_pos < P )
        {
            // take the right subarray
            lo = pivot_pos + 1;
        }
        else if ( pivot_pos > P )
        {
            // take the left subarray
            hi = pivot_pos - 1;
        }
        else  // if pivot_pos == P
        {
            // Now the ptrTable[P] points to the P.th smallest element
            break;  // out of the for loop
        }
    }
}


/*
 * Finds k1.th and optionally K2.th smallest or largest element of
 * the vector 'x' and assigns them/it to the location(s) pointed by
 * 'a1' and optionally 'a2'.
 *
 * If 'a2' equals NULL only K1.th smallest/largest element will be selected
 * and assigned to the location pointed by 'a1'.
 *
 * The function is only called internally, hence it is assumed that 'a1'
 * will never be NULL.
 *
 * @param x - vector of elements
 * @param K1 - first ordinal (between 0 and size(x)-1)
 * @param K2 - second ordinal (between 0 and size(x)-1), ignored if a2==NULL
 * @param a1 - pointer to assign the value of first ordinal
 * @param a2 - pointer to assign the value of the second ordinal, ignored if NULL
 * @param smallest - if TRUE, K1.th and optionally K2.th smallest elements will be returned, K1.th and optionally K2.th largest otherwise
 *
 * @throw SelectionException if 'K1' or 'K2' is invalid or if internal allocation of memory failed
 */
template <typename F>
void __selectMult(
            const std::vector<F>& x,
            const std::size_t K1,
            const std::size_t K2,
            F* const a1,
            F* const a2,
            const bool smallest
          ) throw(math::SelectionException)
{
    if ( NULL == a1 )
    {
        // should never occur, handle it anyway just in case
        return;
    }

    const std::size_t N = x.size();

    // internal vector of indices, necessary to keep 'x' immutable
    std::vector<const F*> ptrTable;

    // sanity check
    if ( ( K1>=N ) || (NULL!=a2 && K2>=N) )
    {
        throw math::SelectionException(math::SelectionException::ARG_OUT_OF_RANGE);
    }


    // Regardless of 'smallest', the algorithm will always select the P.th smallest value
    const std::size_t P1 = ( true==smallest ? K1 : N-K1-1 );
    const std::size_t P2 = ( true==smallest ? K2 : N-K2-1 );

    // Allocates and fills the vector of indices
    math::Selection::__private::__fillPointers(x, ptrTable);


    /*
     * TODO 
     * Consider random shuffling the array of pointers and thus reduce
     * probability of the worst case scenario when the array is already
     * (almost) sorted.
     */


    // Selects the P1.th smallest element, initially in the whole range
    F& ret1 = *a1;
    math::Selection::__private::__selectRange(ptrTable, P1, 0, N-1);
    ret1 = *(ptrTable.at(P1));

    if ( NULL != a2 )
    {
        // If requested, select the P2.th smallest element as well.
        // At this point we are sure that the element is located in
        // the left or right subrange of 'ptrTable'.
        const std::size_t lb = ( P2>P1 ? P1 : 0 );
        const std::size_t ub = ( P2>P1 ? N-1 : P1 );

        F& ret2 = *a2;
        math::Selection::__private::__selectRange(ptrTable, P2, lb, ub);
        ret2 = *(ptrTable.at(P2));
    }
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

    try
    {
        // allocates 'dest' and fills it with initial positions 0..(N-1)
        dest.resize(N);
        math::util::fillVectorWithConsecutiveIndices(N, dest);
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::SelectionException(math::SelectionException::OUT_OF_MEMORY);
    }

    if ( 0 == N )
    {
        return dest;
    }

    math::Selection::__private::IndexCmp<F> idxCmp(x, asc);

    /*
     * Performs stable sort of indices, assuring that order
     * of indices remains unmodified in case of ties.
     */
    std::stable_sort( dest.begin(), dest.end(), idxCmp );

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
        dest.at( idx.at(0) ) = 0;

        /*
         * For all the remaining elements, the algorithm above is slightly
         * modified. x[idx[i]] is compared to the previously ranked element
         * (x[idx[i-1]]). If they are equal, just copy the previously ranked
         * element's rank (rank[idx[i-1]]). Otherwise assign the element x[idx[i]]
         * the value 'i'.
         */
        for ( std::size_t i = 1; i<N; ++i )
        {
            dest.at(idx.at(i)) =
                ( x.at( idx.at(i-1) ) == x.at( idx.at(i) ) ?
                  dest.at( idx.at(i-1) ) :  i );
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


/**
 * Returns the K.th smallest or largest element of the vector 'x'.
 *
 * On average, the complexity of the implemented algorithm is O(n),
 * where 'n' denotes the size of 'x'.
 *
 * @note The function does not modify the input vector 'x'
 *
 * @param x - vector of elements
 * @param K - the ordinal (between 0 and size(x)-1)
 * @param smallest - if TRUE, K.th smallest element will be returned, K.th largest otherwise (default: TRUE)
 *
 * @return K.th smallest or largest element of 'x', depending on 'smallest'
 *
 * @throw SelectionException if 'K' is invalid or if internal allocation of memory failed
 */
template <typename F>
F math::Selection::select(
            const std::vector<F>& x,
            const std::size_t K,
            const bool smallest
          ) throw (math::SelectionException)
{
    F retVal;
    math::Selection::__private::__selectMult(x, K, K, &retVal, static_cast<F* const>(NULL), smallest);
    return retVal;
}


/**
 * Selects the K1.th and K2.th smallest or largest elements of the vector 'x'
 * and assigns them to variables 'val1' and 'val2', respectively.
 *
 * On average, the complexity of the implemented algorithm is O(n),
 * where 'n' denotes the size of 'x'.
 *
 * @note The function does not modify the input vector 'x'
 *
 * @param x - vector of elements
 * @param K1 - the first ordinal (between 0 and size(x)-1)
 * @param K1 - the second ordinal (between 0 and size(x)-1)
 * @param val1 - reference to the variable to assign the K1.th smallest/largest value of 'x'
 * @param val2 - reference to the variable to assign the K2.th smallest/largest value of 'x'
 * @param smallest - if TRUE, K1.th and K2.th smallest element will be returned, K1.th and K2.th largest otherwise (default: TRUE)
 *
 * @throw SelectionException if 'K1' or 'K2' is invalid or if internal allocation of memory failed
 */
template <typename F>
void math::Selection::select2(
            const std::vector<F>& x,
            const std::size_t K1,
            const std::size_t K2,
            F& val1,
            F& val2,
            const bool smallest
          ) throw (math::SelectionException)
{
    math::Selection::__private::__selectMult(x, K1, K2, &val1, &val2, smallest);
}
