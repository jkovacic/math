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
