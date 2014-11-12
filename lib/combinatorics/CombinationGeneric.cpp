/*
Copyright 2013, Jernej Kovacic

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
 * Implementation of the class CombinationGeneric that gradually lists
 * all combinations of a sequence of objects.
 * 
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */

// deliberately there is no #include "CombinationGeneric.hpp" !
#include <new>
#include <cstddef>

#include "exception/CombinatoricsException.hpp"


/**
 * Constructor.
 * 
 * @param el - a vector of elements to be combined
 * 
 * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
 */
template<class T>
math::CombinationGeneric<T>::CombinationGeneric(const std::vector<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
        
        if ( 0==N )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        // copy the vector 'el' into 'elems'
        this->elems = el;
        
        __init();
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
    
}

/**
 * Constructor.
 * 
 * @param el - a list of elements to be combined
 * 
 * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
 */
template<class T>
math::CombinationGeneric<T>::CombinationGeneric(const std::list<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
        if ( 0==N )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        // copy elements of 'el' into 'elems'
        elems.clear();
        elems.reserve(N);
        for ( typename std::list<T>::const_iterator it=el.begin(); it!=el.end(); ++it )
        {
            elems.push_back(*it);
        }
        
        __init();
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}
    
/**
 * Constructor.
 * 
 * @param el - a set of elements to be combined
 * 
 * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
 */
template<class T>
math::CombinationGeneric<T>::CombinationGeneric(const std::set<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
        if ( 0==N )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        // copy elements of 'el' into 'elems'
        elems.clear();
        elems.reserve(N);
        for ( typename std::set<T>::const_iterator it=el.begin(); it!=el.end(); ++it )
        {
            elems.push_back(*it);
        }
        
        __init();
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}

/**
 * Constructor.
 * 
 * @param el - a deque of elements to be combined
 * 
 * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
 */
template<class T>
math::CombinationGeneric<T>::CombinationGeneric(const std::deque<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
        if ( 0==N )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        // copy elements of 'el' into 'elems'
        elems.clear();
        elems.reserve(N);
        for ( typename std::deque<T>::const_iterator it=el.begin(); it!=el.end(); ++it )
        {
            elems.push_back(*it);
        }
        
        __init();
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}

/**
 *  Constructor.
 * 
 * @note Usage of this constructor may be dangerous if 'len' is larger then the actual size of the array.
 *       Use it carefully and double check both input arguments!
 * 
 * @param elarray - array of elements
 * @param len - number of elements in the array
 * 
 * @throw CombinatoricsException if input arguments are invalid or if allocation of memory fails
 */
template<class T>
math::CombinationGeneric<T>::CombinationGeneric(const T* elarray, size_t len) throw (math::CombinatoricsException)
{
    try
    {
        // sanity check
        if ( len<=0 || NULL==elarray )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        if ( len>elems.max_size() )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }
        
        elems.clear();
        elems.reserve(len);
        for ( size_t i=0; i<len; ++i )
        {
            elems.push_back(elarray[i]);
        }
        
        __init();
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}

/*
 * Initialization of class's internal data to initial values
 */
template<class T>
void math::CombinationGeneric<T>::__init()
{
    N_size = elems.size();
    addr.clear();
    K = 0;
    moreCombinations = false;
}

/**
 * @return size of a subset (k)
 */
template<class T>
size_t math::CombinationGeneric<T>::getK() const
{
    return this->K;
}

/**
 * Resets the class's internal states and sets anew value for k (size of a subset).
 * Set of elements (set by any of constructors) remains unaffected.
 * 
 * 'k' must be greater than 0 and less than or equal to number of all elements:
 *     0 < k <= N
 * 
 * @param k - size of a subset
 * 
 * @throw CombinatoricsException if input argument is invalid or if allocation of memory fails
 */
template<class T>
void math::CombinationGeneric<T>::setK(size_t k) throw (math::CombinatoricsException)
{
    try
    {
        const size_t& N = N_size;
        
        // "sanity check"
        if ( k<=0 || k>N )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        __init();
        
        /*
         * Prepare 'addr' for the first k-combination:
         * addr[0]=0, addr[1]=1, ..., addr[k-1]=k-1
         * 
         * For more details of the algorithm, see its detailed description
         * at comments in next().
         */
        addr.reserve(k);
        for ( size_t i=0; i<k; ++i )
        {
            addr.push_back(i);
        }
        
        K = k;
        moreCombinations = true;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}

/**
 * Retrieves the next 'n' (or less) k-combinations.
 * 
 * The class is stateful, i.e. only those k-combinations are returned that have not
 * been returned by previous calls of next() on the same instance of the class.
 * 
 * Order of k-combinations is deterministic. Sequences of indexes are lexicographically 
 * sorted in ascending order, for instance for (5 choose 3):
 * 
 *      a[1]a[2]a[3] -> a[1]a[2]a[4] -> a[1]a[2]a[5] -> a[1]a[3]a[4] -> a[1]a[3]a[5] ->
 *   -> a[1]a[4]a[5] -> a[2]a[3]a[4] -> a[2]a[3]a[5] -> a[2]a[4]a[5] -> a[3]a[4]a[5]
 * 
 * @param n - maximum number of combinations to be retrieved (default: 1)
 * 
 * @return a list with max 'n' sets of k-combinations 
 * 
 * @throw CombinatoriscException if n is too large or if allocation of memory fails
 */
template<class T>
std::list<std::set<T> > math::CombinationGeneric<T>::next(size_t n) throw (math::CombinatoricsException)
{
    try
    {
        // K is often used within this function. As it is supposed to
        // remain constant, protect it from unintentional modifications. 
        const size_t& Klen = this->K;
        const size_t& N = N_size;
        std::list<std::set<T> > retVal;
        retVal.clear();

        // sanity check
        if ( Klen<=0 || Klen>N )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
        }
        
        if ( n>=retVal.max_size() )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }
        
        /* 
         * Description of the algorithm:
         * 
         * Indexes are stored in K elements of 'addr', acting as a sort of a "counter". 
         * As the order of elements is not important, indexes can be always "sorted"
         * in ascending order:
         *   addr[i] < addr[j] <== i<j
         * 
         * Similarly as at an ordinary mechanical counter, the rightmost index is 
         * increased by 1. If it exceeds a limit (N), its left "neighbour" is 
         * increased, etc. As indexes are unique (may appear not more than once in 
         * 'addr') and sorted, elements at each position have different limits:
         * N-1 for the rightmost element, N-2 for its left neighbour, N-3 for
         * its left neighbour, etc.
         */
        
        for ( size_t cnt=0; true==moreCombinations && cnt<n; ++cnt )
        {
            // 'addr' has been set by the previous iteration for the current one.
            // Just create a set wuth elements at addr's indexes and append it to retVal:
            std::set<T> temp;
            temp.clear();
        
            for ( std::vector<size_t>::const_iterator it=addr.begin(); it!=addr.end(); ++it )
            {
                temp.insert(elems.at(*it));
            }
            
            retVal.push_back(temp);
            
            // And update 'addr' for the next iteration as described above:
            for ( size_t i = Klen-1; /* i>=0 */ ; --i )
            {
                if ( addr.at(i) < (N-Klen+i) )
                {
                    // limit for the ith element not reached yet:
                    ++addr.at(i);
                    
                    // update the next elements to preserve uniqueness and ascending order: 
                    const size_t upv = addr.at(i); 
                    for ( size_t j = i+1; j<Klen; ++j )
                    {
                        addr.at(j) = upv + j-i;
                    }
                    
                    /* 
                     * If at least one element was modified, more combinations are 
                     * available to be retrieved by next iterations. In this case it 
                     * is important that the next if statement is skipped by 'break'.
                     */
                    break; // out of for i
                }  // if addr[i] < N-Klen+i
                
                // If no element of 'addr' was modified, all combinations
                // have been retrieved. This condition is flagged:
                if ( 0==i )
                {
                    moreCombinations = false;
                    break;  // out of for i
                }
            } // for i
        }  // for cnt
        
        return retVal;
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}

/**
 * @return whether more combinations are available to retrieve
 */
template<class T>
bool math::CombinationGeneric<T>::hasNext() const
{
    return moreCombinations;
}

/**
 * Destructor
 */
template<class T>
math::CombinationGeneric<T>::~CombinationGeneric()
{
    // probably vector's destructors would clean this automatically...
    elems.clear();
    addr.clear();
}
