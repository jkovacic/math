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
 * An internal header file, it should not be included directly.
 * @headername{Combination.h}
 *
 * Declaration and implementation of the class CombinationGeneric that
 * gradually lists all combinations of a sequence of objects.
 */

#ifndef _MATH_COMBINATIONGENERIC_HPP_
#define _MATH_COMBINATIONGENERIC_HPP_

#include "util/mtcopy.hpp"
#include "exception/CombinatoricsException.hpp"

#include <new>
#include <cstddef>
#include <vector>
#include <list>
#include <set>
#include <deque>


namespace math
{

/**
 * @brief A class that gradually lists all k-combinations of a sequence of objects.
 * 
 * The sequence must be passed to one of the class's constructors. Since the number
 * of all k-combinations can be very high (even larger than theoretically allowed size 
 * of STL containers), possibly resulting in failures of memory allocation, it is a 
 * better idea to instantiate a class once, set the size of a subset (K) and to 
 * gradually retrieve a certain number of k-combination sets by calling next(). The 
 * class is stateful, i.e. once a certain k-combination has been retrieved, it cannot 
 * be retrieved again unless the class is reinstantiated or setK() is called which also
 * resets its internal states.
 */
template <class T>
class CombinationGeneric
{

private:
    // Elements will be extensively accessed by their indexes hence
    // vector is the most optimal data structure to store elements to be combined.
    std::vector<T> elems;
    
    // Current sequence of elements' addresses for the current combination
    std::vector<std::size_t> addr;
    
    // Size of 'elems'
    std::size_t N_size;
    // Size of a subset (k)
    std::size_t K;
    
    // More combinations available?
    bool moreCombinations;
    
    /*
     * Initialization of class's internal data to initial values
     */
    void __init()
    {
        this->N_size = this->elems.size();
        this->addr.clear();
        this->K = 0;
        this->moreCombinations = false;
    }
    
public:

    /**
     * Constructor.
     * 
     * @param el - a vector of elements to be combined
     * 
     * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
     */
    CombinationGeneric(const std::vector<T>& el)
    {
        try
        {
            const std::size_t N = el.size();
        
            if ( 0==N )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            // copy the vector 'el' into 'elems'
            math::mtcopy(el, this->elems);

            this->__init();
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
    CombinationGeneric(const std::list<T>& el)
    {
        try
        {
            const std::size_t N = el.size();
            if ( 0==N )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            // copy elements of 'el' into 'elems'
            this->elems.clear();
            this->elems.reserve(N);
            this->elems.assign(el.begin(), el.end());

            this->__init();
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
    CombinationGeneric(const std::set<T>& el)
    {
        try
        {
            const std::size_t N = el.size();
            if ( 0==N )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            // copy elements of 'el' into 'elems'
            this->elems.clear();
            this->elems.reserve(N);
            this->elems.assign(el.begin(), el.end());

            this->__init();
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
    CombinationGeneric(const std::deque<T>& el)
    {
        try
        {
            const std::size_t N = el.size();
            if ( 0==N )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            // copy elements of 'el' into 'elems'
            this->elems.clear();
            this->elems.reserve(N);
            this->elems.assign(el.begin(), el.end());

            this->__init();
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
        }
    }


    /**
     *  Constructor.
     * 
     * @note Usage of this constructor may be dangerous if 'len' is larger then the actual size
     * of the array. Use it carefully and double check both input arguments!
     * 
     * @param elarray - array of elements
     * @param len - number of elements in the array
     * 
     * @throw CombinatoricsException if input arguments are invalid or if allocation of memory fails
     */
    CombinationGeneric(const T* elarray, const std::size_t len)
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

            math::mtcopy(elarray, len, this->elems);

            this->__init();
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
        }
    }
    
    
    /**
     * @return size of a subset (k)
     */
    std::size_t getK() const
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
    void setK(const std::size_t k)
    {
        try
        {
            const std::size_t& N = this->N_size;
        
            // "sanity check"
            if ( k<=0 || k>N )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            this->__init();
        
            /*
             * Prepare 'addr' for the first k-combination:
             * addr[0]=0, addr[1]=1, ..., addr[k-1]=k-1
             * 
             * For more details of the algorithm, see its detailed description
             * at comments in next().
             */
            this->addr.reserve(k);
            for ( std::size_t i=0; i<k; ++i )
            {
                this->addr.push_back(i);
            }
        
            this->K = k;
            this->moreCombinations = true;
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
     * @param ret - a list to be filled with max. n sets of k-combinations
     * @param n - maximum number of combinations to be retrieved (default: 1)
     * 
     * @throw CombinatoriscException if n is too large or if allocation of memory fails
     */
    void next(std::list<std::set<T> >& ret, const std::size_t n=1)
    {
        try
        {
            // K is often used within this function. As it is supposed to
            // remain constant, protect it from unintentional modifications. 
            const std::size_t& Klen = this->K;
            const std::size_t& N = this->N_size;

            // sanity check
            if ( Klen<=0 || Klen>N )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            if ( n>=ret.max_size() )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
            }

            // Clear the return list
            ret.clear();

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
        
            for ( std::size_t cnt=0; true==this->moreCombinations && cnt<n; ++cnt )
            {
                // append an empty set to 'ret'
                ret.push_back(std::set<T>());
                // obtain a reference to this appended set
                std::set<T>& s = ret.back();
                // and clear it (just in case)
                s.clear();
            
                // 'addr' has been set by the previous iteration for the current one.
                // Just append elements at addr's indexes to ret's set:
            
                for ( std::vector<std::size_t>::const_iterator it=addr.begin(); it!=addr.end(); ++it )
                {
                    const std::size_t& currAddr = *it;

                   s.insert(this->elems.at(currAddr));
                }

                // And update 'addr' for the next iteration as described above:
                for ( std::size_t i = Klen-1; /* i>=0 */ ; --i )
                {
                    if ( this->addr.at(i) < (N-Klen+i) )
                    {
                        // limit for the ith element not reached yet:
                        ++this->addr.at(i);
                    
                        // update the next elements to preserve uniqueness and ascending order: 
                        const std::size_t upv = this->addr.at(i); 
                        for ( std::size_t j = i+1; j<Klen; ++j )
                        {
                            this->addr.at(j) = upv + j-i;
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
                        this->moreCombinations = false;
                        break;  // out of for i
                    }
                } // for i
            }  // for cnt

        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
        }
    }
    
    
    /**
     * @return whether more combinations are available to retrieve
     */
    bool hasNext() const
    {
        return this->moreCombinations;
    }
    
    
    /**
     * Destructor
     */
    ~CombinationGeneric()
    {
        // probably vector's destructors would clean this automatically...
        this->elems.clear();
        this->addr.clear();
    }
    
};

}  // namespace math


#endif  // _MATH_COMBINATIONGENERIC_HPP_
