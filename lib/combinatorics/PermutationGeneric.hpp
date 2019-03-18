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
 * @headername{Permutation.h}
 *
 * Declaration and implementation of the class PermutationGeneric that
 * gradually lists all permutations of a sequence of objects.
 */

#ifndef _MATH_PERMUTATIONGENERIC_HPP_
#define _MATH_PERMUTATIONGENERIC_HPP_

#include "util/mtcopy.hpp"
#include "exception/CombinatoricsException.hpp"

#include <new>
#include <cstddef>
#include <vector>
#include <list>
#include <set>
#include <deque>
#include <algorithm>


namespace math
{

/**
 * @brief A class that gradually lists all permutations of a sequence of objects.
 * 
 * The sequence must be passed to one of the class's constructors. Since the number
 * of all permutations can be very high (even larger than theoretically allowed size 
 * of STL containers), possibly resulting in failures of memory allocation, it is a 
 * better idea to instantiate a class once and to gradually retrieve a certain number 
 * of permutation lists by calling next(). The class is stateful and "irreversible",
 * i.e. once a certain permutation has been retrieved, it cannot be retrieved again
 * unless the class is reinstantiated. 
 */
template <class T>
class PermutationGeneric
{
private:
    // Elements will be extensively accessed by their indexes hence
    // vector is the most optimal data structure to store elements to be permutated.
    std::vector<T> elems;
    
    // Current sequence of elements' addresses for the current permutation  
    std::vector<std::size_t> addr;
    
    // Size of elems
    std::size_t N_len;
    
    // More permutations available?
    bool morePermutations;
    
    // Has listing of permutations already started?
    // (important information for the algorithm)
    bool started;
    
    
    /*
     * Initialization of class's internal data to initial values
     * 
     * @throw CombinatoricsException if allocation of memory fails
    */
    void __init()
    {
        try
        {
            this->morePermutations = true;
            this->started = false;
            this->N_len = this->elems.size();
        
            // Initially 'addr' is filled by consecutive integers from 0 to N-1
            this->addr.clear();
            this->addr.reserve( this->N_len );
            for ( std::size_t i=0; i<this->N_len; ++i )
            {
                this->addr.push_back(i);
            }
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
        }
    }
    
public:

    /**
     * Constructor.
     * 
     * @param el - a vector of elements to be permutated
     * 
     * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
     */
    PermutationGeneric(const std::vector<T>& el)
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
     * @param el - a list of elements to be permutated
     * 
     * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
     */
    PermutationGeneric(const std::list<T>& el)
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
     * @param el - a set of elements to be permutated
     * 
     * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
     */
    PermutationGeneric(const std::set<T>& el)
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
     * @param el - a deque of elements to be permutated
     * 
     * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
     */
    PermutationGeneric(const std::deque<T>& el)
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
     * @note Usage of this constructor may be dangerous if 'len' is larger then the actual size of the array.
     *       Use it carefully and double check both input arguments!
     * 
     * @param elarray - array of elements
     * @param len - number of elements in the array
     * 
     * @throw CombinatoricsException if input arguments are invalid or if allocation of memory fails
     */
    PermutationGeneric(const T* elarray, const std::size_t len)
    {
        try
        {
            // sanity check
            if ( len<=0 || NULL==elarray )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::INVALID_INPUT);
            }
        
            if ( len>this->elems.max_size() )
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
     * Retrieves the next 'n' (or less) permutations and pushes them into the
     * specified list 'ret'.
     * 
     * The class is stateful, i.e. only those permutations are returned that have not
     * been returned by previous calls of next() on the same instance of the class.
     * 
     * Order of permutations is deterministic. Sequences of indexes are lexicographically 
     * sorted in ascending order, for instance:
     * 
     *     a(1)a(2)a(3) -> a(1)a(3)a(2) -> a(2)a(1)a(3) ->
     *  -> a(2)a(3)a(1) -> a(3)a(1)a(2) -> a(3)a(2)a(1)
     * 
     * @param ret - a list to be filled with max. n lists of permutations
     * @param n - maximum number of permutations to be retrieved (default: 1)
     * 
     * @throw CombinatoricsException if 'n' is invalid or if allocation of memory fails
     */
    void next(std::list<std::list<T> >& ret, const std::size_t n=1)
    {
        /*
         * The algorithm is based on code, available at
         * https://github.com/silencedrop/permutation
         */
    
        try
        {
            // N_len is used frequently inside this function. As it is intended to
            // remain constant, a const ref. is used to prevent unintentional modifications.
            const std::size_t& N = this->N_len;

            // sanity check
            if ( n>=ret.max_size() )
            {
                throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
            }

            // clear the return list:
            ret.clear();
        
            // At maximum 'n' permutations will be returned
            for ( std::size_t cnt=0; true==this->morePermutations && cnt<n; ++cnt )
            {
                /*
                 * The first permutation ('elems' in unmodified order) 
                 * is not handled by the algorithm so it
                 * must be handled separately. 
                 */
                if ( false==this->started )
                {
                    // append an empty list to 'ret':
                    ret.push_back(std::list<T>());
                    // and obtain a reference to this appended list
                    std::list<T>& l = ret.back();
                    // clear it (just in case)
                    l.clear();
                
                    // just copy elements of 'elems' into 'l':
                    for ( typename std::vector<T>::const_iterator it=this->elems.begin(); it!=this->elems.end(); ++it )
                    {
                        const T& currElem = *it;
                        l.push_back(currElem);
                    }

                    // set the flag indicating that the initial permutation
                    // has already been returned
                    this->started = true;
                
                    // skip the rest of the loop to increment 'cnt'
                    continue;  // for cnt
                }
            
                /*
                 * Unlike in Java, C++ does not support so called "named loops", i.e.
                 * it is not possible to break out of multiple nested loops 
                 * with a single 'break' statement (an alternative would be 'goto'). 
                 * For that reason, this flag has been introduced, indicating that
                 * additional permutations are still possible. 
                 */ 
                bool pfound = false;
            
                // Permutate from the current 'addr' to the smallest possible left
                for ( std::size_t i=N-1; i>0; --i )
                {
                    // if the 'addr' is greater than the previous one
                    if ( this->addr.at(i) > this->addr.at(i-1) )
                    {
                        // Find the smallest 'addr' larger than current one and 
                        // still behind the current one
                        std::size_t k = N - 1;
                        for ( k=N-1; this->addr.at(i-1)>this->addr.at(k); --k );
                        // swap the addresses
                        std::size_t swap = this->addr.at(i-1);
                        this->addr.at(i-1) = this->addr.at(k);
                        this->addr.at(k) = swap;
                    
                        // revert the order behind i-1 
                        std::size_t start;
                        std::size_t end;
                        for ( start=i, end=N-1; start<end; ++start, --end )
                        {
                            std::swap(this->addr.at(start), this->addr.at(end));
                        }
                    
                        // and generate the next permutation based on 'addr':

                        // append an empty list into 'ret':
                        ret.push_back(std::list<T>());
                        // obtain a reference to this appended list:
                        std::list<T>& l = ret.back();
                        // clear it (just in case):
                        l.clear();
                    
                        for ( std::size_t j=0; j<N; ++j )
                        {
                            l.push_back(this->elems.at(this->addr.at(j)));
                        }
                    
                        // set a flag indicating that a permutation has been found
                        pfound = true;
                        // No need to iterate 'i' further for the current permutation (cnt) 
                        break;  // out of for i
                    
                    }  // if addr[i] > addr[i-1]
        
                } // for i
            
                // false indicates that all possible permutations have already been handled
                if ( false==pfound )
                {
                    // mark the class's flag
                    this->morePermutations = false;
                    // no need to break as "for cnt" checks the value of 'morePermutations'
                    //break; // out of for cnt
                }
            }  // for cnt

        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
        }
    }
    
    
    /**
     * @return whether more permutations are available to retrieve
     */
    bool hasNext() const
    {
        return this->morePermutations;
    }
    
    
    /**
     * Destructor
     */
    virtual ~PermutationGeneric()
    {
        // probably vector's destructors would clean this automatically...
        this->elems.clear();
        this->addr.clear();
    }
    
};

}  // namespace math


#endif  // _MATH_PERMUTATIONGENERIC_HPP_
