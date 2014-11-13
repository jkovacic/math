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
 * Implementation of the class PermutationGeneric that gradually lists
 * all permutations of a sequence of objects.
 * 
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */

// deliberately there is no #include "PermutationGeneric.hpp" !
#include <new>
#include <cstddef>
#include <algorithm>

#include "util/mtcopy.hpp"
#include "exception/CombinatoricsException.hpp"


/**
 * Constructor.
 * 
 * @param el - a vector of elements to be permutated
 * 
 * @throw CombinatoricsException if 'el' is empty or allocation of memory fails
 */
template<class T>
math::PermutationGeneric<T>::PermutationGeneric(const std::vector<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
        
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
template<class T>
math::PermutationGeneric<T>::PermutationGeneric(const std::list<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
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
template<class T>
math::PermutationGeneric<T>::PermutationGeneric(const std::set<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
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
template<class T>
math::PermutationGeneric<T>::PermutationGeneric(const std::deque<T>& el) throw (math::CombinatoricsException)
{
    try
    {
        const size_t N = el.size();
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
template<class T>
math::PermutationGeneric<T>::PermutationGeneric(const T* elarray, size_t len) throw (math::CombinatoricsException)
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
    
/*
 * Initialization of class's internal data to initial values
 * 
 * @throw CombinatoricsException if allocation of memory fails
 */
template<class T>
void math::PermutationGeneric<T>::__init() throw (math::CombinatoricsException)
{
    try
    {
        this->morePermutations = true;
        this->started = false;
        this->N_len = this->elems.size();
        
        // Initially 'addr' is filled by consecutive integers from 0 to N-1
        this->addr.clear();
        this->addr.reserve( this->N_len );
        for ( size_t i=0; i<this->N_len; ++i )
        {
            this->addr.push_back(i);
        }
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
template<class T>
void math::PermutationGeneric<T>::next(std::list<std::list<T> >& ret, size_t n) throw (math::CombinatoricsException)
{
    /*
     * The algorithm is based on code, available at
     * https://github.com/silencedrop/permutation
     */
    
    try
    {
        // N_len is used frequently inside this function. As it is intended to
        // remain constant, a const ref. is used to prevent unintentional modifications.
        const size_t& N = this->N_len;

        // sanity check
        if ( n>=ret.max_size() )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }

        // clear the return list:
        ret.clear();
        
        // At maximum 'n' permutations will be returned
        for ( size_t cnt=0; true==this->morePermutations && cnt<n; ++cnt )
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
                    l.push_back(*it);
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
            for ( size_t i=N-1; i>0; --i )
            {
                // if the 'addr' is greater than the previous one
                if ( this->addr.at(i) > this->addr.at(i-1) )
                {
                    // Find the smallest 'addr' larger than current one and 
                    // still behind the current one
                    size_t k = N - 1;
                    for ( k=N-1; this->addr.at(i-1)>this->addr.at(k); --k );
                    // swap the addresses
                    size_t swap = this->addr.at(i-1);
                    this->addr.at(i-1) = this->addr.at(k);
                    this->addr.at(k) = swap;
                    
                    // revert the order behind i-1 
                    size_t start;
                    size_t end;
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
                    
                    for ( size_t j=0; j<N; ++j )
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
template<class T>
bool math::PermutationGeneric<T>::hasNext() const
{
    return this->morePermutations;
}

/**
 * Destructor
 */
template<class T>
math::PermutationGeneric<T>::~PermutationGeneric()
{
    // probably vector's destructors would clean this automatically...
    this->elems.clear();
    this->addr.clear();
}
