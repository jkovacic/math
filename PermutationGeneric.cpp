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
 * @file PermutationGeneric.cpp
 * 
 * Implementation of the class PermutationGeneric that gradually lists
 * all permutations of a sequence of objects.
 * 
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 * 
 * @author Jernej Kovacic
 */

// deliberately there is no #include "PermutationGeneric.h" !
#include <new>

#ifdef _OPENMP
#    include <omp.h>
#endif


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
        this->elems = el;
        
        init();
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
        elems.clear();
        elems.reserve(N);
        for ( typename std::list<T>::const_iterator it=el.begin(); it!=el.end(); ++it )
        {
            elems.push_back(*it);
        }
        
        init();
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
        elems.clear();
        elems.reserve(N);
        for ( typename std::set<T>::const_iterator it=el.begin(); it!=el.end(); ++it )
        {
            elems.push_back(*it);
        }
        
        init();
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
        elems.clear();
        elems.reserve(N);
        for ( typename std::deque<T>::const_iterator it=el.begin(); it!=el.end(); ++it )
        {
            elems.push_back(*it);
        }
        
        init();
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
        
        init();
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
void math::PermutationGeneric<T>::init() throw (math::CombinatoricsException)
{
    try
    {
        morePermutations = true;
        started = false;
        N_len = elems.size();
        
        // Initially 'addr' is filled by consecutive integers from 0 to N-1
        addr.clear();
        addr.reserve(N_len);
        for ( size_t i=0; i<N_len; ++i )
        {
            addr.push_back(i);
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_MEMORY);
    }
}

/**
 * Retrieves the next 'n' (or less) permutations.
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
 * @param n - maximum number of permutations to be retrieved (default: 1)
 * 
 * @return a list with max. n lists of permutations
 * 
 * @throw CombinatoricsException if 'n' is invalid or if allocation of memory fails
 */
template<class T>
std::list<std::list<T> > math::PermutationGeneric<T>::next(size_t n) throw (math::CombinatoricsException)
{
    /*
     * The algorithm is based on code, available at
     * https://github.com/silencedrop/permutation
     */
    
    try
    {
        // N_len is used frequently inside this function. As it is intended to
        // remain constant, a const ref. is used to prevent unintentional modifications.
        const size_t& N = N_len;
        std::list<std::list<T> > retVal;
        retVal.clear();
        
        // sanity check
        if ( n>=retVal.max_size() )
        {
            throw math::CombinatoricsException(math::CombinatoricsException::OUT_OF_RANGE);
        }
        
        // At maximum 'n' permutations will be returned
        for ( size_t cnt=0; true==morePermutations && cnt<n; ++cnt )
        {
            /*
             * The first permutation ('elems' in unmodified order) 
             * is not handled by the algorithm so it
             * must be handled separately. 
             */
            if ( false==started )
            {
                std::list<T> temp;
                temp.clear();
                
                // just copy elements of 'elems' into temp:
                for ( typename std::vector<T>::const_iterator it=elems.begin(); it!=elems.end(); ++it )
                {
                    temp.push_back(*it);
                }
                // and append it to retVal:
                retVal.push_back(temp);
                
                // set the flag indicating that the initial permutation
                // has already been returned
                started = true;
                
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
                if ( addr.at(i) > addr.at(i-1) )
                {
                    // Find the smallest 'addr' larger than current one and 
                    // still behind the current one
                    size_t k = N - 1;
                    for ( k=N-1; addr.at(i-1)>addr.at(k); --k );
                    // swap the addresses
                    size_t swap = addr.at(i-1);
                    addr.at(i-1) = addr.at(k);
                    addr.at(k) = swap;
                    
                    // revert the order behind i-1 
                    size_t start;
                    size_t end;
                    for ( start=i, end=N-1; start<end; ++start, --end )
                    {
                        size_t swap = addr.at(start);
                        addr.at(start) = addr.at(end);
                        addr.at(end) = swap;
                    }
                    
                    // and generate the next permutation based on 'addr'
                    std::list<T> temp;
                    temp.clear();
                    
                    for ( size_t j=0; j<N; ++j )
                    {
                        temp.push_back(elems.at(addr.at(j)));
                    }
                    
                    // append it to retVal
                    retVal.push_back(temp);
                    
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
                morePermutations = false;
                // no need to break as "for cnt" checks the value of 'morePermutations'
                //break; // out of for cnt
            }
        }  // for cnt
        
        return retVal;
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
    return morePermutations;
}

/**
 * Destructor
 */
template<class T>
math::PermutationGeneric<T>::~PermutationGeneric()
{
    // probably vector's destructors would clean this automatically...
    elems.clear();
    addr.clear();
}
