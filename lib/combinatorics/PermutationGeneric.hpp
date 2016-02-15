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
 * @headername{PermutationGeneric.h}
 *
 * Declaration of the class PermutationGeneric that gradually lists
 * all permutations of a sequence of objects.
 */

#ifndef _MATH_PERMUTATIONGENERIC_HPP_
#define _MATH_PERMUTATIONGENERIC_HPP_

#include "exception/CombinatoricsException.hpp"

#include <cstddef>
#include <vector>
#include <list>
#include <set>
#include <deque>

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
    
    // Initializes internal data
    void __init() throw (CombinatoricsException);
    
public:
    // Constructors
    PermutationGeneric(const std::vector<T>& el) throw (CombinatoricsException);
    PermutationGeneric(const std::list<T>& el) throw (CombinatoricsException);
    PermutationGeneric(const std::set<T>& el) throw (CombinatoricsException);
    PermutationGeneric(const std::deque<T>& el) throw (CombinatoricsException);
    PermutationGeneric(const T* elarray, const std::size_t len) throw (CombinatoricsException);       
    
    // Retrieves next n permutations
    void next(std::list<std::list<T> >& ret, const std::size_t n=1) throw (CombinatoricsException);
    
    // More permutations available to be retrieved?
    bool hasNext() const;
    
    // Destructor
    virtual ~PermutationGeneric();
};

}  // namespace math

// DEFINITION
#include "combinatorics/PermutationGeneric.cpp"


#endif  // _MATH_PERMUTATIONGENERIC_HPP_
