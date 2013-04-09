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
 * @file PermutationGeneric.h
 * 
 * Declaration of the class PermutationGeneric that gradually lists
 * all permutations of a sequence of objects.
 * 
 * @author Jernej Kovacic
 */

#ifndef _MATH_PERMUTATIONGENERIC_H_
#define	_MATH_PERMUTATIONGENERIC_H_

#include "CombinatoricsException.h"

#include <cstddef>
#include <vector>
#include <list>
#include <set>

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
template<class T>
class PermutationGeneric
{
private:
    // Elements will be extensively accessed by their indexes hence
    // vector is the most optimal data structure to store elements to be permutated.
    std::vector<T> elems;
    
    // Current sequence of elements' addresses for the current permutation  
    std::vector<size_t> addr;
    
    // More permutations available?
    bool morePermutations;
    
    // Has listing of permutations already started?
    // (important information for the algorithm)
    bool started;
    
    // Initializes internal data
    void init() throw (CombinatoricsException);
    
public:
    // Constructors
    PermutationGeneric(const std::vector<T>& el) throw (CombinatoricsException);
    PermutationGeneric(const std::list<T>& el) throw (CombinatoricsException);
    PermutationGeneric(const std::set<T>& el) throw (CombinatoricsException);
    // TODO add constructors from other data structures         
    
    // Retrieves next n permutations
    std::list<std::list<T> > next(size_t n=1) throw (CombinatoricsException);
    
    // More permutations available to be retrieved?
    bool hasNext() const;
    
    // Destructor
    virtual ~PermutationGeneric();
};

// Definition could be included into the namespace declaraion, but it
// would cause conflicts with some extra included stdlib header files.

}  // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "PermutationGeneric.cpp"


#endif	// _MATH_PERMUTATIONGENERIC_H_
