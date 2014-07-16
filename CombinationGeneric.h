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
 * @file CombinationGeneric.h
 * 
 * Declaration of the class CombinationGeneric that gradually lists
 * all combinations of a sequence of objects.
 * 
 * @author Jernej Kovacic
 */

#ifndef _MATH_COMBINATIONGENERIC_H
#define	_MATH_COMBINATIONGENERIC_H

#include "CombinatoricsException.h"

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
template<class T>
class CombinationGeneric
{

private:
    // Elements will be extensively accessed by their indexes hence
    // vector is the most optimal data structure to store elements to be combined.
    std::vector<T> elems;
    
    // Current sequence of elements' addresses for the current combination
    std::vector<size_t> addr;
    
    // Size of 'elems'
    size_t N_size;
    // Size of a subset (k)
    size_t K;
    
    // More combinations available?
    bool moreCombinations;
    
    // Initializes internal data
    void init();
    
public:
    // Constructors:
    CombinationGeneric(const std::vector<T>& el) throw(CombinatoricsException);
    CombinationGeneric(const std::list<T>& el) throw(CombinatoricsException);
    CombinationGeneric(const std::set<T>& el) throw(CombinatoricsException);
    CombinationGeneric(const std::deque<T>& el) throw (CombinatoricsException);
    CombinationGeneric(const T* elarray, size_t len) throw (CombinatoricsException);
    
    // Getter and setter for size of a subset (k) 
    size_t getK() const;
    void setK(size_t k) throw (CombinatoricsException);
    
    // Retrieves next n combinations
    std::list<std::set<T> > next(size_t n=1) throw (CombinatoricsException);
    
    // More combinations available to be retrieved?
    bool hasNext() const;
    
    // Destructor
    ~CombinationGeneric();
    
};

// Definition could be included into the namespace declaraion, but it
// would cause conflicts with some extra included stdlib header files.

}  // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "CombinationGeneric.cpp"

#endif	// _MATH_COMBINATIONGENERIC_H