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
 * An internal header file, it should not be included directly.
 *
 * Declaration of auxiliary functions within namespace Pivot that
 * perform partial or full pivoting.
 */

#ifndef _MATH_PIVOTGENERIC_HPP_
#define _MATH_PIVOTGENERIC_HPP_

#include <cstddef>
#include <vector>

#include "matrix/MatrixGeneric.hpp"
#include "exception/MatrixException.hpp"


namespace math
{

namespace Pivot
{


void fillVectorsWithInitialPos(
        const size_t N,
        std::vector<size_t>& v1,
        std::vector<size_t>& v2,
        const bool bothVectors
      ) throw(MatrixException);


template <class T>
bool solveGaussJordan(
        const MatrixGeneric<T>& A,
        const MatrixGeneric<T>& b,
        MatrixGeneric<T>& x,
        const bool fullp
      ) throw(MatrixException);


template <class T>
T getDeterminant(const MatrixGeneric<T>& A, const bool fullp) throw(MatrixException);


template <class T>
size_t getRank(const MatrixGeneric<T>& A) throw(MatrixException);

}  // namespace Pivot

}  // namepace math


// DEFINITION
#include "matrix/PivotGeneric.cpp"


#endif	// _MATH_PIVOTGENERIC_HPP_
