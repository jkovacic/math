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
 * @headername{LinearEquationSolverGeneric.h}
 *
 * Declaration of functions within namespace LinearEquationSolver for
 * solving systems of linear equations.
 */

#ifndef _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_
#define _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_

#include "../settings/lineq_settings.h"
#include "exception/MatrixException.hpp"
#include "matrix/MatrixGeneric.hpp"

namespace math
{

/**
 * @brief A namespace with functions for solving systems of linear equations.
 */
namespace LinearEquationSolver
{

    template <class T>
    bool solveGaussJordan(
              const MatrixGeneric<T>& coef,
              const MatrixGeneric<T>& term,
              MatrixGeneric<T>& sol,
              const bool fullp = true
            ) throw (MatrixException);


    template <class T>
    bool solveSOR(
              const MatrixGeneric<T>& coef,
              const MatrixGeneric<T>& term,
              MatrixGeneric<T>& sol,
              const T& w,
              const T& tol = static_cast<T>(LINEQ_TOL_CONV_NUM) / static_cast<T>(LINEQ_TOL_CONV_DEN),
              const size_t maxiter = LINEQ_MAX_ITER
            ) throw (MatrixException);


    template <class T>
    bool solveGaussSeidel(
              const MatrixGeneric<T>& coef,
              const MatrixGeneric<T>& term,
              MatrixGeneric<T>& sol,
              const T& tol = static_cast<T>(LINEQ_TOL_CONV_NUM) / static_cast<T>(LINEQ_TOL_CONV_DEN),
              const size_t maxiter = LINEQ_MAX_ITER
            ) throw (MatrixException);


}  // namespace LinearEquationSolver

}  // namespace math

// DEFINITION
#include "matrix/LinearEquationSolverGeneric.cpp"


#endif // _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_
