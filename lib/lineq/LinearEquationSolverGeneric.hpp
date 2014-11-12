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
 * Declaration of the class LinearEquationSolverGeneric that solves
 * systems of linear equations.
 */

#ifndef _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_
#define _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_

#include <complex>

#include "exception/LinearEquationSolverException.hpp"
#include "matrix/MatrixGeneric.hpp"
#include "matrix/SqMatrixGeneric.hpp"

namespace math
{

/**
 * @brief A class with functionality for solving systems of linear equation.
 */
template<class T>
class LinearEquationSolverGeneric
{

private:
    /*
     * After both members have been set, solve() will try to find x,
     * satisfying the equation: m_coef * x = m_term
     */
    SqMatrixGeneric<T> m_coef;   /// matrix of coefficients of the system
    MatrixGeneric<T>   m_term;   /// vector (or matrix) of constant terms

public:
    // Constructors
    LinearEquationSolverGeneric();
    LinearEquationSolverGeneric(const SqMatrixGeneric<T>& coef, const MatrixGeneric<T>& term) throw (LinearEquationSolverException);

    // Getters
    SqMatrixGeneric<T>& getCoef() const;
    MatrixGeneric<T>& getTerm() const;

    // Setters
    LinearEquationSolverGeneric<T>& setCoef(const SqMatrixGeneric<T>& coef) throw (LinearEquationSolverException);
    LinearEquationSolverGeneric<T>& setTerm(const MatrixGeneric<T>& term) throw (LinearEquationSolverException);

    // Solves a system of linear equations
    MatrixGeneric<T> solve() const throw (LinearEquationSolverException);

};

// Equations with elements of types float, double, long double and complex
// make most sense, therefore the following types are predefined
typedef LinearEquationSolverGeneric<float>       FLinearEquationSolver;
typedef LinearEquationSolverGeneric<double>      LinearEquationSolver;
typedef LinearEquationSolverGeneric<long double> LDLinearEquationSolver;

typedef LinearEquationSolverGeneric<std::complex<float> >       FCLinearEquationsolver;
typedef LinearEquationSolverGeneric<std::complex<double> >      CLinearEquationsolver;
typedef LinearEquationSolverGeneric<std::complex<long double> > LDCLinearEquationSolver;

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "lineq/LinearEquationSolverGeneric.cpp"


#endif // _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_
