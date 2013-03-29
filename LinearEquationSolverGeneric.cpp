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
 * @file LinearEquationSolverGeneric.cpp
 *
 * Implementation of the class LinearEquationSolverGeneric.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 *
 * @author Jernej Kovacic
 */

// deliberately there is no #include "LinearEquationSolverGeneric.h" !
#include "MatrixException.h"
#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"


/**
 * Constructor.
 * 
 * No coefficients or terms are set. They must be set using setCoef() and
 * SetTerm() before solve() can be called.
 */
template<class T>
math::LinearEquationSolverGeneric<T>::LinearEquationSolverGeneric()
{
    // Nothing to do
}

/**
 * Constructor that assigns values of system's coefficients and constant terms.
 * 
 * @note as it is possible to modify matrices of coefficients and/or terms later,
 *       their dimensions are not checked by this constructor.
 * 
 * @param coef - coefficients of the system of linear equations
 * @param term - constant terms of the system of linear equations
 * 
 * @throw LinearEquationSolverException if allocation of memory fails
 */
template<class T>
math::LinearEquationSolverGeneric<T>::LinearEquationSolverGeneric(const SqMatrixGeneric<T>& coef, const MatrixGeneric<T>& term) throw (math::LinearEquationSolverException)
{
    try
    {
        // just copy both matrices into internal members
        m_coef = coef;
        m_term = term;
    }
    catch (const math::MatrixException& mex )
    {
        // exception is only thrown when memory allocation fails
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::OUT_OF_MEMORY);
    }
}

/**
 * @return reference to the internal matrix of coefficients
 */
template<class T>
math::SqMatrixGeneric<T>& math::LinearEquationSolverGeneric<T>::getCoef() const
{
    return m_coef;
}

/**
 * @return reference to the internal vector or matrix of constant terms
 */
template<class T>
math::MatrixGeneric<T>& math::LinearEquationSolverGeneric<T>::getTerm() const
{
    return m_term;
}

/**
 * Sets a square matrix with coefficients of the system of linear equations.
 * 
 * @note Dimensions of matrices are not checked by this method.
 * 
 * @param coef - square matrix with coefficients of the system of linear equations
 * 
 * @return reference to itself
 * 
 * @throw LinearEquationSolverException if allocation of memory fails
 */
template<class T>
math::LinearEquationSolverGeneric<T>& math::LinearEquationSolverGeneric<T>::setCoef(const math::SqMatrixGeneric<T>& coef) throw (math::LinearEquationSolverException)
{
    try
    {
        // just copy the matrix into the internal member
        m_coef = coef;
    }
    catch ( const math::MatrixException& mex )
    {
        // exception is only thrown when memory allocation fails
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::OUT_OF_MEMORY);
    }

    return *this;
}

/**
 * Sets a vector or a matrix with constant terms of the system of linear equations.
 * 
 * @note Dimensions of matrices are not checked by this method.
 * 
 * @param term - vector or matrix with constant terms of the system of linear equations 
 * 
 * @return reference to itself
 * 
 * @throw LinearEquationSolverException if allocation of memory fails
 */
template<class T>
math::LinearEquationSolverGeneric<T>& math::LinearEquationSolverGeneric<T>::setTerm(const math::MatrixGeneric<T>& term) throw (math::LinearEquationSolverException)
{
    try
    {
        // just copy the matrix into the internal member
        m_term = term;
    }
    catch ( const math::MatrixException& mex )
    {
        // exception is only thrown when memory allocation fails
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::OUT_OF_MEMORY);
    }

    return *this;
}

/**
 * Solves the system of linear equations and returns its unique solution if it exists.
 * Matrices of coefficients and terms must be set beforehand.
 * 
 * Number of coef's columns must be equal to the number of term's rows.
 * 
 * If unique solution does not exist (i.e. determinant of 'coef' is 0), an
 * exception will be thrown.
 * 
 * The method does not modify internal matrices (coefficients and constant terms)
 * so the same instance of the class can be reused as many times as desired.
 * 
 * @return A vector or a matrix 'x' that satisfies the condition: COEF * x = TERM
 * 
 * @throw LinearEquationSolverException if a unique solution cannot be found for any reason
 */
template<class T>
math::MatrixGeneric<T> math::LinearEquationSolverGeneric<T>::solve() const throw (math::LinearEquationSolverException)
{
    /*
     * The Gaussian elimination algorithm is implemented:
     * multiples of coef's and term's lines are added to other lines until
     * "coef" appears as a unit matrix. In this case the modified "term" is a
     * unique solution of a system of linear equations. More details about
     * the algorithm at: http://en.wikipedia.org/wiki/Gaussian_elimination
     */
    const size_t N = m_coef.nrColumns();  // Nr. of unknowns
    const size_t NT = m_term.nrColumns(); // Nr. of terms' columns
    const size_t Nmax = (N>=NT ? N : NT); // max. of both values

    // Check of dimensions
    if ( N!=m_term.nrRows() )
    {
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::INVALID_DIMENSION);
    }

    try
    {
        SqMatrixGeneric<T> temp(m_coef);
        MatrixGeneric<T> retVal(m_term);

        size_t r;
        T el;

        // Try to convert the 'temp' into an identity matrix
        // by appropriate adding multiples of other lines to each line
        // (incl. lines of 'retVal')

        for ( size_t i=0; i<N; i++ )
        {
            // first check if the diagonal element equals 0
            if ( true == math::NumericUtil<T>::isZero(temp.at(i, i)) )
            {
                // if it does, try to find another row r where temp(r,i)!=0
                for ( r=0; r<N; r++ )
                {
                    if ( r==i )
                    {
                        // it is known in advance, that temp(i,i)==0, so skip it
                        continue;  // for r
                    }

                    if ( false == math::NumericUtil<T>::isZero(temp.at(r, i)) )
                    {
                        // found, no need to search further
                        break;  // out of for r
                    }
                }  // for r

                if ( N==r )
                {
                    // No temp(r,i)!=0 was found, the matrix 'temp' is non-invertible.
                    // Throw an exception
                    throw math::LinearEquationSolverException(math::LinearEquationSolverException::NO_UNIQUE_SOLUTION);
                }

                // add the r^th line to the i^th one:
                for ( size_t c=0; c<Nmax; c++ )
                {
                    if ( c<N )
                    {
                        temp.at(i, c) += temp.at(r, c);
                    }

                    if ( c<NT )
                    {
                        retVal.at(i, c) += retVal.at(r, c);
                    }
                }

            }  // if temp(i,i)==0

            // Let the diag element be 1. So divide the whole row by temp(i,i)
            //  (columns smaller than i are already 0)
            el = temp.get(i, i);

            for ( size_t c=i; c<N; c++)
            {
                temp.at(i, c) /= el;
            }

            for ( size_t c=0; c<NT; c++ )
            {
                retVal.at(i, c) /= el;
            }


            // set the i^th column of all other rows (r>i) to 0 by
            // adding the appropriate multiple of the i^th row
            for ( r=i+1; r<N; r++ )
            {
                // Nothing to do if temp(r,i) is already 0.
                if ( true == math::NumericUtil<T>::isZero(temp.at(r, i)) )
                {
                    continue;  // for r
                }

                // Subtract a multiple of the i^th row. Note that temp(i,i) is already 1.
                el = temp.get(r, i);

                for ( size_t c=i; c<N; c++ )
                {
                    temp.at(r, c) -= el*temp.at(i, c);
                }

                for ( size_t c=0; c<NT; c++ )
                {
                    retVal.at(r, c) -= el*retVal.at(i, c);
                }
            }  // for r
        }  // for i

        // Now the lower triangle (below diag excl.) is 0, the diagonal consists of 1,
        // The upper triangle (above the diag) must be set to 0 as well.

        for ( r=0; r<N; r++ )
        {
            for ( size_t c=r+1; c<N; c++ )
            {
                // Nothing to do if already 0
                if ( true == math::NumericUtil<T>::isZero(temp.at(r, c)) )
                {
                    continue;  // for c
                }

                // To set temp(r,c) to 0 it is a good idea to add the c^th row to it.
                // temp(c,i); i<c are already 0 (i.e. will not affect anything left of temp(i,c)
                // and temp(c,c) is already 1.

                el = temp.get(r, c);

                for ( size_t i=c; i<N; i++ )
                {
                    temp.at(r, i) -= el*temp.at(c, i);
                }

                for ( size_t i=0; i<NT; i++ )
                {
                    retVal.at(r, i) -= el*retVal.at(c, i);
                }
            }  // for c
        }  // for r


        return retVal;

    }  // try
    catch ( const math::MatrixException& mex )
    {
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::OUT_OF_MEMORY);
    }

}
