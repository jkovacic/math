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
#include <cstddef>

#include "MatrixException.h"
#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"
#include "omp_settings.h"

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
        this->m_coef = coef;
        this->m_term = term;
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
    return this->m_coef;
}

/**
 * @return reference to the internal vector or matrix of constant terms
 */
template<class T>
math::MatrixGeneric<T>& math::LinearEquationSolverGeneric<T>::getTerm() const
{
    return this->m_term;
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
        this->m_coef = coef;
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
        this->m_term = term;
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
    const size_t N = this->m_coef.nrColumns();  // Nr. of unknowns
    const size_t NT = this->m_term.nrColumns(); // Nr. of terms' columns
    const size_t Nmax = (N>=NT ? N : NT); // max. of both values

    // Check of dimensions
    if ( N != this->m_term.nrRows() )
    {
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::INVALID_DIMENSION);
    }

    try
    {
        SqMatrixGeneric<T> temp(this->m_coef);
        MatrixGeneric<T> retVal(this->m_term);

        // Try to convert the 'temp' into an identity matrix
        // by appropriately adding multiples of other lines to each line
        // (incl. lines of 'retVal')

        /*
         * The first part of the algorithm will (try to) ensure there are
         * no zeros among temp's diagonal elements. Additionally it will
         * subtract multiples of the i.th row from all subsequent rows (r>i)
         * so their i.th column will be equal to 0. This requires plenty of
         * additions/subtractions of individual rows so parallelization of
         * this for loop is not possible due to race conditions. IT will be
         * possible to parallelize certain parts of this lop, though.
         */
        for ( size_t i=0; i<N; ++i )
        {
            // first check if the diagonal element equals 0
            if ( true == math::NumericUtil<T>::isZero(temp.at(i, i)) )
            {
                size_t r;

                // if it does, try to find another row r where temp(r,i)!=0
                for ( r=0; r<N; ++r )
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

                // add the r.th line to the i.th one and thus prevent temp(i,i) from being 0:
                #pragma omp parallel default(none) shared(temp, retVal, i, r)
                for ( size_t c=0; c<Nmax; ++c )
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

            // Set the i.th column of all other rows (r>i) to 0 by
            // adding the appropriate multiple of the i.th row
            #pragma omp parallel for default(none) shared(temp, retVal, i)
            for ( size_t r=i+1; r<N; r++ )
            {
                // Nothing to do if temp(r,i) is already 0.
                if ( false == math::NumericUtil<T>::isZero(temp.at(r, i)) )
                {
                    // Subtract a multiple of the i^th row.
                    const T el = temp.get(r, i) / temp.get(i, i);

                    // temp(r,:) = temp(r,:)-el*temp(i,:)
                    for ( size_t c=i; c<N; ++c )
                    {
                        temp.at(r, c) -= el*temp.at(i, c);
                    }

                    // term(r,:) = term(r,:)-el*term(i,:)
                    for ( size_t c=0; c<NT; ++c )
                    {
                        retVal.at(r, c) -= el*retVal.at(i, c);
                    }
                }  // if (temp(r,i) != 0
            }  // for r
        }  // for i

        // Set the diag elements of 'temp' and 'retVal' to 1 by dividing the
        // whole row by temp(r,r). Columns smaller than 'r' are already equal to 0.

        // Normalizing of each row is independent from other rows so it is
        // perfectly safe to parallelize the task by rows.
        #pragma omp parallel for default(none) shared(temp, retVal)
        for ( size_t r=0; r<N; ++r )
        {
            const T el = temp.get(r, r);

            for ( size_t c=r; c<N; ++c )
            {
                temp.at(r, c) /= el;
            }

            for ( size_t c=0; c<NT; ++c )
            {
                retVal.at(r, c) /= el;
            }
        }

        // Now the lower triangle (below diag excl.) is 0, the diagonal consists of 1,
        // The upper triangle (above the diag) must be set to 0 as well.

        // Column 'c' of each row (for r<c) will be set to 0
        // by adding the appropriate multiple of c.th row

        // Parallelization of this for loop is not possible due to race conditions.
        for ( size_t c=1; c<N; ++c )
        {
            // The current row to apply the operation described above

            // It is possible to parallelize this for loop because a row 'c' (not included
            // into the for loop) will be added independently to each "parallelized" row
            // 'r' (between 0 and c-1 incl.), thus no race condition is possible.
            #pragma omp parallel for default(none) shared(temp, retVal, c) schedule(dynamic)
            for ( size_t r=0; r<c; ++r )
            {
                // Nothing to do if temp(r,c) already equals 0
                if ( false == math::NumericUtil<T>::isZero(temp.at(r, c)) )
                {
                	// To set temp(r,c) to 0 it is a good idea to add the c.th row to it.
                    // temp(c,i); i<c are already 0 (i.e. will not affect anything left of temp(i,c)
                    // and temp(c,c) is already 1.

                    const T el = temp.get(r, c);

                    // temp(r,:) = temp(r,:) - el*temp(c,:)
                    for ( size_t i=c; i<N; ++i )
                    {
                        temp.at(r, i) -= el*temp.at(c, i);
                    }

                    // term(r,:) = term(r,:) - el*term(c,:)
                    for ( size_t i=0; i<NT; ++i )
                    {
                        retVal.at(r, i) -= el*retVal.at(c, i);
                    }
                }  // if temp(r,c) != 0
            }  // for r
        }  // for c

        return retVal;

    }  // try
    catch ( const math::MatrixException& mex )
    {
        throw math::LinearEquationSolverException(math::LinearEquationSolverException::OUT_OF_MEMORY);
    }

}
