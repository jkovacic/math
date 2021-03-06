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
 * @headername{LinearEquationSolver.h}
 *
 * Declaration and implementation of functions within namespace
 * LinearEquationSolver for solving systems of linear equations.
 */

#ifndef _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_
#define _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_


#include <cstddef>
#include <vector>
#include <new>
#include <algorithm>

#include "../settings/lineq_settings.h"
#include "util/PseudoFunctionGeneric.hpp"
#include "util/NumericUtil.hpp"
#include "matrix/MatrixGeneric.hpp"

#include "omp/omp_header.h"
#include "../settings/omp_settings.h"
#include "omp/omp_coarse.h"

#include "exception/MatrixException.hpp"


namespace math
{

/**
 * @brief A namespace with functions for solving systems of linear equations.
 */
namespace LinearEquationSolver
{


namespace __private
{


/*
 * Dot product of the a's currentRow.th row and the x's currentCol.th
 * column, excluding the currentRow.th element in both vectors.
 * Permutations of a's rows and columns are taken into account by
 * vectors 'rows' and 'cols', respectively.
 *
 * @note The function assumes that dimensions of all matrices and vectors
 *       are correct.
 *
 * @param a - a square matrix of coefficients (N*N)
 * @param x - a matrix of solutions (with N rows)
 * @param currentRow - the requested row of 'a'
 * @param currentCol - the requested column of 'x'
 * @param rows - permutation vector of a's rows
 * @param cols - permutation vector of a's columns
 *
 * @return dot product of a's currentRow.th row and x's currentCol.th column
 */
template <class T>
T incSumProd(
    const math::MatrixGeneric<T>& a,
    const math::MatrixGeneric<T>& x,
    const std::size_t currentRow,
    const std::size_t currentCol,
    const std::vector<std::size_t>& rows,
    const std::vector<std::size_t>& cols )
{
    /*
     * Calculates the following sum:
     *
     *                 N
     *               -----  /               \
     *               \      |               |
     *   sumprod  =   >     | a     *  x    |
     *               /      |  i,j      j,c |
     *               -----  \               /
     *                j=1
     *                j!=i
     *
     * where 'i' is actually equals 'currentRow' and
     * 'c' actually equals 'currentCol'.
     */

    const std::size_t N = a.nrRows();
    T sumProd = static_cast<T>(0);

    #pragma omp parallel num_threads(ompIdeal(N)) \
            if(N>OMP_CHUNKS_PER_THREAD) \
            default(none) shared(sumProd, rows, cols, a, x)
    {
        OMP_COARSE_GRAINED_PAR_INIT_VARS(N);

        T tempSum = static_cast<T>(0);
        for ( std::size_t j=istart; j<iend; ++j )
        {

            if ( currentRow == j )
            {
                continue;  // for j
            }

            tempSum += a(rows.at(currentRow), cols.at(j)) * x(j, currentCol);
        }  // for j

        // Note that "omp reduction" does not support complex types...
        #pragma omp critical(lineqgeneric_reduce_sumprod)
        {
            sumProd += tempSum;
        }
    }  // omp parallel

    return sumProd;
}

}  // namespace math::LinearEquationsolver::__private



/**
 * Solves the system of linear equations using the Gauss - Jordan elimination
 * method and returns its unique solution if it exists.
 * 
 * 'coef' must be a square matrix and 'term' must have the same number of
 * rows as 'coef'.
 * 
 * The method performs either partial or full pivoting. While partial pivoting
 * should be sufficient in most cases, full pivoting is usually numerically
 * more stable, however it introduces additional overhead.
 * 
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * @param fullp - should the method perform full pivoting (default: TRUE)
 * @param physSwap - should the internal algorithm perform physical swapping of matrix elements (default: FALSE)
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' are invalid or internal allocation of memory failed
 */
template <class T>
bool solveGaussJordan(
		  const MatrixGeneric<T>& coef,
		  const MatrixGeneric<T>& term,
		  MatrixGeneric<T>& sol,
		  const bool fullp = LINEQ_GAUSS_JORDAN_FULL_PIVOT,
		  const bool physSwap = LINEQ_PHYSSWAP_COEF
		)
{
    // Sanity check
    if ( false == coef.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    const std::size_t N = coef.nrColumns();          // Nr. of unknowns

    // Check the dimensions
    if ( N != term.nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    return math::Pivot::solveGaussJordan(coef, term, sol, fullp, physSwap);
}



/**
 * Solves the system of linear equations using the successive over-relaxation
 * (SOR) iterative method and returns its unique solution if it exists.
 * 
 * 'coef' must be a square matrix and 'term' must have the same number of
 * rows as 'coef'. If 'solInitialized' equals TRUE, 'sol' must have the same
 * dimension as 'term'.
 * 
 * @note The method typically converges when it is possible to permute
 * rows and columns of 'coef' into a diagonally dominant matrix.
 *
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * @param w - relaxation factor (omega)
 * @param solInitialized - is 'sol' prefilled with the initial solution? (default: FALSE)
 * @param tol - infinity norm of the maximum residual of tolerance (default: 1e-6)
 * @param maxiter - maximum number of iterations (default: 10000)
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' (and possibly 'sol') are invalid or internal allocation of memory failed
 */
template <class T>
bool solveSOR(
		  const MatrixGeneric<T>& coef,
		  const MatrixGeneric<T>& term,
		  MatrixGeneric<T>& sol,
		  const T& w,
		  const bool solInitialized = false,
		  const T& tol = static_cast<T>(LINEQ_TOL_CONV_NUM) / static_cast<T>(LINEQ_TOL_CONV_DEN),
		  const std::size_t maxiter = LINEQ_MAX_ITER
		)
{
    /*
     * The algorithm of the SOR method is described in detail at:
     * https://en.wikipedia.org/wiki/Successive_over-relaxation
     */

    // Sanity check
    if ( false == coef.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    if ( true == math::NumericUtil::isZero(w) )
    {
        return false;
    }

    const std::size_t N = coef.nrColumns();
    const std::size_t NC = term.nrColumns();

    // Check the dimensions
    if ( N != term.nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // if 'sol' is not initialized, fill it with the initial solution (zeros)...
    if ( false == solInitialized )
    {
        sol = term;
        sol *= static_cast<T>(0);
    }
    else
    {
        // otherwise just check its dimensions
        if ( N != sol.nrRows() || NC != sol.nrColumns() )
        {
            throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
        }
    }

    // Vectors of permutations of rows and columns, respectively
    std::vector<std::size_t> rows;
    std::vector<std::size_t> cols;

    /* 
     * Try to "convert" 'coef' to a diagonally dominant matrix to ensure
     * convergence of the algorithm.
     * The function will already swap rows of 'sol', however it won't modify
     * the matrix 'coef'. This means that the further procedure will require
     * to access elements of 'coef' and 'term' via the permutation vector 'rows',
     * on the other hand this is not necessary to access elements of 'sol'.
     */
    if ( false == math::Pivot::getDiagonallyDominantMatrix(coef, &sol, rows, cols) )
    {
        return false;
    }

    // counter of iterations
    std::size_t cnt;
    // maximum inf. norm of all columns, initially set to something larger than 'tol'
    T maxInfNorm = static_cast<T>(10) * tol;

    // Iterate until the algorithm converges or 'cnt' exceeds 'maxiter'
    for ( cnt=0; cnt<maxiter && math::PseudoFunction::absgt(maxInfNorm, tol) ; ++cnt )
    {
        maxInfNorm = static_cast<T>(0);

        // for each column of term...
        //   (note that each column is processed independently from the others)
        #pragma omp paralel for default(none) \
                shared(maxInfNorm, sol, coef, term, rows, cols, w)
        for ( std::size_t c=0; c<NC; ++c )
        {
            // These variables will store column's maximum abs. value and increment
            T Xpmax = static_cast<T>(0);
            T Dpmax = static_cast<T>(0);

            /*
             * Each unknown will be updated iteratively...
             * 
             * Note: this for loop must be executed sequentially
             *       and CANNOT be parallelized!
             */
            for ( std::size_t i=0; i<N; ++i )
            {
                // update the column's maximum abs. value if necessary
                const T Xpabs = math::PseudoFunction::pabs(sol(i, c));
                if ( true == math::PseudoFunction::absgt(Xpabs, Xpmax) )
                {
                    Xpmax = Xpabs;
                }

                /*
                 *       -----                  -----
                 *       \             (k+1)    \            (k)
                 *   s =  >    a    * x      +   >   a    * x
                 *       /      i,j    j        /     i,j    j
                 *       -----                  -----
                 *        j<i                    j>i
                 */

                const T s = math::LinearEquationSolver::__private::incSumProd(
                           coef, sol, i, c, rows, cols );

                /*
                 * Update the x_i:
                 * 
                 *                                       /        \
                 *    (k+1)              (k)      w      |        |
                 *   x      = (1 - w) * x    +  ------ * | b  - s |
                 *    i                  i       a       |  i     |
                 *                                i,i    \        /
                 * 
                 * It can be further simplified to reduce the nr. of multiplications:
                 * 
                 *                       /                  \
                 *    (k+1)    (k)       |  b_i - s     (k) |
                 *   x      = x    + w * | --------- - x    |
                 *    i        i         |    a         i   |
                 *                       \     i,i          /
                 * 
                 * However, evaluation of the inf. norm will require the increment
                 * to be calculated first:
                 */
                const T dx = w * ( (term(rows.at(i), c) - s)/coef(rows.at(i), cols.at(i)) - sol(i, c) );

                // then update x_i
                sol(i, c) += dx;

                // update the column's maximum abs. increment if necessary:
                const T Dpabs = math::PseudoFunction::pabs(dx);
                if ( true == math::PseudoFunction::absgt(Dpabs, Dpmax) )
                {
                    Dpmax = Dpabs;
                }
            }  // for i

            // Convert maxima's pseudo absolute values to the actual ones:
            const T Xmax = math::PseudoFunction::pabs2abs(Xpmax);
            const T Dmax = math::PseudoFunction::pabs2abs(Dpmax);

            /*
             * Obtain the column's residual infinite norm:
             * 
             *                 || Dx_i ||_inf
             *   resInfNorm = -----------------
             *                  || x_i ||_inf
             * 
             * Where the infinite norm of a vector is defined as the maximum
             * of vector's elements:
             * 
             *   || x ||_inf = max( |x_1|, |x_2|, ..., |x_n| )
             */

            T resInfNorm;

            /*
             * Divide both abs. maximums and also take care of possible division
             * by zero. When a non-zero value is attempted to be divided by 0,
             * just assign it a value greater than 'tol'
             */            
            if ( true == math::NumericUtil::isZero(Xmax) )
            {
                resInfNorm = ( true==math::NumericUtil::isZero(Dmax) ? 
                               static_cast<T>(0) : static_cast<T>(10) * tol );
            }
            else
            {
                resInfNorm = Dmax / Xmax;
            }

            // Finally update the maximum residual inf. norm across all columns
            // Note that the update must be reduced to prevent possible race conditions
            #pragma omp critical(lineqgeneric_maxinfnorm_sor)
            {
                if ( true == math::PseudoFunction::absgt(resInfNorm, maxInfNorm) )
                {
                    maxInfNorm = resInfNorm;
                }
            }
        }  // for c
    }  // for cnt

    // Check if the algorithm has converged
    if ( cnt >= maxiter )
    {
        return false;
    }

    // rearrange sol's rows according to full pivoting permutations of columns:
    math::Pivot::rearrangeMatrixRows(sol, cols);

    return true;
}



/**
 * Solves the system of linear equations using the Gauss - Seidel iterative
 * method and returns its unique solution if it exists.
 * 
 * 'coef' must be a square matrix and 'term' must have the same number of
 * rows as 'coef'. If 'solInitialized' equals TRUE, 'sol' must have the same
 * dimension as 'term'.
 * 
 * @note The method typically converges when it is possible to permute
 * rows and columns of 'coef' into a diagonally dominant matrix.
 *
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * @param solInitialized - is 'sol' prefilled with the initial solution? (default: FALSE)
 * @param tol - infinity norm of the maximum residual of tolerance (default: 1e-6)
 * @param maxiter - maximum number of iterations (default: 10000)
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' (and possibly 'sol') are invalid or internal allocation of memory failed
 */
template <class T>
bool solveGaussSeidel(
		  const MatrixGeneric<T>& coef,
		  const MatrixGeneric<T>& term,
		  MatrixGeneric<T>& sol,
		  const bool solInitialized = false,
		  const T& tol = static_cast<T>(LINEQ_TOL_CONV_NUM) / static_cast<T>(LINEQ_TOL_CONV_DEN),
		  const std::size_t maxiter = LINEQ_MAX_ITER
		)
{
    /*
     * The Gauss - Seidel method is described in detail at:
     * https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
     * 
     * It is actually equivalent to the SOR method with omega equal
     * to 1.
     */

    return math::LinearEquationSolver::solveSOR(
            coef, term, sol, static_cast<T>(1), solInitialized, tol, maxiter );
}



/**
 * Solves the system of linear equations using the weighted Jacobi
 * iterative method and returns its unique solution if it exists.
 * 
 * 'coef' must be a square matrix and 'term' must have the same number of
 * rows as 'coef'. If 'solInitialized' equals TRUE, 'sol' must have the same
 * dimension as 'term'.
 * 
 * @note The method typically converges when it is possible to permute
 * rows and columns of 'coef' into a diagonally dominant matrix.
 *
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * @param w - weight parameter (a.k.a. omega)
 * @param solInitialized - is 'sol' prefilled with the initial solution? (default: FALSE)
 * @param tol - infinity norm of the maximum residual of tolerance (default: 1e-6)
 * @param maxiter - maximum number of iterations (default: 10000)
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' (and possibly 'sol') are invalid or internal allocation of memory failed
 */
template <class T>
bool solveWeightedJacobi(
		  const MatrixGeneric<T>& coef,
		  const MatrixGeneric<T>& term,
		  MatrixGeneric<T>& sol,
		  const T& w,
		  const bool solInitialized = false,
		  const T& tol = static_cast<T>(LINEQ_TOL_CONV_NUM) / static_cast<T>(LINEQ_TOL_CONV_DEN),
		  const std::size_t maxiter = LINEQ_MAX_ITER
		)
{
    /*
     * The algorithm of the weighted Jacobi method is described in detail at:
     * https://en.wikipedia.org/wiki/Jacobi_method#Weighted_Jacobi_method
     */

    // Sanity check
    if ( false == coef.isSquare() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    if ( true == math::NumericUtil::isZero(w) )
    {
        return false;
    }

    const std::size_t N = coef.nrColumns();
    const std::size_t NC = term.nrColumns();

    // Check the dimensions
    if ( N != term.nrRows() )
    {
        throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
    }

    // if 'sol' is not initialized, fill it with the initial solution (zeros)...
    if ( false == solInitialized )
    {
        sol = term;
        sol *= static_cast<T>(0);
    }
    else
    {
        // otherwise just check its dimensions
        if ( N != sol.nrRows() || NC != sol.nrColumns() )
        {
            throw math::MatrixException(math::MatrixException::INVALID_DIMENSION);
        }
    }

    // Vectors of permutations of rows and columns, respectively
    std::vector<std::size_t> rows;
    std::vector<std::size_t> cols;

    /* 
     * Try to "convert" 'coef' to a diagonally dominant matrix to ensure
     * convergence of the algorithm.
     * The function will already swap rows of 'sol', however it won't modify
     * the matrix 'coef'. This means that the further procedure will require
     * to access elements of 'coef' and 'term' via the permutation vector 'rows',
     * on the other hand this is not necessary to access elements of 'sol'.
     */
    if ( false == math::Pivot::getDiagonallyDominantMatrix(coef, &sol, rows, cols) )
    {
        return false;
    }

    // counter of iterations
    std::size_t cnt;
    // maximum inf. norm of all columns, initially set to something larger than 'tol'
    T maxInfNorm = static_cast<T>(10) * tol;

    /*
     * The algorithm requires an extra matrix to store values
     * of the previous (or next) iteration. The algorithm will
     * switch between 'sol' and this temporary matrix.
     */
    math::MatrixGeneric<T> tempMat(sol);

    // Iterate until the algorithm converges or 'cnt' exceeds 'maxiter'
    for ( cnt=0; cnt<maxiter && math::PseudoFunction::absgt(maxInfNorm, tol); ++cnt )
    {
        /*
         * When 'cnt' is even (divisible by 2), 'sol' stores values of the current
         * iteration x^(k) and 'tempMat' of the next iteration x^(k+1).
         * When 'cnt' is odd, the roles of 'sol' and 'tempMat' are swapped.
         */
        const math::MatrixGeneric<T>& xk = ( cnt%2==0 ? sol : tempMat );
        math::MatrixGeneric<T>& xk_1 = ( cnt%2==1 ? sol : tempMat );

        maxInfNorm = static_cast<T>(0);

        // for each column of term...
        //   (note that each column is processed independently from the others)
        #pragma omp paralel for default(none) \
                shared(maxInfNorm, xk, xk_1, coef, term, rows, cols, w)
        for ( std::size_t c=0; c<NC; ++c )
        {
            // These variables will store column's maximum abs. value and increment
            T Xpmax = static_cast<T>(0);
            T Dpmax = static_cast<T>(0);

            /*
             * Each unknown will be updated iteratively...
             */
            #pragma omp parallel for default(none) \
                    shared(coef, term, rows, cols, xk, xk_1, Xpmax, Dpmax, c, w)
            for ( std::size_t i=0; i<N; ++i )
            {
                /*
                 *           N
                 *         -----
                 *         \             (k)
                 *     s =  >    a    * x
                 *         /      i,j    j
                 *         -----
                 *          j=0
                 *          j!=i
                 */

                const T s = math::LinearEquationSolver::__private::incSumProd(
                           coef, xk, i, c, rows, cols );

                /*
                 * Update the x_i:
                 * 
                 *                                       /        \
                 *    (k+1)              (k)      w      |        |
                 *   x      = (1 - w) * x    +  ------ * | b  - s |
                 *    i                  i       a       |  i     |
                 *                                i,i    \        /
                 * 
                 * It can be further simplified to reduce the nr. of multiplications:
                 * 
                 *                       /                  \
                 *    (k+1)    (k)       |  b_i - s     (k) |
                 *   x      = x    + w * | --------- - x    |
                 *    i        i         |    a         i   |
                 *                       \     i,i          /
                 * 
                 * However, evaluation of the inf. norm will require the increment
                 * to be calculated first:
                 */
                const T dx = w * ( (term(rows.at(i), c) - s)/coef(rows.at(i), cols.at(i)) - xk(i, c) );

                // then update x_i
                xk_1(i, c) = xk(i, c) + dx;

                /*
                 * Update the column's maximum abs. value and abs. increment
                 * if necessary. Note that both updates must be reduced in order
                 * to prevent possible race conditions.
                 */

                const T Xpabs = math::PseudoFunction::pabs(xk(i, c));
                const T Dpabs = math::PseudoFunction::pabs(dx);

                #pragma omp critical(lineqgeneric_dpmaxxpmax_jacobi)
                {
                    if ( true == math::PseudoFunction::absgt(Xpabs, Xpmax) )
                    {
                        Xpmax = Xpabs;
                    }

                    if ( true == math::PseudoFunction::absgt(Dpabs, Dpmax) )
                    {
                        Dpmax = Dpabs;
                    }
                }  // omp critical
            }  // for i

            // Convert maxima's pseudo absolute values to the actual ones:
            const T Xmax = math::PseudoFunction::pabs2abs(Xpmax);
            const T Dmax = math::PseudoFunction::pabs2abs(Dpmax);

            /*
             * Obtain the column's residual infinite norm:
             * 
             *                 || Dx_i ||_inf
             *   resInfNorm = -----------------
             *                  || x_i ||_inf
             * 
             * Where the infinite norm of a vector is defined as the maximum
             * of vector's elements:
             * 
             *   || x ||_inf = max( |x_1|, |x_2|, ..., |x_n| )
             */

            T resInfNorm;

            /*
             * Divide both abs. maximums and also take care of possible division
             * by zero. When a non-zero value is attempted to be divided by 0,
             * just assign it a value greater than 'tol'
             */            
            if ( true == math::NumericUtil::isZero(Xmax) )
            {
                resInfNorm = ( true==math::NumericUtil::isZero(Dmax) ? 
                               static_cast<T>(0) : static_cast<T>(10) * tol );
            }
            else
            {
                resInfNorm = Dmax / Xmax;
            }

            // Finally update the maximum residual inf. norm across all columns
            // Note that the update must be reduced to prevent possible race conditions
            #pragma omp critical(lineqgeneric_maxinfnorm_jacobi)
            {
                if ( true == math::PseudoFunction::absgt(resInfNorm, maxInfNorm) )
                {
                    maxInfNorm = resInfNorm;
                }
            }
        }  // for c
    }  // for cnt

    // Check if the algorithm has converged
    if ( cnt >= maxiter )
    {
        return false;
    }

    /*
     * If the total number of iterations (i.e. updated 'cnt' after the
     * final iteration) is odd (not divisible by 2), the most recent solution
     * is stored in 'tempMat' hence this matrix must be copied to 'sol'.
     */
    if ( cnt%2 == 1 )
    {
        sol = tempMat;
    }

    // rearrange sol's rows according to full pivoting permutations of columns:
    math::Pivot::rearrangeMatrixRows(sol, cols);

    return true;
}



/**
 * Solves the system of linear equations using the Jacobi iterative
 * method and returns its unique solution if it exists.
 * 
 * 'coef' must be a square matrix and 'term' must have the same number of
 * rows as 'coef'. If 'solInitialized' equals TRUE, 'sol' must have the same
 * dimension as 'term'.
 * 
 * @note The method typically converges when it is possible to permute
 * rows and columns of 'coef' into a diagonally dominant matrix.
 *
 * @param coef - a square matrix with coefficients of the system of linear equations
 * @param term - a matrix with constant terms of the system of linear equations
 * @param sol - a reference to a matrix to be assigned the solution of equations
 * @param solInitialized - is 'sol' prefilled with the initial solution? (default: FALSE)
 * @param tol - infinity norm of the maximum residual of tolerance (default: 1e-6)
 * @param maxiter - maximum number of iterations (default: 10000)
 * 
 * @return a logical value indicating whether a unique solution was found
 * 
 * @throw MatrixException if dimensions of 'coef' and 'term' (and possibly 'sol') are invalid or internal allocation of memory failed
 */
template <class T>
bool solveJacobi(
		  const MatrixGeneric<T>& coef,
		  const MatrixGeneric<T>& term,
		  MatrixGeneric<T>& sol,
		  const bool solInitialized = false,
		  const T& tol = static_cast<T>(LINEQ_TOL_CONV_NUM) / static_cast<T>(LINEQ_TOL_CONV_DEN),
		  const std::size_t maxiter = LINEQ_MAX_ITER
		)
{
    /*
     * The Jacobi method is described in detail at:
     * https://en.wikipedia.org/wiki/Jacobi_method
     * 
     * It is actually equivalent to the weighted Jacobi method with
     * weight equal to 1.
     */

    return math::LinearEquationSolver::solveWeightedJacobi(
            coef, term, sol, static_cast<T>(1), solInitialized, tol, maxiter );
}


}  // namespace LinearEquationSolver

}  // namespace math


#endif  // _MATH_LINEAREQUATIONSOLVERGENERIC_HPP_
