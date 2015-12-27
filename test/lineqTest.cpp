/*
Copyright 2014, Jernej Kovacic

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
 * A test module to test functions for solving systems
 * of linear equations (namespace LinearEquationSolver).
 */


/*
 * Note: results are reproduced in 'scripts/test/lineq.m'.
 */


#include <iostream>
#include <complex>

#include "MatrixGeneric.h"
#include "LinearEquationSolverGeneric.h"
#include "MatrixException.h"


using namespace std;
using namespace math;


/*
 * Test of the solver of a system of linear equations
 */
void lineqSolverTest()
{
    try
    {
        MatrixGeneric<complex<double> > a(3);
        MatrixGeneric<complex<double> > b(3, 1);

        a.set(0, 0, complex<double>(1, 1)).set(0, 1, complex<double>(2, -1)).set(0, 2, complex<double>(-1, 0.5));
        a.set(1, 0, complex<double>(0, 1)).set(1, 1, complex<double>(0.25, 3)).set(1, 2, complex<double>(-2, -0.5));
        a.set(2, 0, complex<double>(3, 0)).set(2, 1, complex<double>(2, -3)).set(2, 2, complex<double>(-0.5, 1));

        b.set(0, 0, complex<double>(1, 0.2)).set(1, 0, complex<double>(-2, -1)).set(2, 0, complex<double>(1, 0));

        cout << "Matrix of coefficients:" << endl;
        a.display();
        cout << endl;

        cout << "Vector of constant terms:" << endl;
        b.display();
        cout << endl;

        /*
         * Exact solution of the linear equation system a*x=b:
         * x  = [-0.6631640-0.3626125i, 0.1630189+1.050566i, -0.2240929+0.6002903i]'
         */
        MatrixGeneric<complex<double> > x(b);
        if ( true == LinearEquationSolver::solveGaussJordan(a, b, x, true, true) )
        {
            cout << "Solution with full pivoting:" << endl;
            x.display();
        }
        else
        {
            cerr << "Unique solution does not exist" << endl;
        }

        cout << endl;

        if ( true == LinearEquationSolver::solveGaussJordan(a, b, x, false, true) )
        {
            cout << "Solution with partial pivoting:" << endl;
            x.display();
        }
        else
        {
            cerr << "Unique solution does not exist" << endl;
        }

        cout << endl;

        /*
         * Two systems of linear equations with diagonally dominant
         * coefficient matrices to test iterative methods. One system's
         * coefficient matrix ('ad') is real numbered, the other one ('ac')
         * is complex numbered.
         */
        Matrix ad(3, 3);
        ad.set(0, 0, 4.0).set(0, 1, -1.0).set(0, 2, -1.0);
        ad.set(1, 0, -1.0).set(1, 1, 1.0).set(1, 2, 7.0);
        ad.set(2, 0, -2.0).set(2, 1, 6.0).set(2, 2, 1.0);

        Matrix bd(3, 2);
        bd.set(0, 0, 3.0).set(1, 0, -6.0).set(2, 0, 9.0);
        bd.set(0, 1, 1.0).set(1, 1, 7.0).set(2, 1, -3.0);

        MatrixGeneric<complex<double> > ac(4, 4);
        ac.set(0, 0, complex<double>(0.4, 0.4)).set(0, 1, complex<double>(1.0, -0.7));
        ac.set(0, 2, complex<double>(5.0, 3.0)).set(0, 3, complex<double>(2.0, -1.0));
        ac.set(1, 0, complex<double>(6.0, -1.0)).set(1, 1, complex<double>(0.5, 0.3));
        ac.set(1, 2, complex<double>(-0.1, 2.0)).set(1, 3, complex<double>(1.0, -0.7));
        ac.set(2, 0, complex<double>(0.5, 0.4)).set(2, 1, complex<double>(0.2, -0.3));
        ac.set(2, 2, complex<double>(-1.0, 0.3)).set(2, 3, complex<double>(6.0, -4.0));
        ac.set(3, 0, complex<double>(-0.5, 1.0)).set(3, 1, complex<double>(-7.0, -2.0));
        ac.set(3, 2, complex<double>(0.7, 0.3)).set(3, 3, complex<double>(2.0, -1.0));

        MatrixGeneric<complex<double> > bc(4, 2);
        bc.set(0, 0, complex<double>(2.0, -1.0)).set(0, 1, complex<double>(1.6, 2.4));
        bc.set(1, 0, complex<double>(3.0, 3.0)).set(1, 1, complex<double>(2.7, -1.2));
        bc.set(2, 0, complex<double>(-2.0, -4.0)).set(2, 1, complex<double>(-0.8, -3.1));
        bc.set(3, 0, complex<double>(-1.0, 5.0)).set(3, 1, complex<double>(-1.1, 0.4));

        cout << "Matrix with real coefficients:" << endl;
        ad.display();
        cout << endl;

        cout << "Two vectors of constant terms:" << endl;
        bd.display();
        cout << endl;

        cout << "Matrix with complex coefficients:" << endl;
        ac.display();
        cout << endl;

        cout << "Two vectors of constant terms:" << endl;
        bc.display();
        cout << endl;

        Matrix xd(bd);
        MatrixGeneric<complex<double> > xc(bc);

        if ( true == LinearEquationSolver::solveSOR(ad, bd, xd, 1.5) )
        {
            cout << "Solution of the real-numbered system using the SOR method:" << endl;
            xd.display();
        }
        else
        {
            cerr << "SOR method did not converge." << endl;
        }

        cout << endl;

        if ( true == LinearEquationSolver::solveGaussSeidel(ac, bc, xc) )
        {
            cout << "Solution of the complex-numbered system using the Gauss-Seidel method:" << endl;
            xc.display();
        }
        else
        {
            cerr << "Gauss-Seidel method did not converge." << endl;
        }

        cout << endl;

        xd = bd;
        xc = bc;

        if ( true == LinearEquationSolver::solveWeightedJacobi(ac, bc, xc, complex<double>(2.0/3.0)) )
        {
            cout << "Solution of the complex-numbered system using the weighted Jacobi method:" << endl;
            xc.display();
        }
        else
        {
            cerr << "Weighted Jacobi method did not converge." << endl;
        }

        cout << endl;

        if ( true == LinearEquationSolver::solveJacobi(ad, bd, xd) )
        {
            cout << "Solution of the real-numbered system using the Jacobi method:" << endl;
            xd.display();
        }
        else
        {
            cerr << "Jacobi method did not converge." << endl;
        }

        cout << endl;
    }
    catch ( const MatrixException& mex )
    {
        cerr << "MatrixException caught: '";
        mex.what();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}
