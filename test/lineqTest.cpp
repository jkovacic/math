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

#include <iostream>
#include <complex>

#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"
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
        SqMatrixGeneric<complex<float> > a(3);
        MatrixGeneric<complex<float> > b(3, 1);

        a.set(0, 0, complex<float>(1, 1)).set(0, 1, complex<float>(2, -1)).set(0, 2, complex<float>(-1, 0.5));
        a.set(1, 0, complex<float>(0, 1)).set(1, 1, complex<float>(0.25, 3)).set(1, 2, complex<float>(-2, -0.5));
        a.set(2, 0, complex<float>(3, 0)).set(2, 1, complex<float>(2, -3)).set(2, 2, complex<float>(-0.5, 1));

        b.set(0, 0, complex<float>(1, 0.2)).set(1, 0, complex<float>(-2, -1)).set(2, 0, complex<float>(1, 0));

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
        MatrixGeneric<complex<float> > x;
        if ( true == LinearEquationSolver::solveGaussJordan(a, b, x) )
        {
            cout << "Solution:" << endl;
            x.display();
        }
        else
        {
            cerr << "Unique solution does not exist" << endl;
        }
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
