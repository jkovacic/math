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
 * A test module to test rational matrices
 * (Rational, MatrixGeneric)
 */


/*
 * Note: results are reproduced in 'scripts/test/ratmat.jl'.
 */


#include <iostream>

#include "Rational.h"
#include "Matrix.h"

#include "MatrixException.h"
#include "RationalException.h"


using namespace std;
using namespace math;


/*
 * Test of rational matrices
 */
void rationalMatrixTest()
{
    // Test the remaining matrix operations on rational matrices
    try
    {
        MatrixGeneric<Rational> a(2, 2);
        a.set(0, 0, Rational(-1, 2));       a.set(0, 1, Rational(2, 5));
        a.set(1, 0, Rational(3, 4));        a.set(1, 1, Rational(-1, 3));

        // Copy constructor should also implement cast.
        MatrixGeneric<Rational> inv = a.inverse(false);

        cout << "a:" << endl;
        a.display();
        cout << endl;

        // Rational inverse matrix
        cout << "inverse of a" << endl;
        inv.display();
        cout << endl;

        // check the matrix += operator
        a += inv;
        cout << "a += inv:" << endl;
        a.display();
        cout << endl;

        inv *= a;
        cout << "determinant of inv * (a+inv): " << inv.determinant() << endl << endl;

        // test of overloaded matrix operators:
        Rational r(3, 4);
        MatrixGeneric<Rational> b(a);

        b = a + r;
        cout << "a + " << r << " : " << endl;
        b.display();
        cout << endl;

        b = r + a;
        cout << r << " + a : " << endl;
        b.display();
        cout << endl;

        b = a - r;
        cout << "a - " << r << " : " << endl;
        b.display();
        cout << endl;

        b = r - a;
        cout << r << " - a : " << endl;
        b.display();
        cout << endl;

        r.set(1, 2);
        b = a;
        b += r;
        cout << "a + " << r << " : " << endl;
        b.display();
        cout << endl;

        b = a;
        b -= r;
        cout << "a - " << r << " : " << endl;
        b.display();
        cout << endl;

        // Test of division operators:
        b = a / r;
        cout << " a / (" << r << ") : " << endl;
        b.display();
        cout << endl;

        b /= r;
        cout << "a / (" << r << ")^2 : " << endl;
        b.display();
        cout << endl;

        // Test of element wise multiplication and division
        b = matEwMult(a, (a+r));
        cout << "a .* (a+"<< r << ") : " << endl;
        b.display();
        cout << endl;

        b = a;
        b.ewMult(a);
        cout << " a .* a : " << endl;
        b.display();
        cout << endl;

        b = matEwDiv((a-r), a);
        cout << "(a-" << r << ") ./ a : " << endl;
        b.display();
        cout << endl;

        b = a - r;
        b.ewDiv(a+r);
        cout << "(a-" << r << ") ./ (a+" << r << ") : " << endl;
        b.display();
        cout << endl;
    }
    catch ( const RationalException& rex )
    {
        cerr << "RationalException caught: '";
        rex.what();
        cerr << "'" << endl;
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
