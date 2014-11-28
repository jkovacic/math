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
 * A test module to test matrix classes
 * (MatrixGeneric, SqMatrixGeneric)
 */


#include <iostream>

#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"
#include "MatrixException.h"


using namespace std;
using namespace math;


/*
 * Test of float matrices (float is the most common type used in matrices)
 */
void matrixTest()
{
    // Test of several matrix operations.
    // Results may be checked using MATLAB, Octave, Scilab etc.
    try
    {
        // Create a matrix. f(1,0) should be set to 0 by constructor
        FMatrix f(2,3);
        f.set(0, 0, 1.0f).set(0, 1, 0.5f).set(0, 2, 4.5f);
        /*f(1,0)=0*/    f.set(1, 1, 1.0f).set(1, 2, 0.4f);
        cout << "f:" << endl;
        f.display();

        // Multiplication by a scalar
        FMatrix res(1,1);
        res = 3.0f * f;
        cout << "f multiplied by 3" << endl;
        res.display();
        cout << endl;

        // and scalar * matrix:
        res = f * 0.5f;
        cout << "f multiplied by 0.5" << endl;
        res.display();
        cout << endl;

        // Test of transposing of matrices
        FMatrix t = f.transpose();
        cout << "t = f transposed:" << endl;
        t.display();
        cout << endl;

        // Test of matrix product
        FMatrix prod = t*f;
        cout << "t * f:" << endl;
        prod.display();
        cout << endl;

        // product is not commutative...
        prod = f*t;
        cout << "f * t:" << endl;
        prod.display();
        cout << endl;

        // test creation of a unit (diagonal) matrices
        FSqMatrix a(3);
        a.setUnit();
        cout << "3x3 unit matrix:" << endl;
        a.display();
        cout << endl;

        // row/column removal and insert operations will be tested on this matrix
        FMatrix m1(3,4);
        m1.set(0, 0, 1);    m1.set(0, 1, 2.4f);  m1.set(0, 2, -1.4f); m1.set(0, 3, 0);
        m1.set(1, 0, 4.5f); m1.set(1, 1, 1);     m1.set(1, 2, 0);     m1.set(1, 3, -0.5f);
        m1.set(2, 0, 0);    m1.set(2, 1, 1.75f); m1.set(2, 2, 1);     m1.set(2, 3, 2);
        cout << "m1:" << endl;
        m1.display();
        cout << endl;

        m1.removeColumn(2);
        cout << "removed the 3rd column:" << endl;
        m1.display();
        cout << endl;

        m1.insertColumn(1);
        cout << "inserted a column (zeroes) between the 1st and the 2nd column:" << endl;
        m1.display();
        cout << endl;

        m1.removeRow(1);
        cout << "removed the 2nd row" << endl;
        m1.display();
        cout << endl;

        m1.insertRow(1);
        cout << "inserted a row (zeroes) between the 1st and 2nd row:" << endl;
        m1.display();
        cout << endl;

        cout << "transposed the previous matrix:" << endl;
        m1.transposed().display();
        cout << endl;

        // Test of copy constructor
        // SqMatrix's copy constructor must accept generic matrices where rows == cols
        FMatrix sq(2, 2);
        sq.set(0,0,1).set(0,1,2);
        sq.set(1,0,3).set(1,1,4);
        FSqMatrix kv(sq);
        cout << "kv:" << endl;
        kv.display();
        cout << endl;

        // test of a unary operator -
        FMatrix m3 = -kv;
        cout << "-kv:" << endl;
        m3.display();
        cout << endl;

        // Square matrix to test determinant and inverse
        FSqMatrix a1(3);
        a1.set(0, 0, 1.0f).set(0, 1, 2.0f).set(0, 2, 3.0f);
        a1.set(1, 0, 4.0f).set(1, 1, 5.0f).set(1, 2, 6.0f);
        a1.set(2, 0, 7.0f).set(2, 1, 9.0f).set(2, 2, 8.0f);

        // Test calculation of the determinant
        float d = a1.determinant();
        cout << "a1:" << endl;
        a1.display();
        cout << endl;
        cout << "Determinant of a1: " << d << endl << endl;;

        // and calculation of the inverse matrix
        FSqMatrix inv = a1.inverse();
        cout << "inv = inverse of a1:" <<endl;
        inv.display();
        cout << endl;

        // Test if inverse matrix was calculated properly
        // a1*inv must be a unit matrix
        FMatrix prodUnit;
        prodUnit = a1 * inv;
        cout << "a1 * inv   (must be a unit matrix):" << endl;
        prodUnit.display();
        cout << endl;

        // Test self transpose of a square matrix:
        cout << "inv transposed:" << endl;
        FMatrix* pinv = &inv;
        pinv->transposed().display();
        cout << endl;

        cout << "Add 0.5 to inv(2, 0):" << endl;
        inv.at(2, 0) += 0.5f;
        inv.display();
        cout << endl;

        cout << "3x3 unit matrix:" << endl;
        a1 = NumericUtil::unit(a1);
        a1.display();
        cout << endl;
    }
    catch ( const MatrixException& ex )
    {
        cerr << "MatrixException caught: '";
        ex.what();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}
