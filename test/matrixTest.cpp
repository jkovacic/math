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


/*
 * Note: results are reproduced in 'scripts/test/matrix.m'.
 */


#include <iostream>
#include <complex>

#include "MatrixGeneric.h"
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
        Matrix f(2,3);
        f.set(0, 0, 1.0).set(0, 1, 0.5).set(0, 2, 4.5);
        /*f(1,0)=0*/    f.set(1, 1, 1.0).set(1, 2, 0.4);
        cout << "f:" << endl;
        f.display();
        cout << endl;
        
        cout << "f(0,1) = " << f(0, 1) << endl;
        cout << "f(4)   = " << f(4) << endl << endl;

        // Multiplication by a scalar
        Matrix res(1,1);
        res = 3.0 * f;
        cout << "f multiplied by 3" << endl;
        res.display();
        cout << endl;

        // Conjugation of a real numbered matrix:
        cout << "f conjugated:" << endl;
        res = f.conj();
        res.display();
        cout << endl;

        // and scalar * matrix:
        res = +f * 0.5;
        cout << "f multiplied by 0.5" << endl;
        res.display();
        cout << endl;

        // Test of transposing of matrices
        Matrix t = f.transpose();
        cout << "t = f transposed:" << endl;
        t.display();
        cout << endl;
        cout << "rank: " << t.rank() << endl << endl;

        // Test of matrix product
        Matrix prod = t*f;
        cout << "t * f:" << endl;
        prod.display();
        cout << endl;

        // product is not commutative...
        prod = f*t;
        cout << "f * t:" << endl;
        prod.display();
        cout << endl;

        // test creation of a unit (diagonal) matrices
        Matrix a(3);
        a.setUnit_();
        cout << "3x3 unit matrix:" << endl;
        a.display();
        cout << endl;

        // row/column removal and insert operations will be tested on this matrix
        Matrix m1(3,4);
        m1.set(0, 0, 1).set(0, 1, 2.4).set(0, 2, -1.4).set(0, 3, 0);
        m1.set(1, 0, 4.5).set(1, 1, 1).set(1, 2, 0).set(1, 3, -0.5);
        m1.set(2, 0, 0).set(2, 1, 1.75).set(2, 2, 1).set(2, 3, 2);
        cout << "m1:" << endl;
        m1.display();
        cout << endl;
        cout << "rank(m1): " << m1.rank() << endl << endl;

        m1.removeColumn_(2);
        cout << "removed the 3rd column:" << endl;
        m1.display();
        cout << endl;

        m1.insertColumn_(1);
        cout << "inserted a column (zeroes) between the 1st and the 2nd column:" << endl;
        m1.display();
        cout << endl;

        m1.removeRow_(1);
        cout << "removed the 2nd row" << endl;
        m1.display();
        cout << endl;

        m1.insertRow_(1);
        cout << "inserted a row (zeroes) between the 1st and 2nd row:" << endl;
        m1.display();
        cout << endl;
        cout << "rank: " << m1.rank() << endl << endl;

        cout << "transposed the previous matrix:" << endl;
        m1.transpose_().display();
        cout << endl;

        // Test of copy constructor
        Matrix sq(2, 2);
        sq.set(0,0,1).set(0,1,2);
        sq.set(1,0,3).set(1,1,4);
        Matrix kv(sq);
        cout << "kv:" << endl;
        kv.display();
        cout << endl;

        // test of a unary operator -
        Matrix m3 = -kv;
        cout << "-kv:" << endl;
        m3.display();
        cout << endl;

        // Square matrix to test determinant and inverse
        Matrix a1(3);
        a1.set(0, 0, 1.0).set(0, 1, 2.0).set(0, 2, 3.0);
        a1.set(1, 0, 4.0).set(1, 1, 5.0).set(1, 2, 6.0);
        a1.set(2, 0, 7.0).set(2, 1, 9.0).set(2, 2, 8.0);

        // Test calculation of the determinant
        double d = a1.determinant();
        cout << "a1:" << endl;
        a1.roundSmallElements().display();
        cout << endl;
        cout << "Determinant of a1: " << d << endl;
        cout << "Rank of a1: " << a1.rank() << endl << endl;

        // and calculation of the inverse matrix
        Matrix inv = a1.inverse();
        cout << "inv = inverse of a1:" <<endl;
        inv.display();
        cout << endl;

        // Test if inverse matrix was calculated properly
        // a1*inv must be a unit matrix
        Matrix prodUnit(inv);
        prodUnit = a1 * inv;
        cout << "a1 * inv   (must be a unit matrix):" << endl;
        prodUnit.roundSmallElements(1e-15).display();
        cout << endl;

        cout << "Upper triangular part (incl. diag) of a1:" << endl;
        m1 = a1.upperTriangularPart();
        m1.display();
        cout << endl;

        cout << "Lower triangular part (excl. diag) of a1:" << endl;
        m1 = a1.lowerTriangularPart(false);
        m1.display();
        cout << endl;

        cout << "Diagonal part of a1:" << endl;
        m1 = a1.diagPart();
        m1.display();
        cout << endl;

        // Test self transpose of a square matrix:
        cout << "inv transposed:" << endl;
        Matrix* pinv = &inv;
        pinv->transpose_().display();
        cout << endl;

        cout << "Add 0.5 to inv(2, 0):" << endl;
        inv(2, 0) += 0.5;
        inv.display();
        cout << endl;

        cout << "Swap the 2nd and the 3rd row:" << endl;
        inv.swapRows_(1, 2);
        inv.display();
        cout << endl;

        cout << "Swap the 2nd and the 3rd column:" << endl;
        inv.swapColumns_(1, 2);
        inv.display();
        cout << endl;

        a1 = 5.5;
        cout << "a1:" << endl;
        a1.display();
        cout << endl;

        // Test of complex conjugation:
        MatrixGeneric<complex<float> > c1(2, 2);
        c1.set(0, 0, complex<float>(1.0f, 1.0f)).set(0, 1, complex<float>(1.0f, -2.0f));
        c1.set(1, 0, complex<float>(2.0f, -3.0f)).set(1, 1, complex<float>(2.0f, 4.0f));
        cout << "c1:" << endl;
        c1.display();
        MatrixGeneric<complex<float> > c2 = c1.conj();
        cout << "c1 conjugated:" << endl;
        c2.display();
        cout << endl;

        FMatrix fm(4, 1);
        fm.set(0, 0, 1.1f).set(1, 0, 2.2f).set(2, 0, 3.3f).set(3, 0, 4.4f);
        cout << "fm:" << endl;
        fm.display();
        cout << endl;
        cout << "fm transposed:" << endl;
        FMatrix fmt = fm.transpose();
        fmt.display();
        cout << endl;
        fmt.transpose_();
        cout << "fm transposed transposed:" << endl;
        fmt.display();
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
