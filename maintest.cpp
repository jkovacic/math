/*
Copyright 2011, Jernej Kovacic

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
* @file maintest.cpp
*
* Collection of basic unit tests for [Sq]MatrixGeneric and Rational
*
* @author Jernej Kovacic
*/

#include <iostream>

#include "Rational.h"
#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"

using namespace std;
using namespace math;

template class MatrixGeneric<float>;
template class MatrixGeneric<Rational>;

/*
 * Test of the class Rational
 */
void rationalTest()
{
    try
    {
        // Constructor test:

        Rational a(15, 20);
        // should output '3/4' and '60/80'
        cout << "a = ";
        a.display();
        cout << " = ";
        a.display(20);
        cout << endl;

        // Test of copy constructors, operator= and operator<<:
        Rational b(a);
        Rational c;
        // should output '3/4', '0/1' and '3/4';
        cout << "b = " << b << "  c = " << c;
        c = a;
        cout << "  c = " << c << endl;

        // now test regular inversion
        // should output '3/4', '4/3', '2/5', '5/2'
        b = a.invert();
        c.set(2, 5);
        cout << "a = " << a << "  b = " << b << "  c = " << c;
        c.inverse();
        cout << "  c^(-1) = " << c << endl;

        // Test of operators
        a.set(1,3);
        b.set(1,2);
        c = a + b;
        // should output sth. like "1/3 + 1/2 = 5/6 = 0.833333"
        cout << a << " + " << b << " = " << c << " = " << c.toFloat() << endl;

        b += a;
        a.set(1,4);
        c = b - a;
        // should output sth. like "5/6 - 1/4 = 7/12 = 0.5833333"
        cout << b << " - " << a << " = " << c << " = " << c.toDouble() << endl;
        a -= c;  //  -1/3
        c = a * b; // 1/24
        // should output "-1/3 * 5/6 = -5/18"
        cout << a << " * " << b << " = " << c << endl;

        c /= b;  // -1/3
        a = b / c; // -5/2
        // should output "5/6 / -1/3 = -5/2"
        cout << b << " / " << c << " = " << a << endl;

        // Test isPositive(),isNegative(), isZero();
        c.set(0,4);
        cout << "a = " << a << " isZero: " << a.isZero() <<
                " isPositive: " << a.isPositive() <<
                " isNegative: " << a.isNegative() << endl;
        cout << "b = " << b << " isZero: " << b.isZero() <<
                " isPositive: " << b.isPositive() <<
                " isNegative: " << b.isNegative() << endl;
        cout << "c = " << c << " isZero: " << c.isZero() <<
                " isPositive: " << c.isPositive() <<
                " isNegative: " << c.isNegative() << endl;

       // Test 'automatic' (via constructor) conversion
       // int -> Rational (a -> a/1)
        c = a + 2;
        // should output "-5/2 + 2 = -1/2"
        cout << a << " + " << 2 << " = " << c << endl;
        // if other Rational operators are defined OK,
        // 'hybrid' operations with integers should work as well

        // Test operators == and !=
        c = a;
        cout << "a = " << a << "  b = " << b << "  c = " << c << endl;
        if ( a == c )
        {
            cout << "a == c";
        }
        else
        {
            cout << "a != c";
        }
        cout << "   ";
        if ( a == b )
        {
            cout << "a == b";
        }
        else
        {
            cout << "a != b";
        }
        cout << endl;
        if ( a != c )
        {
            cout << "a != c";
        }
        else
        {
            cout << "a == c";
        }
        cout << "   ";
        if ( a != b )
        {
            cout << "a != b";
        }
        else
        {
            cout << "a == b";
        }
        cout << endl;
    }  // try
    catch ( const RationalException& ex )
    {
        cerr << "RationalException caught: '";
        ex.display();
        cerr << "'" << endl;
    }  // catch RationalException
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }


    // Test of zero related stuff. Now exceptions MUST be thrown and caught!
    try
    {
        Rational a;
        a.set(3, 0);
        // The following line should not be displayed
        cerr << "Something is wrong, definition of a fraction with zero denominator succeeded" << a << endl;
    }
    catch ( const RationalException& zeroEx )
    {
        if ( RationalException::ZERO_DENOMINATOR == zeroEx.error )
        {
            cout << "Successfully detected attempt of a zero denominator" << endl;
        }
        else
        {
            cerr << "Unexpected exception: ";
            zeroEx.display();
            cerr << endl;
        }
    }

    try
    {
        // and also test one argument constuctor
        Rational zero(0);
        cout << "zero = " << zero << endl;
        zero.inverse();
        // the following should not be displayed
        cerr << "Something is wrong, inversion of zero succeeded" << zero << endl;
    }
    catch ( const RationalException& inversion )
    {
        if ( RationalException::UNINVERTIBLE == inversion.error )
        {
            cout << "zero is uninvertible" << endl;
        }
        else
        {
            cerr << "Unexpected exception: ";
            inversion.display();
            cerr << endl;
        }
    } // catch inversion

    // similar test, now invert() should "fail"
    try
    {
        // and also one argument constuctor
        Rational zero(0);
        Rational inv;
        inv = zero.invert();
        // the following should not be displayed
        cerr << "Something is wrong, inversion of zero succeeded: " <<inv << endl;
    }
    catch ( const RationalException& inversion )
    {
        if ( RationalException::UNINVERTIBLE == inversion.error )
        {
            cout << "Inversion of zero failed again. It's OK" << endl;
        }
        else
        {
            cerr << "Unexpected exception: ";
            inversion.display();
            cerr << endl;
        }
    } // catch inversion

    // Detect division by zero
    try
    {
        Rational a(3, 5);
        Rational b(0, 4);
        Rational c = a / b;
        // The following should not be displayed
        cerr << "Something is wrong, division by zero succeeded: " << c << endl;
    }
    catch ( const RationalException& zerEx )
    {
        if ( RationalException::DIVIDE_BY_ZERO == zerEx.error )
        {
            cout << "Attempt of divison by zero successfully detected." << endl;
        }
        else
        {
            cerr << "Unexpected exception: ";
            zerEx.display();
            cerr << endl;
        }
    }
}

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

        // product is not comutative...
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
    }
    catch ( const MatrixException& ex )
    {
        cerr << "MatrixException caught: '";
        ex.display();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}

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
        SqMatrixGeneric<Rational> inv = ((SqMatrixGeneric<Rational>) a).inverse();

        cout << "a:" << endl;
        a.display();
        cout << endl;

        // Rational inverse matrix
        cout << "inverse of a" << endl;
        inv.display();
        cout << endl;

        // finally check the matrix += operator
        a += inv;
        cout << "a + inv:" << endl;
        a.display();
        cout << endl;
    }
    catch ( const RationalException& rex )
    {
        cerr << "RationalException caught: '";
        rex.display();
        cerr << "'" << endl;
    }
    catch ( const MatrixException& mex )
    {
        cerr << "MatrixException caught: '";
        mex.display();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}

/*
 * Main function that starts three groups of unit tests
 */
int main(int argc, char* argv[])
{
    cout << "R A T I O N A L   T E S T" << endl << endl;
    rationalTest();
    cout << endl << "M A T R I X   T E S T" << endl << endl;
    matrixTest();
    cout << endl << "R A T I O N A L   M A T R I X   T E S T" << endl << endl;
    rationalMatrixTest();

    return 0;
}
