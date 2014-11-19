/*
Copyright 2011, 2013 Jernej Kovacic

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
 * Collection of basic unit tests for all mathematical classes
 */

//TODO rearrange tests, possibly into several files

#include "mtcopy.h"
#include "NumericUtil.h"
#include "Rational.h"
#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"
#include "QuaternionGeneric.h"
#include "PolynomialGeneric.h"
#include "LinearEquationSolverGeneric.h"
#include "PolynomialRegressionGeneric.h"
#include "PolynomialInterpolationGeneric.h"
#include "IntFactorization.h"
#include "IntCombinatorics.h"
#include "PermutationGeneric.h"
#include "CombinationGeneric.h"
#include "SampleStatGeneric.h"
#include "SampleQuantileGeneric.h"
#include "IntegGeneric.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <map>
#include <set>
#include <stdexcept>


using namespace std;
using namespace math;

template class MatrixGeneric<float>;
template class MatrixGeneric<Rational>;
template class QuaternionGeneric<float>;
template class PolynomialGeneric<float>;


/*
 * A handy function to print contents of a std::vector
 */
template <class T>
void printVector(const vector<T>& vec)
{
    for ( typename vector<T>::const_iterator it=vec.begin(); it!=vec.end(); ++it )
    {
        cout << *it << "  ";
    }
    cout << endl;
}

/*
 * Unit test of mtcopy functions
 */
void mtcopyTest()
{
    try
    {
    	const size_t N = 10;
        const int arr[] = { 23, 56, 12, 56, 23, -23, 92, -1, 45, -21 };
        vector<int> v;
        vector<int> v2;

        cout << "Copy the entire array: " << endl;
        mtcopy(arr, arr+N, v);
        cout << "Copy:     ";
        printVector(v);
        cout << "Expected: 23  56  12  56  23  -23  92  -1  45  -21" << endl;

        cout << "Copy a part of an array: " << endl;
        mtcopy(arr+1, 7, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: 56  12  56  23  -23  92  -1" << endl;

        cout << "Copy the entire vector: " << endl;
        sort(v. begin(), v.end());
        mtcopy(v, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: -23  -21  -1  12  23  23  45  56  56  92" << endl;

        cout << "Copy of a part of a vector: " << endl;
        mtcopy(v, 1, 5, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: -21  -1  12  23  23 " << endl;

        cout << "Test vector's range: " << endl;
        mtcopy(v, 7, 5, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: 56  56  92 " << endl;
    }
    catch ( const bad_alloc& ba )
    {
        cerr << "Bad allocation exception called" << endl;
    }
    catch ( const out_of_range& oor )
    {
        cerr << "Out of range exception caught" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}


/*
 * Test of the class Rational
 */
void rationalTest()
{
    try
    {
        // Parsing test
        Rational s("-12");
        cout << "s = -12 = " << s << "\texpected -12/1" << endl;
        s.set("34.824");
        cout << "s = 34.824 = " << s << "\texpected 4353/125" << endl;
        s.set("-4.285714", 6);
        cout << "s = -4.|285714 = " << s << "\texpected -30/7" << endl;
        s.set("3,5167", 3);
        cout << "s = 3.5|167 = " << s << "\texpected 17566/4995" << endl;
        cout << endl;
        
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
        ex.what();
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
            zeroEx.what();
            cerr << endl;
        }
    }

    try
    {
        // and also test one argument constructor
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
            inversion.what();
            cerr << endl;
        }
    } // catch inversion

    // similar test, now invert() should "fail"
    try
    {
        // and also one argument constructor
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
            inversion.what();
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
            zerEx.what();
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

        inv *= a;
        cout << "deteminant of inv * (a+inv): " << inv.determinant() << endl;
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

/*
 * Test of float quaternions
 */
void quaternionTest()
{
    try
    {
        FQuaternion zeroq;
        FQuaternion q(1.0f, -0.5f, 1.2f, -2.3f);

        cout << "zero: " << zeroq << endl;
        cout << "q: " << q << endl;
        cout << "q: 1: " << q.getOne() << "   i: " << q.getI() << "   j: " << q.getJ() << "   k: " << q.getK() << endl;

        FQuaternion uq = q.unit();
        cout << "U(q): " << uq << "\tnorm(U(q)): " << uq.norm() << endl;

        FQuaternion o(1.0f);
        FQuaternion i;   i.setI(1.0f);
        FQuaternion j;   j.setJ(1.0f);
        FQuaternion k;   k.setK(1.0f);

        cout << "i*j: " << i*j << "\tj*i: " << j*i << endl;
        cout << "j*k: " << j*k << "\tk*j: " << k*j << endl;
        cout << "k*i: " << k*i << "\ti*k: " << i*k << endl;

        FQuaternion p(1.0f, 2.0f, 3.0f, 4.0f);
        FQuaternion sum = p+q;
        FQuaternion dif = p-q;
        cout << "p: " << p << "\tq: " << q << endl;
        cout << "p+q: " << sum << "\tp-q: " << dif << endl;
        // Must output (7.6-10.2i+6.8j+5.6k) and (7.6+13.2i+1.6j-2.2k), respectively
        cout << "p*q: " << p*q << "\tq*p: " << q*p << endl;
        cout << "p+5: " << p+5.f << "\t3+p: " << 3.0f + p << endl;
        cout << "p-4: " << p-4.f << "\t1-p: " << 1.0f - p << endl;
        p += 2;
        cout << "p+2: " << p;
        p -= 2.f;
        cout << "\tp+2-2: " << p << endl;

        FQuaternion qc = -0.5f*(q+i*q*i+j*q*j+k*q*k);
        cout << "-0.5f*(q+i*q*i+j*q*j+k*q*k)=" << qc << "  q conj: " << qc << endl;
        FQuaternion qrec = q.reciprocal();
        cout << "qrec: " << qrec << "\tq*qrec=" << q*qrec << "\tqrec*q=" << qrec*q << endl;
    } // try
    catch ( const QuaternionException& qex )
    {
        cerr << "QuaternionException caught: '";
        qex.what();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}

/*
 * Test of polynomial algebra
 */
void polynomialTest()
{
    try
    {
        float a[] = { 2.1f, 1.0f, -0.72f, 1.0f, 0.0f };
        FPolynomial t(a, 5);
        cout << "t(x) = " << t << endl;
        cout << "a0=" << t.get(0) << " a1=" << t.get(1) << " a2=" << t.get(2) << " a3=" << t.get(3) << endl;
        cout << "Remove the 2nd coef. from t(x): " << t.remove(1) << endl;
        cout << "Insert 0.2 to the 3rd pos. of t(x): " << t.insert(2, 0.2f) << endl;
        cout << "Insert 5 to the 8th pos. of t(x): " << t.insert(7, 5.0f) << endl;

        cout << "t's coefficients in reverse order: ";
        vector<float> desc;
        t.getDesc(desc);
        for ( typename vector<float>::const_iterator it=desc.begin();
              it!=desc.end(); ++it)
        {
            cout << *it << " ";
        }
        cout << endl;
        t.get(desc);
        FPolynomial trev;
        trev.setDesc(desc);
        cout << "t reversed: " << trev << endl;

        FPolynomial z(true, 4);
        z.set(3, 0.0f);
        cout << "Zero polynomial: " << z << endl;

        FPolynomial p(true, 6), q(true, 4);

        p.set(0, 2.0f).set(1, 3.0f).set(2, -4.0f).set(3, -7.0f).set(4, 2.0f).set(5, 1.0f);
        q.set(0, -1.0f).set(1, 5.0f).set(2, -3.0f).set(3, 1.0f);

        cout << "p(x)=" << p << endl << "q(x)=" << q << endl;
        cout << "p(-1.2) = " << p.value(-1.2f) << endl;
        cout << "q(0.7) = " << q.value(0.7f) << endl;
        cout << "p+q: " << p+q << endl << "q+p: " << q+p <<endl;
        cout << "p-q: " << p-q << endl << "q-p: " << q-p << endl;
        cout << "p*q: " << p*q << endl << "q*p: " << q*p << endl;

        cout << "0.3 * p(x) = " << p*0.3f << endl;
        cout << "0.5 * q(x) = " << 0.5f*q << endl;
        cout << "p(x) + 2 = " << p+2.0f << endl;
        cout << "2 + p(x) = " << 2.0f+p << endl;
        cout << "q(x) - 3 = " << q-3.0f << endl;
        cout << "4 - q(x) = " << 4.0f-q << endl;
        p +=10.0f;
        cout << "p(x) + 10 = " << p << endl;
        p -= 10.0f;
        cout << "p(x) + 10 - 10 = " << p << endl;

        cout << "p(x)' = " << p.deriv() << endl;
        FPolynomial qint = q.integ();
        cout << "I(q(x))dx = " << qint << endl;
        cout << "int(q(x))' = " << qint.deriv() << endl;
        cout << "q / 5 = " << q/5.0f << endl;
        cout << "q mod 5 = " << q%5.f << endl;
        cout << "2 / q = " << 2.0f/q << endl;
        cout << "2 mod q = " << 2.0f%q << endl;
        q /= 2.0f;
        cout << "q / 2 = " << q << endl;
        q %= 2.0f;
        cout << "q mod 2 = " << q << endl;

        p = FPolynomial(true, 6);
        p.set(0, -1.f).set(1, 0.f).set(2, 1.f).set(3, 2.f).set(4, -1.f).set(5, 4.f);
        q = FPolynomial(true, 3);
        q.set(0, 1.f).set(1, 0.f).set(2, 1.f);
        cout << "p = " << p << endl;
        cout << "q = " << q << endl;
        // Expected: "2 -2*x -1*x^2 +4*x^3"
        cout << "p / q = " << p/q << endl;
        // Expected: "-3 +2*x"
        cout << "p mod q = " << p%q << endl;

        p = FPolynomial(true, 4);
        p.set(0, -4.f).set(1, 0.f).set(2, -2.f).set(3, 1.f);
        q = FPolynomial(true, 2);
        q.set(0, -3.f).set(1, 1.f);
        cout << "p = " << p << endl;
        cout << "q = " << q << endl;
        // Expected: "3 +1*x +1*x^2"
        cout << "p / q = " << p/q << endl;
        // Expected: "3 +1*x +1*x^2"
        cout << "p mod q = " << p%q << endl;

        p = FPolynomial(true, 5);
        p.set(0, -5.f).set(1, 3.f).set(2, 0.f).set(3, -6.f).set(4, 4.f);
        q = FPolynomial(true, 2);
        q.set(0, -1.f).set(1, 2.f);
        FPolynomial ptemp = p;
        cout << "p = " << p << endl;
        cout << "q = " << q << endl;
        ptemp /= q;
        // Expected: "1 -1*x -2*x^2 +2*x^3"
        cout << "p / q = " << ptemp << endl;
        ptemp = p;
        ptemp %= q;
        // Expected: "-4"
        cout << "p mod q = " << ptemp << endl;
    }
    catch ( const PolynomialException& pex )
    {
        cerr << "PolynomialException caught: '";
        pex.what();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }

}

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

        LinearEquationSolverGeneric<complex<float> > solver;
        solver.setCoef(a).setTerm(b);

        /* 
         * Exact solution of the linear equation system a*x=b:
         * x  = [-0.6631640-0.3626125i, 0.1630189+1.050566i, -0.2240929+0.6002903i]'
         */
        MatrixGeneric<complex<float> > x;
        solver.solve(x);
        cout << "Solution:" << endl;
        x.display();
    }
    catch ( const LinearEquationSolverException& leqex )
    {
        cerr << "LinearEquationSolverException caught: '";
        leqex.what();
        cerr << "'" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}

/*
 * Test of curve fitting algorithms
 */
void curveFittingTest()
{
    try
    {
        // 1st, 2nd and 3rd degree regression polynomials, respectively,
        // of y=exp(-x) for x=0..5:
        FPolynomialRegression fprexp1;
        FPolynomialRegression fprexp2;
        FPolynomialRegression fprexp3;

        // exp(-x), points in a random order:
        fprexp1.enterPoint(1.0f, 0.367879441f);
        fprexp1.enterPoint(0.0f, 1.0f);
        fprexp1.enterPoint(5.0f, 0.006737947f);
        fprexp1.enterPoint(4.0f, 0.018315639f);
        fprexp1.enterPoint(2.0f, 0.135335283f);
        fprexp1.enterPoint(3.0f, 0.049787068f);

        // copy points to other classes:
        fprexp2 = fprexp1;
        fprexp3 = fprexp1;

        // Checking of bounds
        cout << "Lower bound: " << fprexp2.lowerBound() << endl;
        cout << "Upper bound: " << fprexp2.upperBound() << endl;

        fprexp1.generateCurve(1);
        fprexp2.generateCurve(2);
        fprexp3.generateCurve(3);

        /*
         Correct regression polynomials as calculated by R:

            > x <- c(0,1,2,3,4,5)
            > y <- exp(-x)
            > print(x)
            [1] 0 1 2 3 4 5

            > print(y)
            [1] 1.000000000 0.367879441 0.135335283 0.049787068 0.018315639 0.006737947

            > p1 <- lm(y ~ x)
            > p2 <- lm(y ~ x + I(x^2))
            > p3 <- lm(y ~ x + I(x^2) + I(x^3))
            > print(p1)

            Call:
            lm(formula = y ~ x)

            Coefficients:
            (Intercept)            x
                 0.6988      -0.1743

            > print(p2)

            Call:
            lm(formula = y ~ x + I(x^2))

            Coefficients:
            (Intercept)            x       I(x^2)
                0.93132     -0.52314      0.06977

            > print(p3)

            Call:
            lm(formula = y ~ x + I(x^2) + I(x^3))

            Coefficients:
            (Intercept)            x       I(x^2)       I(x^3)
                0.99180     -0.79932      0.22096     -0.02016
        */

        cout << "1st degree polynomial: ";
        fprexp1.getPolynomial().display();
        cout << endl << "2nd degree polynomial: ";
        fprexp2.getPolynomial().display();
        cout << endl << "3rd degree polynomial: ";
        fprexp3.getPolynomial().display();
        cout << endl;

        /*
             Example from: http://en.wikipedia.org/wiki/Newton_polynomial#Example

             y=tan(x) for x=-1.5, -0.75, 0, 0.75, 1.5

             This time, 3rd degree regression and the same degree interpolation
             polynomial (as this is an odd function, all its even coefficients
             equal 0) are compared.
         */
        FPolynomialInterpolation fpitan;
        FPolynomialRegression fprtan3;

        fpitan.enterPoint(-1.5f, -14.1014f);
        fpitan.enterPoint(-0.75f, -0.931596f);
        fpitan.enterPoint(0.0f, 0.0f);
        fpitan.enterPoint(0.75f, 0.931596);
        fpitan.enterPoint(1.5f, 14.1014f);

        fpitan.generateCurve();
        fprtan3.copy(&fpitan);
        fprtan3.generateCurve(3);

        cout << "Interpolation polynomial: ";
        fpitan.getPolynomial().display();
        cout << endl << "3rd degree regression polynomial: ";
        fprtan3.getPolynomial().display();
        cout << endl;

        /*
             Finally compare 5th degree interpolation and regression polynomials
             for the points from the first test:
         */
        FPolynomialRegression fprexp5;
        FPolynomialInterpolation fpiexp;
        fprexp5.copy(&fprexp1);
        fpiexp.copy(&fprexp5);
        fprexp5.generateCurve(5);
        fpiexp.generateCurve();
        cout << "Exp's regression polynomial: ";
        fprexp5.getPolynomial().display();
        cout << endl << "Exp's interpolation polynomial: ";
        fpiexp.getPolynomial().display();
        cout << endl;
    }
    catch ( const CurveFittingException& cfex )
    {
        cerr << "CurveFittingException caught: ";
        cfex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}

/*
 * Test of the exponentiation algorithm
 */
void intExponentiaorTest()
{
    
    try
    {
        cout << "5^0 = " << NumericUtil<int>::power(5, 0) << endl;
        // 1024
        cout << "2^10 = " << NumericUtil<int>::power(2, 10) << endl;
        // -129140163
        cout << "(-3)^17 = " << NumericUtil<int>::power(-3, 17) << endl;
        // 244140625
        cout << "5^12 = " << NumericUtil<int>::power(5, 12) << endl;
        // 2401
        cout << "(-7)^4 = " << NumericUtil<int>::power(-7, 4) << endl;
        // 32768
        cout << "sqrt(2)^30 = " << NumericUtil<double>::power(sqrt(2.0), 30) << endl;
        cout << endl;
        
        FSqMatrix m(3);
        m.set(0, 0, 0.1f).set(0, 1, -0.2f).set(0, 2, 0.3f);
        m.set(1, 0, -0.4f).set(1, 1, 0.5f).set(1, 2, -0.6f);
        m.set(2, 0, 0.7f).set(2, 1, -0.8f).set(2, 2, 0.9f);
        cout << "m =" << endl;
        m.display();
        cout << "m^0 =" << endl;
        NumericUtil<FSqMatrix>::power(m, 0).display();
        /*
               | 1.21824   -1.49688   1.77552 |
               |-2.75886    3.38985  -4.02084 |
               | 4.29948   -5.28282   6.26616 |
         */
        cout << "m^5 =" << endl;
        NumericUtil<FSqMatrix>::power(m, 5).display();
        cout << endl;
        
        FQuaternion q(1.0f, -0.8f, 1.2f, -1.5f);
        cout << "q = " << q << endl;
        cout << "q^0 = " << NumericUtil<FQuaternion>::power(q, 0) << endl;
        // (991.414+1609.27i-2413.91j+3017.39k)
        cout << "q^10 = " << NumericUtil<FQuaternion>::power(q, 10) << endl;
        cout << endl;
        
        Rational f(-3, 4);
        cout << "f = " << f << endl;
        cout << "f^0 = " << NumericUtil<Rational>::power(f, 0) << endl;
        // -2187/16384
        cout << "f^7 = " << NumericUtil<Rational>::power(f, 7) << endl;
        cout << endl;
        
        FPolynomial p(true, 3);
        p.set(0, 0.2f).set(1, -0.5f);
        cout << "p(x) = " << p << endl;
        cout << "p(x)^0 = " << NumericUtil<FPolynomial>::power(p, 0) << endl;
        /*
           6.4e-05 -0.00096*x +0.00792*x^2 -0.044*x^3 +0.1815*x^4 -0.5775*x^5 
           +1.45062*x^6 -2.8875*x^7 +4.5375*x^8 -5.5*x^9 +4.95*x^10 -3*x^11 +1*x^12
         */
        cout << "p(x)^6 = " << NumericUtil<FPolynomial>::power(p, 6) << endl;
    }
    catch ( IMathException& ex )
    {
        cerr << "Math exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}

/*
 * Test of integer factorization
 */
void intFactorizationTest()
{
    try
    {
        // Find integer square roots of the following numbers:
        const unsigned int NS = 3;
        unsigned long int s[NS] = { 12, 100, 37423 };

        // Expected results: 3, 10 and 193
        for ( unsigned int i=0; i<NS; ++i )
        {
            cout << "Int sqrt of " << s[i] << ": " << IntFactorization::intSqrt(s[i]) << endl;
        }
        cout << endl;

        // Test primality test and the algorithm for the next prime,
        // Reference: http://primes.utm.edu/lists/small/10000.txt
        const unsigned int NP = 10;
        unsigned long int p[NP] = { 3, 6, 15, 37, 257, 703, 907, 101861, 102601, 104597 };
        /*
         * Expected results:
         * (T, 5), (F, 7), (F, 17), (T, 41), (T, 263), (F, 709), (T, 911)
         * (F, 101863), (F, 102607), (T, 104623)
         */
        for ( unsigned int i=0; i<NP; ++i )
        {
            cout << p[i] << " is ";
            cout << ( true==IntFactorization::isPrime(p[i]) ? "" : "not ");
            cout << "a prime\tNext prime: ";
            cout << IntFactorization::nextPrime(p[i]) << endl;
        }
        cout << endl;

        // Test prime factorization and divisors
        const unsigned int NF = 5;
        unsigned long nf[NF] = {245, 6784, 21737, 195327 ,3428543 };
        /*
         * Expected results:
         * 245 = 5 * 7^2      [1, 5, 7, 35, 49, 245]
         * 6784 = 2^7 * 53    [1, 2, 4, 8, 16, 32, 53, 64, 106, 128, 212, 424, 848, 1696, 3392, 6784]
         * 21737 = 21737      [1, 21737]
         * 195327 = 3^2 * 11 * 1973     [1, 3, 9, 11, 33, 99, 1973, 5919, 17757, 21703, 65109, 195327]
         * 3428543 = 17 * 41 * 4919     [1, 17, 41, 697, 4919, 83623, 201679, 3428543]
         */
        for ( unsigned int i=0; i<NF; ++i )
        {
            map<unsigned long int, unsigned int> f;
            IntFactorization::factor(nf[i], f);
            cout << nf[i] << " = ";
            for ( map<unsigned long int, unsigned int>::iterator it=f.begin();
                  it!=f.end(); ++it )
            {
                if ( it!=f.begin() )
                {
                    cout << " * ";
                }

                cout << it->first;

                if ( it->second > 1 )
                {
                    cout << '^' << it->second;
                }
            }

            set<unsigned long int> d;
            IntFactorization::divisors(nf[i], d);
            cout << "\tdivisors: ";
            for ( set<unsigned long int>::iterator it=d.begin(); it!=d.end(); ++it )
            {
                if ( it!=d.begin() )
                {
                    cout << ", ";
                }
                cout << *it;
            }

            cout << endl;
        }
        cout << endl;

        // Finally test the GCD and LCM algorithms:
        const unsigned int NC = 10;
        unsigned long int c[NC] = { 500, 1000, 85, 3428543 , 15, 100, 3, 57, 234, 7643 };
        /*
         * Expected results:
         * n1 = 500, n2 = 1000: gcd = 500, lcm = 1000
         * n1 = 85, n2 = 3428543: gcd = 17, lcm = 17142715
         * n1 = 15, n2 = 100: gcd = 5, lcm = 300
         * n1 = 3, n2 = 57: gcd = 3, lcm = 57
         * n1 = 234, n2 = 7643: gcd = 1, lcm = 1788462
         */
        for ( unsigned int i=0; i<NC/2; ++i )
        {
            const unsigned long int n1 = c[2*i];
            const unsigned long int n2 = c[2*i+1];
            unsigned long int gcd = IntFactorization::greatestCommonDivisor(n1, n2);
            unsigned long int lcm = IntFactorization::leastCommonMultiple(n1, n2);

            cout << "n1 = " << n1 << ", n2 = " << n2 << ": ";
            cout << "gcd = " << gcd << ", lcm = " << lcm << ", ";
            cout << "gcd*lcm = " << gcd*lcm << ", ";
            cout << "n1*n2 = " << n1*n2;
            cout << endl;
        }

        cout << endl;
    }
    catch ( const IntFactorizationException& ex )
    {
        cerr << "IntFactorization exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}

/*
 * Test of integer combinatorics
 */
void intCombinatoricsTest()
{
    try
    {
        cout << "5! = " << IntCombinatorics::factorial(5) << "  expected: 120" << endl;
        cout << "20! = " << IntCombinatorics::factorial(20) << "  expected: 2432902008176640000 " << endl;
        cout << "15!/5! = " << IntCombinatorics::factorial(15, 6) << "  expected: 10897286400" << endl;
        cout << endl;
        
        cout << "falling factorial: (10)_4 = " << IntCombinatorics::fallingFactorial(10, 4) << "  expected:  5040" << endl;
        cout << "falling factorial: (12)_8 = " << IntCombinatorics::fallingFactorial(12, 8) << "  expected:  19958400" << endl;
        cout << "rising factorial: 5^(6) = " << IntCombinatorics::risingFactorial(5, 6) << "  expected :  151200" << endl;
        cout << "rising factorial: 16^(10) = " << IntCombinatorics::risingFactorial(16, 10) << "  expected :  11861676288000" << endl;
        cout << endl;
        
        // Checked with this Maxima function:
        // multif(n,k) := if (n<k) then 1 else n*multif(n-k,k)$
        cout << "18!^(3) = " << IntCombinatorics::multiFactorial(18, 3) << "  expected:  524880" << endl;
        cout << "19!^(3) = " << IntCombinatorics::multiFactorial(19, 3) << "  expected:  1106560" << endl;
        cout << "20!^(3) = " << IntCombinatorics::multiFactorial(20, 3) << "  expected:  2094400" << endl;
        cout << "21!^(3) = " << IntCombinatorics::multiFactorial(21, 3) << "  expected:  11022480" << endl;
        cout << "15!! = " << IntCombinatorics::doubleFactorial(15) << "  expected:  2027025" << endl;
        cout << "22!! = " << IntCombinatorics::doubleFactorial(22) << "  expected:  81749606400" << endl;
        cout << "27!! = " << IntCombinatorics::doubleFactorial(27) << "  expected:  213458046676875" << endl;
        cout << endl;
        
        cout << "14 choose 4 = " << IntCombinatorics::binom(14, 4) << "  expected: 1001" << endl;
        cout << "14 choose 10 = " << IntCombinatorics::binom(14, 10) << "  expected: 1001" << endl;
        cout << "50 choose 41 = " << IntCombinatorics::binom(50, 41) << "  expected:  2505433700" << endl;
    }
    catch ( const CombinatoricsException& ex )
    {
        cerr << "Combinatorics exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}

/*
 * Test of classes that list all permutations and combinations of a sequence
 */
void combinatoricsTest()
{
    try
    {
        cout << "Permutations:" << endl;
        
        vector<char> chvect;
        for ( char ch='a'; ch<='e'; ++ch )
        {
            chvect.push_back(ch);
        }
        
        PermutationGeneric<char> perm(chvect);
        list<list<char> > lp;
        size_t cnt = 1;
        while ( true==perm.hasNext() )
        {
            perm.next(lp, 25);
            for ( list<list<char> >::const_iterator lit=lp.begin(); lit!=lp.end(); ++lit )
            {
                cout << cnt++ << ": ";
                for ( list<char>::const_iterator cit=lit->begin(); cit!=lit->end(); ++cit )
                {
                    cout << *cit;
                }
                cout << endl;
            }
        }
        
        cout << endl << "Combinations:" << endl;
        CombinationGeneric<char> comb(chvect);
        list<set<char> > lc;
        for ( size_t k=1; k<=chvect.size(); ++k )
        {
            cout << endl << "K = " << k << " : " << endl << endl;
            comb.setK(k);
            cnt = 1;
            while ( true==comb.hasNext() )
            {
                comb.next(lc, 5);
                for ( list<set<char> >::const_iterator lit=lc.begin(); lit!=lc.end(); ++lit )
                {
                    cout << cnt++ << ": ";
                    for ( set<char>::const_iterator cit=lit->begin(); cit!=lit->end(); ++cit )
                    {
                        cout << *cit;
                    }
                    cout << endl;
                }
            }
        }
        
        cout << endl;
    }
    catch ( const CombinatoricsException& ex )
    {
        cerr << "Combinatorics exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}


/*
 * Test of numerical integration algorithms
 */
void numIntegTest()
{
    class CFunc : public IIntegFunction
    {
    public:
        double ke, kl, n;

        double func(const double& x) const throw(CalculusException)
        {
            return ke * std::exp(-(x-3)*(x-3)) + kl * x + n;
        }
    };

    try
    {
        CFunc f;
        f.ke = 0.7;   f.kl = 1.0;  f.n = 0.5;

        /*
         * Exact result obtained by Maxima:

            (%i1) float(integrate(0.7*exp(-(x-3)^2) + 1*x + 0.5, x, 0, 5));
            (%o1) 16.23780211731536
         */
        cout << "f(x) = " << f.ke << "*exp(-(x-3)^2) " << showpos << f.kl << "*x " << f.n << noshowpos << endl;
        for ( int method=EIntegAlg::RECTANGLE; method<=EIntegAlg::SIMPSON_3_8; ++method )
        {
            cout << "Method " << method << ": Int(f(x), 0, 5) = " <<
                     Integ::integ(f, 0., 5., 10000, static_cast<EIntegAlg::alg>(method)) << endl;
        }
        cout << "Expected result: 16.23780211731536" << endl;
    }
    catch ( const CalculusException& iex )
    {
        cerr << "Calculus exception caught: ";
        iex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}


/*
 * Test of classes that perform statistical operations
 */
void statisticsTest()
{
    try
    {
        /* Cars' mpg (miles per gallon) from R's data frame 'mtcars' */
        const double ampgs[] =
            { 21.0, 21.0, 22.8, 21.4, 18.7, 18.1, 14.3, 24.4,
              22.8, 19.2, 17.8, 16.4, 17.3, 15.2, 10.4, 10.4,
              14.7, 32.4, 30.4, 33.9, 21.5, 15.5, 15.2, 13.3,
              19.2, 27.3, 26.0, 30.4, 15.8, 19.7, 15.0, 21.4 };

        /* Cars' wt (weight in 1000 lbs) from R's data frame 'mtcars' */
        const double awts[] =
            { 2.620, 2.875, 2.320, 3.215, 3.440, 3.460, 3.570, 3.190,
              3.150, 3.440, 3.440, 4.070, 3.730, 3.780, 5.250, 5.424,
              5.345, 2.200, 1.615, 1.835, 2.465, 3.520, 3.435, 3.840,
              3.845, 1.935, 2.140, 1.513, 3.170, 2.770, 3.570, 2.780 };

        /*
           # Equivalent of the following command in R:
           data(mtcars)
         */
        const size_t LEN = 32;  // Nr. of elements in mtcars
        vector<double> vmpgs;
        vector<double> vwts;

        math::mtcopy(ampgs, LEN, vmpgs);
        math::mtcopy(awts, LEN, vwts);

        /*
           min(mtcars$wt)
           [1] 1.513
           max(mtcars$wt)
           [1] 5.424
         */
        cout << "min(wt): " << SampleStat::min(vwts) << " (expected: 1.513)" << endl;
        cout << "max(wt): " << SampleStat::max(vwts) << " (expected: 5.424)" << endl;
        cout << endl;

        /*
           sum(mtcars$mpg)
           [1] 642.9
           mean(mtcars$mpg)
           [1] 20.09062
           var(mtcars$mpg)
           [1] 36.3241
           sd(mtcars$mpg)
           [1] 6.026948
           n <- length(mtcars$mpg)
           (n-1)/n * var(mtcars$mpg)
           [1] 35.18897
           n <- length(mtcars$mpg)
           sqrt((n-1)/n) * sd(mtcars$mpg)
           [1] 5.93203
         */
        cout << "Sum of all elements: " << SampleStat::sum(vmpgs) << " (expected: 642.9)" << endl;
        cout << "Sample mean: " << SampleStat::mean(vmpgs) << " (expected: 20.09062)" << endl;
        cout << "Sample variance: " << SampleStat::var(vmpgs) << " (expected: 36.3241)" << endl;
        cout << "Sample standard deviation: " << SampleStat::stdev(vmpgs) << " (expected: 6.026948)" << endl;
        cout << "Population variance (w/o Bessel's correction): " << SampleStat::var(vmpgs, false) << " (expected: 35.18897)" << endl;
        cout << "Population standard deviation (w/o Bessel's correction): " << SampleStat::stdev(vmpgs, false) << " (expected: 5.93203)" << endl;

        /*
           cov(mtcars$mpg, mtcars$wt)
           [1] -5.116685
           (n-1)/n * cov(mtcars$mpg, mtcars$wt)
           [1] -4.956788
           cor(mtcars$mpg, mtcars$wt)
           [1] -0.8676594
           cor(mtcars$mpg, mtcars$wt)^2
           [1] 0.7528328
         */
        cout << "Sample covariance: " << SampleStat::cov(vmpgs, vwts) << " (expected: -5.116685)" << endl;
        cout << "Population covariance (w/o B.c.): " << SampleStat::cov(vmpgs, vwts, false) << " (expected: -4.956788)" << endl;
        cout << "Pearson's r: " << SampleStat::cor(vmpgs, vwts) << " (expected: -0.8676594)" << endl;
        cout << "r^2: " << SampleStat::r2(vmpgs, vwts) << " (expected: 0.7528328)" << endl;

        /*
         * R code to perform basic unit test of quantiles:

           data(mtcars)
           median(mtcars$mpg)
           [1] 19.2
           quantile(mtcars$mpg, c(0.25, 0.75))
              25%    75%
           15.425 22.800
           IQR(mtcars$mpg)
           [1] 7.375
           quantile(mtcars$mpg, 0.63)
              63%
           21.212
         */
        SampleQuantile q(vmpgs);
        cout << "Median: " << q.median() << " (expected: 19.2)" << endl;
        cout << "1st quartile: " << q.quantile(1, 4) << " (expected: 15.425)" << endl;
        cout << "3rd quartile: " << q.quantile(3, 4) << " (expected: 22.800)" << endl;
        cout << "IQR: " << q.iqr() << " (expected: 7.375)" << endl;
        cout << "63th percentile: " << q.qntl(0.63) << " (expected: 21.212)" << endl;

        const double probs[] = { 0.01, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 0.99 };
        const size_t N_PROBS = 9;
        const char* exp[] =
        {
            "10.4\t14.3\t15.2\t17.3\t19.2\t21.0\t22.8\t30.4\t33.9",
            "10.40\t14.30\t15.35\t17.55\t19.20\t21.20\t22.80\t30.40\t33.90",
            "10.4\t13.3\t15.2\t17.3\t19.2\t21.0\t22.8\t30.4\t33.9",
            "10.40\t13.50\t15.20\t17.30\t19.20\t21.00\t22.80\t29.78\t33.42",
            "10.40\t14.00\t15.35\t17.55\t19.20\t21.20\t22.80\t30.40\t33.90",
            "10.4000\t13.6000\t15.2750\t17.4875\t19.2000\t21.2500\t22.8000\t30.4000\t33.9000",
            "10.4000\t14.3400\t15.4250\t17.6125\t19.2000\t21.1500\t22.8000\t30.0900\t33.4350",
            "10.40\t13.8667\t15.325\t17.5292\t19.200\t21.2167\t22.80\t30.40\t33.90",
            "10.40\t13.90\t15.3312\t17.5344\t19.20\t21.2125\t22.80\t30.40\t33.90"
        };

        /*
         * R code to test various probabilities and various methods:
         *
             data(mtcars);   x <- mtcars$mpg
             types <- 1:9
             p <- c(0.01, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 0.99)
             for (t in types)
             {
               print(quantile(x, p, type=t, names=FALSE))
             }
         */

        cout << "Test of various quantile methods:" << endl;

        // Not really the best practice, but as long as the enum is contiguous...
        for ( int type=EQntlType::R1; type<=EQntlType::R9; ++type )
        {
            cout << "R" << 1+type << ":\t";
            for ( size_t i=0; i<N_PROBS; ++i )
            {
                cout << q.qntl(probs[i], static_cast<EQntlType::type>(type));
                if ( i < N_PROBS-1 )
                {
                    cout << "\t";
                }
            }
            cout << endl;
            cout << "Exp.:\t" << exp[type] << endl;
        }
        cout << endl;

        /*
         * R code to find outliers:
         *
             q <- quantile(mtcars$mpg, c(0.25, 0.75), name=FALSE)
             d <- IQR(mtcars$mpg)
             mtcars$mpg < (q[1]-1.0*d) | mtcars$mpg > (q[2]+1.0*d)
              [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
             [14] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
             [27] FALSE  TRUE FALSE FALSE FALSE FALSE

             sort(unique(mtcars$mpg[mtcars$mpg<(q[1]-0.5*d) | mtcars$mpg>(q[2]+0.5*d)]))
             [1] 10.4 27.3 30.4 32.4 33.9
         */
        cout << "Min mpg: " << q.min() << " (expected 10.4)" << endl;
        cout << "Max mpg: " << q.max() << " (expected 33.9)" << endl;
        typename vector<double>::const_iterator mpgit;
        for ( mpgit=vmpgs.begin(); mpgit!=vmpgs.end(); ++ mpgit )
        {
            cout << *mpgit << " in range [8.05, 30.175]: ";
            cout << q.isOutlier(*mpgit, 1.0) << endl;
        }
        cout << "Outliers for iqr=0.5: [";
        set<double> oul;
        q.outliers(oul, 0.5);
        typename set<double>::const_iterator oit;
        for ( oit=oul.begin(); oit!=oul.end(); ++oit )
        {
            cout << *oit << " ";
        }
        cout << "] expected [10.4 27.3 30.4 32.4 33.9]" << endl;
    }
    catch ( const StatisticsException& ex )
    {
        cerr << "Statistics exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch ( const bad_alloc& ba )
    {
        cerr << "Could not allocate memory." << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}

/*
 * Main function that starts several groups of unit tests
 */
int main(int argc, const char* argv[])
{

    cout << "M T C O P Y   T E S T" << endl << endl;
    mtcopyTest();

    cout << endl << "Q U A T E R N I O N   T E S T" << endl << endl;
    quaternionTest();

    cout << endl << "R A T I O N A L   T E S T" << endl << endl;
    rationalTest();

    cout << endl << "M A T R I X   T E S T" << endl << endl;
    matrixTest();

    cout << endl << "R A T I O N A L   M A T R I X   T E S T" << endl << endl;
    rationalMatrixTest();

    cout << endl << "L I N E A R   E Q U A T I O N   S O L V E R   T E S T" << endl << endl;
    lineqSolverTest();

    cout << endl << "P O L Y N O M I A L   T E S T" << endl << endl;
    polynomialTest();

    cout << endl << "C U R V E   F I T T I N G   T E S T" << endl << endl;
    curveFittingTest();
    
    cout << endl << "I N T E X P O N E N T I A T O R   T E S T" << endl << endl;
    intExponentiaorTest();
    
    cout << endl << "I N T   F A C T O R I Z A T I O N   T E S T" << endl << endl;
    intFactorizationTest();

    cout << endl << "I N T   C O M B I N A T O R I C S   T E S T" << endl << endl;
    intCombinatoricsTest();
    
    cout << endl << "C O M B I N A T O R I C S   T E S T" << endl << endl;
    combinatoricsTest();

    cout << endl << "N U M   I N T E G   T E S T" << endl << endl;
    numIntegTest();

    cout << endl << "S T A T I S T I C S   T E S T" << endl << endl;
    statisticsTest();

    return 0;

    (void) argc;
    (void) argv;
}
