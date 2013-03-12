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
 * @file maintest.cpp
 *
 * Collection of basic unit tests for all mathematical classes
 *
 * @author Jernej Kovacic
 */

#include "Rational.h"
#include "MatrixGeneric.h"
#include "SqMatrixGeneric.h"
#include "QuaternionGeneric.h"
#include "PolynomialGeneric.h"
#include "LinearEquationSolverGeneric.h"
#include "PolynomialRegressionGeneric.h"
#include "PolynomialInterpolationGeneric.h"

#include <iostream>
#include <complex>

using namespace std;
using namespace math;

template class MatrixGeneric<float>;
template class MatrixGeneric<Rational>;
template class QuaternionGeneric<float>;
template class PolynomialGeneric<float>;

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

        inv *= a;
        cout << "deteminant of inv * (a+inv): " << inv.determinant() << endl;
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
        cout << "p+5: " << p+5 << "\t3+p: " << 3.0f + p << endl;
        cout << "p-4: " << p-4 << "\t1-p: " << 1.0f - p << endl;
        p += 2;
        cout << "p+2: " << p;
        p -= 2;
        cout << "\tp+2-2: " << p << endl;

        FQuaternion qc = -0.5f*(q+i*q*i+j*q*j+k*q*k);
        cout << "-0.5f*(q+i*q*i+j*q*j+k*q*k)=" << qc << "  q conj: " << qc << endl;
        FQuaternion qrec = q.reciprocal();
        cout << "qrec: " << qrec << "\tq*qrec=" << q*qrec << "\tqrec*q=" << qrec*q << endl;
    } // try
    catch ( const QuaternionException& qex )
    {
        cerr << "QuaternionException caught: '";
        qex.display();
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

        FPolynomial z(4);
        z.set(3, 0.0f);
        cout << "Zero polynomial: " << z << endl;

        FPolynomial p(6), q(4);

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
    }
    catch ( const PolynomialException& pex )
    {
        cerr << "PolynomialException caught: '";
        pex.display();
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

        /* Exact solution of the linear equation system a*x=b:
         * x  = [-0.6631640-0.3626125i, 0.1630189+1.050566i, -0.2240929+0.6002903i]'
         */
        MatrixGeneric<complex<float> > x = solver.solve();
        cout << "Solution:" << endl;
        x.display();
    }
    catch ( const LinearEquationSolverException& leqex )
    {
        cerr << "LinearEquationSolverException caught: '";
        leqex.display();
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
        cfex.display();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}

/*
 * Main function that starts three groups of unit tests
 */
int main(int argc, const char* argv[])
{

    cout << "Q U A T E R N I O N   T E S T" << endl << endl;
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

    return 0;
}
