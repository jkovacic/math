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
 * A test module to test rational classes
 * (Rational)
 */

#include <iostream>

#include "RationalGeneric.h"
#include "RationalException.h"


using namespace std;
using namespace math;

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
        c = +a + b;
        // should output sth. like "1/3 + 1/2 = 5/6 = 0.833333"
        cout << a << " + " << b << " = " << c << " = " << c.toNum<float>() << endl;

        b += a;
        a.set(1,4);
        c = b - a;
        // should output sth. like "5/6 - 1/4 = 7/12 = 0.5833333"
        cout << b << " - " << a << " = " << c << " = " << c.toNum<double>() << endl;
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

        // Test of overloaded arithmetic operators:
        cout << a << " + 2 = " << a + 2L << endl;
        cout << "5 + (" << a << ") = " << 5L + a << endl;
        cout << a << " - 3 = " << a - 3L << endl;
        cout << "1 - (" << a << ") = " << 1L - a << endl;
        cout << a << " * (-4) = " << a * (-4L) << endl;
        cout << "2 * (" << a << ") = " << 2L * a << endl;
        cout << a << " / 8 = " << a / 8L << endl;
        cout << "10 / (" << a << ") = " << 10L / a << endl;

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

        // Test of overloaded comparison operators
        cout << a << " == 3 : " << (a==3L ? "TRUE" : "FALSE") << endl;
        cout << "2 == " << b << " : " << (2L==b ? "TRUE" : "FALSE") << endl;
        cout << a << " != 1 : " << (a!=1L ? "TRUE" : "FALSE") << endl;
        cout << "-2 != " << b << " : " << (-2L!=b ? "TRUE" : "FALSE") << endl;
        cout << a << " > 1 : " << (a>1L ? "TRUE" : "FALSE") << endl;
        cout << "3 > " << b << " : " << (3L>b ? "TRUE" : "FALSE") << endl;
        cout << a << " >= 2 : " << (a>=2L ? "TRUE" : "FALSE") << endl;
        cout << "5 >= " << b << " : " << (5L>=b ? "TRUE" : "FALSE") << endl;
        cout << a << " < 6 : " << (a<6L ? "TRUE" : "FALSE") << endl;
        cout << "-2 < " << b << " : " << (-2L<b ? "TRUE" : "FALSE") << endl;
        cout << a << " <= 4 : " << (a<=4L ? "TRUE" : "FALSE") << endl;
        cout << "3 <= " << b << " : " << (3L<=b ? "TRUE" : "FALSE") << endl;

        // Test of overloaded compound assignment operators:
        cout << a << " + 3 = ";   a += 3;  cout << a << endl;
        cout << a << " - 1 = ";   a -= 1;  cout << a << endl;
        cout << a << " * 5 = ";   a *= 5;  cout << a << endl;
        cout << a << " / 2 = ";   a /= 2;  cout << a << endl;
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
            cout << "Attempt of division by zero successfully detected." << endl;
        }
        else
        {
            cerr << "Unexpected exception: ";
            zerEx.what();
            cerr << endl;
        }
    }
}
