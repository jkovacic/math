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
 * A test module to test polynomial classes
 * (PolynomialGeneric)
 */


/*
 * Note: results are reproduced in 'scripts/test/poly.mac'.
 */


#include <iostream>
#include <vector>

#include "Polynomial.h"
#include "PolynomialException.h"

using namespace std;
using namespace math;


/*
 * Test of polynomial algebra
 */
void polynomialTest()
{
    try
    {
        FPolynomial pol1(-4.1f);
        cout << "pol1 = " << pol1 << endl;

        float a[] = { 2.1f, 1.0f, -0.72f, 1.0f, 0.0f };
        FPolynomial t(a, 5);
        cout << "t(x) = " << t << endl;
        cout << "a0=" << t.get(0) << " a1=" << t.get(1) << " a2=" << t.get(2) << " a3=" << t.get(3) << endl;
        cout << "Remove the 2nd coef. from t(x): " << t.remove_(1) << endl;
        cout << "Insert 0.2 to the 3rd pos. of t(x): " << t.insert_(2, 0.2f) << endl;
        cout << "Insert 5 to the 8th pos. of t(x): " << t.insert_(7, 5.0f) << endl;

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

        cout << "Round all t's coefficients below 0.75 to 0: " << t.roundSmallCoefficients(0.75f) << endl;

        t = 78.12f;
        cout << "t = " << t << endl;

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
        ptemp = +p;
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
