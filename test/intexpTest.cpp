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
 * A test module to test int. exponentiation functionality
 * (NumericUtil::power)
 */

#include <iostream>

#include "NumericUtil.h"
#include "SqMatrixGeneric.h"
#include "QuaternionGeneric.h"
#include "Rational.h"
#include "PolynomialGeneric.h"
#include "MatrixException.h"
#include "QuaternionException.h"
#include "RationalException.h"
#include "PolynomialException.h"

using namespace std;
using namespace math;

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
