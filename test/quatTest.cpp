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
 * A test module to test quaternion classes
 * (QuaternionGeneric)
 */


/*
 * Note: results are reproduced in 'scripts/test/quat.m'.
 */


#include <iostream>

#include "QuaternionGeneric.h"
#include "QuaternionException.h"


using namespace std;
using namespace math;


/*
 * Test of quaternions
 */
void quaternionTest()
{
    try
    {
        Quaternion zeroq;
        Quaternion q(1.0, -0.5, 1.2, -2.3);

        cout << "zero: " << zeroq << endl;
        cout << "q: " << q << endl;
        cout << "q: 1: " << q.getOne() << "   i: " << q.getI() << "   j: " << q.getJ() << "   k: " << q.getK() << endl;
        cout << "abs(q): " << abs(q) << " (expected: 2.824889)" << endl;

        Quaternion uq = q.unit();
        cout << "U(q): " << uq << "\tnorm(U(q)): " << uq.norm() << endl;
        uq = 2.3;
        cout << "uq = " << uq << endl;

        Quaternion o(1.0);
        Quaternion i;   i.setI(1.0);
        Quaternion j;   j.setJ(1.0);
        Quaternion k;   k.setK(1.0);

        cout << "i*j: " << i*j << "\tj*i: " << j*i << endl;
        cout << "j*k: " << j*k << "\tk*j: " << k*j << endl;
        cout << "k*i: " << k*i << "\ti*k: " << i*k << endl;

        Quaternion p(1.0, 2.0, 3.0, 4.0);
        Quaternion sum = p+q;
        Quaternion dif = p-q;
        cout << "p: " << +p << "\tq: " << q << endl;
        cout << "p+q: " << sum << "\tp-q: " << dif << endl;
        // Must output (7.6-10.2i+6.8j+5.6k) and (7.6+13.2i+1.6j-2.2k), respectively
        cout << "p*q: " << p*q << "\tq*p: " << q*p << endl;
        cout << "p+5: " << p+5. << "\t3+p: " << 3.0 + p << endl;
        cout << "p-4: " << p-4. << "\t1-p: " << 1.0 - p << endl;
        p += 2;
        cout << "p+2: " << p;
        p -= 2.;
        cout << "\tp+2-2: " << p << endl;

        Quaternion qc = -0.5*(q+i*q*i+j*q*j+k*q*k);
        cout << "-0.5*(q+i*q*i+j*q*j+k*q*k)=" << qc << "  q conj: " << q.conj() << endl;
        Quaternion qrec = q.reciprocal();
        cout << "qrec: " << qrec << "\tq*qrec=" << (q*qrec).roundSmallElements() 
                << "\tqrec*q=" << (qrec*q).roundSmallElements(1e-15) << endl;
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
