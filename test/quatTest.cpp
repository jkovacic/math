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


#include <iostream>

#include "QuaternionGeneric.h"
#include "QuaternionException.h"


using namespace std;
using namespace math;


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
