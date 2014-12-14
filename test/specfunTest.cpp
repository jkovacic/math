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
 * A test module to test special functions (gamma, beta, etc.),
 * implemented in the namespace SpecFun.
 */

#include <iostream>
#include <complex>

#include "SpecFunGeneric.h"
#include "SpecFunException.h"

using namespace std;
using namespace math;


void specfunTest()
{
    try
    {
        /*
         * Unit test results obtained by Maxima:
         (%i1)  float(gamma(3));
         (%o1)  2.0
         (%i2)  gamma(4.2);
         (%o2)  7.756689535793178
         (%i3)  gamma(0.51);
         (%o3)  1.738415068463864
         (%i4)  gamma(0.23);
         (%o4)  3.95980372335778
         (%i5)  gamma(-1.7);
         (%o5)  2.513923519065201
         (%i6)  gamma(-4.56);
         (%o6)  −0.055452106757363
         (%i7)  float(gamma(2-%i));
         (%o7)  0.65296549642016−0.34306583981654*%i
         (%i8)  gamma(0.7+1.2*%i);
         (%o8)  0.31277447328267−0.23707614564003*%i
         (%i9)  gamma(-2.3-2.1*%i);
         (%o9)  0.0011346997482228*%i+0.0064712261389299
         (%i10)  float(beta(4,5));
         (%o10)  0.0035714285714285
         (%i11)  beta(3.2, 1-0.5*%i);
         (%o11)  0.20932777350065*%i+0.16569818762925
         (%i12)  beta(-2+3.7*%i, 0.7-1.3*%i);
         (%o12)  0.0056684836950215*%i+0.005533333385785
         (%i13)  erf(-1.2);
         (%o13)  −0.91031397822963
         (%i14)  erf(0.7);
         (%o14)  0.67780119383741
         (%i15)  erfc(0.2);
         (%o15)  0.77729741078952
         */

        cout << "Gamma(3):    " << SpecFun::gamma(3.0) << " (expected: 2)" << endl;
        cout << "Gamma(4.2):  " << SpecFun::gamma(4.2) << " (expected: 7.756689535793178)" << endl;
        cout << "Gamma(0.51): " << SpecFun::gamma(0.51) << " (expected: 1.738415068463864)" << endl;
        cout << "Gamma(0.23): " << SpecFun::gamma(0.23) << " (expected: 3.95980372335778)" << endl;
        cout << "Gamma(-1.7): " << SpecFun::gamma(-1.7) << " (expected: 2.513923519065201)" << endl;
        cout << "Gamma(-4.56): " << SpecFun::gamma(-4.56) << " (expected: −0.055452106757363)" << endl;
        cout << "Gamma(2-i):  " << SpecFun::gamma(complex<double>(2.0, -1.0)) << " (expected: (0.65296549642016, −0.34306583981654))" << endl;
        cout << "Gamma(0.7+1.2i): " << SpecFun::gamma(complex<double>(0.7, 1.2)) << " (expected: (0.31277447328267, −0.23707614564003))" << endl;
        cout << "Gamma(-2.3-2.1i): " << SpecFun::gamma(complex<double>(-2.3, -2.1)) << " (expected: (0.0064712261389299, 0.0011346997482228))" << endl;

        cout << endl;
        cout << "Beta(4, 5):   " << SpecFun::beta(4.0, 5.0) << " (expected: 0.0035714285714285)" << endl;
        cout << "Beta(3.2, 1-0.5i): " << SpecFun::beta<complex<double> >(3.2, complex<double>(1.0, -0.5)) << " (expected: (0.16569818762925, 0.20932777350065))" << endl;
        cout << "Beta(-2+3.7i, 0.7-1.3i)" << SpecFun::beta(complex<double>(-2.0, 3.7), complex<double>(0.7, -1.3)) << " (expected: (0.005533333385785, 0.0056684836950215))" << endl;

        cout << endl;
        cout << "erf(-1.2):  " << SpecFun::erf(-1.2) << " (expected: −0.91031397822963)" << endl;
        cout << "erf(0.7):   " << SpecFun::erf(0.7) << " (expected: 0.67780119383741)" <<endl;
        cout << "erfc(0.2):  " << SpecFun::erfc(0.2) << " (expected: 0.77729741078952)" << endl;
    }
    catch ( const SpecFunException& ex )
    {
        cerr << "Statistics exception caught: ";
        ex.what();
        cerr << endl;
    }
}
