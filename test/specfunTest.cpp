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
         (%o6)  -0.055452106757363
         (%i7)  float(gamma(2-%i));
         (%o7)  0.65296549642016-0.34306583981654*%i
         (%i8)  gamma(0.7+1.2*%i);
         (%o8)  0.31277447328267-0.23707614564003*%i
         (%i9)  gamma(-2.3-2.1*%i);
         (%o9)  0.0011346997482228*%i+0.0064712261389299
        */

        cout << "Gamma(3):    " << SpecFun::gamma(3.0) << " (expected: 2)" << endl;
        cout << "Gamma(4.2):  " << SpecFun::gamma(4.2) << " (expected: 7.756689535793178)" << endl;
        cout << "Gamma(0.51): " << SpecFun::gamma(0.51) << " (expected: 1.738415068463864)" << endl;
        cout << "Gamma(0.23): " << SpecFun::gamma(0.23) << " (expected: 3.95980372335778)" << endl;
        cout << "Gamma(-1.7): " << SpecFun::gamma(-1.7) << " (expected: 2.513923519065201)" << endl;
        cout << "Gamma(-4.56): " << SpecFun::gamma(-4.56) << " (expected: -0.055452106757363)" << endl;
        cout << "Gamma(2-i):  " << SpecFun::gamma(complex<double>(2.0, -1.0)) << " (expected: (0.65296549642016, -0.34306583981654))" << endl;
        cout << "Gamma(0.7+1.2i): " << SpecFun::gamma(complex<double>(0.7, 1.2)) << " (expected: (0.31277447328267, -0.23707614564003))" << endl;
        cout << "Gamma(-2.3-2.1i): " << SpecFun::gamma(complex<double>(-2.3, -2.1)) << " (expected: (0.0064712261389299, 0.0011346997482228))" << endl;


        /*
         (%i10)  float(beta(4,5));
         (%o10)  0.0035714285714285
         (%i11)  beta(3.2, 1-0.5*%i);
         (%o11)  0.20932777350065*%i+0.16569818762925
         (%i12)  beta(-2+3.7*%i, 0.7-1.3*%i);
         (%o12)  0.0056684836950215*%i+0.005533333385785
        */

        cout << endl;
        cout << "Beta(4, 5):   " << SpecFun::beta(4.0, 5.0) << " (expected: 0.0035714285714285)" << endl;
        cout << "Beta(3.2, 1-0.5i): " << SpecFun::beta<complex<double> >(3.2, complex<double>(1.0, -0.5)) << " (expected: (0.16569818762925, 0.20932777350065))" << endl;
        cout << "Beta(-2+3.7i, 0.7-1.3i)" << SpecFun::beta(complex<double>(-2.0, 3.7), complex<double>(0.7, -1.3)) << " (expected: (0.005533333385785, 0.0056684836950215))" << endl;


        /*
         (%i13)  gamma_incomplete(2, 0.5);
         (%o13)  0.90979598956895
         (%i14)  gamma_incomplete(2, 5.0);
         (%o14)  0.040427681994512
         (%i15)  gamma(3) - gamma_incomplete(3, 2.5);
         (%o15)  0.91237376823334
         (%i16)  gamma(1.5) - gamma_incomplete(1.5, 4);
         (%o16)  0.84545011298495
         (%i17)  gamma_incomplete_regularized(3, 0.5);
         (%o17)  0.98561232203302
         (%i18)  1 - gamma_incomplete_regularized(1.5, 3.2);
         (%o18)  0.90630920959237
         (%i19)  gamma_incomplete(2.0+%i, 1.0-%i);
         (%o19)  1.123545237341651*%i+0.29474062291749
         (%i20)  gamma(1.4-0.7*%i) - gamma_incomplete(1.4-0.7*%i, 0.3+0.1*%i);
         (%o20)  0.13052172604872*%i-0.016291214690452
         (%i21)  rectform(gamma_incomplete(-1.1+0.8*%i, 2.4-1.3*%i) / gamma(-1.1+0.8*%i));
         (%o21)  0.016252537791964*%i+0.0019715754028325
         (%i22)  1 - rectform(gamma_incomplete(-3.2-1.4*%i, -2-0.5*%i) / gamma(-3.2-1.4*%i));
         (%o22)  0.43578355528153*%i+0.23971636701011
        */

        cout << endl;
        cout << "Upper inc. gamma(2, 0.5): " << SpecFun::incGammaUpper(2.0, 0.5) << " (expected: 0.90979598956895)" << endl;
        cout << "Upper inc. gamma(2, 5):   " << SpecFun::incGammaUpper(2.0, 5.0) << " (expected: 0.040427681994512)" << endl;
        cout << "Lower inc. gamma(3, 2.5): " << SpecFun::incGammaLower(3.0, 2.5) << " (expected: 0.91237376823334)" << endl;
        cout << "Lower inc. gamma(1.5, 4): " << SpecFun::incGammaLower(1.5, 4.0) << " (expected: 0.84545011298495)" << endl;
        cout << "Reg. upper inc. gamma(3, 0.5):   " << SpecFun::incGammaUpperReg(3.0, 0.5) << " (expected: 0.98561232203302)" << endl;
        cout << "Reg. lower inc. gamma(1.5, 3.2): " << SpecFun::incGammaLowerReg(1.5, 3.2) << " (expected: 0.90630920959237)" << endl;
        cout << "Upper inc. gamma(2+i, 1-i): " << SpecFun::incGammaUpper(complex<double>(2.0, 1.0), complex<double>(1.0, -1.0)) << " (expected: (0.29474062291749, 1.123545237341651))" << endl;
        cout << "Lower inc. gamma(1.4-0.7i, 0.3+0.1i): " << SpecFun::incGammaLower(complex<double>(1.4, -0.7), complex<double>(0.3, 0.1)) << " (expected: (-0.016291214690452, 0.13052172604872))" << endl;
        cout << "Reg. upper inc. gamma(-1.1+0.8i, 2.4-1.3i): " << SpecFun::incGammaUpperReg(complex<double>(-1.1, 0.8), complex<double>(2.4, -1.3)) << " (expected: (0.0019715754028325, 0.016252537791964))" << endl;
        cout << "Reg. lower inc. gamma(-3.2-1.4i, -2-0.5i): " << SpecFun::incGammaLowerReg(complex<double>(-3.2, -1.4), complex<double>(-2.0, -0.5)) << " (expected: (0.23971636701011, 0.43578355528153))" << endl;


        /*
         (%i23)  beta_incomplete(2, 5, 0.2);
         (%o23)  0.011488
         (%i24)  beta_incomplete(2, 5, 0.7);
         (%o24)  0.032968833333333
         (%i25)  beta(1, 2) - beta_incomplete(1, 2, 0.15);
         (%o25)  0.36125
         (%i26)  beta(1, 2) - beta_incomplete(1, 2, 0.82);
         (%o26)  0.0162
         (%i27)  beta_incomplete_regularized(3, 2, 0.12);
         (%o27)  0.0062899199999999
         (%i28)  1 - beta_incomplete_regularized(4, 1, 0.32);
         (%o28)  0.98951424
        */
        
        cout << endl;
        cout << "Lower inc. beta(0.2, 2, 5):  " << SpecFun::incBetaLower(0.2, 2.0, 5.0) << " (expected: 0.011488)" << endl;
        cout << "Lower inc. beta(0.7, 2, 5):  " << SpecFun::incBetaLower(0.7, 2.0, 5.0) << " (expected: 0.032968833333333)" << endl;
        cout << "Upper inc. beta(0.15, 1, 2): " << SpecFun::incBetaUpper(0.15, 1.0, 2.0) << " (expected: 0.36125)" << endl;
        cout << "Upper inc. beta(0.82, 1, 2): " << SpecFun::incBetaUpper(0.82, 1.0, 2.0) << " (expected: 0.0162)" << endl;
        cout << "Reg. lower inc. beta(0.12, 3, 2): " << SpecFun::incBetaLowerReg(0.12, 3.0, 2.0) << " (expected: 0.0062899199999999)" << endl;
        cout << "Reg. upper inc. beta(0.32, 4, 1): " << SpecFun::incBetaUpperReg(0.32, 4.0, 1.0) << " (expected: 0.98951424)" << endl;


        /*
         (%i29)  erf(-1.2);
         (%o29)  -0.91031397822963
         (%i30)  erf(0.7);
         (%o30)  0.67780119383741
         (%i31)  erfc(0.2);
         (%o31)  0.77729741078952
         */
        
        cout << endl;
        cout << "erf(-1.2):  " << SpecFun::erf(-1.2) << " (expected: -0.91031397822963)" << endl;
        cout << "erf(0.7):   " << SpecFun::erf(0.7) << " (expected: 0.67780119383741)" <<endl;
        cout << "erfc(0.2):  " << SpecFun::erfc(0.2) << " (expected: 0.77729741078952)" << endl;
    }
    catch ( const SpecFunException& ex )
    {
        cerr << "Statistics exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
