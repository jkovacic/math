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
 * A test module to test calculus functions
 * (Integ, Diff)
 */


#include <iostream>
#include <cmath>

#include "IFunctionGeneric.h"
#include "IntegGeneric.h"
#include "DiffGeneric.h"
#include "FunctionException.h"
#include "CalculusException.h"


using namespace std;
using namespace math;



class CFunc : public IFunction
{
public:
    double ke, kl, n;

    double operator()(const double& x) const throw(FunctionException)
    {
        return ke * std::exp(-(x-3)*(x-3)) + kl * x + n;
    }
};


class CNormPdf : public IFunction
{

    double operator()(const double& x) const throw(FunctionException)
    {
        return 0.39894228040143268 * std::exp( -0.5 * x * x );
    }
};



/*
 * Test of numerical integration and derivation algorithms
 */
void calculusTest()
{
    try
    {
        CFunc f;
        f.ke = 0.7;   f.kl = 1.0;  f.n = 0.5;

        /*
         * Exact result obtained by Maxima:
         *
           (%i1) f(x):=0.7*exp(-(x-3)^2)+1*x+0.5$
           (%i2) float(integrate(f(x), x, 0, 5));
           (%o2) 16.23780211731536
         */
        cout << "f(x) = " << f.ke << "*exp(-(x-3)^2) " << showpos << f.kl << "*x " << f.n << noshowpos << endl << endl;

        cout << "Numerical integration:" << endl;
        for ( int method=EIntegAlg::RECTANGLE; method<=EIntegAlg::BOOLE; ++method )
        {
            cout << "Method " << method << ": Int(f(x), 0, 5) = " <<
                     Integ::integ(f, 0., 5., 10000, static_cast<EIntegAlg::alg>(method)) << endl;
        }

        // A different 'n' is required for the Romberg's method:
        cout << "Romberg's method: Int(f(x), 0, 5) = " <<
                     Integ::integ(f, 0., 5., 6, EIntegAlg::ROMBERG) << endl;

        cout << "Expected result: 16.23780211731536" << endl;


        /*
         * Exact results obtained by Maxima:
         *
           (%i3)  fn(x) := 0.39894228040143268 * exp(-x^2/2)$
           (%i4)  float(integrate(fn(x), x, 1, inf));
           (%o4)  0.15865525393145
           (%i5)  float(integrate(fn(x), x, 1.3, inf));
           (%o5)  0.09680048458561
           (%i6)  float(integrate(fn(x), x, -2, inf));
           (%o6)  0.97724986805182
           (%i7)  float(integrate(fn(x), x, -1.8, inf));
           (%o7)  0.96406968088707
           (%i8)  float(integrate(fn(x), x, -inf, -2));
           (%o8)  0.022750131948179
           (%i9)  float(integrate(fn(x), x, -inf, -1.7));
           (%o9)  0.044565462758543
           (%i10) float(integrate(fn(x), x, -inf, 1));
           (%o10) 0.84134474606854
           (%i11) float(integrate(fn(x), x, -inf, 0.5));
           (%o12) 0.69146246127401
           (%i13) float(integrate(fn(x), x, -inf, inf));
           (%o13) 1.0
           (%i14) float(integrate(fn(x), x, -inf, inf));
           (%o14) 1.0
         */
        const CNormPdf fn;

        cout << endl;
        cout << "fn(x) = 0.39894228040143268 * exp(-x^2/2) :" << endl;
        cout << "Int(fn(x), 1, inf)     = " << Integ::integImpPosInf(fn, 1.0) << " (expected: 0.15865525393145)" << endl;
        cout << "Int(fn(x), 1.3, inf)   = " << Integ::integImpPosInfH(fn, 1.3) << " (expected: 0.09680048458561)" << endl;
        cout << "Int(fn(x), -2, inf)    = " << Integ::integImpPosInf(fn, -2.0) << " (expected: 0.97724986805182)" << endl;
        cout << "Int(fn(x), -1.8, inf)  = " << Integ::integImpPosInfH(fn, -1.8) << " (expected: 0.96406968088707)" << endl;
        cout << "Int(fn(x), -inf, -2)   = " << Integ::integImpNegInf(fn, -2.0) << " (expected: 0.022750131948179)" << endl;
        cout << "Int(fn(x), -inf, -1.7) = " << Integ::integImpNegInfH(fn, -1.7) << " (expected: 0.044565462758543)" << endl;
        cout << "Int(fn(x), -inf, 1)    = " << Integ::integImpNegInf(fn, 1.0) << " (expected: 0.84134474606854)" << endl;
        cout << "Int(fn(x), -inf, 0.5)  = " << Integ::integImpNegInfH(fn, 0.5) << " (expected: 0.69146246127401)" << endl;
        cout << "Int(fn(x), -inf, inf)  = " << Integ::integImp(fn) << " (expected: 1.0)" << endl;
        cout << "Int(fn(x), -inf, inf)  = " << Integ::integImpH(fn) << " (expected: 1.0)" << endl;

        /*
         * Exact slope obtained by Maxima:
         *
           (%i15) df(x) := ''(diff(f(x), x))$
           (%i16) float(df(4));
           (%o16) 0.48496878235998
         */
        cout << endl << "Numerical differentiation:" << endl;
        for ( int method=EDiffMethod::FORWARD; method<=EDiffMethod::FIVE_POINT; ++method )
        {
            cout << "Method " << method << ": f'(4) = " <<
                    Diff::diff(f, 4.0, 0.001, static_cast<EDiffMethod::method>(method)) << endl;
        }
        cout << "Expected result: 0.48496878235998" << endl;

        /*
         * Exact 2nd order derivative obtained by Maxima:
         *
           (%i17)  d2f(x) := ''(diff(f(x), x, 2))$
           (%i18)  float(d2f(2));
           (%o18)  0.51503121764002
         */
        cout << endl << "f''(2) = " << Diff::diff2(f, 2.0, 0.001) << " (expected: 0.51503121764002)" << endl;
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
