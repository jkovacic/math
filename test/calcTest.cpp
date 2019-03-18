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


/*
 * Note: results are reproduced in 'scripts/test/calc.mac'.
 */


#include <iostream>
#include <cmath>

#include "IFunction.h"
#include "Integ.h"
#include "Diff.h"
#include "FunctionException.h"
#include "CalculusException.h"


using namespace std;
using namespace math;



class CFunc : public IFunction
{
public:
    double ke, kl, n;

    double operator()(const double& x) const
    {
        return ke * std::exp(-(x-3)*(x-3)) + kl * x + n;
    }
};


class CNormPdf : public IFunction
{

    double operator()(const double& x) const
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

        cout << "f(x) = " << f.ke << "*exp(-(x-3)^2) " << showpos << f.kl << "*x " << f.n << noshowpos << endl << endl;

        cout << "Numerical integration:" << endl;
        for ( int method=Integ::EIntegAlg::RECTANGLE; method<=Integ::EIntegAlg::BOOLE; ++method )
        {
            cout << "Method " << method << ": Int(f(x), 0, 5) = " <<
                     Integ::integ(f, 0., 5., 10000, static_cast<Integ::EIntegAlg::alg>(method)) << endl;
        }

        // A different 'n' is required for the Romberg's method:
        cout << "Romberg's method: Int(f(x), 0, 5) = " <<
                     Integ::integ(f, 0., 5., 6, Integ::EIntegAlg::ROMBERG) << endl;

        cout << "Expected result: 16.23780211731536" << endl;


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


        cout << endl << "Numerical differentiation:" << endl;
        for ( int method=Diff::EDiffMethod::FORWARD; method<=Diff::EDiffMethod::FIVE_POINT; ++method )
        {
            cout << "Method " << method << ": f'(4) = " <<
                    Diff::diff(f, 4.0, 0.001, static_cast<Diff::EDiffMethod::method>(method)) << endl;
        }
        cout << "Expected result: 0.48496878235998" << endl;


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
