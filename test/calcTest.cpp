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

    double func(const double& x) const throw(FunctionException)
    {
        return ke * std::exp(-(x-3)*(x-3)) + kl * x + n;
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
        cout << "Expected result: 16.23780211731536" << endl;

        /*
         * Exact slope obtained by Maxima:
           (%i3) df(x) := ''(diff(f(x), x))$
           (%i4) float(df(4));
           (%o5) 0.48496878235998
         */
        cout << endl << "Numerical differentiation:" << endl;
        for ( int method=EDiffMethod::FORWARD; method<=EDiffMethod::FIVE_POINT; ++method )
        {
            cout << "Method " << method << ": f'(4) = " <<
                    Diff::diff(f, 4.0, 0.001, static_cast<EDiffMethod::method>(method)) << endl;
        }
        cout << "Expected result: 0.48496878235998" << endl;
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
