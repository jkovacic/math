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
 * A test module to test root finding classes
 * (RootFindGeneric)
 */


#include <iostream>
#include <cmath>

#include "NumericUtil.h"
#include "IFunctionGeneric.h"
#include "RootFindGeneric.h"
#include "FunctionException.h"
#include "RootFindException.h"


using namespace std;
using namespace math;


class FFunc : public IFunction
{
public:
    // function: sin(x)/x - 0.5
    double func(const double& x) const throw(FunctionException)
    {
   	    if ( true == NumericUtil::isZero<double>(x) )
        {
            /*
             * When x approaches 0, the division by zero can occur.
             * In this case apply the well known limit:
             *
             *        sin(x)
             *  lim  -------- = 1
             *  x->0    x
             *
             *  Recheck the limit in Maxima:
               (%i1)  f(x) := sin(x)/x$
               (%i2)  limit(sin(x)/x, x, 0);
               (%o2)  1
             */
            return 1.0 - 0.5;
        }
        else
        {
           return sin(x)/x - 0.5;
      	}
    }
};


class DFunc : public IFunction
{
public:
    double func(const double& x) const throw(FunctionException)
    {
        /*
         * Derivation of the function f, defined above.
         * Not that the subtrahend 0.5 is cancelled out by differentiation:
         *
         *  df(x)     cos(x)     sin(x)
         * ------- = -------- - --------
         *   dx         x         x**2
         *
         * Verified in Maxima:
           (%i3)  d(x) := ''(diff(f(x), x));
           (%o3)  d(x):=cos(x)/x−sin(x)/x^2
         */

        if ( true == NumericUtil::isZero<double>(x) )
        {
            /*
             * Prevent division by 0 by applying the limit:
             *
             *       /  cos(x)     sin(x)  \
             *  lim  | -------- - -------- | = 0
             *  x->0 \    x         x**2   /
             *
             * Verified in Maxima:
             (%i4)  limit(d(x), x, 0);
             (%o4)  0
             */

           	return 0.0;
        }
        else
        {
            return  cos(x)/x - sin(x)/(x*x);
        }
    }
};


/*
 * Test of root finding algorithms
 */
void rootFindTest()
{
    try
    {
        const FFunc f;
        const DFunc d;
        double x0;

        /*
         * "Exact" numerical solution found by Maxima:
         (%i5)  find_root(f(x), x, 1, 3);
         (%o5)  1.895494267033981
         */
        x0 = RootFind::bisection(f, 1.0, 3.0, 1e-9, 1e-9);
        cout << "Bisection method:\t\tx0 = " << x0 << "\t";
        cout << "f(x0) = " << f.func(x0) << endl;

        x0 = RootFind::regulaFalsi(f, 1.0, 3.0, 1e-9, 1e-9);
        cout << "Regula falsi method:\tx0 = " << x0 << "\t";
        cout << "f(x0) = " << f.func(x0) << endl;

        x0 = RootFind::secant(f, 1.0, 3.0, 1e-9);
        cout << "Secant method:\t\t\tx0 = " << x0 << "\t";
        cout << "f(x0) = " << f.func(x0) << endl;

        x0 = RootFind::newton(f, d, 1.0, 1e-9);
        cout << "Newton's method:\t\tx0 = " << x0 << "\t";
        cout << "f(x0) = " << f.func(x0) << endl;

        x0 = RootFind::quasiNewton(f, 1.0, 1e-9, 0.001);
        cout << "Quasi Newton's method:\tx0 = " << x0 << "\t";
        cout << "f(x0) = " << f.func(x0) << endl;

        cout << "Correct root: 1.895494267033981" << endl;
    }
    catch ( const RootFindException& rfex )
    {
        cerr << "Root find exception caught: ";
        rfex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
