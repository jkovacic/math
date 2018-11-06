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
 * A test module to test root finding functions
 * (namespace RootFind)
 */


/*
 * Note: results are reproduced in 'scripts/test/rootfind.mac'.
 */


#include <iostream>
#include <cmath>

#include "NumericUtil.h"
#include "IFunction.h"
#include "RootFind.h"
#include "FunctionException.h"
#include "RootFindException.h"


using namespace std;
using namespace math;


class FFunc : public IFunction
{
public:
    // function: sin(x)/x - 0.5
    double operator()(const double& x) const
    {
        const double MND = 0.5;

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
             *  The limit can be rechecked in Maxima.
             */
            return 1.0 - MND;
        }
        else
        {
           return sin(x)/x - MND;
        }
    }
};


class DFunc : public IFunction
{
public:
    double operator()(const double& x) const
    {
        /*
         * Derivation of the function f, defined above.
         * Note that the subtrahend 0.5 is canceled out by differentiation:
         *
         *  df(x)     cos(x)     sin(x)
         * ------- = -------- - --------
         *   dx         x         x**2
         *
         * The derivation can be rechecked in Maxima.
         */

        if ( true == NumericUtil::isZero<double>(x) )
        {
            /*
             * Prevent division by 0 by applying the limit:
             *
             *       /  cos(x)     sin(x)  \
             *  lim  | -------- - -------- | = 0
             *  x->0 \    x          x^2   /
             *
             * The limit can be rechecked in Maxima.
             */

            return 0.0;
        }
        else
        {
            return  cos(x)/x - sin(x)/(x*x);
        }
    }
};


class D2Func : public IFunction
{
public:
    double operator()(const double& x) const
    {
        /*
         * 2nd order derivation of the function f, defined above.
         * Note that the subtrahend 0.5 is canceled out by differentiation:
         *
         *   2
         *  d f(x)       sin(x)     2 * cos(x)     2 * sin(x)
         * -------- = - -------- - ------------ + ------------
         *      2          x           x^2            x^3
         *    dx
         *
         * The derivation can be rechecked in Maxima.
         */

        if ( true == NumericUtil::isZero<double>(x) )
        {
            /*
             * Prevent division by 0 by applying the limit:
             *
             *       /    sin(x)     2 * cos(x)     2 * sin(x)  \      1
             *  lim  | - -------- - ------------ + ------------ | = - ---
             *  x->0 \      x           x^2            x^3      /      3
             *
             * The limit can be rechecked in Maxima.
             */

            return -1.0 / 3.0;
        }
        else
        {
            return -sin(x)/x * (1.0 - 2.0/(x*x)) - 2.0 * cos(x)/(x*x);
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
        const D2Func d2;
        double x0;

        const double EPS = 1e-11;

        x0 = RootFind::bisection(f, 1.0, 3.0, 1e-9, 1e-9);
        cout << "Bisection method:            x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::regulaFalsi(f, 1.0, 3.0, 1e-9, 1e-9);
        cout << "Regula falsi method:         x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::secant(f, 1.0, 3.0, 1e-9);
        cout << "Secant method:               x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::newton(f, d, 1.0, 1e-9);
        cout << "Newton's method:             x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::quasiNewton(f, 1.0, 1e-9, 0.001);
        cout << "Quasi Newton's method:       x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::halley(f, d, d2, 1.0, 1e-9);
        cout << "Halleys's method:            x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::quasiHalley(f, 1.0, 1e-9, 0.001);
        cout << "Quasi Halleys's method:      x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::halleyMod(f, d, d2, 1.0, 1e-9);
        cout << "Mod. Halleys's method:       x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

        x0 = RootFind::quasiHalleyMod(f, 1.0, 1e-9, 0.001);
        cout << "Mod. quasi Halleys's method: x0 = " << x0 << "\t";
        cout << "f(x0) ~= " << NumericUtil::smallValToZero(f(x0), EPS) << endl;

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
