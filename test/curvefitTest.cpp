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
 * A test module to test curve fitting classes
 * (PolynomialInterpolationGeneric, PolynomialRegressionGeneric)
 */


/*
 * Note: results are reproduced in 'scripts/test/curvefit.R'.
 */


#include <iostream>

#include "PolynomialRegression.h"
#include "PolynomialInterpolation.h"
#include "CurveFittingException.h"

using namespace std;
using namespace math;

/*
 * Test of curve fitting algorithms
 */
void curveFittingTest()
{
    try
    {
        // 1st, 2nd and 3rd degree regression polynomials, respectively,
        // of y=exp(-x) for x=0..5:
        PolynomialRegression fprexp1;
        PolynomialRegression fprexp2;
        PolynomialRegression fprexp3;

        // exp(-x), points in a random order:
        fprexp1.enterPoint(1.0, 0.367879441);
        fprexp1.enterPoint(0.0, 1.0);
        fprexp1.enterPoint(5.0, 0.006737947);
        fprexp1.enterPoint(4.0, 0.018315639);
        fprexp1.enterPoint(2.0, 0.135335283);
        fprexp1.enterPoint(3.0, 0.049787068);

        // copy points to other classes:
        fprexp2 = fprexp1;
        fprexp3 = fprexp1;

        // Checking of bounds
        cout << "Lower bound: " << fprexp2.lowerBound() << endl;
        cout << "Upper bound: " << fprexp2.upperBound() << endl;

        fprexp1.generateCurve(1);
        fprexp2.generateCurve(2);
        fprexp3.generateCurve(3);

        // Regression polynomials:
        cout << "1st degree polynomial: ";
        fprexp1.getPolynomial().display();
        cout << endl << "2nd degree polynomial: ";
        fprexp2.getPolynomial().display();
        cout << endl << "3rd degree polynomial: ";
        fprexp3.getPolynomial().display();
        cout << endl;

        /*
             Example from: http://en.wikipedia.org/wiki/Newton_polynomial#Example

             y=tan(x) for x=-1.5, -0.75, 0, 0.75, 1.5

             This time, 3rd degree regression and the same degree interpolation
             polynomial (as this is an odd function, all its even coefficients
             equal 0) are compared.
         */
        PolynomialInterpolation fpitan;
        PolynomialRegression fprtan3;

        fpitan.enterPoint(-1.5, -14.1014);
        fpitan.enterPoint(-0.75, -0.931596);
        fpitan.enterPoint(0.0, 0.0);
        fpitan.enterPoint(0.75, 0.931596);
        fpitan.enterPoint(1.5, 14.1014);

        fpitan.generateCurve();
        fprtan3.copy(&fpitan);
        fprtan3.generateCurve(3);

        cout << "Interpolation polynomial: ";
        fpitan.getPolynomial().display();
        cout << endl << "3rd degree regression polynomial: ";
        fprtan3.getPolynomial().display();
        cout << endl;

        /*
             Finally compare 5th degree interpolation and regression polynomials
             for the points from the first test:
         */
        PolynomialRegression fprexp5;
        PolynomialInterpolation fpiexp;
        fprexp5.copy(&fprexp1);
        fpiexp.copy(&fprexp5);
        fprexp5.generateCurve(5);
        fpiexp.generateCurve();
        cout << "Exp's regression polynomial: ";
        fprexp5.getPolynomial().display();
        cout << endl << "Exp's interpolation polynomial: ";
        fpiexp.getPolynomial().display();
        cout << endl;
    }
    catch ( const CurveFittingException& cfex )
    {
        cerr << "CurveFittingException caught: ";
        cfex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}
