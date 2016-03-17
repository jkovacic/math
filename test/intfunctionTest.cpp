/*
Copyright 2016, Jernej Kovacic

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
 * A test module to test functions in the namespace "IntFunction"
 * that "specialize" some elementary functions for integer types.
 */


/*
 * Note: results are reproduced in 'scripts/test/intfunction.jl'.
 */


#include <iostream>

#include "IntFunction.h"
#include "IntFunctionException.h"


using namespace std;
using namespace math;


/*
 * Test of integer functions
 */
void intFunctionTest()
{

    try
    {
        // Find integer square roots of the following numbers:
        const unsigned int NS = 3;
        const int s[NS] = { 12, 100, 37423 };

        // Expected results: 3, 10 and 193
        for ( unsigned int i=0; i<NS; ++i )
        {
            cout << "floor(sqrt(" << s[i] << ")) = " << IntFunction::intSqrt(s[i]) << endl;
        }
        cout << endl;

        // Expected results: 4, 10 and 194
        for ( unsigned int i=0; i<NS; ++i )
        {
            cout << "ceil(sqrt(" << s[i] << ")) = " << IntFunction::intSqrt(s[i], true) << endl;
        }
        cout << endl;

        // Log2 of signed integers:
        const unsigned int NL2S = 7;
        const int l2s[NL2S] = { 1, 8, 63, 64, 1024, 1025, 65000 };

        // Expected results: 0, 3, 5, 6, 10, 10, 15
        for ( unsigned int i=0; i<NL2S; ++i )
        {
            cout << "floor(log2((" << l2s[i] << ")) = " << IntFunction::intLog2(l2s[i]) << endl;
        }
        cout << endl;

        // Expected results: 0, 3, 6, 6, 10, 11, 16
        for ( unsigned int i=0; i<NL2S; ++i )
        {
            cout << "ceil(log2((" << l2s[i] << ")) = " << IntFunction::intLog2(l2s[i], true) << endl;
        }
        cout << endl;

        // Log2 of unsigned integers:
        const unsigned int NL2U = 3;
        const unsigned int l2u[NL2U] = { 1, 16, 65 };

        // Expected results: 0, 4, 6
        for ( unsigned int i=0; i<NL2U; ++i )
        {
            cout << "floor(log2(" << l2u[i] << ")) = " << IntFunction::intLog2(l2u[i]) << endl;
        }
        cout << endl;

        // Expected results: 0, 4, 7
        for ( unsigned int i=0; i<NL2U; ++i )
        {
            cout << "ceil(log2(" << l2u[i] << ")) = " << IntFunction::intLog2(l2u[i], true) << endl;
        }

    }
    catch ( const IntFunctionException& ex )
    {
        cerr << "IntFunction exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
