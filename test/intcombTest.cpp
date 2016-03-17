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
 * A test module to test integer combinatorics functions
 * in the namespace IntCombinatorics
 */


/*
 * Note: results are reproduced in 'scripts/test/intcomb.mac'.
 */


#include <iostream>

#include "IntCombinatorics.h"
#include "CombinatoricsException.h"


using namespace std;
using namespace math;

/*
 * Test of integer combinatorics
 */
void intCombinatoricsTest()
{
    try
    {
        cout << "5! = " << IntCombinatorics::factorial(5U) << "  expected: 120" << endl;
        cout << "20! = " << IntCombinatorics::factorial(20ULL) << "  expected: 2432902008176640000 " << endl;
        cout << "15!/5! = " << IntCombinatorics::factorial(15ULL, 6ULL) << "  expected: 10897286400" << endl;
        cout << endl;

        cout << "falling factorial: (10)_4 = " << IntCombinatorics::fallingFactorial(10U, 4U) << "  expected:  5040" << endl;
        cout << "falling factorial: (12)_8 = " << IntCombinatorics::fallingFactorial(12UL, 8UL) << "  expected:  19958400" << endl;
        cout << "rising factorial: 5^(6) = " << IntCombinatorics::risingFactorial(5UL, 6UL) << "  expected :  151200" << endl;
        cout << "rising factorial: 16^(10) = " << IntCombinatorics::risingFactorial(16ULL, 10ULL) << "  expected :  11861676288000" << endl;
        cout << endl;

        // Checked with this Maxima function:
        // multif(n,k) := if (n<k) then 1 else n*multif(n-k,k)$
        cout << "18!^(3) = " << IntCombinatorics::multiFactorial(18UL, 3UL) << "  expected:  524880" << endl;
        cout << "19!^(3) = " << IntCombinatorics::multiFactorial(19UL, 3UL) << "  expected:  1106560" << endl;
        cout << "20!^(3) = " << IntCombinatorics::multiFactorial(20UL, 3UL) << "  expected:  2094400" << endl;
        cout << "21!^(3) = " << IntCombinatorics::multiFactorial(21UL, 3UL) << "  expected:  11022480" << endl;
        cout << "15!! = " << IntCombinatorics::doubleFactorial(15UL) << "  expected:  2027025" << endl;
        cout << "22!! = " << IntCombinatorics::doubleFactorial(22ULL) << "  expected:  81749606400" << endl;
        cout << "27!! = " << IntCombinatorics::doubleFactorial(27ULL) << "  expected:  213458046676875" << endl;
        cout << endl;

        cout << "14 choose 4 = " << IntCombinatorics::binom(14UL, 4UL) << "  expected: 1001" << endl;
        cout << "14 choose 10 = " << IntCombinatorics::binom(14UL, 10UL) << "  expected: 1001" << endl;
        cout << "50 choose 41 = " << IntCombinatorics::binom(50ULL, 41ULL) << "  expected:  2505433700" << endl;
    }
    catch ( const CombinatoricsException& ex )
    {
        cerr << "Combinatorics exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
