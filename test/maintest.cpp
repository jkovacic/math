/*
Copyright 2011, 2013 Jernej Kovacic

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
 * A collection of basic smoke tests for all mathematical classes
 * and functions.
 */


#include <iostream>


// Declaration of functions implemented in other files
extern void numutilTest();
extern void mtcopyTest();
extern void sampleOrderTest();
extern void quaternionTest();
extern void rationalTest();
extern void matrixTest();
extern void rationalMatrixTest();
extern void polynomialTest();
extern void lineqSolverTest();
extern void curveFittingTest();
extern void intExponentiaorTest();
extern void intFactorizationTest();
extern void intCombinatoricsTest();
extern void combinatoricsTest();
extern void specfunTest();
extern void calculusTest();
extern void rootFindTest();
extern void statisticsTest();
extern void probDistributionTest();


using namespace std;


/*
 * Main function that starts several smoke test modules
 */
int main(int argc, const char* argv[])
{

    cout << "N U M E R I C   U T I L   T E S T" << endl << endl;
    numutilTest();

    cout << endl << "M T C O P Y   T E S T" << endl << endl;
    mtcopyTest();

    cout << endl << "S A M P L E   O R D E R   T E S T" << endl << endl;
    sampleOrderTest();

    cout << endl << "Q U A T E R N I O N   T E S T" << endl << endl;
    quaternionTest();

    cout << endl << "R A T I O N A L   T E S T" << endl << endl;
    rationalTest();

    cout << endl << "M A T R I X   T E S T" << endl << endl;
    matrixTest();

    cout << endl << "R A T I O N A L   M A T R I X   T E S T" << endl << endl;
    rationalMatrixTest();

    cout << endl << "L I N E A R   E Q U A T I O N   S O L V E R   T E S T" << endl << endl;
    lineqSolverTest();

    cout << endl << "P O L Y N O M I A L   T E S T" << endl << endl;
    polynomialTest();

    cout << endl << "C U R V E   F I T T I N G   T E S T" << endl << endl;
    curveFittingTest();
    
    cout << endl << "I N T E X P O N E N T I A T O R   T E S T" << endl << endl;
    intExponentiaorTest();
    
    cout << endl << "I N T   F A C T O R I Z A T I O N   T E S T" << endl << endl;
    intFactorizationTest();

    cout << endl << "I N T   C O M B I N A T O R I C S   T E S T" << endl << endl;
    intCombinatoricsTest();
    
    cout << endl << "C O M B I N A T O R I C S   T E S T" << endl << endl;
    combinatoricsTest();

    cout << endl << "S P E C I A L   F U N C T I O N   T E S T" << endl << endl;
    specfunTest();

    cout << endl << "C A L C U L U S   T E S T" << endl << endl;
    calculusTest();

    cout << endl << "R O O T   F I N D   T E S T" << endl <<endl;
    rootFindTest();

    cout << endl << "S T A T I S T I C S   T E S T" << endl << endl;
    statisticsTest();

    cout << endl << "P R O B A B I L I T Y   D I S T R I B U T I O N   T E S T" << endl << endl;
    probDistributionTest();

    return 0;

    (void) argc;
    (void) argv;
}
