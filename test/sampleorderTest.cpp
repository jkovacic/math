/*
Copyright 2015, Jernej Kovacic

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
 * A test module to test functionality in namespaces that perform
 * searching of elements' indices in stably sorted vectors (SampleOrder).
 */


/*
 * Note: results are reproduced in 'scripts/test/sampleorder.R'.
 */


#include <cstddef>
#include <iostream>
#include <vector>
#include <new>

#include "mtcopy.h"
#include "SampleOrderGeneric.h"
#include "SampleOrderException.h"

using namespace std;
using namespace math;


// Nr. of entries in mtcars:
#define LEN                    ( 32 )



/*
 * A convenience function that prints a vector of
 * integers as a space separated string.
 */
void printIndices(const vector<size_t>& vi)
{
    for ( vector<size_t>::const_iterator it=vi.begin();
          it!=vi.end(); ++it )
    {
        cout << *it << " ";
    }

    cout << endl;
}


/*
 * Test of functions that search indices of a sorted sample vector
 */
void sampleOrderTest()
{
    try
    {
        /* Cars' mpg (miles per gallon) from R's data frame 'mtcars' */
        const double ampgs[ LEN ] =
            { 21.0, 21.0, 22.8, 21.4, 18.7, 18.1, 14.3, 24.4,
              22.8, 19.2, 17.8, 16.4, 17.3, 15.2, 10.4, 10.4,
              14.7, 32.4, 30.4, 33.9, 21.5, 15.5, 15.2, 13.3,
              19.2, 27.3, 26.0, 30.4, 15.8, 19.7, 15.0, 21.4 };

        /*
           # Equivalent of the following command in R:
           data(mtcars)
         */

        vector<double> vmpgs;
    	mtcopy(ampgs, LEN, vmpgs);

        /*
         * R code to test rearrangement of vector's indices in ascending and descending order:
         *

           order(mtcars$mpg) - 1
            [1] 14 15 23  6 16 30 13 22 21 28 11 12 10  5  4  9 24 29  0  1  3 31 20  2
           [25]  8  7 26 25 18 27 17 19

           order(mtcars$mpg, decreasing=TRUE) - 1
            [1] 19 17 18 27 25 26  7  2  8 20  3 31  0  1 29  9 24  4  5 10 12 11 28 21
           [25] 13 22 30 16  6 23 14 15
         */

        vector<size_t> idx;

        SampleOrder::order(vmpgs, idx);
        cout << "Indices in ascending order:" << endl;
        printIndices(idx);

        cout << endl << "Indices in descending order:" << endl;
        SampleOrder::order(vmpgs, idx, false);
        printIndices(idx);

    }
    catch ( const SampleOrderException& ex )
    {
        cerr << "SampleOrder exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch ( const bad_alloc& ba )
    {
        cerr << "Could not allocate memory." << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
