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


        vector<double> vmpgs;
    	mtcopy(ampgs, LEN, vmpgs);

        vector<size_t> idx;

        SampleOrder::order(vmpgs, idx);
        cout << "Indices in ascending order:" << endl;
        printIndices(idx);

        cout << endl << "Indices in descending order:" << endl;
        SampleOrder::order(vmpgs, idx, false);
        printIndices(idx);
        cout << endl;

        cout << "min(mpg): " << SampleOrder::min(vmpgs) << " (expected: 10.4)" << endl;
        cout << "max(mpg): " << SampleOrder::max(vmpgs) << " (expected: 33.9)" << endl;
        cout << "Index of the smallest mpg: " << SampleOrder::whichMin(vmpgs) << " (expected: 14)" << endl;
        cout << "Index of the largest mpg:  " << SampleOrder::whichMax(vmpgs) << " (expected: 19)" << endl;

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
