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
 * A test module to test multi-threaded vector copy
 * functionality (mtcopy)
 */


#include <cstddef>
#include <vector>
#include <algorithm>
#include <iostream>
#include <new>
#include <stdexcept>

#include "mtcopy.h"


using namespace std;
using namespace math;



/*
 * A handy function to print contents of a std::vector
 */
template <class T>
void printVector(const vector<T>& vec)
{
    for ( typename vector<T>::const_iterator it=vec.begin(); it!=vec.end(); ++it )
    {
        cout << *it << "  ";
    }
    cout << endl;
}


/*
 * Unit test of mtcopy functions
 */
void mtcopyTest()
{
    try
    {
        const size_t N = 10;
        const int arr[] = { 23, 56, 12, 56, 23, -23, 92, -1, 45, -21 };
        vector<int> v;
        vector<int> v2;

        cout << "Copy the entire array: " << endl;
        mtcopy(arr, arr+N, v);
        cout << "Copy:     ";
        printVector(v);
        cout << "Expected: 23  56  12  56  23  -23  92  -1  45  -21" << endl;

        cout << "Copy a part of an array: " << endl;
        mtcopy(arr+1, 7, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: 56  12  56  23  -23  92  -1" << endl;

        cout << "Copy the entire vector: " << endl;
        sort(v.begin(), v.end());
        mtcopy(v, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: -23  -21  -1  12  23  23  45  56  56  92" << endl;

        cout << "Copy of a part of a vector: " << endl;
        mtcopy(v, 1, 5, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: -21  -1  12  23  23 " << endl;

        cout << "Test vector's range: " << endl;
        mtcopy(v, 7, 5, v2);
        cout << "Copy:     ";
        printVector(v2);
        cout << "Expected: 56  56  92 " << endl;
    }
    catch ( const bad_alloc& ba )
    {
        cerr << "Bad allocation exception called" << endl;
    }
    catch ( const out_of_range& oor )
    {
        cerr << "Out of range exception caught" << endl;
    }
    catch (...)
    {
        cerr << "Some other exception caught" << endl;
    }
}
