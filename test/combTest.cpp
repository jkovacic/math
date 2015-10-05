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
 * A test module to test combination classes
 * (PermutationGeneric, CombinationGeneric)
 */


/*
 * Note: results are reproduced in 'scripts/test/comb.py'.
 */


#include <iostream>
#include <cstddef>
#include <vector>
#include <list>
#include <set>

#include "PermutationGeneric.h"
#include "CombinationGeneric.h"
#include "CombinatoricsException.h"


using namespace std;
using namespace math;

/*
 * Test of classes that list all permutations and combinations of a sequence
 */
void combinatoricsTest()
{
    try
    {
        cout << "Permutations:" << endl;

        vector<char> chvect;
        for ( char ch='a'; ch<='e'; ++ch )
        {
            chvect.push_back(ch);
        }

        PermutationGeneric<char> perm(chvect);
        list<list<char> > lp;
        size_t cnt = 1;
        while ( true==perm.hasNext() )
        {
            perm.next(lp, 25);
            for ( list<list<char> >::const_iterator lit=lp.begin(); lit!=lp.end(); ++lit )
            {
                cout << cnt++ << ": ";
                for ( list<char>::const_iterator cit=lit->begin(); cit!=lit->end(); ++cit )
                {
                    cout << *cit;
                }
                cout << endl;
            }
        }

        cout << endl << "Combinations:" << endl;
        CombinationGeneric<char> comb(chvect);
        list<set<char> > lc;
        for ( size_t k=1; k<=chvect.size(); ++k )
        {
            cout << endl << "K = " << k << " : " << endl << endl;
            comb.setK(k);
            cnt = 1;
            while ( true==comb.hasNext() )
            {
                comb.next(lc, 5);
                for ( list<set<char> >::const_iterator lit=lc.begin(); lit!=lc.end(); ++lit )
                {
                    cout << cnt++ << ": ";
                    for ( set<char>::const_iterator cit=lit->begin(); cit!=lit->end(); ++cit )
                    {
                        cout << *cit;
                    }
                    cout << endl;
                }
            }
        }

        cout << endl;
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
