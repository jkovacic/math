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
 * A test module to test int. factorization classes
 * (IntFactorization)
 */

#include <iostream>
#include <map>
#include <set>

#include "IntFactorization.h"
#include "IntFactorizationException.h"


using namespace std;
using namespace math;


/*
 * Test of integer factorization
 */
void intFactorizationTest()
{
    try
    {
        // Find integer square roots of the following numbers:
        const unsigned int NS = 3;
        unsigned long int s[NS] = { 12, 100, 37423 };

        // Expected results: 3, 10 and 193
        for ( unsigned int i=0; i<NS; ++i )
        {
            cout << "Int sqrt of " << s[i] << ": " << IntFactorization::intSqrt(s[i]) << endl;
        }
        cout << endl;

        // Test primality test and the algorithm for the next prime,
        // Reference: http://primes.utm.edu/lists/small/10000.txt
        const unsigned int NP = 10;
        unsigned long int p[NP] = { 3, 6, 15, 37, 257, 703, 907, 101861, 102601, 104597 };
        /*
         * Expected results:
         * (T, 5), (F, 7), (F, 17), (T, 41), (T, 263), (F, 709), (T, 911)
         * (F, 101863), (F, 102607), (T, 104623)
         */
        for ( unsigned int i=0; i<NP; ++i )
        {
            cout << p[i] << " is ";
            cout << ( true==IntFactorization::isPrime(p[i]) ? "" : "not ");
            cout << "a prime\tNext prime: ";
            cout << IntFactorization::nextPrime(p[i]) << endl;
        }
        cout << endl;

        // Test prime factorization and divisors
        const unsigned int NF = 5;
        unsigned long nf[NF] = {245, 6784, 21737, 195327, 3428543 };

        /*
         * Expected results:
         * 245 = 5 * 7^2      [1, 5, 7, 35, 49, 245]
         * 6784 = 2^7 * 53    [1, 2, 4, 8, 16, 32, 53, 64, 106, 128, 212, 424, 848, 1696, 3392, 6784]
         * 21737 = 21737      [1, 21737]
         * 195327 = 3^2 * 11 * 1973     [1, 3, 9, 11, 33, 99, 1973, 5919, 17757, 21703, 65109, 195327]
         * 3428543 = 17 * 41 * 4919     [1, 17, 41, 697, 4919, 83623, 201679, 3428543]
         */
        for ( unsigned int i=0; i<NF; ++i )
        {
            map<unsigned long int, unsigned int> f;
            IntFactorization::factor(nf[i], f);
            cout << nf[i] << " = ";
            for ( map<unsigned long int, unsigned int>::iterator it=f.begin();
                  it!=f.end(); ++it )
            {
                if ( it!=f.begin() )
                {
                    cout << " * ";
                }

                cout << it->first;

                if ( it->second > 1 )
                {
                    cout << '^' << it->second;
                }
            }

            set<unsigned long int> d;
            IntFactorization::divisors(nf[i], d);
            cout << "\tdivisors: ";
            for ( set<unsigned long int>::iterator it=d.begin(); it!=d.end(); ++it )
            {
                if ( it!=d.begin() )
                {
                    cout << ", ";
                }
                cout << *it;
            }

            cout << endl;
        }
        cout << endl;

        // Finally test the GCD and LCM algorithms:
        const unsigned int NC = 10;
        unsigned long int c[NC] = { 500, 1000, 85, 3428543 , 15, 100, 3, 57, 234, 7643 };
        /*
         * Expected results:
         * n1 = 500, n2 = 1000: gcd = 500, lcm = 1000
         * n1 = 85, n2 = 3428543: gcd = 17, lcm = 17142715
         * n1 = 15, n2 = 100: gcd = 5, lcm = 300
         * n1 = 3, n2 = 57: gcd = 3, lcm = 57
         * n1 = 234, n2 = 7643: gcd = 1, lcm = 1788462
         */
        for ( unsigned int i=0; i<NC/2; ++i )
        {
            const unsigned long int n1 = c[2*i];
            const unsigned long int n2 = c[2*i+1];
            unsigned long int gcd = IntFactorization::greatestCommonDivisor(n1, n2);
            unsigned long int lcm = IntFactorization::leastCommonMultiple(n1, n2);

            cout << "n1 = " << n1 << ", n2 = " << n2 << ": ";
            cout << "gcd = " << gcd << ", lcm = " << lcm << ", ";
            cout << "gcd*lcm = " << gcd*lcm << ", ";
            cout << "n1*n2 = " << n1*n2;
            cout << endl;
        }

        cout << endl;
    }
    catch ( const IntFactorizationException& ex )
    {
        cerr << "IntFactorization exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
