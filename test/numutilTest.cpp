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
 * A test module to test functions from the
 * namespace NumericUtil.
 */


#include <iostream>

#include "NumericUtil.h"


using namespace std;
using namespace math;


void numutilTest()
{
    const long double defEps = NumericUtil::getEPS<long double>();
    cout << "Default 'eps' for long double: " << defEps << endl;

    cout << "Setting 'eps' to 1e-9." << endl;
    NumericUtil::setEPS<long double>(1e-9L);
    cout << "'eps' set to: " << NumericUtil::getEPS<long double>() << endl;

    cout << "Setting 'eps' to 100 machine epsilons." << endl;
    NumericUtil::setMultEPS<long double>(100.0L);
    cout << "'eps' set to: " << NumericUtil::getEPS<long double>() << endl;

    cout << "Resetting 'eps' to the system default, i.e. machine epsilon." << endl;
    NumericUtil::setEPS<long double>();
    cout << "'eps' set to: " << NumericUtil::getEPS<long double>() << endl;

    cout << "Setting 'eps' back to the original value." << endl;
    NumericUtil::setEPS<long double>(defEps);
    cout << "'eps' set to: " << NumericUtil::getEPS<long double>() << endl;
}
