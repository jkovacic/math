/*
Copyright 2011, Jernej Kovacic

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
@file NumericUtil.h

Declaration of the class NumericUtil, a collection of some useful numerical utilities

@author Jernej Kovacic
*/

#ifndef NUMERICUTIL_H
#define	NUMERICUTIL_H

template<class T>
class NumericUtil
{
private:
    // Equivalent to eps in MATLAB (and its clones)
    // Its value depends on type T
    static T EPS;
public:
    // Does the value equal zero? Note that in numerical mathematics
    // (where mostly float or double values are in use), "equals" has a
    // different meaning than in discrete mathematics (int etc.)
    static bool isZero(T value);
};

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.

#include "NumericUtil.cpp"

#endif	/* NUMERICUTIL_H */
