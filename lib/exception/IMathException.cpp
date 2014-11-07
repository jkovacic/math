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
 * Implementation of an empty destructor to
 * suppress warnings in derived classes.
 *
 * For more info, see:
 * http://stackoverflow.com/questions/127426/gnu-compiler-warning-class-has-virtual-functions-but-non-virtual-destructor
 */


#include "IMathException.hpp"

/**
 * Destructor
 */
math::IMathException::~IMathException()
{
    // Empty function
}