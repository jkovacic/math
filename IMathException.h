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
 * @file IMathException.h
 *
 * Common interface (pure virtual class) for various mathematical exceptions,
 * e.g. MatrixException, RationalException, etc.
 *
 * @author Jernej Kovacic
 */

#ifndef _MATH_IMATHEXCEPTION_H_
#define	_MATH_IMATHEXCEPTION_H_

#include <iostream>

namespace math
{

/**
 * @brief An "interface" class for all mathematical exceptions
 */
struct IMathException
{
    /**
    * Outputs a short error description to stdout
    *
    * @param str - stream, the error description will be written in (default: cerr) 
    */
    virtual void display(std::ostream& str = std::cerr) const = 0;

    virtual ~IMathException();
};

} // namespace math

#endif	// _MATH_IMATHEXCEPTION_H_

