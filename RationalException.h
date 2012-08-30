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
@file RationalException.h

Declaration of the class RationalException

@author Jernej Kovacic
*/

#ifndef _RATIONALEXCEPTION_H_
#define	_RATIONALEXCEPTION_H_

#include "IMathException.h"

struct RationalException : public IMathException
{
    /// Enum with possible error codes
    enum err_codes {
        ZERO_DENOMINATOR,   /// attempted to set a zero denominator
        UNINVERTIBLE,       /// the fraction cannot be inverted
        DIVIDE_BY_ZERO,     /// attempt of division by zero
        OVERFLOW            /// integer overflow
    };

    err_codes error;    /// type of an error
    // Constructor
    RationalException(err_codes err);
    // Output a short description of the error
    void display(std::ostream& str = std::cerr) const;
};

#endif	/* RATIONALEXCEPTION_H */
