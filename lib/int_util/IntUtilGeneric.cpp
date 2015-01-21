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

#include "IntUtilGeneric.hpp"


/**
 * @file
 * @author Jernej Kovacic
 * 
 * Implementation of auxiliary functions that handle integer variables.
 */

// no #include "IntUtilGeneric.hpp" !!!
#include <limits>


// A namespace with "private" functions and classes
namespace math {  namespace IntUtil {  namespace __private
{

/*
 * A set of clases that check whether an integer is negative.
 * 
 * An integer value can only be negative if its type is signed. This could be
 * quickly checked by comparing (I)(-1) and (I)(0), however, if I is an unsigned type,
 * a compiler would raise a warning that "comparison of unsigned expression < 0 is 
 * always false". One workaround this would be partial specialization. As C++ does
 * allow partial specialization of templated functions, thiswork around using an
 * additional templated static class has been introduced.
 */

template <typename I, bool isSigned>
struct __CheckSign
{
    /*
     * A "general" implementation of the class when isSigned equals 'true'
     */
    static inline bool isNeg(const I& n)
    {
        // As I is a signed type, 'n' must be compared to 0:
        return ( n < static_cast<I>(0) );
    }    
};

// Partial specialization of __checkSignImpl for isSigned = false
template <typename I>
struct __CheckSign<I, false>
{
    /*
     * A specialization for isSigned == false, i.e. I is an unsigned integer type
     */
    static inline bool isNeg(const I& n)
    {
        // Since I is an unsigned type, 'n' can never be negative 
        (void) n;
        return false;
    }    
};

}}}  // namespace math::IntUtil::__private



/*
 * Checks the sign of 'n'.
 * Both, signed and unsigned, integer types are handled correctly.
 * 
 * @param n - integer value to check
 * 
 * @return 'true' if 'n' is strictly negative, 'false' otherwise
 */
template <typename I>
bool math::IntUtil::isNegative(const I& n)
{
    return 
      math::IntUtil::__private::__CheckSign<I, std::numeric_limits<I>::is_signed>::isNeg(n);
}


/*
 * A simple integer implementation of abs
 *
 * @param n
 * @return absolute value of 'n'
 */
template <typename I>
inline I math::IntUtil::absolute(const I& n)
{
    return ( false==math::IntUtil::isNegative<I>(n) ? n : -n );
}
