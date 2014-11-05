/*
Copyright 2013, Jernej Kovacic

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
 * Declaration of the class CurveFittingGenericAb.h. This is an
 * abstract class, its derived classes implement various
 * curve fitting algorithms.
 */

#ifndef _MATH_CURVEFITTINGGENERICAB_H_
#define _MATH_CURVEFITTINGGENERICAB_H_

#include <list>
#include <cstddef>

#include "CurveFittingException.h"


namespace math
{

/**
 * @brief Common functionality for all curve fitting/interpolation algorithms.
 * 
 * The class is abstract and cannot be instantiated. More specialized algorithms
 * are implemented by its derived classes.
 */
template<class T>
class CurveFittingGenericAb
{

protected:
    /*
     * Internal class that stores coordinates of 2D points
     */
    struct CPoint
    {
        T p_x;          // abscissa
        T p_y;          // ordinate
        // required by sort algorithms:
        bool operator<(const CPoint& p) const;
    };

protected:
    /*
     * A container to store entered points. As an arbitrary number of points
     * may be entered and the points are going to be rearranged (sorted),
     * a linked list is the most convenient container.
     */
    std::list<CPoint> points;

    // has the fitting curve been generated yet?
    bool curveGenerated;

    // have any "duplicate" points been entered?
    bool duplicatePoints() const;
    // some necessary checks prior to generation of the curve. It also sorts points
    void curveGenerationCheck() throw (CurveFittingException);

    // sorts entered points
    void sortPoints();

    // some constructor-like functionality
    void init();

public:
    // assignment operator (only applicable for variables of the same type)
    virtual CurveFittingGenericAb<T>& operator=(const CurveFittingGenericAb<T>& orig) throw (CurveFittingException);
    // copies points from 'porig' to this (applicable for all types derived from this one)
    CurveFittingGenericAb<T>& copy(const CurveFittingGenericAb<T>* porig) throw (CurveFittingException);

    // enter a point
    // Note: the function is virtual as derived classes may perform a sort of input control
    virtual CurveFittingGenericAb<T>& enterPoint(const T& x, const T& y) throw (CurveFittingException);

    // Number of points entered
    size_t nrPoints() const;

    // has the curve been generated?
    bool curveReady() const;

    // The lowest and highest abscissa value of all entered points
    T lowerBound() const throw (CurveFittingException);
    T upperBound() const throw (CurveFittingException);

    // generate a curve that fits entered points best
    virtual void generateCurve(size_t degree) throw (CurveFittingException) = 0;

    // value of the curve at the given abscissa
    virtual T valueAt(const T& x, bool strict=true) const throw (CurveFittingException) = 0;

    // Destructor
    virtual ~CurveFittingGenericAb();
};

// Definition could be included into the namespace declaration, but it
// would cause conflicts with some extra included stdlib header files.
}  // namespace math

// DEFINITION

// This is a templated class, so its definition must follow its declaration.
// When building, THIS file must be compiled.
// Alternatively the definition can be included into this file.
#include "CurveFittingGenericAb.cpp"

#endif // _MATH_CURVEFITTINGGENERICAB_H_
