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
 * An internal header file, it should not be included directly.
 *
 * Declaration of the class CurveFittingGenericAb.h. This is an
 * abstract class, its derived classes implement various
 * curve fitting algorithms.
 */

#ifndef _MATH_CURVEFITTINGGENERICAB_HPP_
#define _MATH_CURVEFITTINGGENERICAB_HPP_

#include <list>
#include <cstddef>

#include "exception/CurveFittingException.hpp"


namespace math
{

/**
 * @brief Common functionality for all curve fitting/interpolation algorithms.
 * 
 * The class is abstract and cannot be instantiated. More specialized algorithms
 * are implemented by its derived classes.
 */
template <typename F>
class CurveFittingGenericAb
{

protected:
    /*
     * Internal class that stores coordinates of 2D points
     */
    struct CPoint
    {
        F m_x;          // abscissa
        F m_y;          // ordinate
        // required by sort algorithms:
        bool operator<(const CPoint& p) const;
    };

protected:
    /*
     * A container to store entered points. As an arbitrary number of points
     * may be entered and the points are going to be rearranged (sorted),
     * a linked list is the most convenient container.
     */
    std::list<CPoint> m_points;

    // has the fitting curve been generated yet?
    bool m_curveGenerated;

    // have any "duplicate" points been entered?
    bool _duplicatePoints() const;
    // some necessary checks prior to generation of the curve. It also sorts points
    void _curveGenerationCheck() throw (CurveFittingException);

    // sorts entered points
    void _sortPoints();

    // some constructor-like functionality
    void _init();

public:
    // assignment operator (only applicable for variables of the same type)
    virtual CurveFittingGenericAb<F>& operator=(const CurveFittingGenericAb<F>& orig) throw (CurveFittingException);
    // copies points from 'porig' to this (applicable for all types derived from this one)
    CurveFittingGenericAb<F>& copy(const CurveFittingGenericAb<F>* porig) throw (CurveFittingException);

    // enter a point
    // Note: the function is virtual as derived classes may perform a sort of input control
    virtual CurveFittingGenericAb<F>& enterPoint(const F& x, const F& y) throw (CurveFittingException);

    // Number of points entered
    std::size_t nrPoints() const;

    // has the curve been generated?
    bool curveReady() const;

    // The lowest and highest abscissa value of all entered points
    F lowerBound() const throw (CurveFittingException);
    F upperBound() const throw (CurveFittingException);

    // generate a curve that fits entered points best
    virtual void generateCurve(const std::size_t degree) throw (CurveFittingException) = 0;

    // value of the curve at the given abscissa
    virtual F valueAt(const F& x, const bool strict=true) const throw (CurveFittingException) = 0;

    // Destructor
    virtual ~CurveFittingGenericAb();
};

}  // namespace math

// DEFINITION
#include "curve_fit/CurveFittingGenericAb.cpp"

#endif  // _MATH_CURVEFITTINGGENERICAB_HPP_
