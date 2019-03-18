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
 * Declaration and partial implementation of the class CurveFittingGenericAb.h.
 * This is an abstract class, its derived classes implement various
 * curve fitting algorithms.
 */

#ifndef _MATH_CURVEFITTINGGENERICAB_HPP_
#define _MATH_CURVEFITTINGGENERICAB_HPP_


#include <new>
#include <cstddef>
#include <list>

#include "util/NumericUtil.hpp"
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
 
        /*
         * Implementation of internal class's operator< (required for sort)
         *
         * @param p - instance of point whose abscissa will be compared to this one's
         *
         * @return whether this's abscissa is smaller than p's
         */
        bool operator<(const CPoint& p) const
        {
            // just compare both points' abscissas
            return ( this->m_x < p.m_x );
        }
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


    /*
     * A convenience function to detect if any "duplicate" (with the same abscissa)
     * points have been entered.
     *
     * @note In order for this function to work properly, list of points must be sorted beforehand!
     *
     * @return true/false
     */
    bool _duplicatePoints() const
    {
        // No duplicate points are possible if the list is empty
        if ( true == this->m_points.empty() )
        {
            return false;
        }

        /*
         * Just traverse all points (except the last one) and compare its abscissa
         * to the next point's abscissa. This is why the list must be sorted.
         */
        bool retVal = false;
        for ( 
            typename std::list<typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator it=this->m_points.begin(); 
            it != this->m_points.end(); 
            ++it )
        {
            /*
             * The function also requires a pointer to the next element. As operator+ is not defined
             * for STL list's iterator, another iterator is declared, its initial value is equal to
             * it, then its operator++ (which is defined) is called immediately.
             */
            typename std::list< typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator next = it;
            ++next;

            // the last point, nothing else to compare its abscissa to
            if ( this->m_points.end()==next )
            {
                break;  // out of for
            }

            // abscissas are equal if their difference is zero (or a very small value)
            if ( true==math::NumericUtil::isZero<F>( it->m_x - next->m_x ) )
            {
                // no need to search further when one duplicate abscissa is found
                retVal = true;
                break;  // out of for
            }
        }

        return retVal;
    }
    
    
    /*
     * Performs some necessary checks prior curve generation. If any check
     * fails, the appropriate exception will be thrown.
     *
     * @note The function also sorts entered points
     *
     * @throw appropriate CurveFittingException if any check is not passed
     */
    void _curveGenerationCheck()
    {
        // Curve must not be generated yet
        if ( true==this->m_curveGenerated )
        {
            throw math::CurveFittingException(math::CurveFittingException::CURVE_ALREADY_GENERATED);
        }

        // At least one point must be entered
        if ( true==this->m_points.empty() )
        {
            throw math::CurveFittingException(math::CurveFittingException::NO_POINTS);
        }

        // sort the points by abscissae's values in ascending order
        this->_sortPoints();

        // and detect any "duplicate" points
        if ( true == this->_duplicatePoints() )
        {
            throw math::CurveFittingException(math::CurveFittingException::DUPLICATE_POINTS);
        }
    }


    /*
     * Sort entered points by points' abscissa values in ascending order.
     * CPoint's 'operaator<' method is used as the comparison criteria.
     */
    void _sortPoints()
    {
        // Nothing to do if no points have been entered yet
        if ( false==this->m_points.empty() )
        {
            // apply list's sort method
            this->m_points.sort();
        }
    }


    /*
     * As an abstract class cannot have constructors, this function performs
     * some initialization tasks, common to all derived classes.
     */
    void _init()
    {
        // As an instance has just been created, a curve cannot be generated yet.
        this->m_curveGenerated = false;
        // Init the linked list (delete all points from the list if any exists)
        this->m_points.clear();
    }

public:

    /**
     * Assignment operator (=). It is only applicable for instances of the same class.
     *
     * @param orig - instance of a derived class, whose points will be copied to this one
     *
     * @return reference to itself
     *
     * @throw CurveFittingException if allocation of memory fails
     */
    virtual CurveFittingGenericAb<F>& operator=(const CurveFittingGenericAb<F>& orig)
    {
        // nothing to do if attempting to assign to itself
        if ( &orig==this )
        {
            return *this;
        }

        try
        {
            // initialize
            this->_init();
            // and copy all points
            this->m_points = orig.m_points;
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
        }

        return *this;
    }



    /**
     * A convenience function that copies points from 'porig' to 'this'.
     * Unlike operator=(), 'porig' can be a pointer to any instance that is derived
     * from CurveFittingGenericAb.
     *
     * @param porig - a pointer to an instance whose points will be copied to 'this'
     *
     * @throw CurveFittingException if 'porig' is invalid or if allocation of memory fails
     */
    CurveFittingGenericAb<F>& copy(const CurveFittingGenericAb<F>* porig)
    {
        // Sanity check: porig must not be NULL
        if ( NULL==porig )
        {
            throw math::CurveFittingException(math::CurveFittingException::INVALID_ARGUMENT);
        }

        // ... and it must point to an object that is derived from CurveFittingGenericAb
        const math::CurveFittingGenericAb<F>* const dyncporig = dynamic_cast<const math::CurveFittingGenericAb<F>* >(porig);
        if ( NULL==dyncporig )
        {
            throw math::CurveFittingException(math::CurveFittingException::INVALID_ARGUMENT);
        }

        // if all checks have passed, just reuse functionality of operator=
        return this->operator=(*dyncporig);
     }



    /**
     * Adds a point to be used for generation of the best fitting curve.
     *
     * @note Points can be entered in an arbitrary order.
     * @note No "duplicate" points (with the same abscissa) are allowed. This
     *       check is  not performed at this time but later the generation of
     *       the curve will fail.
     * @note Once a point has been entered, it cannot be removed.
     *
     * @param x - point's abscissa
     * @param y - point's ordinate
     *
     * @return reference to itself
     *
     * @throw CurveFittingException if the curve has already been generated or if allocation of memory fails
     */
    virtual CurveFittingGenericAb<F>& enterPoint(const F& x, const F& y)
    {
        // No additional points are allowed once the curve has been generated
        if ( true == this->m_curveGenerated )
        {
            throw math::CurveFittingException(math::CurveFittingException::ADD_POINT_NOT_ALLOWED);
        }

        // Check that the nr. of points would not exceed the max. allowed list's size
        if ( this->m_points.size() == this->m_points.max_size() )
        {
            throw math::CurveFittingException(math::CurveFittingException::ADD_POINT_NOT_ALLOWED);
        }
    
        try
        {
            CPoint p;
            p.m_x = x;
            p.m_y = y;
            this->m_points.push_back(p);
        }
        catch ( const std::bad_alloc& ba )
        {
            throw math::CurveFittingException(math::CurveFittingException::OUT_OF_MEMORY);
        }

        return *this;
    }



    /**
     * @return Number of all points entered
     */
    std::size_t nrPoints() const
    {
        return m_points.size();
    }



    /**
     * @return has the curve already been generated?
     */
    bool curveReady() const
    {
        return this->m_curveGenerated;
    }



    /**
     * @return the smallest abscissa of all points entered till this moment
     *
     * @throw CurveFittingException if no points have been entered
     */
    F lowerBound() const
    {
        // only applicable when at least one point was entered
        if ( true == this->m_points.empty() )
        {
            throw math::CurveFittingException(math::CurveFittingException::NO_POINTS);
        }

        // If the curve has already been generated, points are already sorted,
        // just return the abscissa of the first point in the list.
        if ( true==this->m_curveGenerated )
        {
            return this->m_points.front().m_x;
        }

        // traverse the list and find the smallest abscissa
        typename std::list< typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator currMin = 
            this->m_points.begin();
        for ( typename std::list<typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator it=m_points.begin();
                   it!=this->m_points.end(); ++it )
        {
            if ( currMin->m_x > it->m_x )
            {
                currMin = it;
            }
        }

        return currMin->m_x;
    }

    
    
    /**
     * @return the highest abscissa of all points entered till this moment
     *
     * @throw CurveFittingException if no points have been entered
     */
    F upperBound() const
    {
        // only applicable when at least one point was entered
        if ( true == this->m_points.empty() )
        {
            throw math::CurveFittingException(math::CurveFittingException::NO_POINTS);
        }

        // If the curve has already been generated, points are already sorted,
        // just return the abscissa of the last point in the list.
        if ( true==this->m_curveGenerated )
        {
            return this->m_points.back().m_x;
        }

        // traverse the list and find the highest abscissa
        typename std::list<typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator currMax = 
            this->m_points.begin();
        for ( 
            typename std::list<typename math::CurveFittingGenericAb<F>::CPoint>::const_iterator it=this->m_points.begin();
            it != this->m_points.end(); 
            ++it )
        {
            if ( currMax->m_x < it->m_x )
            {
                currMax = it;
            }
        }

        return currMax->m_x;
    }
    

    // generate a curve that fits entered points best
    virtual void generateCurve(const std::size_t degree) = 0;


    // value of the curve at the given abscissa
    virtual F valueAt(const F& x, const bool strict=true) const = 0;


    /**
     * Destructor
     */
    virtual ~CurveFittingGenericAb()
    {
        // List's destructors would probably clean up this automatically.
        // Anyway, let us clear the list, just to be aware of allocated resources.
        this->m_points.clear();

        // Other dynamically allocated memory (via malloc or new) should be freed here.
        // There are no other resources to release.
    }
 
};

}  // namespace math


#endif  // _MATH_CURVEFITTINGGENERICAB_HPP_
