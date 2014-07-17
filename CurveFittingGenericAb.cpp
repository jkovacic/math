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
 * @file CurveFittingGenericAb.cpp
 *
 * Implementation of the class CurveFittingGenericAb.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 *
 * @author Jernej Kovacic
 */

// Deliberately there is no #include "CurveFittingGenericAb.h" !
#include <new>

#include "NumericUtil.h"
#include "CurveFittingException.h"


/*
 * Implementation of internal class's operator< (required for sort)
 *
 * @param p - instance of point whose abscissa will be compared to this one's
 *
 * @return whether this's abscissa is smaller than p's
 */
template<class T>
bool math::CurveFittingGenericAb<T>::CPoint::operator<(const math::CurveFittingGenericAb<T>::CPoint& p) const
{
    // just compare both points' abscissas
    return ( p_x < p.p_x );
}

/*
 * As an abstract class cannot have constructors, this function performs
 * some initialization tasks, common to all derived classes.
 */
template<class T>
void math::CurveFittingGenericAb<T>::init()
{
    // As an instance has just been created, a curve cannot be generated yet.
    curveGenerated = false;
    // Init the linked list (delete all points from the list if any exists)
    points.clear();
}

/**
 * Assignment operator (=). It is only applicable for instances of the same class.
 *
 * @param orig - instance of a derived class, whose points will be copied to this one
 *
 * @return reference to itself
 *
 * @throw CurveFittingException if allocation of memory fails
 */
template<class T>
math::CurveFittingGenericAb<T>& math::CurveFittingGenericAb<T>::operator=(const math::CurveFittingGenericAb<T>& orig) throw (math::CurveFittingException)
{
    // nothing to do if attempting to assign to itself
    if ( &orig==this )
    {
        return *this;
    }

    try
    {
        // initialize
        init();
        // and copy all points
        points = orig.points;
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
template<class T>
math::CurveFittingGenericAb<T>& math::CurveFittingGenericAb<T>::copy(const math::CurveFittingGenericAb<T>* porig) throw (math::CurveFittingException)
{
    // Sanity check: porig must not be NULL
    if ( NULL==porig )
    {
        throw math::CurveFittingException(math::CurveFittingException::INVALID_ARGUMENT);
    }

    // ... and it must point to an object that is derived from CurveFittingGenericAb
    const math::CurveFittingGenericAb<T>* dyncporig = dynamic_cast<const math::CurveFittingGenericAb<T>* >(porig);
    if ( NULL==dyncporig )
    {
        throw math::CurveFittingException(math::CurveFittingException::INVALID_ARGUMENT);
    }

    // if all checks have passed, just reuse functionality of operator=
    return operator=(*dyncporig);
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
template<class T>
math::CurveFittingGenericAb<T>& math::CurveFittingGenericAb<T>::enterPoint(const T& x, const T& y) throw (math::CurveFittingException)
{
    // No additional points are allowed once the curve has been generated
    if ( true==curveGenerated )
    {
        throw math::CurveFittingException(math::CurveFittingException::ADD_POINT_NOT_ALLOWED);
    }

    // Check that the nr. of points would not exceed the max. allowed list's size
    if ( points.size()==points.max_size() )
    {
        throw math::CurveFittingException(math::CurveFittingException::ADD_POINT_NOT_ALLOWED);
    }
    
    try
    {
        CPoint p;
        p.p_x = x;
        p.p_y = y;
        points.push_back(p);
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
template<class T>
size_t math::CurveFittingGenericAb<T>::nrPoints() const
{
    return points.size();
}

/**
 * @return has the curve already been generated?
 */
template<class T>
bool math::CurveFittingGenericAb<T>::curveReady() const
{
    return curveGenerated;
}

/**
 * @return the smallest abscissa of all points entered till this moment
 *
 * @throw CurveFittingException if no points have been entered
 */
template<class T>
T math::CurveFittingGenericAb<T>::lowerBound() const throw (math::CurveFittingException)
{
    // only applicable when at least one point was entered
    if ( true==points.empty() )
    {
        throw math::CurveFittingException(math::CurveFittingException::NO_POINTS);
    }

    // If the curve has already been generated, points are already sorted,
    // just return the abscissa of the first point in the list.
    if ( true==curveGenerated )
    {
        return points.front().p_x;
    }

    // traverse the list and find the smallest abscissa
    typename std::list< typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator currMin = points.begin();
    for ( typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator it=points.begin();
               it!=points.end(); ++it )
    {
        if ( currMin->p_x > it->p_x )
        {
            currMin = it;
        }
    }

    return currMin->p_x;
}

/**
 * @return the highest abscissa of all points entered till this moment
 *
 * @throw CurveFittingException if no points have been entered
 */
template<class T>
T math::CurveFittingGenericAb<T>::upperBound() const throw (math::CurveFittingException)
{
    // only applicable when at least one point was entered
    if ( true==points.empty() )
    {
        throw math::CurveFittingException(math::CurveFittingException::NO_POINTS);
    }

    // If the curve has already been generated, points are already sorted,
    // just return the abscissa of the last point in the list.
    if ( true==curveGenerated )
    {
        return points.back().p_x;
    }

    // traverse the list and find the highest abscissa
    typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator currMax = points.begin();
    for ( typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator it=points.begin(); it!=points.end(); ++it )
    {
        if ( currMax->p_x < it->p_x )
        {
            currMax = it;
        }
    }

    return currMax->p_x;
}

/*
 * A convenience function to detect if any "duplicate" (with the same abscissa)
 * points have been entered.
 *
 * @note In order for this function to work properly, list of points must be sorted beforehand!
 *
 * @return true/false
 */
template<class T>
bool math::CurveFittingGenericAb<T>::duplicatePoints() const
{
    // No duplicate points are possible if the list is empty
    if ( true==points.empty() )
    {
        return false;
    }

    /*
        Just traverse all points (except the last one) and compare its abscissa
        to the next point's abscissa. This is why the list must be sorted.
     */
    bool retVal = false;
    for ( typename std::list<typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator it=points.begin(); it!=points.end(); ++it )
    {
        /*
            The function also requires a pointer to the next element. As operator+ is not defined
            for STL list's iterator, another iterator is declared, its initial value is equal to
            it, then its operator++ (which is defined) is called immediately.
        */
        typename std::list< typename math::CurveFittingGenericAb<T>::CPoint>::const_iterator next = it;
        ++next;

        // the last point, nothing else to compare its abscissa to
        if ( points.end()==next )
        {
            break;  // out of for
        }

        // abscissas are equal if their difference is zero (or a very small value)
        if ( true==math::NumericUtil<T>::isZero( it->p_x - next->p_x ) )
        {
            // no need to search further when one duplicate abscissa is found
            retVal = true;
            break;  // out of for
        }
    }

    return retVal;
}

/*
 * Sort entered points by points' abscissa values in ascending order.
 * CPoint's 'operaator<' method is used as the comparison criteria.
 */
template<class T>
void math::CurveFittingGenericAb<T>::sortPoints()
{
    // Nothing to do if no points have been entered yet
    if ( false==this->points.empty() )
    {
        // apply list's sort method
        this->points.sort();
    }
}

/*
 * Performs some necessary checks prior curve generation. If any check
 * fails, the appropriate exception will be thrown.
 *
 * @note The function also sorts entered points
 *
 * @throw appropriate CurveFittingException if any check is not passed
 */
template<class T>
void math::CurveFittingGenericAb<T>::curveGenerationCheck() throw (math::CurveFittingException)
{
    // Curve must not be generated yet
    if ( true==this->curveGenerated )
    {
        throw math::CurveFittingException(math::CurveFittingException::CURVE_ALREADY_GENERATED);
    }

    // At least one point must be entered
    if ( true==this->points.empty() )
    {
        throw math::CurveFittingException(math::CurveFittingException::NO_POINTS);
    }

    // sort the points by abscissae's values in ascending order
    sortPoints();

    // and detect any "duplicate" points
    if ( true==duplicatePoints() )
    {
        throw math::CurveFittingException(math::CurveFittingException::DUPLICATE_POINTS);
    }
}

/**
 * Destructor
 */
 template<class T>
 math::CurveFittingGenericAb<T>::~CurveFittingGenericAb()
 {
    // List's destructors would probably clean up this automatically.
    // Anyway, let us clear the list, just to be aware of allocated resources.
     points.clear();

    // Other dynamically allocated memory (via malloc or new) should be freed here.
    // There are no other resources to release.
 }
