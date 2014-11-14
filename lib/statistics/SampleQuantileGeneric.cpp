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
 * Implementation of the class SampleQuantileGeneric that estimates
 * quantiles of a sample.
 *
 * As the class is templated, this file must not be compiled.
 * Instead it must be included after the class declaration in the .h file
 */


// deliberately there is no #include "SampleQuantileGeneric.hpp" !

#include <vector>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <set>

#include "util/mtcopy.hpp"
#include "exception/StatisticsException.hpp"
#include "util/NumericUtil.hpp"
#include "SampleQuantileGeneric.hpp"


/*
 * Definitions of the most commonly used numerical constants:
 */
template <class T>
const T math::SampleQuantileGeneric<T>::ONE ( math::NumericUtil<T>::ONE );

template <class T>
const T math::SampleQuantileGeneric<T>::TWO ( static_cast<T>(2) );

template <class T>
const T math::SampleQuantileGeneric<T>::THREE ( static_cast<T>(3) );

template <class T>
const T math::SampleQuantileGeneric<T>::FOUR ( static_cast<T>(4) );

template <class T>
const T math::SampleQuantileGeneric<T>::FIVE ( static_cast<T>(5) );

template <class T>
const T math::SampleQuantileGeneric<T>::EIGHT ( static_cast<T>(8) );

template <class T>
const T math::SampleQuantileGeneric<T>::HALF ( ONE / static_cast<T>(2) );

template <class T>
const T math::SampleQuantileGeneric<T>::QUARTER ( ONE / static_cast<T>(4) );


/**
 * Constructor.
 * Creates its own copy of the sample vector.
 *
 * @param sample - a vector of samples
 *
 * @throw StatisticsException if 'sample' is empty or allocation of memory for its copy failed
 */
template <class T>
math::SampleQuantileGeneric<T>::SampleQuantileGeneric(const std::vector<T>& sample) throw(math::StatisticsException)
{
    try
    {
        if ( 0==sample.size() )
        {
            throw math::StatisticsException(math::StatisticsException::SAMPLE_EMPTY);
        }

        // copy the 'sample' into an internal vector:
        math::mtcopy(sample, this->m_v);
        this->m_N = this->m_v.size();

        // and sort it in ascending order
        std::sort(this->m_v.begin(), this->m_v.end());
    }
    catch ( const std::bad_alloc& ex )
    {
        throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
    }
}


/**
 * @return size of the sample
 */
template <class T>
size_t math::SampleQuantileGeneric<T>::sampleSize() const
{
    return this->m_N;
}


/**
 * Estimates a quantile for a probability expressed as a ratio
 * of two integers.
 *
 * Note that 'den' must be greater or equal as 2 and that 'num'
 * must be strictly greater than 0 and strictly less than 'den':
 *   den >= 2
 *   1 <= num <= den-1
 *
 * @param num - probability's numerator
 * @param den - probabilty's denominator
 * @param method - one of the supported methods to estimate the quantile (default: R7)
 *
 * @return quantile for the probabilty 'num'/'den'
 *
 * @throw StatisticsException if any input argument is invalid
 */
template <class T>
T math::SampleQuantileGeneric<T>::quantile(size_t num, size_t den, math::EQntlType::type method) const throw(math::StatisticsException)
{
    if ( 0==num || den<2 || num>=den )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    return this->qntl( static_cast<T>(num) / static_cast<T>(den), method);
}


/*
 * Internal (and private) method that estimates the quantile from the given value 'h'.
 * The function was introduced because several methods actually use the same formula
 * to estimate a quantile from 'h' (calculated differently, depending on a method).
 *
 * @param h - h value, calculated depending on a method
 *
 * @return estimated quantile as a function of 'h'
 */
template <class T>
T math::SampleQuantileGeneric<T>::linIntrp(const T& h) const
{
	// "rename" the vector as it is referred in statistical publications:
    const std::vector<T>& x = this->m_v;
    const T hf = std::floor(h);

    /*
     * if floor(h) is referred as 'hf', the quantile is estimated as:
     *
     *   qp = x[hf] + (h - hf) * (x[hf+1] - x[hf])
     *
     * Note that the expression assumes one-based arrays.
     * In C++ all indices must be additionally decremented by 1
     * prior accessing elements of x.
     *
     * Hence:
     *
     *   qp = x[hf-1] + (h - hf) * (x[hf] - x[hf-1])
     */

    return ( x.at(toInt(hf)-1) + (h - hf) * (x.at(toInt(hf)) - x.at(toInt(hf)-1)) );
}


/**
 * Estimates a quantile for the given probability 'p'.
 *
 * See https://en.wikipedia.org/wiki/Quantile for more
 * details about supported methods.
 *
 * Note that 'p' must be greater than 0 and less than 1:
 *   0 <= p <= 1
 *
 * @param p - probabilty
 * @param method - one of the supported methods to estimate the quantile (default: R7)
 *
 * @return quantile for the probabilty 'p'
 *
 * @throw StatisticsException if any input argument is invalid
 */
template <class T>
T math::SampleQuantileGeneric<T>::qntl(const T& p, math::EQntlType::type method) const throw(math::StatisticsException)
{
    T retVal;

    // Use references to "rename" the variables into names,
    // typically used in statistical publications.
    const size_t& N = this->m_N;
    const std::vector<T>& x = this->m_v;

    // N casted to , used often in nonint expressions
    const T NT = static_cast<T>(N);

    if ( p<math::NumericUtil<T>::ZERO || p>ONE )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }


    /*
     * Implementation of all methods described in:
     * https://en.wikipedia.org/wiki/Quantile
     * Note that the Wikipedia article assumes one-based arrays.
     * In a language with zero-based arrays (e.g. C++), the index must
     * be decremented by 1 before an array element is accessed.
     */

    switch (method)
    {
    case math::EQntlType::R1 :
    case math::EQntlType::SAS3 :
    {
        // Inverse of empirical distribution function
        if ( true == math::NumericUtil<T>::isZero(p) )
        {
            retVal = x.at(0);
        }
        else
        {
            const T h = NT * p + HALF;
            retVal = x.at(toInt(std::ceil(h-HALF)) - 1);
        }

        break;
    }

    case math::EQntlType::R2 :
    case math::EQntlType::SAS5 :
    {
        // Inverse of empirical distribution function
        // with averaging at discontinuities
        if ( true == math::NumericUtil<T>::isZero(p) )
        {
            retVal = x.at(0);
        }
        else if ( true == math::NumericUtil<T>::isZero(ONE-p) )
        {
            retVal = x.at(N - 1);
        }
        else
        {
            const T h = NT * p + HALF;
            retVal = ( x.at(toInt(std::ceil(h-HALF)) - 1) + 
                       x.at(toInt(std::floor(h+HALF)) - 1) ) * HALF;
        }

        break;
    }

    case math::EQntlType::R3 :
    case math::EQntlType::SAS2 :
    {
        // Observation closest to N*p
        if ( p < HALF/NT )
        {
            retVal = x.at(0);
        }
        else
        {
            const T h = NT * p;
            const T frac = h - std::floor(h);
            size_t idx = toInt(std::floor(h));

            /*
             * The index is h, rounded to the nearest integer.
             * In case of a tie, round to the even integer,
             * assuming one-based indexing.
             */

            if ( true == math::NumericUtil<T>::isZero(frac-HALF) )
            {
                // in this case round to the nearest even (divisible by 2) integer idx
                if ( 0 != idx%2 )
                {
                    ++idx;
                }
            }
            else if ( frac > HALF )
            {
                ++idx;
            }
            // no need to to do anything with 'idx' if frac<0.5

            retVal = x.at(idx-1);
        }

        break;
    }

    case math::EQntlType::R4 :
    case math::EQntlType::SAS1 :
    case math::EQntlType::SCIPY_0_1 :
    {
        // Linear interpolation of the empirical distribution function
        if ( p < ONE/NT )
        {
            retVal = x.at(0);
        }
        else if ( true == math::NumericUtil<T>::isZero(ONE-p) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp(NT * p);
        }

        break;
    }

    case math::EQntlType::R5 :
    case math::EQntlType::SCIPY_05_05 :
    {
        // Piecewise linear function where the knots are the values midway
        // through the steps of the empirical distribution function
        if ( p < HALF/NT )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT-HALF)/NT )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp(NT * p + HALF);
        }

        break;
    }

    case math::EQntlType::R6 :
    case math::EQntlType::SAS4 :
    case math::EQntlType::SCIPY_0_0 :
    {
        // Linear interpolation of the expectations for the order statistics
        // for the uniform distribution on [0,1]
        if ( p < ONE/(NT+ONE) )
        {
            retVal = x.at(0);
        }
        else if ( p >= NT/(NT+ONE) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + ONE) * p);
        }

        break;
    }

    case math::EQntlType::R7 :
    case math::EQntlType::SCIPY_1_1 :
    case math::EQntlType::EXCEL :
    {
        // Linear interpolation of the modes for the order statistics for
        // the uniform distribution on [0,1]
        if ( true == math::NumericUtil<T>::isZero(ONE-p) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT - ONE) * p + ONE);
        }

        break;
    }

    case math::EQntlType::R8 :
    case math::EQntlType::SCIPY_13_13 :
    {
        const T third = ONE / THREE;

        // Linear interpolation of the approximate medians for order statistics
        if ( p < (TWO*third) / (NT+third) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT-third) / (NT+third) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + third) * p + third);
        }

        break;
    }

    case math::EQntlType::R9 :
    case math::EQntlType::SCIPY_38_38 :
    {
        const T F_3_8 = THREE / EIGHT;
        const T F_5_8 = FIVE / EIGHT;

        // The resulting quantile estimates are approximately unbiased for the
        // expected order statistics if x is normally distributed
        if ( p < F_5_8/(NT+QUARTER) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT-F_3_8)/(NT+QUARTER) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + QUARTER) * p + F_3_8);
        }

        break;
    }

    case math::EQntlType::SCIPY_N05_N05 :
    {
        // If h were rounded, this would give the order statistic with the least
        // expected square deviation relative to p
        if ( p < (THREE/TWO) / (NT+TWO) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (NT+HALF)/(NT+TWO) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((NT + TWO) * p - HALF);
        }

        break;
    }

    default:
    	throw math::StatisticsException(math::StatisticsException::UNSUPPORTED_QUANTILE_METHOD);
    }  // switch

    return retVal;
}


/**
 * Median of the sample.
 *
 * If the sample's number of elements is odd, the middle element
 * is returned:
 *
 *   median = sorted_vector[(N-1)/2]
 *
 * If the sample's number of elements is even, the mean of the middle
 * two elements is returned:
 *
 *             sorted_vector[N/2 - 1] + sorted_vector[N/2]
 *   median = ---------------------------------------------
 *                                  2
 *
 *
 * @return median of the sample
 */
template <class T>
T math::SampleQuantileGeneric<T>::median() const
{
    T retVal;
    const size_t& N = this->m_N;
    const std::vector<T>& x = this->m_v;

    if ( 0 == N % 2 )
    {
        // even number of elements
    	size_t h = N / 2;
        retVal = (x.at(h-1) + x.at(h)) * HALF;
    }
    else
    {
        // odd number of elements
        retVal = x.at((N-1) / 2);
    }

    return retVal;
}


/**
 * Estimation of a quartile of the sample.
 *
 * @param q - desired quartile (between 0 and 4 )
 * @param method - one of the supported methods to estimate the quartile (default: R7)
 *
 * @return estimation of the desired quartile
 *
 * @throw StatisticsException  if any input argument is invalid
 */
template <class T>
T math::SampleQuantileGeneric<T>::quartile(size_t q, math::EQntlType::type method) const throw(math::StatisticsException)
{
	T retVal;

    if ( 0 == q )
    {
        retVal = this->m_v.at(0);
    }
    else if ( 4 == q )
    {
        retVal = this->m_v.at(this->m_N - 1);
    }
    else if ( q>=1 && q<=3 )
    {
        retVal = this->quantile(q, 4, method);
    }
    else
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    return retVal;
}


/**
 * @param method - one of the supported methods to estimate both quartiles (default: R7)
 *
 * @return sample's interquartile range (difference between the 3rd and the 1st quartile)
 */
template <class T>
T math::SampleQuantileGeneric<T>::iqr(math::EQntlType::type method) const
{
    return ( this->quartile(3, method) - this->quartile(1, method) );
}


/**
 * @return minimum value of the sample
 */
template <class T>
T math::SampleQuantileGeneric<T>::min() const
{
    // Note that this->m_v is already sorted in ascending order
    return this->m_v.at(0);
}


/**
 * @return maximum value of the sample
 */
template <class T>
T math::SampleQuantileGeneric<T>::max() const
{
    // Note that this->m_v is already sorted in ascending order
    return this->m_v.at(this->m_N-1);
}


/**
 * Test if a value is outlier regarding the sample.
 * 
 * An element is an outlier if it lies either below the first quartile - iqrs * IQR
 * or above the third quartile + iqrs * IQR.
 * 
 * The exact value of 'iqrs' can be passed as an argument and is typically
 * equal to 1.5.
 *
 * Note that quartiles are estimated from the original sample that was passed to
 * the constructor. In other words, if 'val' is not a member of the sample,
 * quartiles will not be reestimated to include this sample.
 * 
 * @param val - value to be checked
 * @param iqrs - number of interquartile ranges below/above Q1 and Q2 (default 1.5)
 * @param method - method to estimate quartiles (default R7)
 * 
 * @return a logical value indicating whether 'val' is an outlier or not 
 */
template <class T>
bool math::SampleQuantileGeneric<T>::isOutlier(
                const T& val,
                const T& iqrs, 
                math::EQntlType::type method) const
{
    const T q1 = this->quartile(1, method);
    const T q3 = this->quartile(3, method);
    const T diff = q3 - q1;   // IQR

    return ( ( val < (q1-iqrs * diff) ) ||
             ( val > (q3+iqrs * diff) ) );
}


/**
 * Fills all sample's outliers into the set 'outl'.
 * 
 * An element is an outlier if it lies either below the first quartile - iqrs * IQR
 * or above the third quartile + iqrs * IQR.
 * 
 * The exact value of 'iqrs' can be passed as an argument and is typically
 * equal to 1.5.
 * 
 * @param outl - a reference to a set to be filled with outliers
 * @param iqrs - number of interquartile ranges below/above Q1 and Q2 (default 1.5)
 * @param method - method to estimate quartiles (default R7)
 * 
 * @throw StatisticsException if allocation of memory fails
 */
template <class T>
void math::SampleQuantileGeneric<T>::outliers(
                std::set<T>& outl, 
                const T& iqrs, 
                math::EQntlType::type method) const
            throw (math::StatisticsException)
{
    try
    {
        const T q1 = this->quartile(1, method);
        const T q3 = this->quartile(3, method);
        const T diff = q3 - q1;   // IQR
        const T lowerBound = q1 - iqrs * diff;
        const T upperBound = q3 + iqrs * diff;

        /*
         * Vector this->m_v is already sorted in ascending order. That said,
         * it is sufficient to start iterating from the vector's start as long
         * as the elements are below the lower bound. Then another iteration is
         * performed from the vector's end as long as the elements are above
         * the upper bound. 
         */

        typename std::vector<T>::const_iterator it;
        for ( it=this->m_v.begin(); 
              it!=this->m_v.end() && *it<lowerBound;
              ++it )
        {
            outl.insert(*it);
        }

        typename std::vector<T>::const_reverse_iterator rit;
        for ( rit = this->m_v.rbegin();
              rit!=this->m_v.rend() && *rit>upperBound;
              ++rit )
        {
            outl.insert(*rit);
        }
    }
    catch ( const std::bad_alloc& ba )
    {
        throw math::StatisticsException(math::StatisticsException::OUT_OF_MEMORY);
    }
}
    
    
/**
 * Destructor
 */
template<class T>
math::SampleQuantileGeneric<T>::~SampleQuantileGeneric()
{
    // Vector's destructors would probably clean up this automatically.
    // Anyway let us clear the vector, just to be aware of allocated resources.
    this->m_v.clear();

    // No other resources to release.
}
