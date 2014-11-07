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

#include "exception/StatisticsException.hpp"
#include "util/NumericUtil.hpp"


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
        this->m_v = sample;
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
size_t math::SampleQuantileGeneric<T>::sampleSize()
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
 * @param method - one of supported methods to estimate the quantile
 *
 * @return quantile for the probabilty 'num'/'den'
 *
 * @throw StatisticsException if any input argument is invalid
 */
template <class T>
T math::SampleQuantileGeneric<T>::quantile(size_t num, size_t den, quant_types method) throw(math::StatisticsException)
{
    if ( 0==num || den<2 || num>=den )
    {
        throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
    }

    return this->qntl( static_cast<double>(num) / static_cast<double>(den), method);
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
T math::SampleQuantileGeneric<T>::linIntrp(double h)
{
	// "rename" the vector as it is referred in statistical publications:
    const std::vector<T>& x = this->m_v;
    const double hf = std::floor(h);

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

    return ( x.at(dbl2int(hf)-1) + (h - hf) * (x.at(dbl2int(hf)) - x.at(dbl2int(hf)-1)) );
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
 * @param method - one of the supported methods to estimate the quantile
 *
 * @return quantile for the probabilty 'p'
 *
 * @throw StatisticsException if any input argument is invalid
 */
template <class T>
T math::SampleQuantileGeneric<T>::qntl(double p, quant_types method) throw(math::StatisticsException)
{
    T retVal;

    // Use references to "rename" the variables into names,
    // typically used in statistical publications.
    const size_t& N = this->m_N;
    const std::vector<T>& x = this->m_v;

    if ( p<0.0 || p>1.0 )
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
    case R1 :
    case SAS3 :
    {
        // Inverse of empirical distribution function
        if ( true == math::NumericUtil<double>::isZero(p) )
        {
            retVal = x.at(0);
        }
        else
        {
            const double h = N * p + 0.5;
            retVal = x.at(dbl2int(std::ceil(h-0.5)) - 1);
        }

        break;
    }

    case R2 :
    case SAS5 :
    {
        // Inverse of empirical distribution function
        // with averaging at discontinuities
        if ( true == math::NumericUtil<double>::isZero(p) )
        {
            retVal = x.at(0);
        }
        else if ( true == math::NumericUtil<double>::isZero(1.0-p) )
        {
            retVal = x.at(N - 1);
        }
        else
        {
            const double h = N * p + 0.5;
            retVal = ( x.at(dbl2int(std::ceil(h-0.5)) - 1) + x.at(dbl2int(std::floor(h+0.5)) - 1) )
                / static_cast<T>(2);

        }

        break;
    }

    case R3 :
    case SAS2 :
    {
        // Observation closest to N*p
        if ( p < 0.5/N )
        {
            retVal = x.at(0);
        }
        else
        {
            const double h = N * p;
            const double frac = h - std::floor(h);
            size_t idx = dbl2int(std::floor(h));

            /*
             * The index is h, rounded to the nearest integer.
             * In case of a tie, round to the even integer,
             * assuming one-based indexing.
             */

            if ( true == math::NumericUtil<double>::isZero(frac-0.5) )
            {
                // in this case round to the nearest even (divisible by 2) integer idx
                if ( 0 != idx%2 )
                {
                    ++idx;
                }
            }
            else if ( frac > 0.5 )
            {
                ++idx;
            }
            // no need to to do anything with 'idx' if frac<0.5

            retVal = x.at(idx-1);
        }

        break;
    }

    case R4 :
    case SAS1 :
    case SCIPY_0_1 :
    {
        // Linear interpolation of the empirical distribution function
        if ( p < 1.0/N )
        {
            retVal = x.at(0);
        }
        else if ( true == math::NumericUtil<double>::isZero(1.0-p) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp(N * p);
        }

        break;
    }

    case R5 :
    case SCIPY_05_05 :
    {
        // Piecewise linear function where the knots are the values midway
        // through the steps of the empirical distribution function
        if ( p < 0.5/N )
        {
            retVal = x.at(0);
        }
        else if ( p >= (N-0.5)/N )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp(N * p + 0.5);
        }

        break;
    }

    case R6 :
    case SAS4:
    case SCIPY_0_0 :
    {
        // Linear interpolation of the expectations for the order statistics
        // for the uniform distribution on [0,1]
        if ( p < 1.0/(N+1.0) )
        {
            retVal = x.at(0);
        }
        else if ( p >= N/(N+1.0) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((N + 1.0) * p);
        }

        break;
    }

    case R7 :
    case SCIPY_1_1:
    case EXCEL :
    {
        // Linear interpolation of the modes for the order statistics for
        // the uniform distribution on [0,1]
        if ( true == math::NumericUtil<T>::isZero(1.0-p) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((N - 1.0) * p + 1.0);
        }

        break;
    }

    case R8 :
    case SCIPY_13_13 :
    {
        // Linear interpolation of the approximate medians for order statistics
        if ( p < (2.0/3.0)/(N+1.0/3.0) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (N-1.0/3.0)/(N+1.0/3.0) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((N + 1.0/3.0) * p + 1.0/3.0);
        }

        break;
    }

    case R9:
    case SCIPY_38_38 :
    {
        // The resulting quantile estimates are approximately unbiased for the
        // expected order statistics if x is normally distributed
        if ( p < 0.625/(N+0.25) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (N-0.375)/(N+0.25) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((N + 0.25) * p + 0.375);
        }

        break;
    }

    case SCIPY_N05_N05 :
    {
        // If h were rounded, this would give the order statistic with the least
        // expected square deviation relative to p
        if ( p < 1.5/(N+2.0) )
        {
            retVal = x.at(0);
        }
        else if ( p >= (N+0.5)/(N+2.0) )
        {
            retVal = x.at(N-1);
        }
        else
        {
            retVal = this->linIntrp((N + 2.0) * p - 0.25);
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
 *   median = sorted_vecor[(N-1)/2]
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
T math::SampleQuantileGeneric<T>::median()
{
    T retVal;
    const size_t& N = this->m_N;
    const std::vector<T>& x = this->m_v;

    if ( 0 == N % 2 )
    {
        // even number of elements
    	size_t h = N / 2;
        retVal = (x.at(h-1) + x.at(h)) / static_cast<T>(2);
    }
    else
    {
        // odd number of elements
        retVal = x.at((N-1) / 2);
    }

    return retVal;
}


/**
 * @return sample's interquartile range (difference between the 3rd and the 1st quartile)
 */
template <class T>
T math::SampleQuantileGeneric<T>::iqr(quant_types method)
{
    return this->qntl(0.75, method) - this->qntl(0.25, method);
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
