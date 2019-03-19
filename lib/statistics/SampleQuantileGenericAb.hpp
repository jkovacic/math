/*
Copyright 2016, Jernej Kovacic

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
 * Declaration and implmentation of the class SampleQuantileGenericAb,
 * a superclass of all classes that estimate quantiles of a sample.
 */


#ifndef _MATH_SAMPLEQUANTILEGENERICAB_HPP_
#define _MATH_SAMPLEQUANTILEGENERICAB_HPP_


#include <cstddef>
#include <vector>
#include <cmath>
#include <set>

#include "int_util/IntUtilGeneric.hpp"
#include "util/NumericUtil.hpp"
#include "../settings/stat_settings.h"

#include "exception/StatisticsException.hpp"


namespace math
{

/**
 * @brief Supported methods to estimate a sample quantile.
 */
struct EQntlType
{
    enum type
    {
        R1,                       /// R quantile, type 1
        R2,                       /// R quantile, type 2
        R3,                       /// R quantile, type 3
        R4,                       /// R quantile, type 4
        R5,                       /// R quantile, type 5
        R6,                       /// R quantile, type 6
        R7,                       /// R quantile, type 7
        R8,                       /// R quantile, type 8
        R9,                       /// R quantile, type 9
        SAS1,                     /// SAS QNTL, method 1
        SAS2,                     /// SAS QNTL, method 2
        SAS3,                     /// SAS QNTL, method 3
        SAS4,                     /// SAS QNTL, method 4
        SAS5,                     /// SAS QNTL, method 5
        EXCEL,                    /// MS Excel, function PERCENTILE
        SCIPY_0_0,                /// SciPy scipy.stats.mstats.mquantiles, alphap=0, betap=0
        SCIPY_0_1,                /// SciPy scipy.stats.mstats.mquantiles, alphap=0, betap=1
        SCIPY_05_05,              /// SciPy scipy.stats.mstats.mquantiles, alphap=0.5, betap=0.5
        SCIPY_1_1,                /// SciPy scipy.stats.mstats.mquantiles, alphap=1, betap=1
        SCIPY_13_13,              /// SciPy scipy.stats.mstats.mquantiles, alphap=1/3, betap=1/3
        SCIPY_38_38,              /// SciPy scipy.stats.mstats.mquantiles, alphap=3/8, betap=3/8
        SCIPY_N05_N05,            /// SciPy scipy.stats.mstats.mquantiles, alphap=-1/2, betap=-1/2
    };
};


/**
 * @brief A base abstract class for all classes that estimate sample's quantile
 * for any probability within a valid range. It supports several estimation
 * methods, however, it relies on several pure virtual methods that are implemented
 * in derived classes.
 *
 * As this class is abstract, it cannot be instantiated directly.
 */
template <typename F>
class SampleQuantileGenericAb
{

protected:
    const std::size_t m_N;

protected:

   /*
    * Constructor.
    * It assigns 'm_N' the sample size.
    * Since this a constructor of an abstract class it can only be
    * called by derived classes' constructors.
    *
    * @param N - size of the sample
    */
   SampleQuantileGenericAb(const std::size_t N) : m_N(N)
   {
       // nothing else to do
   }


private:

    /*
     * Internal (and private) method that estimates the quantile from the given value 'h'.
     * The function was introduced because several methods actually use the same formula
     * to estimate a quantile from 'h' (calculated differently, depending on a method).
     *
     * @param h - h value, calculated depending on a method
     *
     * @return estimated quantile as a function of 'h'
     */
    F __linIntrp(const F& h) const
    {
        const F hf = std::floor(h);

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

        const std::size_t ihf = math::NumericUtil::intRound<F, std::size_t>(hf);
        F xhf, xhf_1;
        this->_select2(ihf, ihf-1, xhf, xhf_1 );

        return ( xhf_1 + (h - hf) * (xhf - xhf_1));
    }


private:
    // pure virtual methods that select (depending on implementation) n.th smallest element(s)
    virtual F _select(const std::size_t n) const = 0;

    virtual void _select2(
            const std::size_t n1,
            const std::size_t n2,
            F& val1,
            F& val2
          ) const = 0;

protected:

    /*
     * Obtains the range that is used to detect outliers. Observation 'x' is
     * an outlier if:
     *   x < (q1 - iqrs * IQR)  OR  x > (q2 + iqrs * IQR)
     * 
     * @param lower - reference to the variable to assign the lower boundary of the range
     * @param upper - reference to the variable to assign the lower boundary of the range
     * @param iqrs - number of interquartile ranges to subtract from q1 and add to q3
     * @param method - method to estimate quartiles
     */
    void _outlierBounds(
            F& lower,
            F& upper,
            const F& iqrs,
            const EQntlType::type method
          ) const
    {
        const F q1 = this->quartile(1, method);
        const F q3 = this->quartile(3, method);
        const F IQR = q3 - q1;

        lower = q1 - iqrs * IQR;
        upper = q3 + iqrs * IQR;
    }


public:

    /**
     * @return size of the sample
     */
    std::size_t sampleSize() const
    {
        return this->m_N;
    }



    /**
     * Estimates a quantile for a probability expressed as a ratio
     * of two integers.
     *
     * Note that 'den' must be greater than or equal to 2 and that 'num'
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
    template <typename I>
    F quantile(
           const I& num,
           const I& den, 
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG 
         ) const
    {
        if ( true == math::IntUtil::isNegative<I>(num) ||
             true == math::IntUtil::isNegative<I>(den) ||
             static_cast<I>(0) == num || 
             den < static_cast<I>(2) || num >= den )
        {
            throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
        }

        return this->qntl( static_cast<F>(num) / static_cast<F>(den), method);
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
    F qntl(
           const F& p,
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG 
         ) const
    {
        F retVal;

        // Use references to "rename" the variables into names,
        // typically used in statistical publications.
        const std::size_t& N = this->m_N;

        // N casted to T, used often in nonint expressions
        const F NT = static_cast<F>(N);

        if ( p<static_cast<F>(0) || p>static_cast<F>(1) )
        {
            throw math::StatisticsException(math::StatisticsException::INVALID_PROBABILTY);
        }

        // Frequently used numerical constants:
        const F ONE = static_cast<F>(1);
        const F TWO = static_cast<F>(2);
        const F HALF = ONE / TWO;
        const F QUARTER = ONE / static_cast<F>(4);

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
            if ( true == math::NumericUtil::isZero<F>(p) )
            {
                retVal = this->_select(0);
            }
            else
            {
                const F h = NT * p + HALF;
                retVal = this->_select(math::NumericUtil::intCeil<F, std::size_t>(h-HALF) - 1);
            }

            break;
        }

        case math::EQntlType::R2 :
        case math::EQntlType::SAS5 :
        {
            // Inverse of empirical distribution function
            // with averaging at discontinuities
            if ( true == math::NumericUtil::isZero<F>(p) )
            {
                retVal = this->_select(0);
            }
            else if ( true == math::NumericUtil::isZero<F>(ONE-p) )
            {
                retVal = this->_select(N - 1);
            }
            else
            {
                const F h = NT * p + HALF;
                F xhn, xhp;
                this->_select2(
                        math::NumericUtil::intCeil<F, std::size_t>(h-HALF) - 1,
                        math::NumericUtil::intFloor<F, std::size_t>(h+HALF) - 1,
                        xhn, xhp );
                retVal = ( xhn + xhp ) * HALF;
            }

            break;
        }

        case math::EQntlType::R3 :
        case math::EQntlType::SAS2 :
        {
            // Observation closest to N*p
            if ( p < HALF/NT )
            {
                retVal = this->_select(0);
            }
            else
            {
                const F h = NT * p;
                const F frac = h - std::floor(h);
                std::size_t idx = math::NumericUtil::intFloor<F, std::size_t>(h);

                /*
                 * The index is h, rounded to the nearest integer.
                 * In case of a tie, round to the even integer,
                 * assuming one-based indexing.
                 */

                if ( true == math::NumericUtil::isZero<F>(frac-HALF) )
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

                retVal = this->_select(idx-1);
            }

            break;
        }

        case math::EQntlType::R4 :
        case math::EQntlType::SAS1 :
        case math::EQntlType::SCIPY_0_1 :
        {
            // Linear interpolation of the empirical distribution function
            if ( p < static_cast<F>(1)/NT )
            {
                retVal = this->_select(0);
            }
            else if ( true == math::NumericUtil::isZero<F>(ONE-p) )
            {
                retVal = this->_select(N-1);
            }
            else
            {
                retVal = this->__linIntrp(NT * p);
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
                retVal = this->_select(0);
            }
            else if ( p >= (NT-HALF)/NT )
            {
                retVal = this->_select(N-1);
            }
            else
            {
                retVal = this->__linIntrp(NT * p + HALF);
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
                retVal = this->_select(0);
            }
            else if ( p >= NT/(NT+ONE) )
            {
                retVal = this->_select(N-1);
            }
            else
            {
                retVal = this->__linIntrp((NT + ONE) * p);
            }

            break;
        }

        case math::EQntlType::R7 :
        case math::EQntlType::SCIPY_1_1 :
        case math::EQntlType::EXCEL :
        {
            // Linear interpolation of the modes for the order statistics for
            // the uniform distribution on [0,1]
            if ( true == math::NumericUtil::isZero<F>(ONE-p) )
            {
                retVal = this->_select(N-1);
            }
            else
            {
                retVal = this->__linIntrp((NT - ONE) * p + ONE);
            }

            break;
        }

        case math::EQntlType::R8 :
        case math::EQntlType::SCIPY_13_13 :
        {
            const F third = ONE / static_cast<F>(3);

            // Linear interpolation of the approximate medians for order statistics
            if ( p < (TWO * third) / (NT+third) )
            {
                 retVal = this->_select(0);
            }
            else if ( p >= (NT-third) / (NT+third) )
            {
                 retVal = this->_select(N-1);
            }
            else
            {
                 retVal = this->__linIntrp((NT + third) * p + third);
            }

            break;
        }

        case math::EQntlType::R9 :
        case math::EQntlType::SCIPY_38_38 :
        {
            const F F_3_8 = static_cast<F>(3) / static_cast<F>(8);
            const F F_5_8 = static_cast<F>(5) / static_cast<F>(8);

            // The resulting quantile estimates are approximately unbiased for the
            // expected order statistics if x is normally distributed
            if ( p < F_5_8/(NT+QUARTER) )
            {
                retVal = this->_select(0);
            }
            else if ( p >= (NT-F_3_8)/(NT+QUARTER) )
            {
                retVal = this->_select(N-1);
            }
            else
            {
                retVal = this->__linIntrp((NT + QUARTER) * p + F_3_8);
            }

            break;
        }

        case math::EQntlType::SCIPY_N05_N05 :
        {
            // If h were rounded, this would give the order statistic with the least
            // expected square deviation relative to p
            if ( p < (static_cast<F>(3) / TWO) / (NT + TWO) )
            {
                retVal = this->_select(0);
            }
            else if ( p >= (NT + HALF) / (NT + TWO) )
            {
                retVal = this->_select(N-1);
            }
            else
            {
                retVal = this->__linIntrp((NT + TWO) * p - HALF);
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
     * If the number of observations is odd, the middle element
     * is returned:
     *
     *   median = sorted_vector[(N-1)/2]
     *
     * If the number of observations is even and 'approx' equals FALSE, the
     * exact median, i.e. the mean of the middle two elements, is returned:
     *
     *             sorted_vector[N/2 - 1] + sorted_vector[N/2]
     *   median = ---------------------------------------------
     *                                  2
     *
     * If the number of observations is even and 'approx' equals TRUE, the
     * approximate median is returned:
     *
     *   median = sorted_vector[N/2]
     *
     * @param approx - what to return in case of even number of observations, see above (default: FALSE)
     *
     * @return median of the sample
     */
    F median(const bool approx = false) const
    {
        F retVal;
        const std::size_t& N = this->m_N;

        if ( 0 == N % 2 )
        {
            // even number of observations
            if ( true == approx )
            {
                /*
                 * A quote from the Numerical Recipes, 3rd Edition, page 432
                 * (http://numerical.recipes/):
                 *
                 * "
                 *   For N > 100 we usually use k = N/2 as the median
                 *   element, formalists be damned.
                 * "
                 */

                retVal = this->_select(N/2);
            }
            else
            {
                std::size_t h = N / 2;
                F elem1, elem2;
                this->_select2(h-1, h, elem1, elem2);
                retVal = (elem1 + elem2) / static_cast<F>(2);
            }
        }
        else
        {
            // odd number of observations
            retVal = this->_select((N-1) / 2);
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
    F quartile(
           const unsigned short int q,
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG 
         ) const
    {
        F retVal;

        if ( 0 == q )
        {
            retVal = this->_select(0);
        }
        else if ( 4 == q )
        {
            retVal = this->_select(this->m_N - 1);
        }
        else if ( q>=1 && q<=3 )
        {
            retVal = this->quantile<short int>(q, 4, method);
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
    F iqr(const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG) const
    {
        return ( this->quartile(3, method) - this->quartile(1, method) );
    }



    /**
     * Test if an observation is an outlier regarding the sample.
     *
     * An observation is an outlier if it lies either below the first quartile - iqrs * IQR
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
     * @param iqrs - number of interquartile ranges below/above Q1 and Q3 (default 1.5)
     * @param method - method to estimate quartiles (default R7)
     *
     * @return a logical value indicating whether 'val' is an outlier or not
     */
    bool isOutlier(
               const F& val,
               const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
               const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
            )  const
    {
        F lower, upper;

        this->_outlierBounds(lower, upper, iqrs, method);

        return ( ( val < lower ) || ( val > upper ) );
    }


    // The following methods are pure virtual, i.e. implemented in derived classes
    virtual F ecdf(const F& t) const = 0;

    virtual F min() const = 0;

    virtual F max() const = 0;

    virtual F elem(
           const std::size_t n, 
           const bool largest = true,
           const bool zerobase = STAT_DEFAULT_ZERO_BASE
         ) const = 0;

    virtual void outliers(
           std::set<F>& outl,
           const F& iqrs = static_cast<F>(STAT_OUTLIER_IQRS_NUM) / static_cast<F>(STAT_OUTLIER_IQRS_DEN),
           const EQntlType::type method = STAT_DEFAULT_QUANTILE_ALG
         ) const = 0;

    /**
     * Destructor
     */
    virtual ~SampleQuantileGenericAb() = 0;


};  // class SampleQuantileGenericAb



/*
 * Pure virtual destructor must be "implemented" separately
 */
template <class F>
SampleQuantileGenericAb<F>::~SampleQuantileGenericAb()
{
    // Actually nothing to do
}


}  // namespace math


#endif  // _MATH_SAMPLEQUANTILEGENERICAB_HPP_
