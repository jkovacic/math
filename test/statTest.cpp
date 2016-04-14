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
 * A test module to test functionality in statistics related classes 
 * and namespaces (SampleStat, SampleQuantileGeneric)
 */


/*
 * Note: results are reproduced in 'scripts/test/stat.R'.
 */


#include <iostream>
#include <cstddef>
#include <vector>
#include <set>
#include <new>

#include "mtcopy.h"
#include "SampleStat.h"
#include "SampleQuantileSortedArray.h"
#include "SampleQuantileSelection.h"
#include "StatisticsException.h"

using namespace std;
using namespace math;


// Nr. of entries in mtcars:
#define LEN      ( 32 )

// Nr. of probabilities in quantile test:
#define N_PROBS  ( 9 )


/*
 * An auxiliary function that tests functionality of classes derived
 * from SampleQuantileGenericAb. Those classes obtain quantiles and
 * (check for) outliers.
 */
static void testQuantilesOutliers(
       const SampleQuantileGenericAb<double>& q,
       const vector<double>& vmpgs )
{
    // Quantiles:
    cout << endl;
    cout << "Median: " << q.median() << " (expected: 19.2)" << endl;
    cout << "Approx. median: " << q.median(true) << " (expected: 19.2)" << endl;
    cout << "1st quartile: " << q.quantile(1, 4) << " (expected: 15.425)" << endl;
    cout << "3rd quartile: " << q.quantile(3, 4) << " (expected: 22.800)" << endl;
    cout << "IQR: " << q.iqr() << " (expected: 7.375)" << endl;
    cout << "63th percentile: " << q.qntl(0.63) << " (expected: 21.212)" << endl;
    cout << "ecdf(25):   " << q.ecdf(25.0) << " (expected: 0.8125)" << endl;
    cout << "ecdf(14.5): " << q.ecdf(14.5) << " (expected: 0.125)" << endl;

    const double probs[] = { 0.01, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 0.99 };
    const char* exp[ N_PROBS ] =
    {
        "10.4\t14.3\t15.2\t17.3\t19.2\t21.0\t22.8\t30.4\t33.9",
        "10.40\t14.30\t15.35\t17.55\t19.20\t21.20\t22.80\t30.40\t33.90",
        "10.4\t13.3\t15.2\t17.3\t19.2\t21.0\t22.8\t30.4\t33.9",
        "10.40\t13.50\t15.20\t17.30\t19.20\t21.00\t22.80\t29.78\t33.42",
        "10.40\t14.00\t15.35\t17.55\t19.20\t21.20\t22.80\t30.40\t33.90",
        "10.4000\t13.6000\t15.2750\t17.4875\t19.2000\t21.2500\t22.8000\t30.4000\t33.9000",
        "10.4000\t14.3400\t15.4250\t17.6125\t19.2000\t21.1500\t22.8000\t30.0900\t33.4350",
        "10.40\t13.8667\t15.325\t17.5292\t19.200\t21.2167\t22.80\t30.40\t33.90",
        "10.40\t13.90\t15.3312\t17.5344\t19.20\t21.2125\t22.80\t30.40\t33.90"
    };


    // Quantile methods:
    cout << endl << "Test of various quantile methods:" << endl;

    // Not really the best practice, but as long as the enum is contiguous...
    for ( int type=EQntlType::R1; type<=EQntlType::R9; ++type )
    {
        cout << "R" << 1+type << ":\t";
        for ( size_t i=0; i<N_PROBS; ++i )
        {
            cout << q.qntl(probs[i], static_cast<EQntlType::type>(type));
            if ( i < N_PROBS-1 )
            {
                cout << "\t";
            }
        }
        cout << endl;
        cout << "Exp.:\t" << exp[type] << endl;
    }
    cout << endl;


    // n.th largest/smallest element:
    cout << "Min mpg: " << q.min() << " (expected 10.4)" << endl;
    cout << "Max mpg: " << q.max() << " (expected 33.9)" << endl;
    cout << "5th smallest mpg: " << q.elem(5-1, false) << " (expected: 14.7)" << endl;
    cout << "5th largest mpg:  " << q.elem(5-1) << " (expected: 27.3)" << endl;
    cout << "3rd smallest mpg: " << q.elem(3, false, false) << " (expected: 13.3)" << endl;
    cout << "6th largest mpg:  " << q.elem(6, true, false) << " (expected: 26)" << endl;
    cout << endl;


    // Outliers:
    typename vector<double>::const_iterator mpgit;
    for ( mpgit=vmpgs.begin(); mpgit!=vmpgs.end(); ++ mpgit )
    {
        cout << *mpgit << " in range [8.05, 30.175]: ";
        cout << q.isOutlier(*mpgit, 1.0) << endl;
    }
    cout << "Outliers for iqr=0.5: [";
    set<double> oul;
    q.outliers(oul, 0.5);
    typename set<double>::const_iterator oit;
    for ( oit=oul.begin(); oit!=oul.end(); ++oit )
    {
        cout << *oit << " ";
    }
    cout << "] expected [10.4 27.3 30.4 32.4 33.9]" << endl;
}

/*
 * Test of classes that perform statistical operations
 */
void statisticsTest()
{
    try
    {
        /* Cars' mpg (miles per gallon) from R's data frame 'mtcars' */
        const double ampgs[ LEN ] =
            { 21.0, 21.0, 22.8, 21.4, 18.7, 18.1, 14.3, 24.4,
              22.8, 19.2, 17.8, 16.4, 17.3, 15.2, 10.4, 10.4,
              14.7, 32.4, 30.4, 33.9, 21.5, 15.5, 15.2, 13.3,
              19.2, 27.3, 26.0, 30.4, 15.8, 19.7, 15.0, 21.4 };

        /* Cars' wt (weight in 1000 lbs) from R's data frame 'mtcars' */
        const double awts[ LEN ] =
            { 2.620, 2.875, 2.320, 3.215, 3.440, 3.460, 3.570, 3.190,
              3.150, 3.440, 3.440, 4.070, 3.730, 3.780, 5.250, 5.424,
              5.345, 2.200, 1.615, 1.835, 2.465, 3.520, 3.435, 3.840,
              3.845, 1.935, 2.140, 1.513, 3.170, 2.770, 3.570, 2.780 };


        vector<double> vmpgs;
        vector<double> vwts;

        mtcopy(ampgs, LEN, vmpgs);
        mtcopy(awts, LEN, vwts);


        cout << "Sum of all elements: " << SampleStat::sum(vmpgs) << " (expected: 642.9)" << endl;
        cout << "Sample mean: " << SampleStat::mean(vmpgs) << " (expected: 20.09062)" << endl;
        cout << "Sample variance: " << SampleStat::var(vmpgs) << " (expected: 36.3241)" << endl;
        cout << "Sample standard deviation: " << SampleStat::stdev(vmpgs) << " (expected: 6.026948)" << endl;
        cout << "Population variance (w/o Bessel's correction): " << SampleStat::var(vmpgs, false) << " (expected: 35.18897)" << endl;
        cout << "Population standard deviation (w/o Bessel's correction): " << SampleStat::stdev(vmpgs, false) << " (expected: 5.93203)" << endl;
        cout << "4th central moment about the mean: " << SampleStat::centralMoment(vmpgs, 4) << " (expected: 3466.479)" << endl;
        cout << "4th moment about the origin: " << SampleStat::moment(vmpgs, 4) << " (expected: 262350.3)" << endl;
        cout << "Sample skewness: " << SampleStat::skewness(vmpgs) << " (expected: 0.6404399)" << endl;
        cout << "Sample skewness (with sample sd): " << SampleStat::skewness(vmpgs, true) << " (expected: 0.610655)" << endl;
        cout << "Sample kurtosis: " << SampleStat::kurtosis(vmpgs) << " (expected: 2.799467)" << endl;
        cout << "Sample excess kurtosis: " << SampleStat::kurtosis(vmpgs, true) << " (expected: -0.2005332)" << endl;
        cout << "Sample ecdf(<=26): " << SampleStat::ecdf(vmpgs, 26.0) << " (expected: 0.84375)" << endl;
        cout << "Sample ecdf(<26):  " << SampleStat::ecdf(vmpgs, 26.0, false) << " (expected: 0.8125)" << endl;


        cout << "Sample covariance: " << SampleStat::cov(vmpgs, vwts) << " (expected: -5.116685)" << endl;
        cout << "Population covariance (w/o B.c.): " << SampleStat::cov(vmpgs, vwts, false) << " (expected: -4.956788)" << endl;
        cout << "Pearson's r: " << SampleStat::cor(vmpgs, vwts) << " (expected: -0.8676594)" << endl;
        cout << "r^2: " << SampleStat::r2(vmpgs, vwts) << " (expected: 0.7528328)" << endl;

        // Classes that obtain quantiles and outliers are tested by an auxiliary
        // function that accepts all classes derived from SampleQuantileGenericAb.

        cout << endl << "= Test of quantile algorithms based on sorted arrays: =" << endl;
        SampleQuantileSortedArray qsa(vmpgs);
        testQuantilesOutliers(qsa, vmpgs);

        cout << endl << "= Test of quantile algorithms based on selection (copy of the original array): =" << endl;
        SampleQuantileSortedArray qselc(vmpgs);
        testQuantilesOutliers(qselc, vmpgs);

        cout << endl << "= Test of quantile algorithms based on selection (the original array): =" << endl;
        SampleQuantileSortedArray qselorig(vmpgs);
        testQuantilesOutliers(qselorig, vmpgs);
    }
    catch ( const StatisticsException& ex )
    {
        cerr << "Statistics exception caught: ";
        ex.what();
        cerr << endl;
    }
    catch ( const bad_alloc& ba )
    {
        cerr << "Could not allocate memory." << endl;
    }
    catch (...)
    {
        cerr << "Other exception caught.";
    }
}
