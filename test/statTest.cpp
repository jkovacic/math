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
 * and namespaces (SampleStat, NormalDist, StudentDist, ChiSquareDist,
 * FDist, SampleQuantileGeneric)
 */


#include <iostream>
#include <cstddef>
#include <vector>
#include <set>
#include <new>

#include "mtcopy.h"
#include "SampleStatGeneric.h"
#include "SampleQuantileGeneric.h"
#include "NormalDistGeneric.h"
#include "StudentDistGeneric.h"
#include "ChiSquareDistGeneric.h"
#include "FDistGeneric.h"
#include "StatisticsException.h"

using namespace std;
using namespace math;


// Nr. of entries in mtcars:
#define LEN      ( 32 )

// Nr. of probabilities in quantile test:
#define N_PROBS  ( 9 )


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


        /*
           # Equivalent of the following command in R:
           data(mtcars)
         */

        vector<double> vmpgs;
        vector<double> vwts;

        mtcopy(ampgs, LEN, vmpgs);
        mtcopy(awts, LEN, vwts);


        /*
           min(mtcars$wt)
           [1] 1.513
           max(mtcars$wt)
           [1] 5.424
         */

        cout << "min(wt): " << SampleStat::min(vwts) << " (expected: 1.513)" << endl;
        cout << "max(wt): " << SampleStat::max(vwts) << " (expected: 5.424)" << endl;
        cout << endl;


        /*
           sum(mtcars$mpg)
           [1] 642.9
           mean(mtcars$mpg)
           [1] 20.09062
           var(mtcars$mpg)
           [1] 36.3241
           sd(mtcars$mpg)
           [1] 6.026948
           n <- length(mtcars$mpg)
           (n-1)/n * var(mtcars$mpg)
           [1] 35.18897
           n <- length(mtcars$mpg)
           sqrt((n-1)/n) * sd(mtcars$mpg)
           [1] 5.93203
         */

        cout << "Sum of all elements: " << SampleStat::sum(vmpgs) << " (expected: 642.9)" << endl;
        cout << "Sample mean: " << SampleStat::mean(vmpgs) << " (expected: 20.09062)" << endl;
        cout << "Sample variance: " << SampleStat::var(vmpgs) << " (expected: 36.3241)" << endl;
        cout << "Sample standard deviation: " << SampleStat::stdev(vmpgs) << " (expected: 6.026948)" << endl;
        cout << "Population variance (w/o Bessel's correction): " << SampleStat::var(vmpgs, false) << " (expected: 35.18897)" << endl;
        cout << "Population standard deviation (w/o Bessel's correction): " << SampleStat::stdev(vmpgs, false) << " (expected: 5.93203)" << endl;


        /*
           cov(mtcars$mpg, mtcars$wt)
           [1] -5.116685
           (n-1)/n * cov(mtcars$mpg, mtcars$wt)
           [1] -4.956788
           cor(mtcars$mpg, mtcars$wt)
           [1] -0.8676594
           cor(mtcars$mpg, mtcars$wt)^2
           [1] 0.7528328
         */

        cout << "Sample covariance: " << SampleStat::cov(vmpgs, vwts) << " (expected: -5.116685)" << endl;
        cout << "Population covariance (w/o B.c.): " << SampleStat::cov(vmpgs, vwts, false) << " (expected: -4.956788)" << endl;
        cout << "Pearson's r: " << SampleStat::cor(vmpgs, vwts) << " (expected: -0.8676594)" << endl;
        cout << "r^2: " << SampleStat::r2(vmpgs, vwts) << " (expected: 0.7528328)" << endl;


        /*
         * R code to test normal distribution related functions:
         *
           dnorm(4.5, mean=2, sd=3)
           [1] 0.09397063

           (6.7 - 2) / 3
           [1] 1.566667

           2 - 1.3 * 3
           [1] -1.9

           pnorm(1.72, mean=2, sd=3)
           [1] 0.4628194

           pnorm(2.48, mean=2, sd=3, lower.tail=FALSE)
           [1] 0.4364405

           pnorm(3, mean=2, sd=3) - pnorm(1, mean=2, sd=3)
           [1] 0.2611173

           qnorm(0.75, mean=2, sd=3)
           [1] 4.023469

           qnorm(0.52, mean=2, sd=3, lower.tail=FALSE)
           [1] 1.849539
         */

        cout << endl;
        cout << "N(2,3): z for x = 6.7:  " << NormalDist::getZ(6.7, 2.0, 3.0) << " (expected: 1.566667)" << endl;
        cout << "N(2,3): x for z = -1.3: " << NormalDist::getX(-1.3, 2.0, 3.0) << " (expected: -1.9)" << endl;
        cout << "N(2,3) at x=4.5: " << NormalDist::pdf(4.5, 2.0, 3.0) << " (expected: 0.09397063)" << endl;
        cout << "N(2,3): P(X<1.72): " << NormalDist::prob(1.72, 2.0, 3.0) << " (expected: 0.4628194)" << endl;
        cout << "N(2,3): P(X>2.48): " << NormalDist::prob(2.48, 2.0, 3.0, false) << " (expected: 0.4364405)" << endl;
        cout << "N(2,3): P(1<X<3): " << NormalDist::probInt(3.0, 1.0, 2.0, 3.0) << " (expected: 0.2611173)" << endl;
        cout << "N(2,3): q(p>0.75): " << NormalDist::quant(0.75, 2.0, 3.0) << " (expected: 4.023469)" << endl;
        cout << "N(2,3): q(p<0.52): " << NormalDist::quant(0.52, 2.0, 3.0, false) << " (expected: 1.849539)" << endl;
        cout << endl;


        /*
         * R code to test Student's distribution related functions:
         *
           (3-2) / (1.5/sqrt(10))
           [1] 2.108185

           2 - 1.2 * 1.5/sqrt(10)
           [1] 1.43079

           dt(2, df=5)
           [1] 0.06509031

           pt(2, df=12)
           [1] 0.9656725

           pt(1.1, df=12, lower.tail=FALSE)
           [1] 0.1464549

           pt(1, df=12) - pt(-0.5, df=12)
           [1] 0.5184167

           qt(0.75, df=12)
           [1] 0.6954829

           qt(0.52, df=12, lower.tail=FALSE)
           [1] -0.05121096
         */

        cout << "T(n=10, mu=2, s=1.5): t for x=3:    " << StudentDist::getT(3.0, 10, 2.0, 1.5) << " (expected: 2.108185)" << endl;
        cout << "T(n=10, mu=2, s=1.5): x for t=-1.2: " << StudentDist::getX(-1.2, 10, 2.0, 1.5) << " (expected: 1.43079)" << endl;
        cout << "T(df=5):  pdf at x=2:      " << StudentDist::pdf(2.0, 5.0) << " (expected: 0.06509031)" << endl;
        cout << "T(df=12): P(X<2):   " << StudentDist::prob(2.0, 12.0) << " (expected: 0.9656725)" << endl;
        cout << "T(df=12): P(X>1.1): " << StudentDist::prob(1.1, 12.0, false) << " (expected: 0.1464549)" << endl;
        cout << "T(df=12): P(-0.5<X<1): " << StudentDist::probInt(-0.5, 1.0, 12.0) << " (expected: 0.5184167)" << endl;
        cout << "T(df=12): q(p>0.75): " << StudentDist::quant(0.75, 12.0) << " (expected: 0.6954829)" << endl;
        cout << "T(df=12): q(p<0.52): " << StudentDist::quant(0.52, 12.0, false) << " (expected: -0.05121096)" << endl;
        cout << endl;


        /*
         * R code to test chi-squared distribution related functions:
         *
           dchisq(1.2, df=2)
           [1] 0.2744058 

           dchisq(3.1, df=7)
           [1] 0.0955139

           pchisq(2.7, df=1)
           [1] 0.8996518

           pchisq(1.8, df=4)
           [1] 0.2275176

           pchisq(3.4, df=0.3, lower.tail=FALSE)
           [1] 0.01365495

           pchisq(1.7, df=5, lower.tail=FALSE)
           [1] 0.8888998

           pchisq(3, df=1.3) - pchisq(2, df=1.3)
           [1] 0.0975555

           pchisq(3, df=4.2) - pchisq(2, df=4.2)
           [1] 0.1737052

           qchisq(0.25, df=0.75)
           [1] 0.03672361

           qchisq(0.25, df=3.8)
           [1] 1.776557

           qchisq(0.25, df=0.8, lower.tail=FALSE)
           [1] 1.009612

           qchisq(0.25, df=6, lower.tail=FALSE)
           [1] 7.840804
         */

        cout << "ChiSq(df=2) : pdf at x=1.2: " << ChiSquareDist::pdf(1.2, 2.0) << " (expected: 0.2744058)" << endl;
        cout << "ChiSq(df=7) : pdf at x=3.1: " << ChiSquareDist::pdf(3.1, 7.0) << " (expected: 0.0955139)" << endl;
        cout << "ChiSq(df=1): P(X<2.7): " <<  ChiSquareDist::prob(2.7, 1.0) << " (expected: 0.8996518)" << endl;
        cout << "ChiSq(df=4): P(X<1.8): " << ChiSquareDist::prob(1.8, 4.0) << " (expected: 0.2275176)" << endl;
        cout << "ChiSq(df=0.3): P(X>3.4): " << ChiSquareDist::prob(3.4, 0.3, false) << " (expected: 0.01365495)" << endl;
        cout << "ChiSq(df=5):   P(X>1.7): " << ChiSquareDist::prob(1.7, 5.0, false) << " (expected: 0.8888998)" << endl;
        cout << "ChiSq(df=1.3): P(2<X<3): " << ChiSquareDist::probInt(2.0, 3.0, 1.3) << " (expected: 0.0975555)" << endl;
        cout << "ChiSq(df=4.2): P(2<X<3): " << ChiSquareDist::probInt(3.0, 2.0, 4.2) << " (expected: 0.1737052)" << endl;
        cout << "ChiSq(df=0.75): q(p>0.25): " << ChiSquareDist::quant(0.25, 0.75) << " (expected: 0.03672361)" << endl;
        cout << "ChiSq(df=3.8):  q(p>0.25): " << ChiSquareDist::quant(0.25, 3.8) << " (expected: 1.776557)" << endl;
        cout << "ChiSq(df=0.8):  q(p<0.25): " << ChiSquareDist::quant(0.25, 0.8, false) << " (expected: 1.009612)" << endl;
        cout << "ChiSq(df=6):    q(p<0.25): " << ChiSquareDist::quant(0.25, 6.0, false) << " (expected: 7.840804)" << endl;
        cout << endl;


        /*
         * R code to test F-distribution related functions:
         *
           df(2.1, df1=1, df2=3)
           [1] 0.08776311

           df(3.5, df1=4, df2=3)
           [1] 0.05386789

           pf(4, df1=0.7, df2=2.5)
           [1] 0.8499816

           pf(4, df1=2.5, df2=0.7)
           [1] 0.5759108

           pf(3, df1=0.8, df2=3.5, lower.tail=FALSE)
           [1] 0.1645458

           pf(3, df1=3.5, df2=1.5, lower.tail=FALSE)
           [1] 0.3174175

           pf(3, df1=0.5, df2=0.5) - pf(1, df1=0.5, df2=0.5)
           [1] 0.1022432

           pf(3, df1=4, df2=0.2) - pf(1, df1=4, df2=0.2)
           [1] 0.07963281

           qf(0.63, df1=0.7, df2=0.3)
           [1] 45.799

           qf(0.63, df1=5, df2=6)
           [1] 1.313811

           qf(0.72, df1=0.3, df2=0.7, lower.tail=FALSE)
           [1] 0.003393905

           qf(0.72, df1=6, df2=5, lower.tail=FALSE)
           [1] 0.6081648
         */

        cout << "F(1,3): pdf at x=2.1: " << FDist::pdf(2.1, 1.0, 3.0) << " (expected: 0.08776311)" << endl;
        cout << "F(4,3): pdf at x=3.5: " << FDist::pdf(3.5, 4.0, 3.0) << " (expected: 0.05386789)" << endl;
        cout << "F(0.7, 2.5): P(X<4): " << FDist::prob(4.0, 0.7, 2.5) << " (expected: 0.8499816)" << endl;
        cout << "F(2.5, 0.7): P(X<4): " << FDist::prob(4.0, 2.5, 0.7) << " (expected: 0.5759108)" << endl;
        cout << "F(0.8, 3.5): P(X>3): " << FDist::prob(3.0, 0.8, 3.5, false) << " (expected: 0.1645458)" << endl;
        cout << "F(3.5, 1.5): P(X>3): " << FDist::prob(3.0, 3.5, 1.5, false) << " (expected: 0.3174175)" << endl;
        cout << "F(0.5, 0.5): P(1<X<3): " << FDist::probInt(1.0, 3.0, 0.5, 0.5) << " (expected: 0.1022432)" << endl;
        cout << "F(4, 0.2):   P(1<X<3): " << FDist::probInt(1.0, 3.0, 4.0, 0.2) << " (expected: 0.07963281)" << endl;
        cout << "F(0.7, 0.3): q(p>0.63): " << FDist::quant(0.63, 0.7, 0.3) << " (expected: 45.799)" << endl;
        cout << "F(5, 6):     q(p>0.63): " << FDist::quant(0.63, 5.0, 6.0) << " (expected: 1.313811)" << endl;
        cout << "F(0.3, 0.7): q(p<0.72): " << FDist::quant(0.72, 0.3, 0.7, false) << " (expected: 0.003393905)" << endl;
        cout << "F(6, 5):     q(p<0.72): " << FDist::quant(0.72, 6.0, 5.0, false) << " (expected: 0.6081648)" << endl;
        cout << endl;


        /*
         * R code to perform basic unit test of quantiles:

           data(mtcars)
           median(mtcars$mpg)
           [1] 19.2
           quantile(mtcars$mpg, c(0.25, 0.75))
              25%    75%
           15.425 22.800
           IQR(mtcars$mpg)
           [1] 7.375
           quantile(mtcars$mpg, 0.63)
              63%
           21.212
         */

        SampleQuantile q(vmpgs);
        cout << "Median: " << q.median() << " (expected: 19.2)" << endl;
        cout << "1st quartile: " << q.quantile(1, 4) << " (expected: 15.425)" << endl;
        cout << "3rd quartile: " << q.quantile(3, 4) << " (expected: 22.800)" << endl;
        cout << "IQR: " << q.iqr() << " (expected: 7.375)" << endl;
        cout << "63th percentile: " << q.qntl(0.63) << " (expected: 21.212)" << endl;


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


        /*
         * R code to test various probabilities and various methods:
         *
             data(mtcars);   x <- mtcars$mpg
             types <- 1:9
             p <- c(0.01, 0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 0.9, 0.99)
             for (t in types)
             {
               print(quantile(x, p, type=t, names=FALSE))
             }
         */

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


        /*
         * R code to find outliers:
         *
             q <- quantile(mtcars$mpg, c(0.25, 0.75), name=FALSE)
             d <- IQR(mtcars$mpg)
             mtcars$mpg < (q[1]-1.0*d) | mtcars$mpg > (q[2]+1.0*d)
              [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
             [14] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
             [27] FALSE  TRUE FALSE FALSE FALSE FALSE

             sort(unique(mtcars$mpg[mtcars$mpg<(q[1]-0.5*d) | mtcars$mpg>(q[2]+0.5*d)]))
             [1] 10.4 27.3 30.4 32.4 33.9
         */

        cout << "Min mpg: " << q.min() << " (expected 10.4)" << endl;
        cout << "Max mpg: " << q.max() << " (expected 33.9)" << endl;
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
