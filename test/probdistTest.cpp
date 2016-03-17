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
 * A test module to test functionality in namespaces that handle various
 * probabilty distributions (NormalDist, StudentDist, ChiSquareDist,
 * FDist, ContUniformDist, BinomDist, PoissonDist)
 */


/*
 * Note: results are reproduced in 'scripts/test/probdist.R'.
 */


#include <iostream>

#include "NormalDist.h"
#include "StudentDist.h"
#include "ChiSquareDist.h"
#include "FDist.h"
#include "ContUniformDist.h"
#include "BinomDist.h"
#include "PoissonDist.h"

#include "StatisticsException.h"


using namespace std;
using namespace math;


/*
 * Test of functions that handle probability distributions
 */

void probDistributionTest()
{
    try
    {
        // Normal distribution:
        cout << "N(2,3): z for x = 6.7:  " << NormalDist::getZ(6.7, 2.0, 3.0) << " (expected: 1.566667)" << endl;
        cout << "N(2,3): x for z = -1.3: " << NormalDist::getX(-1.3, 2.0, 3.0) << " (expected: -1.9)" << endl;
        cout << "N(2,3) at x=4.5: " << NormalDist::pdf(4.5, 2.0, 3.0) << " (expected: 0.09397063)" << endl;
        cout << "N(2,3): P(X<1.72): " << NormalDist::prob(1.72, 2.0, 3.0) << " (expected: 0.4628194)" << endl;
        cout << "N(2,3): P(X>2.48): " << NormalDist::prob(2.48, 2.0, 3.0, false) << " (expected: 0.4364405)" << endl;
        cout << "N(2,3): P(1<X<3): " << NormalDist::probInt(3.0, 1.0, 2.0, 3.0) << " (expected: 0.2611173)" << endl;
        cout << "N(2,3): q(p>0.75): " << NormalDist::quant(0.75, 2.0, 3.0) << " (expected: 4.023469)" << endl;
        cout << "N(2,3): q(p<0.52): " << NormalDist::quant(0.52, 2.0, 3.0, false) << " (expected: 1.849539)" << endl;
        cout << endl;


        // Student's distribution:
        cout << "T(n=10, mu=2, s=1.5): t for x=3:    " << StudentDist::getT(3.0, 10, 2.0, 1.5) << " (expected: 2.108185)" << endl;
        cout << "T(n=10, mu=2, s=1.5): x for t=-1.2: " << StudentDist::getX(-1.2, 10, 2.0, 1.5) << " (expected: 1.43079)" << endl;
        cout << "T(df=5):  pdf at x=2:      " << StudentDist::pdf(2.0, 5.0) << " (expected: 0.06509031)" << endl;
        cout << "T(df=12): P(X<2):   " << StudentDist::prob(2.0, 12.0) << " (expected: 0.9656725)" << endl;
        cout << "T(df=12): P(X>1.1): " << StudentDist::prob(1.1, 12.0, false) << " (expected: 0.1464549)" << endl;
        cout << "T(df=12): P(-0.5<X<1): " << StudentDist::probInt(-0.5, 1.0, 12.0) << " (expected: 0.5184167)" << endl;
        cout << "T(df=12): q(p>0.75): " << StudentDist::quant(0.75, 12.0) << " (expected: 0.6954829)" << endl;
        cout << "T(df=12): q(p<0.52): " << StudentDist::quant(0.52, 12.0, false) << " (expected: -0.05121096)" << endl;
        cout << endl;


        // Chi-squared distribution:
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


        // F-distribution:
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


        // Continuous uniform distribution:
        cout << "U(1,3): pdf at x=0:   " << ContUniformDist::pdf(0.0, 1.0, 3.0) << " (expected: 0)" << endl;
        cout << "U(1,3): pdf at x=2:   " << ContUniformDist::pdf(2.0, 1.0, 3.0) << " (expected: 0.5)" << endl;
        cout << "U(1,3): pdf at x=4.5: " << ContUniformDist::pdf(4.5, 1.0, 3.0) << " (expected: 0)" << endl;
        cout << "U(1,3): P(X<0):   " << ContUniformDist::prob(0.0, 1.0, 3.0) << " (expected: 0)" << endl;
        cout << "U(1,3): P(X<1.5): " << ContUniformDist::prob(1.5, 1.0, 3.0) << " (expected: 0.25)" << endl;
        cout << "U(1,3): P(X<3.5): " << ContUniformDist::prob(3.5, 1.0, 3.0) << " (expected: 1)" << endl;
        cout << "U(1,3): P(X>-2):  " << ContUniformDist::prob(-2.0, 1.0, 3.0, false) << " (expected: 1)" << endl;
        cout << "U(1,3): P(X>2.7): " << ContUniformDist::prob(2.7, 1.0, 3.0, false) << " (expected: 0.15)" << endl;
        cout << "U(1,3): P(X>4):   " << ContUniformDist::prob(4.0, 1.0, 3.0, false) << " (expected: 0)" << endl;
        cout << "U(1,3): P(0.25<X<1.4): " << ContUniformDist::probInt(0.25, 1.4, 1.0, 3.0) << " (expected: 0.2)" << endl;
        cout << "U(1,3): P(1.9<X<3.7):  " << ContUniformDist::probInt(3.7, 1.9, 1.0, 3.0) << " (expected: 0.55)" << endl;
        cout << "U(1,3): q(p>0.42): " << ContUniformDist::quant(0.42, 1.0, 3.0) << " (expected: 1.84)" << endl;
        cout << "U(1,3): q(p>0.57): " << ContUniformDist::quant(0.57, 1.0, 3.0) << " (expected: 2.14)" << endl;
        cout << "U(1,3): q(p<0.34): " << ContUniformDist::quant(0.34, 1.0, 3.0, false) << " (expected: 2.32)" << endl;
        cout << "U(1,3): q(p<0.81): " << ContUniformDist::quant(0.81, 1.0, 3.0, false) << " (expected: 1.38)" << endl;
        cout << endl;


        // Binomial distribution:
        cout << "Binom(5, 0.6): pmf at k=2:   " << BinomDist::pmf(2, 5, 0.6) << " (expected: 0.2304)" << endl;
        cout << "Binom(5, 0.6): exp. value:   " << BinomDist::mean(5, 0.6) << " (expected: 3)" << endl;
        cout << "Binom(5, 0.6): variance:     " << BinomDist::var(5, 0.6) << " (expected: 1.2)" << endl;
        cout << "Binom(5, 0.6): std. dev.:    " << BinomDist::stdev(5, 0.6) << " (expected: 1.095445)" << endl;
        cout << "Binom(20, 0.6): normal approx: " << (BinomDist::normalApprox(20, 0.6) ? "true" : "false") << " (expected: FALSE)" << endl;
        cout << "Binom(30, 0.6): normal approx: " << (BinomDist::normalApprox(30, 0.6) ? "true" : "false") << " (expected: TRUE)" << endl;
        cout << "Binom(10, 0.6): P(X<=7):     " << BinomDist::prob(7, 10, 0.6) << " (expected: 0.8327102)" << endl;
        cout << "Binom(10, 0.6): P(X<7):      " << BinomDist::prob(7, 10, 0.6, false) << " (expected: 0.6177194)" << endl;
        cout << "Binom(10, 0.6): P(X>=6):     " << BinomDist::prob(6, 10, 0.6, true, false) << " (expected: 0.6331033)" << endl;
        cout << "Binom(10, 0.6): P(X>6):      " << BinomDist::prob(6, 10, 0.6, false, false) << " (expected: 0.3822806)" << endl;
        cout << "Binom(10, 0.6): P(X<0):      " << BinomDist::prob(0, 10, 0.6, false) << " (expected: 0)" << endl;
        cout << "Binom(10, 0.6): P(X<=0):     " << BinomDist::prob(0, 10, 0.6) << " (expected: 0.0001048576)" << endl;
        cout << "Binom(10, 0.6): P(X>0):      " << BinomDist::prob(0, 10, 0.6, false, false) << " (expected: 0.9998951)" << endl;
        cout << "Binom(10, 0.6): P(X>=0):     " << BinomDist::prob(0, 10, 0.6, true, false) << " (expected: 1)" << endl;
        cout << "Binom(10, 0.6): P(X<10):     " << BinomDist::prob(10, 10, 0.6, false) << " (expected: 0.9939534)" << endl;
        cout << "Binom(10, 0.6): P(X<=10):    " << BinomDist::prob(10, 10, 0.6) << " (expected: 1)" << endl;
        cout << "Binom(10, 0.6): P(X>10):     " << BinomDist::prob(10, 10, 0.6, false, false) << " (expected: 0)" << endl;
        cout << "Binom(10, 0.6): P(X>=10):    " << BinomDist::prob(10, 10, 0.6, true, false) << " (expected: 0.006046618)" << endl;
        cout << "Binom(10, 0.6): P(5<=X<=7):  " << BinomDist::probInt(5, 7, 10, 0.6) << " (expected: 0.6664716)" << endl;
        cout << "Binom(10, 0.6): P(5<X<=7):   " << BinomDist::probInt(5, 7, 10, 0.6, false) << " (expected: 0.4658135)" << endl;
        cout << "Binom(10, 0.6): P(4<=X<9):   " << BinomDist::probInt(4, 9, 10, 0.6, true, false) << " ( expected: 0.8988807)" << endl;
        cout << "Binom(10, 0.6): P(3<X<8):    " << BinomDist::probInt(3, 8, 10, 0.6, false, false) << " (expected: 0.7779484)" << endl;
        cout << "Binom(10, 0.6): q(p<0.4):    " << BinomDist::quant(0.4, 10, 0.6) << " (expected: 6)" << endl;
        cout << "Binom(10, 0.6): q(p<=0.15):  " << BinomDist::quant(0.15, 10, 0.6, false) << " (expected: 3)" << endl;
        cout << "Binom(10, 0.6): q(p>0.3):    " << BinomDist::quant(0.3, 10, 0.6, true, false) << " (expected: 8)" << endl;
        cout << "Binom(10, 0.6): q(p>=0.3):   " << BinomDist::quant(0.3, 10, 0.6, false, false) << " (expected: 7)" << endl;
        cout << endl;


        // Poisson distribution:
        cout << "Poisson(4): pmf at k=3: " << PoissonDist::pmf(3, 4.0) << " (expected: 0.1953668)" << endl;
        cout << "Poisson(4): P(X<=3):    " << PoissonDist::prob(3, 4.0) << " (expected: 0.4334701)" << endl;
        cout << "Poisson(4): P(X<5):     " << PoissonDist::prob(5, 4.0, false) << " (expected: 0.6288369)" << endl;
        cout << "Poisson(4): P(X>=6):    " << PoissonDist::prob(6, 4.0, true, false) << " (expected: 0.2148696)" << endl;
        cout << "Poisson(4): P(X>2):     " << PoissonDist::prob(2, 4.0, false, false) << " (expected: 0.7618967)" << endl;
        cout << "Poisson(4): P(X<=0):    " << PoissonDist::prob(0, 4.0) << " (expected: 0.01831564)" << endl;
        cout << "Poisson(4): P(X<0):     " << PoissonDist::prob(0, 4.0, false) << " (expected: 0)" << endl;
        cout << "Poisson(4): P(X>=0):    " << PoissonDist::prob(0, 4.0, true, false) << " (expected: 1)" << endl;
        cout << "Poisson(4): P(X>0):     " << PoissonDist::prob(0, 4.0, false, false) << " (expected: 0.9816844)" << endl;
        cout << "Poisson(4): P(5<=X<=7): " << PoissonDist::probInt(5, 7, 4.0) << " (expected: 0.3200294)" << endl;
        cout << "Poisson(4): P(5<X<=7):  " << PoissonDist::probInt(5, 7, 4.0, false) << " (expected: 0.163736)" << endl;
        cout << "Poisson(4): P(5<=X<7):  " << PoissonDist::probInt(5, 7, 4.0, true, false) << " (expected: 0.2604891)" << endl;
        cout << "Poisson(4): P(5<X<7):   " << PoissonDist::probInt(5, 7, 4.0, false, false) << " (expected: 0.1041956)" << endl;
        cout << "Poisson(4): q(p<0.3):   " << PoissonDist::quant<double, int>(0.3, 4.0) << " (expected: 3)" << endl;
        cout << "Poisson(4): q(p<=0.7):  " << PoissonDist::quant<double, int>(0.7, 4.0, false) << " (expected: 4)" << endl;
        cout << "Poisson(4): q(p>0.4):   " << PoissonDist::quant<double, int>(0.4, 4.0, true, false) << " (expected: 5)" << endl;
        cout << "Poisson(4): q(p>0.6):   " << PoissonDist::quant<double, int>(0.6, 4.0, false, false) << " (expected: 3)" << endl;

    }
    catch ( const StatisticsException& ex )
    {
        cerr << "Statistics exception caught: ";
        ex.what();
        cerr << endl;
    }
}
