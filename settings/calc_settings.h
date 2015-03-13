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
 * Default algorithms and step sizes for implemented numerical
 * integration and differentiation algorithms. 
 */


#ifndef _MATH_CALC_SETTINGS_H_
#define _MATH_CALC_SETTINGS_H_


/**
 * Default step size for numerical integration.
 * 
 * For better flexibility with generic (templated) functions
 * the step size is given by its numerator and denominator.
 */
#define INTEG_STEP_NUM             ( 1 )
#define INTEG_STEP_DEN             ( 10000 )

/**
 * Default number of integration steps
 */
#define INTEG_DEFAULT_STEPS        ( 10000 )


/**
 * Default integration method
 */
#define INTEG_DEFAULT_METHOD       math::EIntegAlg::SIMPSON


/**
 * Breakpoint for evaluation of improper integrals
 */
#define INTEG_IMP_INT_BREAKPOINT_NUM     ( 5 )
#define INTEG_IMP_INT_BREAKPOINT_DEN     ( 1 )



/**
 * Default step size for numerical differentiation.
 * 
 * For better flexibility with generic (templated) functions
 * the step size is given by its numerator and denominator.
 */
#define DIFF_STEP_NUM              ( 1 )
#define DIFF_STEP_DEN              ( 10000 )

/**
 * Default 1st order differentiation method
 */
#define DIFF_DEFAULT_METHOD        math::EDiffMethod::CENTRAL


#endif  /* _MATH_CALC_SETTINGS_H_ */
