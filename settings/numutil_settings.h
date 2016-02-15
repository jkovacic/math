/*
Copyright 2015, Jernej Kovacic

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
 * Default values for epsilons, used to determine whether a value
 * of a floating point type is "close enough" to zero.
 */


#ifndef _MATH_NUMUTIL_SETTINGS_H_
#define _MATH_NUMUTIL_SETTINGS_H_


/**
 * Default epsilons for each floating point type (float, double and
 * long double).
 * 
 * If *_DIR equals 'true', the default epsilon for the specified type
 * will be set directly to *_VAL. Otherwise it will be set to *_VAL,
 * multiplied by the system specific machine epsilon (the difference
 * between 1 and the least value greater than 1 that is representable). 
 */


// Settings for 'float':
#define NUMUTIL_FLOAT_DIR                          ( false )
#define NUMUTIL_FLOAT_VAL                          ( 5.0f )

// Settings for 'double':
#define NUMUTIL_DOUBLE_DIR                         ( false )
#define NUMUTIL_DOUBLE_VAL                         ( 10.0 )

// Settings for 'long double':
#define NUMUTIL_LONGDOUBLE_DIR                     ( false )
#define NUMUTIL_LONGDOUBLE_VAL                     ( 10.0L )


#endif  /* _MATH_NUMUTIL_SETTINGS_H_ */
