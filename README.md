##About
A simple and robust C++ library for quaternion, matrix, rational and
polynomial arithmetics, suitable for scientific applications. See 
_maintest.cpp_ as an example of its usage. Or use Doxygen to generate 
documentation about the API.

##Usage
Only a few classes can be compiled into a static or dynamic library and
linked to an application. All other classes are templated and cannot be
compiled separately. Instead you should just include their corresponding
header files (*.h) into an application. See _Makefile_ for more info about
templated classes.

##Build
An example _Makefile_ with instructions to build a simple demo application
is provided. It assumes you have GCC installed and in your path. You are
free to modify a few variables to use another C++ compiler or toolchain.
For more details about official _make_ targets, please type:

`make help`

##OpenMP support
Mathematical algorithms in classes are currently being rewritten to support
OpenMP (parallelization using multiple threads on multicore systems).
Please note that OpenMP is currently disabled by default. To enable it,
build the application as:

`make openmp`

For more details, type

`make help`

Build instructions assume that GCC is used as a compiler. If you wish to
use any other compiler, you should edit _Makefile_ and modify appropriate
variables and/or build rules.

To build a parallelized version of the application, your compiler must
support __OpenMP version 3.0__ or newer.

##License
The library is licenced under the
[Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0).
See LICENSE.txt for more details.

##Author
The author of the library is Jernej Kova&#x010d;i&#x010d;.
