About
-----
A simple and robust C++ library for mathematical and scientific applications. 
It supports matrix, rational, polynomial and quaternion arithmetics,
as well as various calculus, root-finding, combinatorics, statistical, curve
fitting algorithms, etc. See `test/maintest.cpp` as a basic demonstration. Or use
Doxygen to generate documentation about the API.

Usage
-----
Only headers (*.h files) from the directory `include` should be included
into applications. The included headers may additionally include several
templated classes that cannot be compiled separately. See _Makefile_ for
more info about the templated classes.

Additionally a simple demo application `test/maintest.cpp` is included
that performs a few basic unit tests. Its expected output can be found
in the file `test/test_output.out`.

Build
-----
An example _Makefile_ with instructions to build a simple demo application
is provided. It assumes you have GCC installed and in your path. You are
free to modify a few variables to use another C++ compiler or toolchain.
For more details about official _make_ targets, please type:

`make help`

OpenMP support
--------------
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
support __OpenMP version 3.0__ or higher.

License
-------
The library is licensed under the
[Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0).
See LICENSE.txt for more details.

Author
------
The author of the library is Jernej Kova&#x010d;i&#x010d;.
