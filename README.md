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

##License
The library is licenced under the Apache 2.0 license. See LICENSE.txt and
[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)
for more details.

##Author
The author of the library is Jernej Kova&#x010d;i&#x010d;.
