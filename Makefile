# Copyright 2014, Jernej Kovacic
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# If necessary, add a path to compiler commands.
# Useful for cross compiling.
TOOLCHAIN =

# Compiler commands.
# The commands are set for the GCC compiler.
# If any other compiler is used, the commands must be modified appropriately.
CC = $(TOOLCHAIN)gcc
CPP = $(TOOLCHAIN)g++
LINKER = $(TOOLCHAIN)g++
AS = $(TOOLCHAIN)as
OBJCOPY = $(TOOLCHAIN)objcopy
AR = $(TOOLCHAIN)ar

# Optional compiler flags
CPPFLAGS =

# Typical file sufixes
OBJSUFFIX = .o
CPPSUFFIX = .cpp
HEADERSUFFIX = .h

# Directory for all *.o and other intermediate files:
OBJDIR = obj/

# Directory for the target binary
BUILDDIR = build/

# Target binary, without any prefixes and suffixes
TARGETROOT = maintest


# Exception class names. Any nonapplicable classes may be commented out.
# Note that appropriate file sufixes will be appended later.
EXCEPTIONCLASS =
EXCEPTIONCLASS += MatrixException
EXCEPTIONCLASS += PolynomialException
EXCEPTIONCLASS += RationalException
EXCEPTIONCLASS += QuaternionException
EXCEPTIONCLASS += LinearEquationSolverException
EXCEPTIONCLASS += CurveFittingException
EXCEPTIONCLASS += CombinatoricsException
EXCEPTIONCLASS += IntFactorizationException

# Nontemplated classes, i.e. their source files will be compiled.
# Note that appropriate file sufixes will be appended later.
COMPILECLASS =
COMPILECLASS += Rational
COMPILECLASS += IntCombinatorics
COMPILECLASS += IntFactorization

# Templated classes. These classes are not compiled directly, instead 
# their source code will be included into files that need it.
# Note that appropriate file sufixes will be appended later.
GENERICCLASS =
GENERICCLASS += NumericUtil
GENERICCLASS += MatrixGeneric
GENERICCLASS += SqMatrixGeneric
GENERICCLASS += PolynomialGeneric
GENERICCLASS += QuaternionGeneric
GENERICCLASS += CurveFittingGenericAb
GENERICCLASS += PolynomialFittingGenericAb
GENERICCLASS += PolynomialInterpolationGeneric
GENERICCLASS += PolynomialRegressionGeneric
GENERICCLASS += LinearEquationSolverGeneric
GENERICCLASS += CombinationGeneric
GENERICCLASS += PermutationGeneric


# Append file name extensions to exception classes
EXCEPTIONHEADER = $(addsuffix $(HEADERSUFFIX), $(EXCEPTIONCLASS) )
EXCEPTIONSRC = $(addsuffix $(CPPSUFFIX), $(EXCEPTIONCLASS) )
EXCEPTIONOBJ = $(addsuffix $(OBJSUFFIX), $(EXCEPTIONCLASS) )

# Append file name extensions to nontemplated (compiled) classes
COMPILEHEADER = $(addsuffix $(HEADERSUFFIX), $(COMPILECLASS) )
COMPILESRC = $(addsuffix $(CPPSUFFIX), $(COMPILECLASS) )
COMPILEOBJ = $(addsuffix $(OBJSUFFIX), $(COMPILECLASS) )

# Append file name extensions to templated classes.
# Note that these classes are not compiled so no *.o filenames are created
GENERICHEADER = $(addsuffix $(HEADERSUFFIX), $(GENERICCLASS) )
GENERICSRC = $(addsuffix $(CPPSUFFIX), $(GENERICCLASS) )

# Append file name extensions and a prefix to the final binary
TARGETSRC = $(addsuffix $(CPPSUFFIX), $(TARGETROOT) )
TARGETOBJ = $(addsuffix $(OBJSUFFIX), $(TARGETROOT) )
TARGET = $(addprefix $(BUILDDIR), $(TARGETROOT) )


#Join all desired *.o files into a single variable... 
OBJS = $(EXCEPTIONOBJ) $(COMPILEOBJ) $(TARGETOBJ)
# and append OBJDIR to each one.
OBJS := $(addprefix $(OBJDIR), $(OBJS) )

# Join all source files to be compiled into a single variable.
SRC = $(EXCEPTIONSRC) $(COMPILESRC) $(TARGETSRC)


#
# Make rules:
#

all : $(TARGET)

rebuild : clean all

$(OBJDIR) :
	mkdir -p $@

$(BUILDDIR) :
	mkdir -p $@

#Build rules for exception classes
$(OBJDIR)MatrixException$(OBJSUFFIX) : MatrixException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)PolynomialException$(OBJSUFFIX) : PolynomialException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)RationalException$(OBJSUFFIX) : RationalException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)QuaternionException$(OBJSUFFIX) : QuaternionException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)LinearEquationSolverException$(OBJSUFFIX) : LinearEquationSolverException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)CurveFittingException$(OBJSUFFIX) : CurveFittingException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)CombinatoricsException$(OBJSUFFIX) : CombinatoricsException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)IntFactorizationException$(OBJSUFFIX) : IntFactorizationException.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@


# Build rules for nontemplated classes
$(OBJDIR)Rational$(OBJSUFFIX) : Rational.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)IntCombinatorics$(OBJSUFFIX) : IntCombinatorics.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(OBJDIR)IntFactorization$(OBJSUFFIX) : IntFactorization.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@


# Build rule for the application that uses the library
$(OBJDIR)maintest$(OBJSUFFIX) : maintest.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@


# Build rule for the final binary
$(TARGET) : $(OBJDIR) $(BUILDDIR) $(OBJS) $(GENERICHEADER) $(GENERICSRC)
	$(LINKER) $(CPPFLAGS) $(OBJS) -o $@ 

# Cleanup directives:

clean_intermediate :
	rm -rf $(OBJDIR)

clean : clean_intermediate
	rm -rf $(BUILDDIR)

.PHONY : all rebuild clean clean_intermediate
