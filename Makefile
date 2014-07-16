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


#
# Type "make help" for more details.
#

# If necessary, add a path to compiler commands.
# Useful for cross compiling.
TOOLCHAIN =

# Compiler commands.
# The commands are set for the GCC compiler.
# If any other compiler is used, the commands must be modified appropriately.
CC = $(TOOLCHAIN)gcc
CPP = $(TOOLCHAIN)g++
FC = $(TOOLCHAIN)gfortran
LINKER = $(TOOLCHAIN)g++
AS = $(TOOLCHAIN)as
OBJCOPY = $(TOOLCHAIN)objcopy
AR = $(TOOLCHAIN)ar

# Compiler flags to produce deugging symbols and
# support OpenMP, respectively.
#
# Note: you should edit these variables if you
# use any other compiler than gcc.
DEBUG_FLAG = -g
OPENMP_FLAG = -fopenmp

# These preprocessor macros are predefined to build
# targets 'debug' and/or 'openmp'.
#
# Note: you should edit the macros if your compiler
# does not use -D to define macros.
DEBUG_MACRO = -DDEBUG
OPENMP_MACRO = -DOPENMP

# Optional compiler flags
CPPFLAGS =

# Optional preprocesor macros
MACROS =

# Optional linker flags
LDFLAGS =

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

debug : _debug_flags all

debug_rebuild : _debug_flags rebuild

openmp : _openmp_flags all

openmp_rebuild : _openmp_flags rebuild

debug_openmp : _debug_flags _openmp_flags all

debug_openmp_rebuild : _debug_flags _openmp_flags rebuild 

_debug_flags :
	$(eval CPPFLAGS += $(DEBUG_FLAG))
	$(eval MACROS += $(DEBUG_MACRO))

_openmp_flags :
	$(eval CPPFLAGS += $(OPENMP_FLAG))
	$(eval MACROS += $(OPENMP_MACRO))


#Build rules for exception classes
$(OBJDIR)MatrixException$(OBJSUFFIX) : MatrixException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)PolynomialException$(OBJSUFFIX) : PolynomialException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)RationalException$(OBJSUFFIX) : RationalException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)QuaternionException$(OBJSUFFIX) : QuaternionException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)LinearEquationSolverException$(OBJSUFFIX) : LinearEquationSolverException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)CurveFittingException$(OBJSUFFIX) : CurveFittingException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)CombinatoricsException$(OBJSUFFIX) : CombinatoricsException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntFactorizationException$(OBJSUFFIX) : IntFactorizationException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@


# Build rules for nontemplated classes
$(OBJDIR)Rational$(OBJSUFFIX) : Rational.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntCombinatorics$(OBJSUFFIX) : IntCombinatorics.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntFactorization$(OBJSUFFIX) : IntFactorization.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@


# Build rule for the application that uses the library
$(OBJDIR)maintest$(OBJSUFFIX) : maintest.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@


# Build rule for linking the final binary
$(TARGET) : $(OBJDIR) $(BUILDDIR) $(OBJS) $(GENERICHEADER) $(GENERICSRC)
	$(LINKER) $(LDFLAGS) $(OBJS) -o $@ 

# Cleanup directives:

clean_intermediate :
	$(RM) -r $(OBJDIR)

clean : clean_intermediate
	$(RM) -r $(BUILDDIR)

# Short help instructions:

help :
	@echo
	@echo Valid targets:
	@echo - all: builds missing dependencies and creates the target executable \'$(TARGET)\'.
	@echo - rebuild: rebuilds all dependencies and creates the target executable \'$(TARGET)\'.
	@echo - debug: same as \'all\', also includes debugging symbols to \'$(TARGET)\'.
	@echo - debug_rebuild: same as \'rebuild\', also includes debugging symbols to \'$(TARGET)\'.
	@echo - openmp: same as \'all\', also enables support for OpenMP.
	@echo - openmp_rebuild: same as \'rebuild\', also enables support for OpenMP.
	@echo - debug_openmp: same as \'all\', also includes debuging symbols to \'$(TARGET)\' and enables support for OpenMP.
	@echo - debug_openmp_rebuild: same as \'rebuild\', also includes debuging symbols to \'$(TARGET)\' and enables support for OpenMP.
	@echo - clean_intermediate: deletes all intermediate binaries, only keeps the target executable \'$(TARGET)\'.
	@echo - clean: deletes all intermediate binaries, incl. the target executable \'$(TARGET)\'.
	@echo - help: displays these help instructions.
	@echo


.PHONY : all rebuild debug debug_rebuild openmp openmp_rebuild \
         debug_openmp debug_openmp_rebuild clean clean_intermediate \
         _debug_flags _openmp_flags help
