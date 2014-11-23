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

# Typical file suffixes
OBJSUFFIX = .o
CPPSUFFIX = .cpp
HEADERSUFFIX = .hpp

# Top directory with the library's source code
LIBDIR = lib/

# Applications should include header files from this directory
APPINCDIR = include/

# Subdirectories of LIBDIR:
LIBEXCPDIR = $(LIBDIR)exception/
LIBUTILDIR = $(LIBDIR)util/
LIBMATRIXDIR = $(LIBDIR)matrix/
LIBLINEQDIR = $(LIBDIR)lineq/
LIBPOLYDIR = $(LIBDIR)polynomial/
LIBCOMBDIR = $(LIBDIR)combinatorics/
LIBSTATDIR = $(LIBDIR)statistics/
LIBCURVEFITDIR = $(LIBDIR)curve_fit/
LIBCALCULUSDIR = $(LIBDIR)calculus/
LIBROOTFINDDIR = $(LIBDIR)root_find/
LIBQUATDIR = $(LIBDIR)quaternion/
LIBRATIONALDIR = $(LIBDIR)rational/
LIBINTUTILDIR = $(LIBDIR)int_util/
LIBOMPDIR = $(LIBDIR)omp/


# Directory with unit testing files
TESTDIR = test/

# Directory with settings file(s)
SETTINGDIR = settings/

# Directory for all *.o and other intermediate files:
OBJDIR = obj/

# Directory for the target binary
BUILDDIR = build/

# Target binary, without any prefixes and suffixes
TARGETROOT = maintest


# Exception class names. Any non-applicable classes may be commented out.
# Note that appropriate file suffixes will be appended later.
EXCEPTIONCLASS = $(LIBEXCPDIR)IMathException
EXCEPTIONCLASS += $(LIBEXCPDIR)MatrixException
EXCEPTIONCLASS += $(LIBEXCPDIR)PolynomialException
EXCEPTIONCLASS += $(LIBEXCPDIR)RationalException
EXCEPTIONCLASS += $(LIBEXCPDIR)QuaternionException
EXCEPTIONCLASS += $(LIBEXCPDIR)LinearEquationSolverException
EXCEPTIONCLASS += $(LIBEXCPDIR)CurveFittingException
EXCEPTIONCLASS += $(LIBEXCPDIR)CombinatoricsException
EXCEPTIONCLASS += $(LIBEXCPDIR)IntFactorizationException
EXCEPTIONCLASS += $(LIBEXCPDIR)StatisticsException
EXCEPTIONCLASS += $(LIBEXCPDIR)FunctionException
EXCEPTIONCLASS += $(LIBEXCPDIR)CalculusException
EXCEPTIONCLASS += $(LIBEXCPDIR)RootFindException


# Nontemplated classes, i.e. their source files will be compiled.
# Note that appropriate file suffixes will be appended later.
COMPILECLASS =
COMPILECLASS += $(LIBRATIONALDIR)Rational
COMPILECLASS += $(LIBCOMBDIR)IntCombinatorics
COMPILECLASS += $(LIBINTUTILDIR)IntFactorization

# Templated classes. These classes are not compiled directly, instead 
# their source code will be included into files that need it.
# Note that appropriate file suffixes will be appended later.
GENERICCLASS =
GENERICCLASS += $(LIBUTILDIR)mtcopy
GENERICCLASS += $(LIBUTILDIR)NumericUtil
GENERICCLASS += $(LIBMATRIXDIR)MatrixGeneric
GENERICCLASS += $(LIBMATRIXDIR)SqMatrixGeneric
GENERICCLASS += $(LIBPOLYDIR)PolynomialGeneric
GENERICCLASS += $(LIBQUATDIR)QuaternionGeneric
GENERICCLASS += $(LIBCURVEFITDIR)CurveFittingGenericAb
GENERICCLASS += $(LIBCURVEFITDIR)PolynomialFittingGenericAb
GENERICCLASS += $(LIBCURVEFITDIR)PolynomialInterpolationGeneric
GENERICCLASS += $(LIBCURVEFITDIR)PolynomialRegressionGeneric
GENERICCLASS += $(LIBLINEQDIR)LinearEquationSolverGeneric
GENERICCLASS += $(LIBCOMBDIR)CombinationGeneric
GENERICCLASS += $(LIBCOMBDIR)PermutationGeneric
GENERICCLASS += $(LIBSTATDIR)SampleStatGeneric
GENERICCLASS += $(LIBSTATDIR)SampleQuantileGeneric
GENERICCLASS += $(LIBUTILDIR)IFunctionGeneric
GENERICCLASS += $(LIBCALCULUSDIR)IntegGeneric
GENERICCLASS += $(LIBCALCULUSDIR)DiffGeneric
GENERICCLASS += $(LIBROOTFINDDIR)RootFindGeneric


# Append file name extensions to exception classes
EXCEPTIONHEADER = $(addsuffix $(HEADERSUFFIX), $(EXCEPTIONCLASS) )
EXCEPTIONSRC = $(addsuffix $(CPPSUFFIX), $(EXCEPTIONCLASS) )
EXCEPTIONOBJ = $(addsuffix $(OBJSUFFIX), $(notdir $(EXCEPTIONCLASS)) )

# Append file name extensions to nontemplated (compiled) classes
COMPILEHEADER = $(addsuffix $(HEADERSUFFIX), $(COMPILECLASS) )
COMPILESRC = $(addsuffix $(CPPSUFFIX), $(COMPILECLASS) )
COMPILEOBJ = $(addsuffix $(OBJSUFFIX), $(notdir $(COMPILECLASS)) )

# Append file name extensions to templated classes.
# Note that these classes are not compiled so no *.o filenames are created
GENERICHEADER = $(addsuffix $(HEADERSUFFIX), $(GENERICCLASS) )
GENERICSRC = $(addsuffix $(CPPSUFFIX), $(GENERICCLASS) )

# Append file name extensions and a prefix to the final binary
TARGETSRC = $(addsuffix $(CPPSUFFIX), $(TARGETROOT) )
TARGETOBJ = $(addsuffix $(OBJSUFFIX), $(TARGETROOT) )
TARGET = $(addprefix $(BUILDDIR), $(TARGETROOT) )


# Join all desired *.o files into a single variable... 
OBJS = $(EXCEPTIONOBJ) $(COMPILEOBJ) $(TARGETOBJ)
# and prepend OBJDIR to each one.
OBJS := $(addprefix $(OBJDIR), $(OBJS) )

# Join all compilable source files into a single variable.
SRC = $(EXCEPTIONSRC) $(COMPILESRC) $(TARGETSRC)


# Compiler flags to produce debugging symbols and
# support OpenMP, respectively.
#
# Note: you should edit these variables if you
# use any other compiler than gcc.
DEBUG_FLAG = -g
OPENMP_FLAG = -fopenmp

# Compiler flags to produce debugging symbols and
# support OpenMP, respectively.
#
# Note: you should edit these variables if you
# use any other compiler than gcc.
DEBUG_MACRO = -D_DEBUG
OPENMP_MACRO = -D_OPENMP

# GCC flag and paths to include directories
INCLUDEFLAG = -I
INCLIB = $(LIBDIR)
LIBINCFLAG = $(INCLUDEFLAG)$(INCLIB)
APPINCFLAG = $(INCLUDEFLAG)$(APPINCDIR)

# Dependencies of OpenMP related files
OMPSETTINGDEP = $(SETTINGDIR)omp_settings.h
OMPLIBDEP = $(addprefix $(LIBOMPDIR), omp_header.h omp_coarse.h)

# Optional compiler flags
CPPFLAGS = -Wall -Wextra -Wno-unknown-pragmas $(LIBINCFLAG)

# Optional preprocessor macros
MACROS =

# Optional linker flags
LDFLAGS =

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
	$(eval LDFLAGS += $(OPENMP_FLAG))
#	
#	Typically a C++ compiler should automatically
#	predefine the _OPENMP macro if OpenMP is enabled.
#	If this is not the case, the line below must be uncommented:
#	
#	$(eval MACROS += $(OPENMP_MACRO))


# Build rules for exception classes
$(OBJDIR)IMathException$(OBJSUFFIX) : $(LIBEXCPDIR)IMathException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)MatrixException$(OBJSUFFIX) : $(LIBEXCPDIR)MatrixException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)PolynomialException$(OBJSUFFIX) : $(LIBEXCPDIR)PolynomialException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)RationalException$(OBJSUFFIX) : $(LIBEXCPDIR)RationalException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)QuaternionException$(OBJSUFFIX) : $(LIBEXCPDIR)QuaternionException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)LinearEquationSolverException$(OBJSUFFIX) : $(LIBEXCPDIR)LinearEquationSolverException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)CurveFittingException$(OBJSUFFIX) : $(LIBEXCPDIR)CurveFittingException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)CombinatoricsException$(OBJSUFFIX) : $(LIBEXCPDIR)CombinatoricsException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntFactorizationException$(OBJSUFFIX) : $(LIBEXCPDIR)IntFactorizationException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)StatisticsException$(OBJSUFFIX) : $(LIBEXCPDIR)StatisticsException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)FunctionException$(OBJSUFFIX) : $(LIBEXCPDIR)FunctionException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)CalculusException$(OBJSUFFIX) : $(LIBEXCPDIR)CalculusException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)RootFindException$(OBJSUFFIX) : $(LIBEXCPDIR)RootFindException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@


# Build rules for nontemplated classes
$(OBJDIR)Rational$(OBJSUFFIX) : $(LIBRATIONALDIR)Rational.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntCombinatorics$(OBJSUFFIX) : $(LIBCOMBDIR)IntCombinatorics.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntFactorization$(OBJSUFFIX) : $(LIBINTUTILDIR)IntFactorization.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@


# Build rule for the application that uses the library
$(OBJDIR)maintest$(OBJSUFFIX) : $(TESTDIR)maintest.cpp $(GENERICHEADER) $(GENERICSRC) \
                                $(OMPLIBDEP) $(OMPSETTINGDEP)
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@


# Build rule for linking the final binary
$(TARGET) : $(OBJDIR) $(BUILDDIR) $(OBJS) 
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
