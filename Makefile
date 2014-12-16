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
LIBSPECFUNDIR = $(LIBDIR)specfun/
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

# These objects should probably always be compiled and 
# "implicitly" linked into any application:
IMPL_OBJS = IMathException

# Test modules (source files in $(TESTDIR)) that will be linked 
# to the final test application.
# Unnecessary modules may be commented out
TESTFILES =
TESTFILES += calcTest
TESTFILES += combTest
TESTFILES += curvefitTest
TESTFILES += intcombTest
TESTFILES += intexpTest
TESTFILES += intfactorTest
TESTFILES += lineqTest
TESTFILES += matrixTest
TESTFILES += mtcopyTest
TESTFILES += polyTest
TESTFILES += quatTest
TESTFILES += specfunTest
TESTFILES += rationalTest
TESTFILES += ratmatTest
TESTFILES += rootfindTest
TESTFILES += statTest

# Prepend path and append suffixes to the selected test modules
TESTOBJS = $(addprefix $(OBJDIR), $(addsuffix $(OBJSUFFIX), $(TESTFILES) ))


# Dependencies of OpenMP headers

DEP_OMPSETTINGS = $(SETTINGDIR)omp_settings.h
DEP_OMPHEADER = $(LIBOMPDIR)omp_header.h
DEP_OMPCOARSE = $(LIBOMPDIR)omp_coarse.h

# Dependencies of other settings headers
DEP_STAT_SETTINGS = $(SETTINGDIR)stat_settings.h


# Object dependencies of templated classes
# Note #1: Some dependencies may repeat several times.
#          This is not a problem as all duplicated dependencies 
#          will be removed from the list later
# Note #2: No Suffixes (.hpp and .cpp) must be appended to file names
#          as it will be handled later.
#          The only exception are OMP dependencies (*.h files) that are
#          actually not objects
DEP_NUMUTIL = $(LIBUTILDIR)NumericUtil
DEP_MATHCONST = $(LIBUTILDIR)math_constant.h
DEP_MTCOPY = $(LIBUTILDIR)mtcopy $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_MTVECTOP = $(LIBUTILDIR)mtvectop $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_IFUNCTION = $(LIBUTILDIR)IFunctionGeneric

DEP_QUATERNION = $(LIBQUATDIR)QuaternionGeneric $(DEP_NUMUTIL) $(DEP_OMPSETTINGS)

DEP_POLYNOMIAL = $(LIBPOLYDIR)PolynomialGeneric $(DEP_NUMUTIL) $(DEP_MTCOPY) $(DEP_VECTOP) $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)

DEP_COMBINATION = $(LIBCOMBDIR)CombinationGeneric $(DEP_MTCOPY)
DEP_PERMUTATION =$(LIBCOMBDIR)PermutationGeneric $(DEP_MTCOPY)

DEP_INTEXP =

DEP_SPECFUN = $(LIBSPECFUNDIR)SpecFunGeneric $(DEP_NUMUTIL) $(DEP_MATHCONST)

DEP_MATRIX = $(LIBMATRIXDIR)MatrixGeneric $(DEP_MTCOPY) $(DEP_MTVECTOP) $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_SQMATRIX = $(LIBMATRIXDIR)SqMatrixGeneric $(DEP_NUMUTIL) $(DEP_MTCOPY) $(LIBLINEQDIR)LinearEquationSolverGeneric $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_LINEQ = $(LIBLINEQDIR)LinearEquationSolverGeneric $(DEP_NUMUTIL) $(DEP_MATRIX) $(DEP_SQMATRIX) $(DEP_OMPSETTINGS)

DEP_CURVEFITAB = $(LIBCURVEFITDIR)CurveFittingGenericAb $(DEP_NUMUTIL)
DEP_CURVEFITPOLY = $(LIBCURVEFITDIR)PolynomialFittingGenericAb $(DEP_CURVEFITAB) $(DEP_POLYNOMIAL)
DEP_POLYINT = $(LIBCURVEFITDIR)PolynomialInterpolationGeneric $(DEP_CURVEFITPOLY) $(DEP_POLYNOMIAL) $(DEP_OMPSETTINGS)
DEP_POLYREG = $(LIBCURVEFITDIR)PolynomialRegressionGeneric $(DEP_CURVEFITPOLY) $(DEP_POLYNOMIAL) $(DEP_LINEQ) $(DEP_SQMATRIX) $(DEP_MATRIX)

DEP_INTEG = $(LIBCALCULUSDIR)IntegGeneric $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE) 
DEP_DIFF = $(LIBCALCULUSDIR)DiffGeneric $(DEP_NUMUTIL) $(DEP_IFUNCTION)

DEP_ROOTFIND = $(LIBROOTFINDDIR)RootFindGeneric $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_DIFF)

DEP_SAMPLESTAT = $(LIBSTATDIR)SampleStatGeneric $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_SAMPLEQUANT = $(LIBSTATDIR)SampleQuantileGeneric $(DEP_NUMUTIL) $(DEP_MTCOPY)
DEP_NORMDIST = $(LIBSTATDIR)NormalDistGeneric $(DEP_NUMUTIL) $(DEP_MATHCONST) $(DEP_IFUNCTION) $(DEP_STAT_SETTINGS) $(DEP_ROOTFIND) $(DEP_SPECFUN)
DEP_STUDDIST = $(LIBSTATDIR)StudentDistGeneric $(DEP_NUMUTIL) $(DEP_STAT_SETTINGS) $(DEP_IFUNCTION) $(DEP_SPECFUN) $(DEP_INTEG)



# Object files (.o) that must be built and linked to
# the final test application if it includes the selected
# test module.
# Note #1: No paths and file suffixes must be appended to file names
#          as it will be handled later.
# Note #2: no problem if the same object file repeats among several test module
#          dependencies as it will be handled later
TEST_MTCOPY_OBJDEP =
TEST_QUAT_OBJDEP = QuaternionException
TEST_RAT_OBJDEP = Rational RationalException
TEST_MATRIX_OBJDEP = MatrixException
TEST_RATMAT_OBJDEP = Rational MatrixException RationalException
TEST_POLY_OBJDEP = PolynomialException
TEST_LINEQ_OBJDEP = LinearEquationSolverException
TEST_CURVEFIT_OBJDEP = CurveFittingException
TEST_INTEXP_OBJDEP = Rational MatrixException QuaternionException RationalException PolynomialException
TEST_INTFACTOR_OBJDEP = IntFactorization IntFactorizationException
TEST_INTCOMB_OBJDEP = IntCombinatorics CombinatoricsException
TEST_SPECFUN_OBJDEP = SpecFunException
TEST_COMB_OBJDEP = CombinatoricsException
TEST_CALC_OBJDEP = FunctionException CalculusException
TEST_ROOTFIND_OBJDEP = FunctionException RootFindException
TEST_STAT_OBJDEP = StatisticsException


# Templated classes included into test modules.
# The purpose of these variables is to check whether
# any templated class file has changed before the test
# module is built.
# Note: only include DEP_* variables declared above,
#       any duplicates will be removed later
TEST_MTCOPY_GENDEP = $(DEP_MTCOPY)
TEST_QUAT_GENDEP = $(DEP_QUATERNION)
TEST_RAT_GENDEP = 
TEST_MATRIX_GENDEP = $(DEP_MATRIX) $(DEP_SQMATRIX)
TEST_RATMAT_GENDEP = $(DEP_MATRIX)
TEST_POLY_GENDEP = $(DEP_POLYNOMIAL)
TEST_LINEQ_GENDEP = $(DEP_MATRIX) $(DEP_SQMATRIX) $(DEP_LINEQ)
TEST_CURVEFIT_GENDEP = $(DEP_POLYREG) $(DEP_POLYINT)
TEST_INTEXP_GENDEP = $(DEP_INTEXP) $(DEP_SQMATRIX) $(DEP_QUATERNION) $(DEP_POLYNOMIAL)
TEST_INTFACTOR_GENDEP = 
TEST_INTCOMB_GENDEP = 
TEST_SPECFUN_GENDEP = $(DEP_SPECFUN)
TEST_COMB_GENDEP = $(DEP_PERMUTATION) $(DEP_COMBINATION)
TEST_CALC_GENDEP = $(DEP_IFUNCTION) $(DEP_INTEG) $(DEP_DIFF)
TEST_ROOTFIND_GENDEP = $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_ROOTFIND)
TEST_STAT_GENDEP = $(DEP_MTCOPY) $(DEP_SAMPLESTAT) $(DEP_SAMPLEQUANT) $(DEP_NORMDIST) $(DEP_STUDDIST)


# Join object file dependencies for selected test modules.
# Lines for unnecessary test modules may be commented out.
# All duplicates will be removed in the next step.
# Note: the first line should include implicit object files
# that should be linked into any application.
TEST_LINKOBJ = $(IMPL_OBJS)
TEST_LINKOBJ += $(TEST_MTCOPY_OBJDEP)
TEST_LINKOBJ += $(TEST_QUAT_OBJDEP)
TEST_LINKOBJ += $(TEST_RAT_OBJDEP)
TEST_LINKOBJ += $(TEST_MATRIX_OBJDEP)
TEST_LINKOBJ += $(TEST_RATMAT_OBJDEP)
TEST_LINKOBJ += $(TEST_POLY_OBJDEP)
TEST_LINKOBJ += $(TEST_LINEQ_OBJDEP)
TEST_LINKOBJ += $(TEST_CURVEFIT_OBJDEP)
TEST_LINKOBJ += $(TEST_INTEXP_OBJDEP)
TEST_LINKOBJ += $(TEST_INTFACTOR_OBJDEP)
TEST_LINKOBJ += $(TEST_INTCOMB_OBJDEP)
TEST_LINKOBJ += $(TEST_COMB_OBJDEP)
TEST_LINKOBJ += $(TEST_SPECFUN_OBJDEP)
TEST_LINKOBJ += $(TEST_CALC_OBJDEP)
TEST_LINKOBJ += $(TEST_ROOTFIND_OBJDEP)
TEST_LINKOBJ += $(TEST_STAT_OBJDEP)

# Prepend a path and append $(OBJSUFFIX) to dependencies for
# selected test modules, remove dependencies using the Make's sort command
TEST_LINKOBJS = $(addprefix $(OBJDIR), $(addsuffix $(OBJSUFFIX), $(sort $(TEST_LINKOBJ)) ) )


# A convenience "function" that prepares list of all dependencies.
# Its input (referred as $1) is a list of files that declare and implement
# templated classes (a *_GENDEP variable), the function will first filter out
# OMP dependencies (*h. files), the remaining names will be appended
# *.hpp and *.cpp suffixes, finally the unmodified list of *.h files
# will be joined to the list.
#
# Usage of the "function":
#    $(call gen_deps,<list_of_files>)
#
# Note: any white spaces after the comma may be joined to the input list 
gen_deps = $(sort $(addsuffix $(HEADERSUFFIX),$(filter-out %.h,$1))) \
           $(sort $(addsuffix $(CPPSUFFIX),$(filter-out %.h,$1))) \
           $(sort $(filter %.h,$1))


# Append file name extensions and a prefix to the final binary
TARGETSRC = $(addsuffix $(CPPSUFFIX), $(TARGETROOT) )
TARGETOBJ = $(addprefix $(OBJDIR), $(addsuffix $(OBJSUFFIX), $(TARGETROOT) ))
TARGET = $(addprefix $(BUILDDIR), $(TARGETROOT) )


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

# Optional compiler flags
CPPFLAGS = -Wall -Wextra -Wno-unknown-pragmas $(LIBINCFLAG)

# Optional preprocessor macros
MACROS =

# Optional linker flags
LDFLAGS = -Wl,--allow-multiple-definition

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

$(OBJDIR)SpecFunException$(OBJSUFFIX) : $(LIBEXCPDIR)SpecFunException.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@


# Build rules for nontemplated classes
$(OBJDIR)Rational$(OBJSUFFIX) : $(LIBRATIONALDIR)Rational.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntCombinatorics$(OBJSUFFIX) : $(LIBCOMBDIR)IntCombinatorics.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@

$(OBJDIR)IntFactorization$(OBJSUFFIX) : $(LIBINTUTILDIR)IntFactorization.cpp
	$(CPP) -c $(CPPFLAGS) $(MACROS) $< -o $@



# Build rules for test modules
$(OBJDIR)mtcopyTest$(OBJSUFFIX) : $(TESTDIR)mtcopyTest.cpp $(call gen_deps,$(TEST_MTCOPY_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)quatTest$(OBJSUFFIX) : $(TESTDIR)quatTest.cpp $(call gen_deps,$(TEST_QUAT_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)rationalTest$(OBJSUFFIX) : $(TESTDIR)rationalTest.cpp $(call gen_deps,$(TEST_RAT_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)matrixTest$(OBJSUFFIX) : $(TESTDIR)matrixTest.cpp $(call gen_deps,$(TEST_MATRIX_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)ratmatTest$(OBJSUFFIX) : $(TESTDIR)ratmatTest.cpp $(call gen_deps,$(TEST_RATMAT_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)polyTest$(OBJSUFFIX) : $(TESTDIR)polyTest.cpp $(call gen_deps,$(TEST_POLY_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)lineqTest$(OBJSUFFIX) : $(TESTDIR)lineqTest.cpp $(call gen_deps,$(TEST_LINEQ_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)curvefitTest$(OBJSUFFIX) : $(TESTDIR)curvefitTest.cpp $(call gen_deps,$(TEST_CURVEFIT_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)intexpTest$(OBJSUFFIX) : $(TESTDIR)intexpTest.cpp $(call gen_deps,$(TEST_INTEXP_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)intfactorTest$(OBJSUFFIX) : $(TESTDIR)intfactorTest.cpp $(call gen_deps,$(TEST_INTFACTOR_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)intcombTest$(OBJSUFFIX) : $(TESTDIR)intcombTest.cpp $(call gen_deps,$(TEST_INTCOMB_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)specfunTest$(OBJSUFFIX) : $(TESTDIR)specfunTest.cpp $(call gen_deps,$(TEST_SPECFUN_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)combTest$(OBJSUFFIX) : $(TESTDIR)combTest.cpp $(call gen_deps,$(TEST_COMB_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)calcTest$(OBJSUFFIX) : $(TESTDIR)calcTest.cpp $(call gen_deps,$(TEST_CALC_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)rootfindTest$(OBJSUFFIX) : $(TESTDIR)rootfindTest.cpp $(call gen_deps,$(TEST_ROOTFIND_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@

$(OBJDIR)statTest$(OBJSUFFIX) : $(TESTDIR)statTest.cpp $(call gen_deps,$(TEST_STAT_GENDEP))
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@



# Build rule for the main test module
$(OBJDIR)maintest$(OBJSUFFIX) : $(TESTDIR)maintest.cpp
	$(CPP) -c $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< -o $@


# Build (link) rule for the final test application
$(TARGET) : $(OBJDIR) $(BUILDDIR) $(TARGETOBJ) $(TESTOBJS) $(TEST_LINKOBJS) 
	$(LINKER) $(LDFLAGS) $(TARGETOBJ) $(TESTOBJS) $(TEST_LINKOBJS) -o $@ 


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
