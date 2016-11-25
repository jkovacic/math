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
# If any other compiler is used, the commands must be modified accordingly.
#CC = $(TOOLCHAIN)gcc
CPP = $(TOOLCHAIN)g++
#FC = $(TOOLCHAIN)gfortran
LINKER = $(TOOLCHAIN)g++
#AS = $(TOOLCHAIN)as
#OBJCOPY = $(TOOLCHAIN)objcopy
#AR = $(TOOLCHAIN)ar

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
LIBPROBDISTDIR = $(LIBSTATDIR)dist/

# Directory with smoke test source files
TESTDIR = test/

# Directory with settings file(s)
SETTINGDIR = settings/

# Directory for all *.o and other intermediate files:
OBJDIR = obj/

# Directory for the target binary
BUILDDIR = build/
 
# Target binary, without any prefixes and suffixes
TARGETROOT = maintest


# Dependencies of OpenMP headers

DEP_OMPSETTINGS = $(SETTINGDIR)omp_settings.h
DEP_OMPHEADER = $(LIBOMPDIR)omp_header.h
DEP_OMPCOARSE = $(LIBOMPDIR)omp_coarse.h

# Dependencies of other settings' headers
DEP_NUMUTIL_SETTINGS = $(SETTINGDIR)numutil_settings.h
DEP_STAT_SETTINGS = $(SETTINGDIR)stat_settings.h
DEP_PROBDIST_SETTINGS = $(SETTINGDIR)probdist_settings.h
DEP_SPECFUN_SETTINGS = $(SETTINGDIR)specfun_settings.h
DEP_CALC_SETTINGS = $(SETTINGDIR)calc_settings.h
DEP_ROOTFIND_SETTINGS = $(SETTINGDIR)rootfind_settings.h
DEP_MATRIX_SETTINGS = $(SETTINGDIR)matrix_settings.h
DEP_LINEQ_SETTINGS = $(SETTINGDIR)lineq_settings.h

# Header dependencies of templated classes
# Note #1: Dependencies may be listed several times and all
#          duplicated dependencies will be removed from the list later
# Note #2: No Suffixes (.hpp and .cpp) must be appended to file names
#          as it will be handled later.
#          The only exception are OMP and settings dependencies (*.h files) that are
#          actually not objects.
DEP_NUMUTIL = $(LIBUTILDIR)NumericUtil $(DEP_NUMUTIL_SETTINGS)
DEP_MATHCONST = $(LIBUTILDIR)math_constant.h
DEP_FILLVECTORS = $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_PSEUDOFUNC = $(LIBUTILDIR)PseudoFunctionGeneric
DEP_MTCOPY = $(LIBUTILDIR)mtcopy $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_MTVECTOP = $(LIBUTILDIR)mtvectop $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE) $(DEP_NUMUTIL)
DEP_MTSWAP = $(LIBUTILDIR)mtswap $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_IFUNCTION = $(LIBUTILDIR)IFunctionGeneric

DEP_SELECTION = $(LIBUTILDIR)SelectionGeneric $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)

DEP_QUATERNION = $(LIBQUATDIR)QuaternionGeneric $(DEP_NUMUTIL) $(DEP_OMPSETTINGS)

DEP_POLYNOMIAL = $(LIBPOLYDIR)PolynomialGeneric $(DEP_NUMUTIL) $(DEP_MTCOPY) $(DEP_VECTOP) $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)

DEP_COMBINATION = $(LIBCOMBDIR)CombinationGeneric $(DEP_MTCOPY)
DEP_PERMUTATION = $(LIBCOMBDIR)PermutationGeneric $(DEP_MTCOPY)

DEP_INTUTIL = $(LIBINTUTILDIR)IntUtilGeneric
DEP_INTFUNCTION = $(LIBINTUTILDIR)IntFunctionGeneric $(DEP_INTUTIL)
DEP_INTFACT = $(LIBINTUTILDIR)IntFactorizationGeneric $(DEP_INTUTIL) $(DEP_INTFUNCTION)
DEP_INTEXP = $(LIBINTUTILDIR)IntExponentiatorGeneric $(DEP_INTUTIL)
DEP_INTCOMB = $(LIBCOMBDIR)IntCombinatoricsGeneric $(DEP_INTUTIL)

DEP_RAT = $(LIBRATIONALDIR)RationalGeneric $(DEP_INTUTIL) $(DEP_INTFACT)

DEP_SPECFUN = $(LIBSPECFUNDIR)SpecFunGeneric $(DEP_NUMUTIL) $(DEP_MATHCONST) $(DEP_SPECFUN_SETTINGS) $(LIBSPECFUNDIR)lanczos_coef.h
DEP_SPECFUN += $(LIBSPECFUNDIR)CtdFracGeneric $(DEP_NUMUTIL) $(DEP_SPECFUN_SETTINGS)

DEP_PIVOT = $(LIBMATRIXDIR)PivotGeneric $(DEP_NUMUTIL) $(DEP_PSEUDOFUNC) $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_MATRIX = $(DEP_MATRIX_SETTINGS) $(LIBMATRIXDIR)MatrixGeneric $(DEP_PIVOT) \
             $(DEP_MTCOPY) $(DEP_MTVECTOP) $(DEP_MTSWAP) $(DEP_OMPSETTINGS) \
             $(DEP_OMPHEADER) $(DEP_OMPCOARSE) $(DEP_NUMUTIL)
DEP_LINEQ = $(LIBMATRIXDIR)LinearEquationSolverGeneric $(DEP_NUMUTIL) $(DEP_PSEUDOFUNC) $(DEP_MATRIX) $(DEP_PIVOT) $(DEP_LINEQ_SETTINGS) \
            $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)

DEP_CURVEFITAB = $(LIBCURVEFITDIR)CurveFittingGenericAb $(DEP_NUMUTIL)
DEP_CURVEFITPOLY = $(LIBCURVEFITDIR)PolynomialFittingGenericAb $(DEP_CURVEFITAB) $(DEP_POLYNOMIAL)
DEP_POLYINT = $(LIBCURVEFITDIR)PolynomialInterpolationGeneric $(DEP_CURVEFITPOLY) $(DEP_POLYNOMIAL) $(DEP_OMPSETTINGS)
DEP_POLYREG = $(LIBCURVEFITDIR)PolynomialRegressionGeneric $(DEP_CURVEFITPOLY) $(DEP_POLYNOMIAL) $(DEP_LINEQ) $(DEP_MATRIX)

DEP_INTEG = $(LIBCALCULUSDIR)IntegGeneric $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_CALC_SETTINGS) $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE) 
DEP_DIFF = $(LIBCALCULUSDIR)DiffGeneric $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_CALC_SETTINGS)

DEP_ROOTFIND = $(LIBROOTFINDDIR)RootFindGeneric $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_DIFF) $(DEP_ROOTFIND_SETTINGS)

DEP_SAMPLESTAT = $(LIBSTATDIR)SampleStatGeneric $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE) $(DEP_INTUTIL) $(DEP_INTEXP) $(DEP_NUMUTIL)
DEP_SAMPLEQUANTAB = $(LIBSTATDIR)SampleQuantileGenericAb $(DEP_STAT_SETTINGS) $(DEP_NUMUTIL) $(DEP_INTUTIL)
DEP_SAMPLEQUANTSORTEDARRAY = $(LIBSTATDIR)SampleQuantileSortedArrayGeneric $(DEP_SAMPLEQUANTAB) $(DEP_STAT_SETTINGS) $(DEP_MTCOPY)
DEP_SAMPLEQUANTSELECTION = $(LIBSTATDIR)SampleQuantileSelectionGeneric $(DEP_SAMPLEQUANTAB) $(DEP_STAT_SETTINGS) $(DEP_MTCOPY) $(DEP_SELECTION) $(DEP_SAMPLESTAT) \
                           $(DEP_OMPSETTINGS) $(DEP_OMPHEADER) $(DEP_OMPCOARSE)
DEP_NORMDIST = $(LIBPROBDISTDIR)NormalDistGeneric $(DEP_NUMUTIL) $(DEP_PROBDIST_SETTINGS) $(DEP_MATHCONST) $(DEP_SPECFUN)
DEP_STUDDIST = $(LIBPROBDISTDIR)StudentDistGeneric $(DEP_NUMUTIL) $(DEP_PROBDIST_SETTINGS) $(DEP_MATHCONST) $(DEP_SPECFUN)
DEP_CHISQDIST = $(LIBPROBDISTDIR)ChiSquareDistGeneric $(DEP_NUM_UTIL) $(DEP_PROBDIST_SETTINGS) $(DEP_SPECFUN)
DEP_FDIST = $(LIBPROBDISTDIR)FDistGeneric $(DEP_NUM_UTIL) $(DEP_PROBDIST_SETTINGS) $(DEP_SPECFUN)
DEP_CONTUNIFDIST = $(LIBPROBDISTDIR)ContUniformDistGeneric $(DEP_NUM_UTIL)
DEP_TRIANGULARDIST = $(LIBPROBDISTDIR)TriangularDistGeneric $(DEP_NUM_UTIL)
DEP_BINOMDIST = $(LIBPROBDISTDIR)BinomDistGeneric $(DEP_PROBDIST_SETTINGS) $(DEP_NUM_UTIL) $(DEP_SPECFUN) $(DEP_INTUTIL) $(DEP_INTEXP) $(DEP_INTCOMB)
DEP_POISSONDIST = $(LIBPROBDISTDIR)PoissonDistGeneric $(DEP_PROBDIST_SETTINGS) $(DEP_NUM_UTIL) $(DEP_SPECFUN) $(DEP_INTUTIL) $(DEP_INTEXP)



# Object dependencies, i.e. a list of other files (typically exception files)
# that must be compiled and linked to the application if certain header file
# is included into any module. Dependencies may be listed several
# times as all duplicates will be removed later.
# Note: object dependencies must not contain their directories and suffixes
#       as this will be handled later.
OBJDEP_NUMUTIL = 
OBJDEP_PSEUDOFUNC = 
OBJDEP_MTCOPY = 
OBJDEP_MTVECTOP = 
OBJDEP_MTSWAP = 
OBJDEP_FILLVECTORS = FillVectors
OBJDEP_IFUNCTION = FunctionException
OBJDEP_SELECTION = SelectionException $(OBJDEP_FILLVECTORS)
OBJDEP_QUATERNION = QuaternionException $(OBJDEP_NUMUTIL)
OBJDEP_POLYNOMIAL = PolynomialException $(OBJDEP_NUMUTIL) $(OBJDEP_MTCOPY) $(OBJDEP_MTVECTOP)
OBJDEP_COMBINATION = CombinatoricsException $(OBJDEP_MTCOPY)
OBJDEP_PERMUTATION = CombinatoricsException $(OBJDEP_MTCOPY)
OBJDEP_INTUTIL = 
OBJDEP_INTFUNCTION = IntFunctionException $(OBJDEP_INTUTIL)
OBJDEP_INTFACT = IntFactorizationException $(OBJDEP_INTUTIL) $(OBJDEP_INTFUNCTION)
OBJDEP_INTEXP = IntFactorizationException $(OBJDEP_INTUTIL)
OBJDEP_INTCOMB = CombinatoricsException $(OBJDEP_MTCOPY)
OBJDEP_RAT = RationalException $(OBJDEP_INTUTIL) $(OBJDEP_INTFACT)
OBJDEP_SPECFUN = SpecFunException $(OBJDEP_IFUNCTION) $(OBJDEP_NUMUTIL)
OBJDEP_PIVOT = MatrixException $(OBJDEP_FILLVECTORS) $(OBJDEP_PSEUDOFUNC) $(OBJDEP_NUMUTIL)
OBJDEP_MATRIX = MatrixException $(OBJDEP_PIVOT) $(OBJDEP_NUMUTIL) \
                $(OBJDEP_MTCOPY) $(OBJDEP_MTVECTOP) $(OBJDEP_MTSWAP) 
OBJDEP_LINEQ = MatrixException $(OBJDEP_MATRIX) $(OBJDEP_PSEUDOFUNC) \
               $(OBJDEP_NUMUTIL) $(OBJDEP_PSEUDOFUNC)
OBJDEP_CURVEFITAB = CurveFittingException $(OBJDEP_NUMUTIL)
OBJDEP_CURVEFITPOLY = $(OBJDEP_CURVEFITAB) $(OBJDEP_POLYNOMIAL)
OBJDEP_POLYINT = $(OBJDEP_CURVEFITPOLY)
OBJDEP_POLYREG = $(OBJDEP_CURVEFITPOLY) $(OBJDEP_MATRIX) $(OBJDEP_LINEQ)
OBJDEP_INTEG = CalculusException $(OBJDEP_IFUNCTION) $(OBJDEP_NUMUTIL)
OBJDEP_DIFF = CalculusException $(OBJDEP_IFUNCTION) $(OBJDEP_NUMUTIL)
OBJDEP_ROOTFIND = RootFindException $(OBJDEP_IFUNCTION) $(OBJDEP_NUMUTIL) $(OBJDEP_DIFF)
OBJDEP_SAMPLESTAT = StatisticsException $(OBJDEP_NUMUTIL) $(OBJDEP_INTUTIL) $(OBJDEP_INTEXP)
OBJDEP_SAMPLEQUANTAB = StatisticsException $(OBJDEP_INTUTIL) $(OBJDEP_NUMUTIL)
OBJDEP_SAMPLEQUANTSORTEDARRAY = $(OBJDEP_SAMPLEQUANTAB) $(OBJDEP_MTCOPY)
OBJDEP_SAMPLEQUANTSELECTION = $(OBJDEP_SAMPLEQUANTAB) $(OBJDEP_SELECTION) \
                              $(OBJDEP_SAMPLESTAT) $(OBJDEP_MTCOPY)
OBJDEP_NORMDIST = StatisticsException $(OBJDEP_SPECFUN) $(OBJDEP_NUMUTIL)
OBJDEP_STUDDIST = StatisticsException $(OBJDEP_SPECFUN) $(OBJDEP_NUMUTIL)
OBJDEP_CHISQDIST = StatisticsException $(OBJDEP_SPECFUN) $(OBJDEP_NUMUTIL)
OBJDEP_FDIST = StatisticsException $(OBJDEP_SPECFUN) $(OBJDEP_NUMUTIL)
OBJDEP_CONTUNIFDIST = StatisticsException $(OBJDEP_NUMUTIL)
OBJDEP_TRIANGULARDIST = StatisticsException $(OBJDEP_NUMUTIL)
OBJDEP_BINOMDIST = StatisticsException $(OBJDEP_SPECFUN) $(OBJDEP_NUMUTIL) \
                   $(OBJDEP_INTUTIL) $(OBJDEP_INTCOMB) $(OBJDEP_INTEXP)
OBJDEP_POISSONDIST = StatisticsException $(OBJDEP_SPECFUN) $(OBJDEP_NUMUTIL) \
                     $(OBJDEP_INTUTIL) $(OBJDEP_INTEXP)


# These objects should probably always be compiled and 
# "implicitly" linked into each application:
IMPL_OBJS = IMathException


# Object files (.o) that must be built and linked to
# the final test application if it includes the selected
# test module. All lines with unnecessary modules may be commented out.
# Note #1: No paths and file suffixes must be appended to file names
#          as it will be handled later.
# Note #2: the same object file may be listed several times
#          as duplicates will be handled later.
# Note #3: the first line should include implicit object files
#          that should be linked into each application.
TESTLINKOBJ = $(IMPL_OBJS)
TESTLINKOBJ += numutilTest $(OBJDEP_NUMUTIL)
TESTLINKOBJ += selectionTest $(OBJDEP_SELECTION) $(OBJDEP_MTCOPY)
TESTLINKOBJ += calcTest $(OBJDEP_IFUNCTION) $(OBJDEP_INTEG) $(OBJDEP_DIFF)
TESTLINKOBJ += combTest $(OBJDEP_PERMUTATION) $(OBJDEP_COMBINATION)
TESTLINKOBJ += curvefitTest $(OBJDEP_POLYINT) $(OBJDEP_POLYREG)
TESTLINKOBJ += intcombTest $(OBJDEP_INTCOMB)
TESTLINKOBJ += intfunctionTest $(OBJDEP_INTFUNCTION)
TESTLINKOBJ += intexpTest $(OBJDEP_INTEXP) $(OBJDEP_MATRIX) $(OBJDEP_QUATERNION) $(OBJDEP_RAT) $(OBJDEP_POLYNOMIAL)
TESTLINKOBJ += intfactorTest $(OBJDEP_INTFACT)
TESTLINKOBJ += lineqTest $(OBJDEP_LINEQ) $(OBJDEP_MATRIX)
TESTLINKOBJ += matrixTest $(OBJDEP_MATRIX)
TESTLINKOBJ += mtcopyTest $(OBJDEP_MTCOPY)
TESTLINKOBJ += polyTest $(OBJDEP_POLYNOMIAL)
TESTLINKOBJ += quatTest $(OBJDEP_QUATERNION)
TESTLINKOBJ += specfunTest $(OBJDEP_SPECFUN)
TESTLINKOBJ += rationalTest $(OBJDEP_RAT)
TESTLINKOBJ += ratmatTest $(OBJDEP_RAT) $(OBJDEP_MATRIX)
TESTLINKOBJ += rootfindTest $(OBJDEP_NUMUTIL) $(OBJDEP_IFUNCTION) $(OBJDEP_ROOTFIND)
TESTLINKOBJ += statTest $(OBJDEP_MTCOPY) $(OBJDEP_SAMPLESTAT) $(OBJDEP_SAMPLEQUANTSORTEDARRAY) $(OBJDEP_SAMPLEQUANTSELECTION)
TESTLINKOBJ += probdistTest $(OBJDEP_NORMDIST) $(OBJDEP_STUDDIST) $(OBJDEP_CHISQDIST) $(OBJDEP_FDIST) $(OBJDEP_CONTUNIFDIST) \
                            $(OBJDEP_TRIANGULARDIST) $(OBJDEP_BINOMDIST) $(OBJDEP_POISSONDIST)


# Prepend a path and append $(OBJSUFFIX) to object dependencies for
# selected test modules, remove duplicates using the Make's sort command.
# Note: even though this is not necessary, the list will contain all test
# objects (*Test files) first, followed by all other object dependencies. 
TEST_LINKOBJS_TEST = $(filter %Test, $(sort $(TESTLINKOBJ)) )
TEST_LINKOBJS_OTHER = $(filter-out %Test, $(sort $(TESTLINKOBJ)) )
TEST_LINKOBJS = $(addprefix $(OBJDIR), $(addsuffix $(OBJSUFFIX), $(TEST_LINKOBJS_TEST) $(TEST_LINKOBJS_OTHER) ))


# Templated classes included into test modules.
# The purpose of these variables is to check whether
# any templated class file has changed before the test
# module is built.
# Note: only include DEP_* variables declared above,
#       all duplicates will be removed later
TEST_NUMUTIL_GENDEP = $(DEP_NUMUTIL)
TEST_MTCOPY_GENDEP = $(DEP_MTCOPY)
TEST_SELECTION_GENDEP = $(DEP_MTCOPY) $(DEP_SELECTION)
TEST_QUAT_GENDEP = $(DEP_QUATERNION)
TEST_RAT_GENDEP = $(DEP_RAT)
TEST_MATRIX_GENDEP = $(DEP_MATRIX)
TEST_RATMAT_GENDEP = $(DEP_RAT) $(DEP_MATRIX)
TEST_POLY_GENDEP = $(DEP_POLYNOMIAL)
TEST_LINEQ_GENDEP = $(DEP_MATRIX) $(DEP_LINEQ)
TEST_CURVEFIT_GENDEP = $(DEP_POLYREG) $(DEP_POLYINT)
TEST_INTFUNCTION_GENDEP = $(DEP_INTFUNCTION)
TEST_INTEXP_GENDEP = $(DEP_INTEXP) $(DEP_MATRIX) $(DEP_QUATERNION) $(DEP_POLYNOMIAL) $(DEP_RAT)
TEST_INTFACTOR_GENDEP = $(DEP_INTFACT)
TEST_INTCOMB_GENDEP = $(DEP_INTCOMB) 
TEST_SPECFUN_GENDEP = $(DEP_SPECFUN)
TEST_COMB_GENDEP = $(DEP_PERMUTATION) $(DEP_COMBINATION)
TEST_CALC_GENDEP = $(DEP_IFUNCTION) $(DEP_INTEG) $(DEP_DIFF)
TEST_ROOTFIND_GENDEP = $(DEP_NUMUTIL) $(DEP_IFUNCTION) $(DEP_ROOTFIND)
TEST_STAT_GENDEP = $(DEP_MTCOPY) $(DEP_SAMPLESTAT) $(DEP_SAMPLEQUANTSORTEDARRAY) $(DEP_SAMPLEQUANTSELECTION)
TEST_PROBDIST_GENDEP = $(DEP_NORMDIST) $(DEP_STUDDIST) $(DEP_CHISQDIST) $(DEP_FDIST) \
                       $(DEP_CONTUNIFDIST) $(DEP_TRIANGULARDIST) $(DEP_BINOMDIST) $(DEP_POISSONDIST)


# A convenience "function" that prepares list of all dependencies.
# Its input (referred as $1) is a list of files that declare and implement
# templated classes (a *_GENDEP variable), the function will first filter out
# OMP and settings dependencies (*h. files), the remaining names will be appended
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

# GCC compile and output flag
CFLAG = -c
OFLAG = -o

# GCC flag and paths to include directories
INCLUDEFLAG = -I
INCLIB = $(LIBDIR)
LIBINCFLAG = $(INCLUDEFLAG)$(INCLIB)
APPINCFLAG = $(INCLUDEFLAG)$(APPINCDIR)

# Optional compiler flags
CPPFLAGS = -Wall -Wextra -Wno-unknown-pragmas -Werror $(LIBINCFLAG)

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
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)SelectionException$(OBJSUFFIX) : $(LIBEXCPDIR)SelectionException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)MatrixException$(OBJSUFFIX) : $(LIBEXCPDIR)MatrixException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)PolynomialException$(OBJSUFFIX) : $(LIBEXCPDIR)PolynomialException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)RationalException$(OBJSUFFIX) : $(LIBEXCPDIR)RationalException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)QuaternionException$(OBJSUFFIX) : $(LIBEXCPDIR)QuaternionException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)CurveFittingException$(OBJSUFFIX) : $(LIBEXCPDIR)CurveFittingException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)CombinatoricsException$(OBJSUFFIX) : $(LIBEXCPDIR)CombinatoricsException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)IntFunctionException$(OBJSUFFIX) : $(LIBEXCPDIR)IntFunctionException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)IntFactorizationException$(OBJSUFFIX) : $(LIBEXCPDIR)IntFactorizationException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)StatisticsException$(OBJSUFFIX) : $(LIBEXCPDIR)StatisticsException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)FunctionException$(OBJSUFFIX) : $(LIBEXCPDIR)FunctionException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)CalculusException$(OBJSUFFIX) : $(LIBEXCPDIR)CalculusException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)RootFindException$(OBJSUFFIX) : $(LIBEXCPDIR)RootFindException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)SpecFunException$(OBJSUFFIX) : $(LIBEXCPDIR)SpecFunException.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@


# Build rules for object dependencies
$(OBJDIR)FillVectors$(OBJSUFFIX) : $(LIBUTILDIR)FillVectors.cpp $(LIBUTILDIR)FillVectors.h $(call gen_deps,$(DEP_FILLVECTORS))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(MACROS) $< $(OFLAG) $@


# Build rules for test modules
$(OBJDIR)numutilTest$(OBJSUFFIX) : $(TESTDIR)numutilTest.cpp $(call gen_deps,$(TEST_NUMUTIL_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)mtcopyTest$(OBJSUFFIX) : $(TESTDIR)mtcopyTest.cpp $(call gen_deps,$(TEST_MTCOPY_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)selectionTest$(OBJSUFFIX) : $(TESTDIR)selectionTest.cpp $(call gen_deps,$(TEST_SELECTION_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)quatTest$(OBJSUFFIX) : $(TESTDIR)quatTest.cpp $(call gen_deps,$(TEST_QUAT_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)rationalTest$(OBJSUFFIX) : $(TESTDIR)rationalTest.cpp $(call gen_deps,$(TEST_RAT_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)matrixTest$(OBJSUFFIX) : $(TESTDIR)matrixTest.cpp $(call gen_deps,$(TEST_MATRIX_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)ratmatTest$(OBJSUFFIX) : $(TESTDIR)ratmatTest.cpp $(call gen_deps,$(TEST_RATMAT_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)polyTest$(OBJSUFFIX) : $(TESTDIR)polyTest.cpp $(call gen_deps,$(TEST_POLY_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)lineqTest$(OBJSUFFIX) : $(TESTDIR)lineqTest.cpp $(call gen_deps,$(TEST_LINEQ_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)curvefitTest$(OBJSUFFIX) : $(TESTDIR)curvefitTest.cpp $(call gen_deps,$(TEST_CURVEFIT_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)intfunctionTest$(OBJSUFFIX) : $(TESTDIR)intfunctionTest.cpp $(call gen_deps,$(TEST_INTFUNCTION_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@
	
$(OBJDIR)intexpTest$(OBJSUFFIX) : $(TESTDIR)intexpTest.cpp $(call gen_deps,$(TEST_INTEXP_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)intfactorTest$(OBJSUFFIX) : $(TESTDIR)intfactorTest.cpp $(call gen_deps,$(TEST_INTFACTOR_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)intcombTest$(OBJSUFFIX) : $(TESTDIR)intcombTest.cpp $(call gen_deps,$(TEST_INTCOMB_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)specfunTest$(OBJSUFFIX) : $(TESTDIR)specfunTest.cpp $(call gen_deps,$(TEST_SPECFUN_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)combTest$(OBJSUFFIX) : $(TESTDIR)combTest.cpp $(call gen_deps,$(TEST_COMB_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)calcTest$(OBJSUFFIX) : $(TESTDIR)calcTest.cpp $(call gen_deps,$(TEST_CALC_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)rootfindTest$(OBJSUFFIX) : $(TESTDIR)rootfindTest.cpp $(call gen_deps,$(TEST_ROOTFIND_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)statTest$(OBJSUFFIX) : $(TESTDIR)statTest.cpp $(call gen_deps,$(TEST_STAT_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@

$(OBJDIR)probdistTest$(OBJSUFFIX) : $(TESTDIR)probdistTest.cpp $(call gen_deps,$(TEST_PROBDIST_GENDEP))
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@


# Build rule for the main test module
$(OBJDIR)maintest$(OBJSUFFIX) : $(TESTDIR)maintest.cpp
	$(CPP) $(CFLAG) $(CPPFLAGS) $(APPINCFLAG) $(MACROS) $< $(OFLAG) $@


# Build (link) rule for the final test application
$(TARGET) : $(OBJDIR) $(BUILDDIR) $(TARGETOBJ) $(TEST_LINKOBJS)
	$(LINKER) $(LDFLAGS) $(TARGETOBJ) $(TEST_LINKOBJS) $(OFLAG) $@ 


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
