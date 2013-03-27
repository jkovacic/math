#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=Cygwin_4.x-Windows
CND_DLIB_EXT=dll
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/IntFactorizationException.o \
	${OBJECTDIR}/PolynomialException.o \
	${OBJECTDIR}/QuaternionException.o \
	${OBJECTDIR}/LinearEquationSolverException.o \
	${OBJECTDIR}/IntFactorization.o \
	${OBJECTDIR}/maintest.o \
	${OBJECTDIR}/Rational.o \
	${OBJECTDIR}/IntCombinatorics.o \
	${OBJECTDIR}/MatrixException.o \
	${OBJECTDIR}/CombinatoricsException.o \
	${OBJECTDIR}/CurveFittingException.o \
	${OBJECTDIR}/RationalException.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/math.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/math.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/math ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/IntFactorizationException.o: IntFactorizationException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/IntFactorizationException.o IntFactorizationException.cpp

${OBJECTDIR}/PolynomialException.o: PolynomialException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/PolynomialException.o PolynomialException.cpp

${OBJECTDIR}/QuaternionException.o: QuaternionException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/QuaternionException.o QuaternionException.cpp

${OBJECTDIR}/LinearEquationSolverException.o: LinearEquationSolverException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/LinearEquationSolverException.o LinearEquationSolverException.cpp

${OBJECTDIR}/IntFactorization.o: IntFactorization.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/IntFactorization.o IntFactorization.cpp

${OBJECTDIR}/maintest.o: maintest.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/maintest.o maintest.cpp

${OBJECTDIR}/Rational.o: Rational.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/Rational.o Rational.cpp

${OBJECTDIR}/IntCombinatorics.o: IntCombinatorics.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/IntCombinatorics.o IntCombinatorics.cpp

${OBJECTDIR}/MatrixException.o: MatrixException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixException.o MatrixException.cpp

${OBJECTDIR}/CombinatoricsException.o: CombinatoricsException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/CombinatoricsException.o CombinatoricsException.cpp

${OBJECTDIR}/CurveFittingException.o: CurveFittingException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/CurveFittingException.o CurveFittingException.cpp

${OBJECTDIR}/RationalException.o: RationalException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/RationalException.o RationalException.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/math.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
