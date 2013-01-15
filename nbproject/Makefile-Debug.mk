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
FC=
AS=as

# Macros
CND_PLATFORM=Cygwin_4.x-Windows
CND_CONF=Debug
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/QuaternionException.o \
	${OBJECTDIR}/maintest.o \
	${OBJECTDIR}/Rational.o \
	${OBJECTDIR}/MatrixException.o \
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
	"${MAKE}"  -f nbproject/Makefile-Debug.mk dist/Debug/Cygwin_4.x-Windows/math.exe

dist/Debug/Cygwin_4.x-Windows/math.exe: ${OBJECTFILES}
	${MKDIR} -p dist/Debug/Cygwin_4.x-Windows
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/math ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/QuaternionException.o: QuaternionException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Werror -MMD -MP -MF $@.d -o ${OBJECTDIR}/QuaternionException.o QuaternionException.cpp

${OBJECTDIR}/maintest.o: maintest.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Werror -MMD -MP -MF $@.d -o ${OBJECTDIR}/maintest.o maintest.cpp

${OBJECTDIR}/Rational.o: Rational.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Werror -MMD -MP -MF $@.d -o ${OBJECTDIR}/Rational.o Rational.cpp

${OBJECTDIR}/MatrixException.o: MatrixException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Werror -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixException.o MatrixException.cpp

${OBJECTDIR}/RationalException.o: RationalException.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Werror -MMD -MP -MF $@.d -o ${OBJECTDIR}/RationalException.o RationalException.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Debug
	${RM} dist/Debug/Cygwin_4.x-Windows/math.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
