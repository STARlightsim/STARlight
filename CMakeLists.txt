###########################################################################
#
#    Copyright 2010
#
#    This file is part of Starlight.
#
#    Starlight is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#	  
#    Starlight is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#	  
#    You should have received a copy of the GNU General Public License
#    along with Starlight. If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################
#
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      Starlight build file
#
#
###########################################################################


# check if cmake has the required version
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)


# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
# the same effect can be achieved by calling 'make VERBOSE=1'


# define project
project(starlight)


# load some common cmake macros
# set path, where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
message(STATUS "Using cmake module path '${CMAKE_MODULE_PATH}'")
include(CommonMacros)


# force out-of-source builds.
enforce_out_of_source_build()


# warn user if system is not UNIX
if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()
message(STATUS "Detected host system '${CMAKE_HOST_SYSTEM_NAME}' version '${CMAKE_HOST_SYSTEM_VERSION}' architecture '${CMAKE_HOST_SYSTEM_PROCESSOR}'")
message(STATUS "Compiling for system '${CMAKE_SYSTEM_NAME}' version '${CMAKE_SYSTEM_VERSION}' architecture '${CMAKE_SYSTEM_PROCESSOR}'")


# define build types
# set a default build type for single-configuration CMake generators, if no build type is set.
set(CMAKE_BUILD_TYPE Debug)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type was specified. Setting build type to 'Release'.")
  set(CMAKE_BUILD_TYPE Release)
endif()
# common compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")
# flags for specific build types
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_LDFLAGS_DEBUG "-g")
# report global build settings
message(STATUS "Using CXX compiler '${CMAKE_CXX_COMPILER}'")
message(STATUS "Using CXX general compiler flags '${CMAKE_CXX_FLAGS}'")
foreach(_BUILD_TYPE "DEBUG" "MINSIZEREL" "RELEASE" "RELWITHDEBINFO")
  message(STATUS "Using CXX compiler flags '${CMAKE_CXX_FLAGS_${_BUILD_TYPE}}' for build type ${_BUILD_TYPE}")
endforeach()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'")


# redirect output files
#set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib")
message(STATUS "Using library output path '${LIBRARY_OUTPUT_PATH}'")
#set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")
message(STATUS "Using executable output path '${EXECUTABLE_OUTPUT_PATH}'")


# make CMAKE_SOURCE_DIR accessible in source code via predefined macro CMAKE_SOURCE_DIR
if(CMAKE_SOURCE_DIR)
  add_definitions(-D'CMAKE_SOURCE_DIR=\"${CMAKE_SOURCE_DIR}\"')
else()
  add_definitions(-D'CMAKE_SOURCE_DIR=\"\"')
endif()


# make SVN version string accessible in source code via predefined macro SVN_VERSION
find_package(Subversion)
if(Subversion_FOUND)
  # unfortunately CMAKE only parses 'svn info'
  find_program(SVNVERSION_EXECUTABLE
    svnversion
		)
  if(NOT SVNVERSION_EXECUTABLE)
    message(STATUS "Could not find subversion command 'svnversion'. Repository version unknown.")
  else()
    execute_process(
      COMMAND ${SVNVERSION_EXECUTABLE} "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE SVN_VERSION
      RESULT_VARIABLE _SVNVERSION_RETURN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT ${_SVNVERSION_RETURN})
      message(STATUS "Subversion repository revision is '${SVN_VERSION}'")
    else()
      message(STATUS "Error running 'svnversion'. Repository version unknown.")
      set(SVN_VERSION "")
    endif()
  endif()
else()
  message(STATUS "Could not find subversion installation. Repository version unknown.")
endif()
if(SVN_VERSION)
  add_definitions(-D'SVN_VERSION=\"${SVN_VERSION}\"')
else()
  add_definitions(-D'SVN_VERSION=\"\"')
endif()


# setup doxygen
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  message(WARNING "Cannot find Doxygen. No HTML documentation will be generated.")
else()
  set(DOXYGEN_TARGET  "doxygen")
  set(DOXYGEN_DOC_DIR "${CMAKE_SOURCE_DIR}/doxygen")
  set(DOXYGEN_CONF    "${CMAKE_SOURCE_DIR}/starlightDoxyfile.conf")
  message(STATUS "Run 'make ${DOXYGEN_TARGET}' to create Doxygen documentation files in '${DOXYGEN_DOC_DIR}'")
  add_custom_target(${DOXYGEN_TARGET}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF}
    DEPENDS ${DOXYGEN_CONF}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		)
endif()


# setup Pythia 8
option (ENABLE_PYTHIA  "Enable compilation against pythia (necessary for certain processes)" OFF) 
if(ENABLE_PYTHIA)
  find_package(Pythia8)
  if(PYTHIA8_FOUND)
    find_package(LHAPDF REQUIRED) # implemented for dummy version in Pythia8 
    option(ENABLE_PYTHIA  "Should we use the Pythia8 library" ON) 
  else() 
    option(ENABLE_PYTHIA  "Should we use the Pythia8 library" OFF) 
  endif()
endif()


# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_SOURCE_DIR}/include
	#	${PROJECT_BINARY_DIR}
	${PYTHIA8_INCLUDE_DIR}
	)
include_directories(${INCLUDE_DIRECTORIES})

# Set our source files, include the generated dictionary
set(SOURCES
  src/bessel.cpp
  src/beam.cpp
  src/inputParameters.cpp
  src/beambeamsystem.cpp
  src/starlightparticle.cpp
  src/gammaaluminosity.cpp
  src/randomgenerator.cpp
  src/nucleus.cpp
  src/eventchannel.cpp
  src/gammaavm.cpp
  src/psifamily.cpp
  src/gammagammasingle.cpp
  src/photonNucleusCrossSection.cpp
  src/wideResonanceCrossSection.cpp
  src/narrowResonanceCrossSection.cpp
  src/readinluminosity.cpp
  src/twophotonluminosity.cpp
  src/gammagammaleptonpair.cpp
  src/starlight.cpp
  src/upcevent.cpp
  src/vector3.cpp
  src/lorentzvector.cpp
  src/filewriter.cpp
  src/eventfilewriter.cpp
  src/starlightparticlecodes.cpp
  src/starlightStandalone.cpp
  src/nBodyPhaseSpaceGen.cpp
  src/inputParser.cpp
  )
if(ENABLE_PYTHIA) 
  set (SOURCES
		${SOURCES}
		src/PythiaStarlight.cpp
		)
endif()


# add Starlight library to the build system
set(THIS_LIB "Starlib")
add_library(${THIS_LIB} STATIC ${SOURCES})
#make_shared_library("${THIS_LIB}" "${SOURCES}"
#	"${PYTHIA8_LIBRARY}"
#	"${LHAPDF_LIBRARIES}"
#)


# add starlight executable to the build system
add_executable(starlight src/main.cpp)
target_link_libraries(starlight Starlib ${PYTHIA8_LIBRARY} ${LHAPDF_LIBRARIES})
#make_executable(starlight src/main.cpp ${THIS_LIB})


message(STATUS "Cmake did not find any errors. run 'make' to build the project.")
message(STATUS "On multi-core machines 'make -j#', where # is the number of parallel jobs, can speedup compilation considerably.")
