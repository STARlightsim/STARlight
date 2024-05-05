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
# $Rev:: 247                         $: revision of last commit
# $Author:: butter                   $: author of last commit
# $Date:: 2016-03-05 01:59:24 +0100 #$: date of last commit
#
# Description:
#      Starlight build file
#
#
###########################################################################


# check if cmake has the required version
# note: Compatibility with CMake < 3.5 will be removed from a future version of CMake.
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# set verbosity
set(CMAKE_VERBOSE_MAKEFILE 0)  # if set to 1 compile and link commands are displayed during build
# the same effect can be achieved by calling 'make VERBOSE=1'


# The version number. 9999 indicates trunk
set (Starlight_VERSION_MAJOR 9999) 
set (Starlight_VERSION_MINOR 1)
set (Starlight_VERSION_MINOR_MINOR 0)

# define project
project(starlight)
#find_package (Threads)

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

# Set the C++ standard
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    message(STATUS "Compiler ${CMAKE_CXX_COMPILER} with C++11 support - OK.")
    set(CMAKE_CXX_STANDARD 11)
else()
    message(WARNING "Compiler ${CMAKE_CXX_COMPILER} w/o C++11 support. Trying with C++98 standard.")
    set(CMAKE_CXX_STANDARD 98)
endif()

# common compiler flags
if (CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wno-error=unused-but-set-variable -Wno-error=unused-but-set-parameter")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")
endif()

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
    set(optionalLibs ${optionalLibs} ${PYTHIA8_LIBRARY})
#    find_package(LHAPDF REQUIRED) # implemented for dummy version in Pythia8
#    set(optionalLibs ${optionalLibs} ${LHAPDF_LIBRARIES})#liblhapdfdummy
#    removed from v8.2, if you want to add your own lhapdf set, uncomment
#    find_package + set  - jb 05192015
    option(ENABLE_PYTHIA  "Should we use the Pythia8 library" ON) 
  else() 
    option(ENABLE_PYTHIA  "Should we use the Pythia8 library" OFF) 
  endif()
endif()

# setup Pythia 6
option (ENABLE_PYTHIA6  "Enable compilation against pythia 6 (necessary for certain processes)" OFF) 
if(ENABLE_PYTHIA6)
   find_package(Pythia6 REQUIRED)
   if(PYTHIA6_FOUND)
      set(optionalLibs ${optionalLibs} ${PYTHIA6_LIBRARY})
      option (ENABLE_PYTHIA6  "Enable compilation against pythia 6 (necessary for certain processes)" ON) 
      include_directories(pythia6)
   else(PYTHIA6_FOUND) 
      option (ENABLE_PYTHIA6  "Enable compilation against pythia 6 (necessary for certain processes)" OFF) 
   endif(PYTHIA6_FOUND)
endif()

# setup DPMJET
option (ENABLE_DPMJET  "Enable compilation against DPMJet" OFF)
if(ENABLE_DPMJET)
   find_package(DPMJet REQUIRED)
   if(DPMJET_FOUND)
      option (ENABLE_DPMJET  "Enable compilation against DPMJet" ON)       
   else(DPMJET_FOUND)
      option (ENABLE_DPMJET  "Enable compilation against DPMJet" OFF)
   endif(DPMJET_FOUND)
 endif(ENABLE_DPMJET)

 option (ENABLE_HEPMC3  "Enable compilation against hepmc3 (necessary for hepmc3 output)" OFF) 
if(ENABLE_HEPMC3)
  find_package(HepMC3 HINTS ${HepMC3_DIR} ${HepMC3} )
  if(HepMC3_FOUND)
    message ( "-- HepMC3 library found at ${HEPMC3_LIB} --> Adding HepMC sources and headers" )
    set(optionalLibs ${optionalLibs} ${HEPMC3_LIB})
    option(ENABLE_HEPMC3  "Should we use the HepMC3 library" ON) 
  else()
    message ( WARNING "Hepmc3 not found.")
    option(ENABLE_HEPMC3  "Should we use the HepMC3 library" OFF) 
  endif()
endif()



# set include directories
set(INCLUDE_DIRECTORIES
	${CMAKE_SOURCE_DIR}/include
	${PROJECT_BINARY_DIR}
	${PYTHIA8_INCLUDE_DIR}#uncommented 05192015
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
  src/gammagammasingle.cpp
  src/photonNucleusCrossSection.cpp
  src/wideResonanceCrossSection.cpp
  src/narrowResonanceCrossSection.cpp
  src/readinluminosity.cpp
  src/twophotonluminosity.cpp
  src/gammagammaleptonpair.cpp
  src/starlight.cpp
  src/upcevent.cpp
  src/upcXevent.cpp
  src/vector3.cpp
  src/lorentzvector.cpp
  src/filewriter.cpp
  src/eventfilewriter.cpp
  src/starlightparticlecodes.cpp
  src/starlightStandalone.cpp
  src/nBodyPhaseSpaceGen.cpp
  src/inputParser.cpp
  src/incoherentPhotonNucleusLuminosity.cpp
  src/incoherentVMCrossSection.cpp
  )
if(ENABLE_PYTHIA) 
  set (SOURCES
		${SOURCES}
		#src/PythiaStarlight.cpp
                src/pythiadecayer.cpp
		)
  include_directories(${PYTHIA8_INCLUDE_DIR})
endif()
if(ENABLE_PYTHIA6)
  set (SOURCES 
		${SOURCES}
		src/starlightpythia.cpp
		src/spectrum.cpp
		src/spectrumprotonnucleus.cpp
		)
endif()

if (ENABLE_HEPMC3)
  set (SOURCES
		${SOURCES}
                src/hepmc3writer.cpp
		)
  include_directories(${HEPMC3_INCLUDE_DIR})
  add_definitions( -DHEPMC3_ON )
endif()

if(ENABLE_DPMJET)
  set (SOURCES 
		${SOURCES}
		src/starlightdpmjet.cpp
		src/spectrum.cpp
		src/spectrumprotonnucleus.cpp
		)
endif()

# add Starlight library to the build system
set(THIS_LIB "Starlib")
add_library(${THIS_LIB} STATIC ${SOURCES})
#make_shared_library("${THIS_LIB}" "${SOURCES}"
#	"${PYTHIA8_LIBRARY}"
#	"${LHAPDF_LIBRARIES}"
#)

if(ENABLE_DPMJET)
  enable_language(Fortran)
  set(DPMJET_LIB "DpmJetLib")
  
  # ###
  # old:
  #
  # message(STATUS "DPMJet objects: ${DPMJET_OBJECTS}")
  # add_library(${DPMJET_LIB} STATIC dpmjet/dpmjetint.f ${DPMJET_OBJECTS})
  # ###

  # ###
  # new (Luan Arbeletche):
  #
  message(STATUS "DPMJet external library: ${DPMJET_EXTERNAL_LIB}")
  add_library(${DPMJET_LIB} STATIC dpmjet/dpmjetint.f)
  #
  # ###

  target_include_directories(${DPMJET_LIB} PRIVATE $ENV{DPMJETDIR}/include/dpmjet)
  set(optionalLibs ${optionalLibs} ${DPMJET_LIB})
  set(optionalLibs ${optionalLibs} "gfortran")
  set(optionalLibs ${optionalLibs} ${DPMJET_EXTERNAL_LIB})
endif()

if(ENABLE_PYTHIA6)
  enable_language(Fortran)
endif()

# add starlight executable to the build system
add_executable(starlight src/main.cpp)
#target_link_libraries(starlight Starlib ${PYTHIA8_LIBRARY} ${LHAPDF_LIBRARIES} ${PYTHIA6_LIBRARY} ${DPMJET_LIB} )
target_link_libraries(starlight Starlib ${optionalLibs})# ${CMAKE_THREAD_LIBS_INIT}) 

#make_executable(starlight src/main.cpp ${THIS_LIB})

configure_file (
  "${PROJECT_SOURCE_DIR}/starlightconfig.h.in"
  "${PROJECT_BINARY_DIR}/starlightconfig.h"
  )
# Erase xsec values in case changes in code affects the xsec, executed during make process
add_custom_command (TARGET Starlib POST_BUILD COMMAND touch ARGS slight.txt)
add_custom_command (TARGET Starlib POST_BUILD COMMAND cp ARGS slight.txt slight.txt.bak)
add_custom_command (TARGET Starlib POST_BUILD COMMAND echo ARGS '' > slight.txt )

# place things in the right place - add an install target
# use with 'make install' to configure use "cmake <src_dir_where_this_file> -DCMAKE_INSTALL_PREFIX=<install_dir>  -DCMAKE_BUILD_TYPE=Release"
include(GNUInstallDirs)
install(TARGETS starlight Starlib
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

message(STATUS "Cmake did not find any errors. run 'make' to build the project.")
message(STATUS "On multi-core machines 'make -j#', where # is the number of parallel jobs, can speedup compilation considerably.")
