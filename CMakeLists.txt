cmake_minimum_required(VERSION 2.6)

project(starlight)

# The version number.
# moved to starlight.h / the same is true for enable Pythia 
#set (Starlight_VERSION_MAJOR 4)
#set (Starlight_VERSION_MINOR 0)

# The options we have (to use: cmake <pathtosourcetree> -D<OPTION>=<VALUE>)
option (ENABLE_PYTHIA  "Enable compilation against pythia (necessary for certain processes)" OFF) 

# configure a header file to pass some of the CMake settings
# to the source code

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")

set(PYTHIAINCDIR $ENV{PYTHIADIR}/include)
set(PYTHIALIBDIR $ENV{PYTHIADIR}/lib/archive)

set(optionalLibs)
set(optionalIncludes)

if(ENABLE_PYTHIA)
   find_package(Pythia8)
   if(PYTHIA8_FOUND)
      set(optionalLibs ${optionalLibs} ${PYTHIA8_LIBRARY})
      set(optionalIncludes ${optionalIncludes} ${PYTHIA8_INCLUDE_DIR})
      find_package(LHAPDF REQUIRED) # implemented for dummy version in Pythia8 
      set(optionalLibs ${optionalLibs} ${LHAPDF_LIBRARIES})
      option (ENABLE_PYTHIA  "Should we use the Pythia8 library" ON) 
   else(PYTHIA8_FOUND) 
      option (ENABLE_PYTHIA  "Should we use the Pythia8 library" OFF) 
   endif(PYTHIA8_FOUND)
endif(ENABLE_PYTHIA)

#configure_file (
#  "${PROJECT_SOURCE_DIR}/StarlightConfig.h.in"
#  "${PROJECT_BINARY_DIR}/StarlightConfig.h"
#  )

# Set the path to header files for which we want to create a ROOT dictionary
#set(ROOTDEPINCLUDEPATH -I${CMAKE_CURRENT_SOURCE_DIR}/include)
#Set the header files themselves
#set(ROOTDEPHEADERS Particle.h)
   
# Command for generating the dictionary
#add_custom_command(
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/StarlightDICT.cxx
#   COMMAND cd ${CMAKE_CURRENT_BINARY_DIR} && rootcint -f ${CMAKE_CURRENT_BINARY_DIR}/StarlightDICT.cxx -c ${ROOTDEPINCLUDEPATH} ${ROOTDEPHEADERS}
#)


# Set the include directories
# add the include directory in the source tree
# add the binary tree to the search path for include files so that we will find TutorialConfig.h
# add the ROOT include directory
# add the optional (Pythia8, LHAPDF) include directories
include_directories(include/ ${PROJECT_BINARY_DIR} ${optionalIncludes})

# Set our source files, include the generated dictionary
set(SRCS
   src/bessel.cpp
   src/beam.cpp
   src/inputparameters.cpp
   src/beambeamsystem.cpp
   src/starlightparticle.cpp
   src/gammaaluminosity.cpp
   src/randomgenerator.cpp
   src/nucleus.cpp
   src/eventchannel.cpp
   src/gammaavm.cpp
   src/psifamily.cpp
   src/gammagammasingle.cpp
   src/gammaacrosssection.cpp
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
   src/starlightstandalone.cpp
   src/nBodyPhaseSpaceGen.cpp
   )

if(ENABLE_PYTHIA) 
   set (SRCS ${SRCS} src/PythiaStarlight.cpp)
endif(ENABLE_PYTHIA)

# Add the Starlight library to the build system
add_library(Starlight STATIC ${SRCS})

add_executable(starlight src/main.cpp)
target_link_libraries(starlight Starlight ${LHAPDF_LIBRARIES} ${optionalLibs})

