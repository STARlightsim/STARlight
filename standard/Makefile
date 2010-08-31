
#
# The PROGRAM macro defines the name of the program or project.  It
# allows the program name to be changed by editing in only one
# location
#

PROGRAM = starlight

# The INCLUDEDIRS macro contains a list of include directories
# to pass to the compiler so it can find necessary header files.
#
# The LIBDIRS macro contains a list of library directories
# to pass to the linker so it can find necessary libraries.
#
# The LIBS macro contains a list of libraries that the the
# executable must be linked against.
#

INCLUDEDIRS = \
	-I/home/gorbunov/pythia/pythia8120/include \
	-I/usr/local/root/include \
	-I../include 	


LIBDIRS = \
	-L/usr/local/root/lib \
	-L/home/gorbunov/pythia/pythia8120/lib/archive


LIBS =  -lCore -lCint -lRIO  -lGpad  -lRint -lGui -lPhysics -lMathCore\
        -lGraf -lGraf3d -lHist -lHtml -lMatrix  -lMinuit -lNet \
        -lPostscript -lProof -lTree -lThread  -lVMC -lGeom \
	-lEG -ldl -lpthread -lTreePlayer -lProofPlayer -lpythia8 -llhapdfdummy


#
# The CXXSOURCES macro contains a list of source files.
#
# The CXXOBJECTS macro converts the CXXSOURCES macro into a list
# of object files.
#
# The CXXFLAGS macro contains a list of options to be passed to
# the compiler.  Adding "-g" to this line will cause the compiler
# to add debugging information to the executable.
#
# The CXX macro defines the C++ compiler.
#
# The LDFLAGS macro contains all of the library and library
# directory information to be passed to the linker.
#
OTHERDIR= ../src
CXXSOURCES = $(OTHERDIR)/main.cpp $(OTHERDIR)/inputparameters.cpp $(OTHERDIR)/nucleus.cpp $(OTHERDIR)/bessel.cpp $(OTHERDIR)/beam.cpp $(OTHERDIR)/beambeamsystem.cpp $(OTHERDIR)/twophotonluminosity.cpp $(OTHERDIR)/upcevent.cpp $(OTHERDIR)/vector3.cpp $(OTHERDIR)/lorentzvector.cpp $(OTHERDIR)/eventfilewriter.cpp $(OTHERDIR)/filewriter.cpp $(OTHERDIR)/eventchannel.cpp $(OTHERDIR)/gammagammaleptonpair.cpp $(OTHERDIR)/readinluminosity.cpp   $(OTHERDIR)/randomgenerator.cpp   $(OTHERDIR)/gammagammasingle.cpp $(OTHERDIR)/gammaavm.cpp $(OTHERDIR)/gammaaluminosity.cpp $(OTHERDIR)/gammaacrosssection.cpp  $(OTHERDIR)/psifamily.cpp $(OTHERDIR)/starlightparticle.cpp $(OTHERDIR)/starlightparticlecodes.cpp $(OTHERDIR)/starlight.cpp $(OTHERDIR)/starlightstandalone.cpp # list of source files 
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)  # expands to list of object files
CXXFLAGS =  $(INCLUDEDIRS)
CXX = g++

LDFLAGS = $(LIBDIRS) $(LIBS)

#
# Default target: the first target is the default target.
# Just type "make -f Makefile.Linux" to build it.
#

all: $(PROGRAM)

#
# Link target: automatically builds its object dependencies before
# executing its link command.
#

$(PROGRAM): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS) $(LDFLAGS)

#
# Object targets: rules that define objects, their dependencies, and
# a list of commands for compilation.
#
#
#basic_sample.o: basic_sample.cpp basic_sample.h
#	$(CXX) $(CXXFLAGS) -c -o basic_sample.o basic_sample.cpp
#
#
# Clean target: "make -f Makefile.Linux clean" to remove unwanted objects and executables.
#

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM) 

