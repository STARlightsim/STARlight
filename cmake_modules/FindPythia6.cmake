
FIND_LIBRARY(PYTHIA6_LIBRARY NAMES pythia6 PATHS /usr/lib /usr/lib64 /usr/local/lib $ENV{PYTHIADIR})

IF (PYTHIA6_LIBRARY)
   SET(PYTHIA6_FOUND TRUE)
ENDIF (PYTHIA6_LIBRARY)

IF (PYTHIA6_FOUND)
   IF (NOT Pythia6_FIND_QUIETLY)
      MESSAGE(STATUS "Found Pythia6: ${PYTHIA6_LIBRARY}")
   ENDIF (NOT Pythia6_FIND_QUIETLY)
ELSE (PYTHIA6_FOUND)
   IF (Pythia6_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Pythia6. We search first in the normal library paths, then in $PYTHIADIR")
   ELSE(Pythia6_FIND_REQUIRED)
      IF(NOT Pythia6_FIND_QUIETLY)
	 MESSAGE(STATUS "Could not find Pythia6. We search first in the normal library paths, then in $PYTHIADIR")
      ENDIF(NOT Pythia6_FIND_QUIETLY)
   ENDIF (Pythia6_FIND_REQUIRED)
   
ENDIF (PYTHIA6_FOUND)

