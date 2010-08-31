
#ifndef STARLIGHT_H
#define STARLIGHT_H
#include <string>
#include "upcevent.h"

class TDatabasePDG;
class Beambeamsystem;
class Inputparameters;
class Beam;
class Eventchannel;

class Starlight
{
   public:
      
      Starlight();
      
      ~Starlight();
      
      int Init();

      UPCEvent ProduceEvent();
      
      std::string GetConfigFileName() const { return fConfigFileName; }
   
      void SetInputParameters(Inputparameters *inputParams) { fInputParameters = inputParams; }
   
   private:
      
      bool CheckForLuminosityTable();

      int CreateEventChannel();
      
      Inputparameters *fInputParameters;
      
      Beam *fBeam0;
      
      Beam *fBeam1;
      
      Beambeamsystem *fBeamSystem;
      
      Eventchannel *fEventChannel;
      
      std::string fConfigFileName;
      
      int fNumberOfEventsPerFile;
      
      unsigned long long fNumberOfEventsToGenerate;
      
      std::string fStandardFilename;
      
      bool fIsInitialised;

};

#define Starlight_VERSION_MAJOR 1
#define Starlight_VERSION_MINOR 0


#endif // STARLIGHT_H
