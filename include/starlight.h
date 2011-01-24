///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef STARLIGHT_H
#define STARLIGHT_H
#define Starlight_VERSION_MAJOR 1
#define Starlight_VERSION_MINOR 0


#include <string>

#include "upcevent.h"


class TDatabasePDG;
class beamBeamSystem;
class inputParameters;
class beam;
class eventChannel;


class starlight
{
   public:
      
      starlight();
      
      ~starlight();
      
      int init();

      upcEvent produceEvent();
      
      std::string getConfigFileName() const { return _configFileName; }
   
      void setInputParameters(inputParameters *inputParams) { _inputParameters = inputParams; }
      int getNTries(); 
      int getNSuccess(); 
   private:
      
      bool checkForLuminosityTable();

      int createEventChannel();
      
      inputParameters* _inputParameters;
      beam* _beam0;
      beam* _beam1;
      beamBeamSystem* _beamSystem;
      eventChannel* _eventChannel;
      std::string _configFileName;
      int _numberOfEventsPerFile;
      unsigned long long _numberOfEventsToGenerate;
      std::string _standardFilename;
      bool _isInitialised;
};


#endif  // STARLIGHT_H
