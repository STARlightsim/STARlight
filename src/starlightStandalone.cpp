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


#include <iostream>

#include "reportingUtils.h"
#include "starlight.h"
#include "inputParameters.h"
#include "eventfilewriter.h"
#include "starlightStandalone.h"


using namespace std;


starlightStandalone::starlightStandalone()
	:	_configFileName   ("slight.in"),
//		_eventDataFileName("slight.out"),
		_starlight        (0),
		_inputParameters  (0),
		_nmbEventsTot     (1),
		_nmbEventsPerFile (_nmbEventsTot)
{ }


starlightStandalone::~starlightStandalone()
{ }


bool
starlightStandalone::init()
{
	_inputParameters = new inputParameters();
	// read input parameters from config file
	_inputParameters->configureFromFile(_configFileName);
	if (!_inputParameters->init()) {
		printWarn << "problems initializing input parameters. cannot initialize starlight." << endl;
		return false;
	}

	// copy input file to one with baseFileName naming scheme
        std::string inputCopyName, _baseFileName;
        _baseFileName = _inputParameters->baseFileName();
         inputCopyName = _baseFileName +".in";

        ofstream inputCopyFile;
	inputCopyFile.open(inputCopyName.c_str());

        std::ifstream infile(_configFileName.c_str());
        if ((!infile) || (!infile.good()))
        {
          return -1;
        }

        int lineSize = 256;
        char tmp[lineSize];
        while (!infile.getline(tmp, lineSize).eof())
        {
     	 cout << tmp << endl;
         inputCopyFile << tmp << endl;
        }
        inputCopyFile.close();


	// get the number of events
	// for now we write everything to one file
	_nmbEventsTot     = _inputParameters->nmbEvents();
	_nmbEventsPerFile = _nmbEventsTot;

	// create the starlight object
	_starlight = new starlight();

        // give starlight the input parameters
        _starlight->setInputParameters(_inputParameters);
	
	// initialize starlight
	return _starlight->init();
}


bool
starlightStandalone::run()
{
	if (!_starlight) {
		printWarn << "null pointer to starlight object. make sure that init() was called. "
		          << "cannot generate events." << endl;
		return false;
	}

	// open output file
	eventFileWriter fileWriter;
	fileWriter.writeFullPythiaInfo(_inputParameters->pythiaFullEventRecord());
        _baseFileName = _inputParameters->baseFileName();
        _eventDataFileName = _baseFileName +".out";
	fileWriter.open(_eventDataFileName);

	printInfo << "generating events:" << endl;
	unsigned int nmbEvents = 0;
	while (nmbEvents < _nmbEventsTot) {
		for (unsigned int iEvent = 0; (iEvent < _nmbEventsPerFile) && (nmbEvents < _nmbEventsTot);
		     ++iEvent, ++nmbEvents) {
			progressIndicator(iEvent, _nmbEventsTot, true, 4);
			upcEvent event = _starlight->produceEvent();
			// Boost event from CMS system to lab system
			boostEvent(event);
			fileWriter.writeEvent(event, iEvent);
		}
	}
	printInfo << "number of attempts = " << _starlight->nmbAttempts() << ", "
	          << "number of accepted events = " << _starlight->nmbAccepted() << endl;
	fileWriter.close();

	return true;
}
void starlightStandalone::boostEvent(upcEvent &event)
{
  
  // Should probably move this calculation to inputParameters (and remove from bbs)
   // Calculate CMS boost 
   double rap1 = acosh(_inputParameters->beam1LorentzGamma());
   double rap2 = -acosh(_inputParameters->beam2LorentzGamma());
   double boost = (rap1+rap2)/2.;

   event.boost(boost);
}

