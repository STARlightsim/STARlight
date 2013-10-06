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
// $Rev:: 102                         $: revision of last commit
// $Author:: odjuvsla                 $: author of last commit
// $Date:: 2012-10-22 16:25:54 -0500 #$: date of last commit
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
		_eventDataFileName("slight.out"),
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
	// read input parameters from config file
	inputParametersInstance.configureFromFile(_configFileName);
	if (!inputParametersInstance.init()) {
		printWarn << "problems initializing input parameters. cannot initialize starlight." << endl;
		return false;
	}

	// get the number of events
	// for now we write everything to one file
	_nmbEventsTot     = inputParametersInstance.nmbEvents();
	_nmbEventsPerFile = _nmbEventsTot;

	// create the starlight object
	_starlight = new starlight();
	
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
	fileWriter.writeFullPythiaInfo(inputParametersInstance.pythiaFullEventRecord());
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
   double rap1 = acosh(inputParametersInstance.beam1LorentzGamma());
   double rap2 = -acosh(inputParametersInstance.beam2LorentzGamma());
   double boost = (rap1+rap2)/2.;

   event.boost(boost);
}

