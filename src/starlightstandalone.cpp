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
#include "inputparameters.h"
#include "eventfilewriter.h"
#include "starlightstandalone.h"


using namespace std;


starlightStandalone::starlightStandalone()
	:	_configFileName       ("slight.in"),
		_dataFileName         ("slight.out"),
		_starlight            (0),
		_inputParameters      (0),
		_nmbEventsTot       (1),
		_nmbEventsPerFile(_nmbEventsTot)
{ }


starlightStandalone::~starlightStandalone()
{ }


bool
starlightStandalone::init()
{
	// read input parameters from config file
	_inputParameters = new inputParameters();
	if (!_inputParameters->init(_configFileName)) {
		printWarn << "problems initializing input parameters. cannot initialize starlight.";
		return false;
	}

	// get the number of events
	// for now we write everything to one file
	_nmbEventsTot     = _inputParameters->nmbEvents();
	_nmbEventsPerFile = _nmbEventsTot;

	// create the starlight object
	_starlight = new starlight();
	// give starlight the input parameters
	_starlight->setInputParameters(_inputParameters);
	// initialize starlight
	_starlight->init();
    
	return true;
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
	fileWriter.open(_dataFileName);

	printInfo << "generating events:" << endl;
	unsigned int nmbEvents = 0;
	while (nmbEvents < _nmbEventsTot) {
		for (unsigned int iEvent = 0; (iEvent < _nmbEventsPerFile) && (nmbEvents < _nmbEventsTot);
		     ++iEvent, ++nmbEvents) {
			progressIndicator(iEvent, _nmbEventsTot, true, 4);
			upcEvent event = _starlight->produceEvent();
			fileWriter.writeEvent(event, iEvent);
		}
	}
	fileWriter.close();

	return true;
}
