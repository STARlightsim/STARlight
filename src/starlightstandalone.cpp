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

#include "starlightstandalone.h"
#include "starlight.h"
#include "inputparameters.h"
#include "eventfilewriter.h"


starlightStandalone::starlightStandalone() :
        _starlight(0)
        ,_inputParameters(0)
        ,_configFileName(std::string("slight.in"))
        ,_numberOfEvents(1)
        ,_numberOfEventsPerFile(_numberOfEvents)
        ,_fileName(std::string("slight.out"))
{ }


starlightStandalone::~starlightStandalone()
{ }


int starlightStandalone::init()
{
    // Reading input parameters from config file
    _inputParameters = new inputParameters();
    _inputParameters->init(_configFileName);

    // Get the number of events, for now we write everything to one file
    _numberOfEvents = _inputParameters->getNumberOfEvents();
    _numberOfEventsPerFile = _numberOfEvents;

    // Creating the starlight object
    _starlight = new starlight();

    // Give starlight the input parameters
    _starlight->setInputParameters(_inputParameters);

    // Initialising starlight
    _starlight->init();
    
    return 0;
}


int starlightStandalone::run()
{
    int ntotev = 0;

    eventFileWriter fw;

    fw.open(_fileName);

    while (ntotev < _numberOfEvents)
    {
        for (int nev = 0; nev < _numberOfEventsPerFile && ntotev < _numberOfEvents; nev++, ntotev++)
        {
	   if(ntotev%10000 == 0)  std::cout << "PRODUCING EVENT #: " << nev << std::endl;
            upcEvent event = _starlight->produceEvent();
            //std::cout << "Number of tracks in event: " << event.getParticles()->size() << std::endl;
            fw.writeEvent(event, nev);
        }
    }
    
    fw.close();

    return 0;
}
