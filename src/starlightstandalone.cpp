// starlightstandalone.cpp
/*
 * $Id: starlightstandalone.cpp,v 1.0 2010/07/04  $
 *
 *
 * /author 
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * $Log: $
 */

#include "starlightstandalone.h"
#include "starlight.h"
#include "inputparameters.h"
#include "eventfilewriter.h"
#include <iostream>

StarlightStandalone::StarlightStandalone() :
        fStarlight(0)
        ,fInputParameters(0)
        ,fConfigFileName(std::string("slight.in"))
        ,fNumberOfEvents(1)
        ,fNumberOfEventsPerFile(fNumberOfEvents)
        ,fFileName(std::string("slight.out"))
{

}

StarlightStandalone::~StarlightStandalone()
{

}

int StarlightStandalone::Init()
{

    // Reading input parameters from config file
    fInputParameters = new Inputparameters();
    fInputParameters->Init(fConfigFileName);

    // Get the number of events, for now we write everything to one file
    fNumberOfEvents = fInputParameters->getnumberofevents();
    fNumberOfEventsPerFile = fNumberOfEvents;

    // Creating the starlight object
    fStarlight = new Starlight();

    // Give starlight the input parameters
    fStarlight->SetInputParameters(fInputParameters);

    // Initialising starlight
    fStarlight->Init();
    
    return 0;
}

int StarlightStandalone::Run()
{
    int ntotev = 0;

    EventFileWriter fw;

    fw.Open(fFileName);

    while (ntotev < fNumberOfEvents)
    {
        for (int nev = 0; nev < fNumberOfEventsPerFile && ntotev < fNumberOfEvents; nev++, ntotev++)
        {
	   if(ntotev%10000 == 0)  std::cout << "PRODUCING EVENT #: " << nev << std::endl;
            UPCEvent event = fStarlight->ProduceEvent();
            //std::cout << "Number of tracks in event: " << event.GetParticles()->size() << std::endl;
            fw.WriteEvent(event, nev);
        }
    }
    
    fw.Close();

    return 0;

}
