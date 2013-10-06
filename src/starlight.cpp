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
#include <fstream>
#include <cstdlib>

#include "starlightconfig.h"

#ifdef ENABLE_PYTHIA
#include "PythiaStarlight.h"
#endif

#ifdef ENABLE_DPMJET
#include "starlightdpmjet.h"
#endif

#ifdef ENABLE_PYTHIA6
#include "starlightpythia.h"
#endif

#include "reportingUtils.h"
#include "inputParameters.h"
#include "eventchannel.h"
#include "gammagammaleptonpair.h"
#include "gammagammasingle.h"
#include "gammaavm.h"
#include "psifamily.h"
#include "twophotonluminosity.h"
#include "gammaaluminosity.h"
#include "incoherentPhotonNucleusLuminosity.h"
#include "upcevent.h"
#include "eventfilewriter.h"
#include "starlight.h"


using namespace std;
using namespace starlightConstants;


starlight::starlight() :
		_beam0                 (0),
		_beam1                 (0),
		_beamSystem            (0),
		_eventChannel          (0),
		_nmbEventsPerFile      (100),
		_nmbEventsToGenerate   (10),
		_configFileName        ("slight.in"),
		_eventDataFileName     ("slight.out"),
		_lumLookUpTableFileName("slight.txt"),
		_isInitialised         (false)
{ }


starlight::~starlight()
{ }


bool
starlight::init()
{
	if(Starlight_VERSION_MAJOR == 9999)
	{
		cout << "##################################" << endl
	     	<< " Initialising Starlight version: trunk..." << endl
	     	<< "##################################" << endl;
	}
	else
	{
		cout << "##################################" << endl
	     	<< " Initialising Starlight version: " << Starlight_VERSION_MAJOR << "."
	     	<< Starlight_VERSION_MINOR << "." << Starlight_VERSION_MINOR_MINOR << "..." << endl
	     	<< "##################################" << endl;
	}

	_nmbEventsPerFile    = inputParametersInstance.nmbEvents();  // for now we write only one file...

	_beamSystem = new beamBeamSystem;
	
// 	cout << "Created beam system with beam lorentz gamma: " << _beamSystem->beamLorentzGamma() << endl;

	// streamsize precision(15);
	cout.setf(ios_base::fixed,ios_base::floatfield);
	cout.precision(15);
	const bool lumTableIsValid = luminosityTableIsValid();
	bool createChannel = true;
	switch (inputParametersInstance.interactionType())	{
	case PHOTONPHOTON:
		if (!lumTableIsValid) {
			printInfo << "creating luminosity table for photon-photon channel" << endl;
			twoPhotonLuminosity(_beamSystem->beam1(), _beamSystem->beam2());
		}
		break;		
	case PHOTONPOMERONNARROW:  // narrow and wide resonances use
	case PHOTONPOMERONWIDE:    // the same luminosity function
		if (!lumTableIsValid) {
			printInfo << "creating luminosity table for coherent photon-Pomeron channel" <<endl;
			photonNucleusLuminosity lum(*_beamSystem);
		}
		break;
        case PHOTONPOMERONINCOHERENT:  // narrow and wide resonances use
                if (!lumTableIsValid) {
                        printInfo << "creating luminosity table for incoherent photon-Pomeron channel" << endl;
                        incoherentPhotonNucleusLuminosity lum(*_beamSystem);
                }
                break;
#ifdef ENABLE_DPMJET
	case PHOTONUCLEARSINGLE:
		createChannel = false;
		_eventChannel = new starlightDpmJet(*_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/DPMJET SINGLE" << std::endl;
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setSingleMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMinGammaEnergy(inputParametersInstance.minGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMaxGammaEnergy(inputParametersInstance.maxGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->init();
		break;
	case PHOTONUCLEARDOUBLE:
		createChannel = false;
		_eventChannel = new starlightDpmJet(*_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/DPMJET DOUBLE" << std::endl;
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setDoubleMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMinGammaEnergy(inputParametersInstance.minGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMaxGammaEnergy(inputParametersInstance.maxGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->init();
		break;
	case PHOTONUCLEARSINGLEPA:
		createChannel = false;
		_eventChannel = new starlightDpmJet(*_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/DPMJET SINGLE" << std::endl;
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setSingleMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setProtonMode();
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMinGammaEnergy(inputParametersInstance.minGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->setMaxGammaEnergy(inputParametersInstance.maxGammaEnergy());
		dynamic_cast<starlightDpmJet*>(_eventChannel)->init();
		break;
#endif
#ifdef ENABLE_PYTHIA6
	case PHOTONUCLEARSINGLEPAPY:
		createChannel = false;
		_eventChannel = new starlightPythia(*_beamSystem);
		std::cout << "CREATING PHOTONUCLEAR/PYTHIA SINGLE" << std::endl;
		dynamic_cast<starlightPythia*>(_eventChannel)->setSingleMode();
		dynamic_cast<starlightPythia*>(_eventChannel)->setMinGammaEnergy(inputParametersInstance.minGammaEnergy());
		dynamic_cast<starlightPythia*>(_eventChannel)->setMaxGammaEnergy(inputParametersInstance.maxGammaEnergy());
		dynamic_cast<starlightPythia*>(_eventChannel)->init(inputParametersInstance.pythiaParams(), inputParametersInstance.pythiaFullEventRecord());
		break;
#endif
	default:
		{
			printWarn << "unknown interaction type '" << inputParametersInstance.interactionType() << "'."
			          << " cannot initialize starlight." << endl;
			return false;
		}
	}
	
	if(createChannel)
	{
	  if (!createEventChannel())
		  return false;
	}

	_isInitialised = true;
	return true;
}


upcEvent
starlight::produceEvent()
{
	if (!_isInitialised) {
		printErr << "trying to generate event but Starlight is not initialised. aborting." << endl;
		exit(-1);
	}
	return _eventChannel->produceEvent();
}


bool
starlight::luminosityTableIsValid() const
{
	printInfo << "using random seed = " << inputParametersInstance.randomSeed() << endl;

	ifstream lumLookUpTableFile(_lumLookUpTableFileName.c_str());
	lumLookUpTableFile.precision(15);
	if ((!lumLookUpTableFile) || (!lumLookUpTableFile.good())) {
		// printWarn << "cannot open file '" << _lumLookUpTableFileName << "'" << endl;
		return false;
	}

	unsigned int beam1Z, beam1A, beam2Z, beam2A;
	double       beamLorentzGamma = 0;
	double       maxW = 0, minW = 0;
	unsigned int nmbWBins;
	double       maxRapidity = 0;
	unsigned int nmbRapidityBins;
	int          productionMode, beamBreakupMode;
	bool         interferenceEnabled = false;
	double       interferenceStrength = 0;
	bool         coherentProduction = false;
	double       incoherentFactor = 0, deuteronSlopePar = 0, maxPtInterference = 0;
	int          nmbPtBinsInterference;
	std::string  validationKey;
	if (!(lumLookUpTableFile
	      >> validationKey
	      >> beam1Z >> beam1A
	      >> beam2Z >> beam2A
	      >> beamLorentzGamma
	      >> maxW >> minW >> nmbWBins
	      >> maxRapidity >> nmbRapidityBins
	      >> productionMode
	      >> beamBreakupMode
	      >> interferenceEnabled >> interferenceStrength
	      >> coherentProduction >> incoherentFactor
	      >> deuteronSlopePar
	      >> maxPtInterference
	      >> nmbPtBinsInterference))
		// cannot read parameters from lookup table file
		return false;
		
	std::string validationKeyEnd;
	while(!lumLookUpTableFile.eof())
	{
	  lumLookUpTableFile >> validationKeyEnd; 
	}
	lumLookUpTableFile.close();
	return (validationKey == inputParametersInstance.parameterValueKey() && validationKeyEnd == validationKey);
	return true;
}


bool
starlight::createEventChannel()
{
	switch (inputParametersInstance.prodParticleType()) {
	case ELECTRON:
	case MUON:
	case TAUON:
		{
			_eventChannel = new Gammagammaleptonpair(*_beamSystem);
			if (_eventChannel)
				return true;
			else {
				printWarn << "cannot construct Gammagammaleptonpair event channel." << endl;
				return false;
			}
		}
	case A2:        // jetset
	case ETA:       // jetset
	case ETAPRIME:  // jetset
	case ETAC:      // jetset
	case F0:        // jetset
		{
#ifdef ENABLE_PYTHIA
			// PythiaOutput = true;
 		        cout<<"Pythia is enabled!"<<endl;
// 			return true;
#else
			printWarn << "Starlight is not compiled against Pythia8; "
			          << "jetset event channel cannot be used." << endl;
 			return false;
#endif
		}
	case F2:
	case F2PRIME:
	case ZOVERZ03:
		{
		  //  #ifdef ENABLE_PYTHIA
	 	        cout<<" This is f2, f2prim, zoverz03 "<<endl; 
			_eventChannel= new Gammagammasingle(*_beamSystem);
			if (_eventChannel)
				return true;
			else {
				printWarn << "cannot construct Gammagammasingle event channel." << endl;
				return false;
			}
			// #endif
			//			printWarn << "Starlight is not compiled against Pythia8; "
			//          << "Gammagammasingle event channel cannot be used." << endl;
			// return false;
		}
	case RHO:
	case RHOZEUS:
	case FOURPRONG:
	case OMEGA:  
	case PHI:
	case JPSI:
	case JPSI2S:
	case JPSI2S_ee:
	case JPSI2S_mumu:
	case JPSI_ee:
	case JPSI_mumu:
	case UPSILON:
	case UPSILON_ee:
	case UPSILON_mumu:
	case UPSILON2S:
	case UPSILON2S_ee:
	case UPSILON2S_mumu:
	case UPSILON3S:
	case UPSILON3S_ee:
	case UPSILON3S_mumu:
		{
			if (inputParametersInstance.interactionType() == PHOTONPOMERONNARROW) {
				_eventChannel = new Gammaanarrowvm(*_beamSystem);
				if (_eventChannel)
					return true;
				else {
					printWarn << "cannot construct Gammaanarrowvm event channel." << endl;
					return false;
				}
			}

			if (inputParametersInstance.interactionType() == PHOTONPOMERONWIDE) {
				_eventChannel = new Gammaawidevm(*_beamSystem);
				if (_eventChannel)
					return true;
				else {
					printWarn << "cannot construct Gammaawidevm event channel." << endl;
					return false;
				}
			}

                        if (inputParametersInstance.interactionType() == PHOTONPOMERONINCOHERENT) {
                                _eventChannel = new Gammaaincoherentvm(*_beamSystem);
                                if (_eventChannel)
                                        return true;
                                else {
                                        printWarn << "cannot construct Gammaanarrowvm event channel." << endl;
                                        return false;
                                }
                        }

			printWarn << "interaction type '" << inputParametersInstance.interactionType() << "' "
			          << "cannot be used with particle type '" << inputParametersInstance.prodParticleType() << "'. "
			          << "cannot create event channel." << endl;
			return false;
		}
	default:
		{
			printWarn << "unknown event channel '" << inputParametersInstance.prodParticleType() << "'."
			          << " cannot create event channel." << endl;
			return false;
		}
	}
}
