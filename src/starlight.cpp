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

#ifdef ENABLE_PYTHIA
#include "PythiaStarlight.h"
#endif

#include "starlight.h"
#include "inputparameters.h"
#include "eventchannel.h"
#include "gammagammaleptonpair.h"
#include "gammagammasingle.h"
#include "gammaavm.h"
#include "psifamily.h"
#include "twophotonluminosity.h"
#include "gammaaluminosity.h"
#include "upcevent.h"
#include "eventfilewriter.h"


starlight::starlight() :
	_inputParameters(0)
	,_beam0(0)
	,_beam1(0)
	,_beamSystem(0)
	,_eventChannel(0)
	,_configFileName("slight.in")
	,_numberOfEventsPerFile(100)
	,_numberOfEventsToGenerate(10)
	,_standardFilename("slight.out")
	,_isInitialised(false)
{ }


starlight::~starlight()
{ }


int starlight::init()
{
	std::cout << "##################################" << std::endl;
	std::cout << " Initialising Starlight v" << Starlight_VERSION_MAJOR << "." << Starlight_VERSION_MINOR << "..." << std::endl;
	std::cout << "##################################" << std::endl;

	_numberOfEventsToGenerate = _inputParameters->nmbEvents();
	_numberOfEventsPerFile = _numberOfEventsToGenerate; // For now we write only one file...

	_beamSystem = new beamBeamSystem(*_inputParameters);

	// std::streamsize precision(15);
	std::cout.setf(std::ios_base::fixed,std::ios_base::floatfield);
	std::cout.precision(15);

	bool flag = checkForLuminosityTable();

	switch (_inputParameters->getInteractionTest())
		{
		case starlightConstants::PHOTONPHOTON:
			if (flag==true) {
				std::cout << "CREATING LUMINOSITY TABLE FOR PHOTONPHOTON" << std::endl;
				twoPhotonLuminosity(_beamSystem->getBeam1(), _beamSystem->getBeam2(), _inputParameters->beamBreakupMode(), *_inputParameters);
			}
			break;
		case starlightConstants::PHOTONPOMERONNARROW://Both narrow and wide use the same luminosity function.
		case starlightConstants::PHOTONPOMERONWIDE:
			if (flag==true) {
				std::cout << "CREATING LUMINOSITY TABLE FOR GAMMA" << std::endl;
				//Luminosity function
				photonNucleusLuminosity(*_inputParameters, *_beamSystem);
			}
			break;
		default :
			std::cout << "Please go back and define an appropriate interaction type. Thank you."<<std::endl;
		}

	int res = createEventChannel();

	if (res)
		{
			return -1;
		}
    
	_isInitialised = true;
    
	return 0;
}


upcEvent starlight::produceEvent()
{
	if(!_isInitialised)
		{
			std::cerr << "Trying to produce event but Starlight is not initialised, exiting..." << std::endl;
			exit(-1);
		}
	return _eventChannel->produceEvent();
}

bool starlight::checkForLuminosityTable()
{
	std::cout<<"ISEED: "<<_inputParameters->randomSeed()<<std::endl;
	std::ifstream wylumfile;
	wylumfile.precision(15);
	wylumfile.open("slight.txt");
	unsigned int Z1test = 0, A1test = 0, Z2test = 0, A2test = 0, numwtest = 0, numytest = 0;
	int gg_or_gPtest=0,ibreakuptest=0,iinterferetest=0,NPTtest=0,in_or_cotest=0;
	double Gammatest=0.,Wmaxtest=0.,Wmintest=0.,Ymaxtest=0.,xinterferetest=0.,ptmaxtest=0.,bfordtest=0.,incoherentfactortest=0.;
	bool flag = false;
	// bool b;
	wylumfile >> Z1test;
	wylumfile >> A1test;
	wylumfile >> Z2test;
	wylumfile >> A2test;
	wylumfile >> Gammatest;
	wylumfile >> Wmaxtest;
	wylumfile >> Wmintest;
	wylumfile >> numwtest;
	wylumfile >> Ymaxtest;
	wylumfile >> numytest;
	wylumfile >> gg_or_gPtest;
	wylumfile >> ibreakuptest;
	wylumfile >> iinterferetest;
	wylumfile >> xinterferetest;
	wylumfile >> in_or_cotest;
	wylumfile >> incoherentfactortest;
	wylumfile >> bfordtest;
	wylumfile >> ptmaxtest;
	wylumfile >> NPTtest;
	wylumfile.close();

	if ( !(
	       _inputParameters->beam1Z() == Z1test
	       && _inputParameters->beam1A() == A1test
	       && _inputParameters->beam2Z() == Z2test
	       && _inputParameters->beam2A() == A2test
	       && _inputParameters->beamLorentzGamma() == Gammatest
	       && _inputParameters->numWBins() == numwtest
	       && _inputParameters->minW() == Wmintest
	       && _inputParameters->maxRapidity() == Ymaxtest
	       && _inputParameters->nmbRapidityBins() == numytest
	       && _inputParameters->productionMode() == gg_or_gPtest
	       && _inputParameters->beamBreakupMode() == ibreakuptest
	       && _inputParameters->interferenceEnabled() == iinterferetest
	       && _inputParameters->interferenceStrength() == xinterferetest
	       && _inputParameters->getbford() == bfordtest
	       && _inputParameters->coherentProduction() == in_or_cotest
	       && _inputParameters->incoherentFactor() == incoherentfactortest
	       && _inputParameters->maxPtInterference() == ptmaxtest
	       && _inputParameters->nmbPtBinsInterference() == NPTtest )
	     )
		{
			//okay, if we are in this loop, it means the input parameters are different than the one's used to create the last set of luminosity tables
			//Now lets create a new set
			flag=true;

		}

	return flag;

}


int starlight::createEventChannel()
{
	switch (_inputParameters->getPidTest()) {
	case starlightConstants::ELECTRON:
	case starlightConstants::MUON:
	case starlightConstants::TAUON:
		{
			_eventChannel = new Gammagammaleptonpair(*_inputParameters, *_beamSystem);
			if (_eventChannel) return 0;
			else return -1;
		}
	case starlightConstants::A2://jetset
	case starlightConstants::ETA://jetset
	case starlightConstants::ETAPRIME://jetset
	case starlightConstants::ETAC://jetset
	case starlightConstants::F0://jetset
		{
#ifdef ENABLE_PYTHIA
			//	    PythiaOutput=true;
			return 0;
#endif
			std::cout << "Starlight is not compiled against Pythia8, jetset cannot be used" << std::endl;
			return -1;
			//This way we can output mother and daughter listings.
		}
	case starlightConstants::F2:
	case starlightConstants::F2PRIME:
	case starlightConstants::ZOVERZ03:
		{
#ifdef ENABLE_PYTHIA
			_eventChannel= new Gammagammasingle(*_inputParameters, *_beamSystem);
			if (_eventChannel) return 0;
			else return -1;
#endif
			std::cout << "Starlight is not compiled against Pythia8, gamma-gamma single cannot be used" << std::endl;
			return -1;

		}
	case starlightConstants::RHO:
	case starlightConstants::RHOZEUS:
	case starlightConstants::FOURPRONG:
	case starlightConstants::OMEGA://Will probably be three body
	case starlightConstants::PHI:
	case starlightConstants::JPSI:
	case starlightConstants::JPSI2S:
	case starlightConstants::JPSI2S_ee:
	case starlightConstants::JPSI2S_mumu:
	case starlightConstants::JPSI_ee:
	case starlightConstants::JPSI_mumu:
	case starlightConstants::UPSILON:
	case starlightConstants::UPSILON_ee:
	case starlightConstants::UPSILON_mumu:
	case starlightConstants::UPSILON2S:
	case starlightConstants::UPSILON2S_ee:
	case starlightConstants::UPSILON2S_mumu:
	case starlightConstants::UPSILON3S:
	case starlightConstants::UPSILON3S_ee:
	case starlightConstants::UPSILON3S_mumu:
		{
			if (_inputParameters->getInteractionTest()==2) {
				_eventChannel = new Gammaanarrowvm(*_inputParameters, *_beamSystem);
				if (_eventChannel) return 0;
				else return -1;
			}

			if (_inputParameters->getInteractionTest()==3) {
				_eventChannel = new Gammaawidevm(*_inputParameters, *_beamSystem);
				if (_eventChannel) return 0;
				else return -1;
			}
			std::cout<<"Please go back and adjust gg_or_gp to 2 or 3 for a VM, main.cpp"<<std::endl;
			return -1;
		}
		//    case starlightConstants::JPSI:
		//    case starlightConstants::JPSI2S:
		//    {
		//        _eventChannel = new psiFamily(*_inputParameters, *_beamSystem);
		//        if (_eventChannel) return 0;
		//        else return -1;
		//    }
		//rhoprime
	default:
		std::cout<<"Hi and welcome to default event channel(null), main::CreateEventChannel"<<std::endl;
		return -1;//Maybe return empty eventChannel object?
	}
}
