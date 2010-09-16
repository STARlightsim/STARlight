// starlight.cpp
/*
 * $Id: starlight.cpp,v 1.0 2010/07/04  $
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
 

#include "starlight.h"
#include <iostream>
#include <fstream>
#include "inputparameters.h"
#include "eventchannel.h"
#include "gammagammaleptonpair.h"
#include "gammagammasingle.h"
#include "gammaavm.h"
#include "psifamily.h"
#include "twophotonluminosity.h"
#include "gammaaluminosity.h"
#include <cstdlib>

#ifdef ENABLE_PYTHIA
#include "PythiaStarlight.h"
#endif
#include <upcevent.h>
#include <eventfilewriter.h>



Starlight::Starlight() :
        fInputParameters(0)
        ,fBeam0(0)
        ,fBeam1(0)
        ,fBeamSystem(0)
        ,fEventChannel(0)
        ,fConfigFileName("slight.in")
        ,fNumberOfEventsPerFile(100)
        ,fNumberOfEventsToGenerate(10)
        ,fStandardFilename("slight.out")
        ,fIsInitialised(false)
{
}

Starlight::~Starlight()
{
}

int Starlight::Init()
{

    std::cout << "##################################" << std::endl;
    std::cout << " Initialising Starlight v" << Starlight_VERSION_MAJOR << "." << Starlight_VERSION_MINOR << "..." << std::endl;
    std::cout << "##################################" << std::endl;

    fNumberOfEventsToGenerate = fInputParameters->getnumberofevents();
    fNumberOfEventsPerFile = fNumberOfEventsToGenerate; // For now we write only one file...

    fBeamSystem = new Beambeamsystem(*fInputParameters);

    std::streamsize precision(15);
    std::cout.setf(std::ios_base::fixed,std::ios_base::floatfield);
    std::cout.precision(15);

    bool flag = CheckForLuminosityTable();

    switch (fInputParameters->getinteractiontest())
    {
    case StarlightConstants::PHOTONPHOTON:
        if (flag==true) {
            std::cout << "CREATING LUMINOSITY TABLE FOR PHOTONPHOTON" << std::endl;
            Twophotonluminosity(fBeamSystem->getBeam1(), fBeamSystem->getBeam2(), fInputParameters->getbreakupmode(), *fInputParameters);
        }
        break;
    case StarlightConstants::PHOTONPOMERONNARROW://Both narrow and wide use the same luminosity function.
    case StarlightConstants::PHOTONPOMERONWIDE:
        if (flag==true) {
            std::cout << "CREATING LUMINOSITY TABLE FOR GAMMA" << std::endl;
            //Luminosity function
            Gammaaluminosity(*fInputParameters, *fBeamSystem);
        }
        break;
    default :
        std::cout << "Please go back and define an appropriate interaction type. Thank you."<<std::endl;
    }

    int res = CreateEventChannel();

    if (res)
    {
        return -1;
    }
    
    fIsInitialised = true;
    
    return 0;
}

UPCEvent Starlight::ProduceEvent()
{
   if(!fIsInitialised)
   {
      std::cerr << "Trying to produce event but Starlight is not initialised, exiting..." << std::endl;
      exit(-1);
   }
   return fEventChannel->ProduceEvent();
}

bool Starlight::CheckForLuminosityTable()
{
    std::cout<<"ISEED: "<<fInputParameters->getseed()<<std::endl;
    std::ifstream wylumfile;
    wylumfile.precision(15);
    wylumfile.open("slight.txt");
    int Z1test=0,A1test=0,Z2test=0,A2test=0,numwtest=0,numytest=0,gg_or_gPtest=0,ibreakuptest=0,iinterferetest=0,NPTtest=0,in_or_cotest=0;
    double Gammatest=0.,Wmaxtest=0.,Wmintest=0.,Ymaxtest=0.,xinterferetest=0.,ptmaxtest=0.,bfordtest=0.,incoherentfactortest=0.;
    bool flag = false;
    bool b;
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
                fInputParameters->getZ1() == Z1test
                && fInputParameters->getA1() == A1test
                && fInputParameters->getZ2() == Z2test
                && fInputParameters->getA2() == A2test
                && fInputParameters->getgamma_em() == Gammatest
                && fInputParameters->getnumw() == numwtest
                && fInputParameters->getWmin() == Wmintest
                && fInputParameters->getYmax() == Ymaxtest
                && fInputParameters->getnumy() == numytest
                && fInputParameters->getgg_or_gP() == gg_or_gPtest
                && fInputParameters->getbreakupmode() == ibreakuptest
                && fInputParameters->getinterferencemode() == iinterferetest
                && fInputParameters->getinterferencepercent() == xinterferetest
                && fInputParameters->getbford() == bfordtest
                && fInputParameters->getincoherentorcoherent() == in_or_cotest
		&& fInputParameters->getincoherentfactor() == incoherentfactortest
                && fInputParameters->getmaximuminterpt() == ptmaxtest
                && fInputParameters->getNPT() == NPTtest )
       )
    {
        //okay, if we are in this loop, it means the input parameters are different than the one's used to create the last set of luminosity tables
        //Now lets create a new set
        flag=true;

    }

    return flag;

}

int Starlight::CreateEventChannel()
{
    switch (fInputParameters->getpidtest()) {
    case StarlightConstants::ELECTRON:
    case StarlightConstants::MUON:
    case StarlightConstants::TAUON:
    {
        fEventChannel = new Gammagammaleptonpair(*fInputParameters, *fBeamSystem);
        if (fEventChannel) return 0;
        else return -1;
    }
    case StarlightConstants::A2://jetset
    case StarlightConstants::ETA://jetset
    case StarlightConstants::ETAPRIME://jetset
    case StarlightConstants::ETAC://jetset
    case StarlightConstants::F0://jetset
    {
#ifdef ENABLE_PYTHIA
//	    PythiaOutput=true;
        return 0;
#endif
        std::cout << "Starlight is not compiled against Pythia8, jetset cannot be used" << std::endl;
        return -1;
        //This way we can output mother and daughter listings.
    }
    case StarlightConstants::F2:
    case StarlightConstants::F2PRIME:
    case StarlightConstants::ZOVERZ03:
    {
#ifdef ENABLE_PYTHIA
        fEventChannel= new Gammagammasingle(*fInputParameters, *fBeamSystem);
        if (fEventChannel) return 0;
        else return -1;
#endif
        std::cout << "Starlight is not compiled against Pythia8, gamma-gamma single cannot be used" << std::endl;
        return -1;

    }
    case StarlightConstants::RHO:
    case StarlightConstants::RHOZEUS:
    case StarlightConstants::OMEGA://Will probably be three body
    case StarlightConstants::PHI:
    case StarlightConstants::JPSI:
    case StarlightConstants::JPSI2S:
    case StarlightConstants::JPSI_ee:
    case StarlightConstants::JPSI_mumu:
    case StarlightConstants::UPSILON:
    case StarlightConstants::UPSILON2S:
    case StarlightConstants::UPSILON3S:
    {
        if (fInputParameters->getinteractiontest()==2) {
            fEventChannel = new Gammaanarrowvm(*fInputParameters, *fBeamSystem);
            if (fEventChannel) return 0;
            else return -1;
        }

        if (fInputParameters->getinteractiontest()==3) {
            fEventChannel = new Gammaawidevm(*fInputParameters, *fBeamSystem);
            if (fEventChannel) return 0;
            else return -1;
        }
        std::cout<<"Please go back and adjust gg_or_gp to 2 or 3 for a VM, main.cpp"<<std::endl;
        return -1;
    }
    //    case StarlightConstants::JPSI:
    //    case StarlightConstants::JPSI2S:
    //    {
    //        fEventChannel = new Psifamily(*fInputParameters, *fBeamSystem);
    //        if (fEventChannel) return 0;
    //        else return -1;
    //    }
    //rhoprime
    default:
        std::cout<<"Hi and welcome to default event channel(null), main::CreateEventChannel"<<std::endl;
        return -1;//Maybe return empty Eventchannel object?
    }
}
