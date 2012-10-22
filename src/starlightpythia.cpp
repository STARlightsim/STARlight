/*
    <one line to give the library's name and an idea of what it does.>
    Copyright (C) 2011  Oystein Djuvsland <oystein.djuvsland@gmail.com>

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <cstdio>
#include <iostream>
#include "starlightpythia.h"
#include "pythiaInterface.h"
#include "spectrumprotonnucleus.h"
#include <cmath>
#include <sstream>

starlightPythia::starlightPythia(inputParameters &input, beamBeamSystem &bbs) : eventChannel(input, bbs)
        ,_spectrum(0)
        ,_doDoubleEvent(false)
        ,_minGammaEnergy(input.minGammaEnergy())
        ,_maxGammaEnergy(input.maxGammaEnergy())
	,_fullEventRecord(false)
{
}

starlightPythia::~starlightPythia()
{

}

int starlightPythia::init(std::string pythiaParams, bool fullEventRecord)
{
   _fullEventRecord = fullEventRecord;
   _spectrum = new spectrumProtonNucleus(&_bbs);
   //_spectrum = new Spectrum(&bbs);
   _spectrum->setRandomGenerator(&_randy);

   _spectrum->setMinGammaEnergy(_minGammaEnergy);
   _spectrum->setMaxGammaEnergy(_maxGammaEnergy);
   
   if(!_doDoubleEvent)
   {
    _spectrum->generateKsingle();
   }
   else 
   {
    _spectrum->generateKdouble();
   }

    // Set to run with varying energies
    pythiaInterface::pygive("mstp(171)=1"); // Varying energies
    pythiaInterface::pygive("mstp(172)=1"); // Set the energy before generation

    std::stringstream ss(pythiaParams);
    std::string p;
    while(std::getline(ss, p, ';')) 
    {
      if(p.size()>1)
      {
        pythiaInterface::pygive(p.c_str());
      }
    }
    //pythiaInterface::pygive("mstp(12)=0");
    //pythiaInterface::pygive("mstj(21)=0"); // Disable decays of particles 
    
    //pythiaInterface::pygive("parp(2)=1.0"); // Cut off c.m. energy (GeV)
  
    pythiaInterface::pyinit("FIXT", "gamma", "p", _maxGammaEnergy); // Fixed target, beam, target, beam momentum (GeV/c)

    return 0;
}

upcEvent starlightPythia::produceEvent()
{
  upcEvent event;
  
  
  
    if (!_doDoubleEvent)
    {
      //int zdirection = (Randy.Rndom()) < 0.5 ? -1 : 1;
      double gammaE = 0;
      do
      {
	gammaE = _spectrum->drawKsingle();
      } while(isnan(gammaE));
      event.addGamma(gammaE);

      char opt[32];
      std::sprintf(opt, "parp(171)=%f", gammaE/_maxGammaEnergy);
      pythiaInterface::pygive(opt); // Set the energy of the photon beam (gammaE/1000 * 1000.0);
      pythiaInterface::pyevnt(); // Generate event
//      pythiaInterface::pyfram(2); // go to CMS 
      double boost = _bbs.cmsBoost();
      vector3 boostVector(0, 0, tanh(boost));
      for(int idx = 0; idx < pyjets_.n; idx++)
      {
//	if(std::abs(pyjets_.k[1][idx]) <= 6) std::cout << "Quark: " << pyjets_.k[1][idx] << ", status: " <<  pyjets_.k[0][idx] << std::endl;
	if(pyjets_.k[0][idx] > 10 && _fullEventRecord==false) continue;
	int pdgCode = pyjets_.k[1][idx];
	int charge = 0;
	if( pdgCode == 12 || pdgCode == 14 || pdgCode == 16 ||  pdgCode == 22 || pdgCode == 111 || pdgCode == 130 || pdgCode == 321 || pdgCode == 2112)
	{
	  charge = 0;
	}
	else
	{
	  charge = (pdgCode > 0) - (pdgCode < 0);
	}
	
	starlightParticle particle(pyjets_.p[0][idx], pyjets_.p[1][idx], pyjets_.p[2][idx], pyjets_.p[3][idx], pyjets_.p[4][idx], pyjets_.k[1][idx], charge);
	if(_fullEventRecord)
	{
	  particle.setParent(pyjets_.k[2][idx]);
	  particle.setFirstDaughter(pyjets_.k[3][idx]);
	  particle.setLastDaughter(pyjets_.k[4][idx]);
	  particle.setStatus(pyjets_.k[0][idx]);
	}
	particle.Boost(boostVector);
        event.addParticle(particle);
      }
      
    }
    return event;
}