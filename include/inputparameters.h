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


#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H


#include "starlightconstants.h"
//This is where we read in our input values.


class inputParameters
{
 public:
  inputParameters();
  ~inputParameters();

  int init(const std::string& filename = "config/slight.in");

  int getZ1();
  int getA1();
  int getZ2();
  int getA2();
  double getgamma_em();
  double getWmax();
  double getWmin();
  int getnumw();
  double getYmax();
  int getnumy();
  int getgg_or_gP();
  int getNumberOfEvents();
  int getParticleId();
  int getSeed();
  int getOutputMode();
  int getBreakupMode();
  int getInterferenceMode();
  double getInterferencePercent();
  int getIncoherentOrCoherent();
  double getIncoherentFactor();
  double getbford();
  double getMaximumInterPt();
  int getNPT();
  double getdpt();
  starlightConstants::particle getPidTest();
  starlightConstants::decaytype getDecayTest();
  starlightConstants::interactiontype getInteractionTest();
  double getf2o4pi();
  double getbslope();
  double getProtonEnergy();
  
 private:

  /** Name of the configuration file (default: slight.in) */
  std::string _configFileName; 

  int    _Z1;          // Atomic Number
  int    _A1;          //Atomic Mass
  int    _Z2;
  int    _A2;
  double _gamma_em;   //gamma for the colliding ions
  double _Wmax;       //gamma-gamma center of mass Energy, max
  double _Wmin;       //min for Energy, -1 is default
  int    _numw;          //number of bins for energy in the calculations
  double _Ymax;       //max value for the rapidity
  int    _numy;       //number of bins used
  int    _gg_or_gP;   //1=2 photon channels,2=vector meson with narrow resonance,3=vector meson with wide resonance(B-W)
  int    _ievents;          //number of events to generate
  int    _pid;           //pdg ID --channel choice //ip
  int    _iseed;        //random number seed
  int    _iout;             //output format,1=text,2=GSTARtext,3=ntuple
  int    _ibreakup;             //1=hard sphere nuclei (b>2),2=both nuclei break up (XnXn),3=a single neutron from each nucleus (1n1n),
  //4=require that neither nucleon break up (with b>2R),5=require that there be no hadronic break up (This is similar to option one, but with the actual hadronic interaction)
  int    _iinterfere;             //0= no interference, 1= interference
  double _xinterfere;            //When there is interference, this gives the %age of interference, 0=none to 1=full
  int    _in_or_co;//0 for incoherent 1 for coherent
  double _incoherentFactor; //Variable allows the user to scale the incoherent contribution in VM production
  double _bford;
  double _ptMax;          //When there is interference, this is the max pt considered
  int    _NPT;           //When there is interference, this is the number of pt bins
  double _dpt;
  starlightConstants::particle _pidTest;   //Testing the starlightconstants particle enumeration here...
  starlightConstants::decaytype _decayTest;
  starlightConstants::interactiontype _interactionTest;
  double _protonEnergy;
};


#endif  // INPUTPARAMETERS_H
