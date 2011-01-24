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
//    Class needed for root output
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "eventchannel.h"


using namespace std;


//______________________________________________________________________________
eventChannel::eventChannel(inputParameters& input, beamBeamSystem& bbsystem)
  : readLuminosity(input), _bbs(bbsystem), _nTries(0), _nSuccess(0)/*, _accCutPt(false), _accCutEta(false), _ptMin(0.1), _ptMax(2.0), _etaMin(-0.9), _etaMax(0.9)*/
{
  _randy.SetSeed(input.randomSeed());
  _accCutPt    = input.getCutPt();
  _accCutEta   = input.getCutEta();
  _ptMin       = input.getMinPt();
  _ptMax       = input.getMaxPt();
  _etaMin      = input.getMinEta();
  _etaMax      = input.getMaxEta();
}


//______________________________________________________________________________
eventChannel::~eventChannel()
{ }


//______________________________________________________________________________
void eventChannel::transform(double betax,double betay,double betaz,double &E,
                             double &px,double &py,double &pz,int &iFbadevent)
{
  //     carries out a lorentz transform of the frame.  (Not a
  //     boost!)
  double beta,gamma,gob;
  double E0,px0,py0,pz0;
                                                                                                                                                          
  E0 = E;
  px0 = px;
  py0 = py;
  pz0 = pz;

  beta = sqrt(betax*betax + betay*betay + betaz*betaz);

  if(beta >= 1.0)  iFbadevent = 1;
  gamma = 1./sqrt(1. - beta*beta);

  gob = (gamma-1)/(beta*beta);
                                                                                                                                                          
  E = gamma*(E0 - betax*px0 - betay*py0 - betaz*pz0);
                                                                                                                                                          
  px = -gamma*betax*E0 + (1. + gob*betax*betax)*px0  +  gob*betax*betay*py0 + gob*betax*betaz*pz0;
                                                                                                                                                          
  py = -gamma*betay*E0 + gob*betay*betax*px0+  (1. + gob*betay*betay)*py0 + gob*betay*betaz*pz0;
                                                                                                                                                          
  pz = -gamma*betaz*E0 +  gob*betaz*betax*px0 +  gob*betaz*betay*py0 + (1. + gob*betaz*betaz)*pz0;
}

double eventChannel::pseudoRapidity(double px, double py, double pz)
{
  double pT = sqrt(px*px + py*py);
  double p = sqrt(pz*pz + pT*pT);
  double eta = -99.9; if((p-pz) != 0){eta = 0.5*log((p+pz)/(p-pz));}
  return eta;
}

