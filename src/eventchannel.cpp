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
: readLuminosity(input), _bbs(bbsystem)
{
  _randy.SetSeed(input.randomSeed());
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
