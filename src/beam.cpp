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
#include <cmath>

#include "beam.h"
#include "starlightconstants.h"
#include "inputparameters.h"
#include "bessel.h"


using namespace std;


//______________________________________________________________________________
beam::beam(int Zin, int Ain, double bdeuteron, int in_or_co, inputParameters& input)
	: nucleus(Zin, Ain, bdeuteron, in_or_co)
{
  //Setting needed inputparameters to protected variables
  _beamInputGamma_em=input.beamLorentzGamma();
}


//______________________________________________________________________________
beam::~beam()
{ }


//______________________________________________________________________________
double beam::nofe(double impactparameter)
{
  //Function for the calculation of the "photon density".
  //nofe=numberofgammas/(energy*area)
  //Assume beta=1.0 and gamma>>1, i.e. neglect the (1/gamma^2)*K0(x) term

  double X=0.,nofex=0.,factor1=0.,factor2=0.,factor3=0.;
  
  X = (impactparameter*_photonEnergy)/(_beamInputGamma_em*starlightConstants::hbarc);
  
  if(X <= 0.0) 
    cout<<"In nofe, X= "<<X<<endl;
  
  factor1 = (double(getZin()*getZin())*starlightConstants::alpha)/
    (starlightConstants::pi*starlightConstants::pi);
  
  factor2 = 1./(_photonEnergy*impactparameter*impactparameter);
  factor3 = X*X*(bessel::dbesk1(X))*(bessel::dbesk1(X));
  nofex   = factor1*factor2*factor3;

  return nofex;
}
