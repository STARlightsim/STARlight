// beam.cpp
/*
 * $Id: beam.cpp,v 1.0 2010/07/04   $
 *
 * /author Joseph Butterwoth
 *
 * $Log: $
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
 */
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#include "beam.h"
#include "starlightconstants.h"
#include "inputparameters.h"
#include "bessel.h"


//______________________________________________________________________________
Beam::Beam(int Zin, int Ain, double bdeuteron, int in_or_co, Inputparameters& input): Nucleus(Zin,Ain,bdeuteron,in_or_co)
{
  
  //Setting needed inputparameters to protected variables
  BeamInputGamma_em=input.getgamma_em();
  
}
//______________________________________________________________________________
double Beam::nofe(double impactparameter)
{
  //Function for the calculation of the "photon density".
  //nofe=numberofgammas/(energy*area)
  //Assume beta=1.0 and gamma>>1, i.e. neglect the (1/gamma^2)*K0(x) term


  double X=0.,nofex=0.,factor1=0.,factor2=0.,factor3=0.;

  
  X = (impactparameter*Egamma)/(BeamInputGamma_em*StarlightConstants::hbarc);
  
  if(X <= 0.0) 
    cout<<"In nofe, X= "<<X<<endl;
  
  factor1 = (double(getZin()*getZin())*StarlightConstants::alpha)/
    (StarlightConstants::pi*StarlightConstants::pi);
  
  factor2 = 1./(Egamma*impactparameter*impactparameter);
  factor3 = X*X*(bessel::dbesk1(X))*(bessel::dbesk1(X));
  nofex   = factor1*factor2*factor3;

  return nofex;
   
}
//______________________________________________________________________________
Beam::~Beam()
{

}
