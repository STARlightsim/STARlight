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
// $Rev:: 45                          $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2011-02-27 20:52:35 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "inputParameters.h"
#include "beambeamsystem.h"
#include "beam.h"
#include "starlightconstants.h"
#include "nucleus.h"
#include "bessel.h"
#include "incoherentPhotonNucleusLuminosity.h"


using namespace std;


//______________________________________________________________________________
incoherentPhotonNucleusLuminosity::incoherentPhotonNucleusLuminosity(inputParameters& input, beamBeamSystem& bbsystem)
  : photonNucleusCrossSection(input, bbsystem), _inputgammaa(input)
{
  cout <<"Creating Luminosity Tables for incoherent vector meson production."<<endl;
  incoherentPhotonNucleusDifferentialLuminosity();
  cout <<"Luminosity Tables created."<<endl;
}


//______________________________________________________________________________
incoherentPhotonNucleusLuminosity::~incoherentPhotonNucleusLuminosity()
{ }


//______________________________________________________________________________
void incoherentPhotonNucleusLuminosity::incoherentPhotonNucleusDifferentialLuminosity()
{
	// double Av,Wgp,cs,cvma;
  double W,dW,dY;
  double Egamma,Y;
  // double t,tmin,tmax;
  double testint,dndWdY;
  // double ax,bx;
  double C;  

  ofstream wylumfile;
  wylumfile.precision(15);
  
  double  bwnorm,Eth;

  dW = (_inputgammaa.maxW()-_inputgammaa.minW())/_inputgammaa.nmbWBins();
  dY  = (_inputgammaa.maxRapidity()-(-1.0)*_inputgammaa.maxRapidity())/_inputgammaa.nmbRapidityBins();
    
  // Write the values of W used in the calculation to slight.txt.  
  wylumfile.open("slight.txt");
  wylumfile << getbbs().beam1().Z() <<endl;
  wylumfile << getbbs().beam1().A() <<endl;
  wylumfile << getbbs().beam2().Z() <<endl;
  wylumfile << getbbs().beam2().A() <<endl;
  wylumfile << _inputgammaa.beamLorentzGamma() <<endl;
  wylumfile << _inputgammaa.maxW() <<endl;
  wylumfile << _inputgammaa.minW() <<endl;
  wylumfile << _inputgammaa.nmbWBins() <<endl;
  wylumfile << _inputgammaa.maxRapidity() <<endl;
  wylumfile << _inputgammaa.nmbRapidityBins() <<endl;
  wylumfile << _inputgammaa.productionMode() <<endl;
  wylumfile << _inputgammaa.beamBreakupMode() <<endl;
  wylumfile << _inputgammaa.interferenceEnabled() <<endl;
  wylumfile << _inputgammaa.interferenceStrength() <<endl;
  wylumfile << _inputgammaa.coherentProduction() <<endl;
  wylumfile << _inputgammaa.incoherentFactor() <<endl;
  wylumfile << _inputgammaa.deuteronSlopePar() <<endl;
  wylumfile << _inputgammaa.maxPtInterference() <<endl;
  wylumfile << _inputgammaa.nmbPtBinsInterference() <<endl;
  
  //     Normalize the Breit-Wigner Distribution and write values of W to slight.txt
  testint=0.0;
  //Grabbing default value for C in the breit-wigner calculation
  C=getDefaultC();
  for(unsigned int i = 0; i <= _inputgammaa.nmbWBins() - 1; ++i) {
    W = _inputgammaa.minW() + double(i)*dW + 0.5*dW;
    testint = testint + breitWigner(W,C)*dW;
    wylumfile << W << endl;
  }
  bwnorm = 1./testint;
  
  //     Write the values of Y used in the calculation to slight.txt.
  for(unsigned int i = 0; i <= _inputgammaa.nmbRapidityBins() - 1; ++i) {
    Y = -1.0*_inputgammaa.maxRapidity() + double(i)*dY + 0.5*dY;
    wylumfile << Y << endl;
  }
 
  //  Eth=0.5*(((_inputgammaa.minW()+starlightConstants::protonMass)*(_inputgammaa.minW()
  //							    +starlightConstants::protonMass)-starlightConstants::protonMass*starlightConstants::protonMass)/
  //	   (Ep + sqrt(Ep*Ep-starlightConstants::protonMass*starlightConstants::protonMass)));
  
  for(unsigned int i = 0; i <= _inputgammaa.nmbWBins() - 1; ++i) {

    W = _inputgammaa.minW() + double(i)*dW + 0.5*dW;

    double Ep = _inputgammaa.getProtonEnergy();

    Eth=0.5*(((W+starlightConstants::protonMass)*(W+starlightConstants::protonMass)-starlightConstants::protonMass*starlightConstants::protonMass)/
	   (Ep + sqrt(Ep*Ep-starlightConstants::protonMass*starlightConstants::protonMass)));
    
    for(unsigned int j = 0; j <= _inputgammaa.nmbRapidityBins() - 1; ++j) {

      Y = -1.0*_inputgammaa.maxRapidity() + double(j)*dY + 0.5*dY;

      int A_1 = getbbs().beam1().A(); 
      int A_2 = getbbs().beam2().A();
      if( A_2 == 1 && A_1 != 1 ){
        // pA, first beam is the nucleus 
        Egamma = 0.5*W*exp(Y);
      } else if( A_1 ==1 && A_2 != 1){
        // pA, second beam is the nucleus 
        Egamma = 0.5*W*exp(-Y); 
      } else {
        Egamma = 0.5*W*exp(Y);        
      }
      
      dndWdY = 0.; 

      if(Egamma > Eth){
	if(Egamma > maxPhotonEnergy())Egamma = maxPhotonEnergy();
        double Wgp = sqrt(2.*Egamma*(Ep+sqrt(Ep*Ep-starlightConstants::protonMass*
                                 starlightConstants::protonMass))+starlightConstants::protonMass*starlightConstants::protonMass);

        double localsig = sigmagp(Wgp); 
        // int localz = 0; 
        // double localbmin = 0; 
        if( A_1 == 1 && A_2 != 1 ){
          // localbmin = getbbs().beam2().nuclearRadius() + 0.7; 
          // localz = getbbs().beam2().Z(); 
	  //   dndWdY = Egamma*localz*localz*nepoint(Egamma,localbmin)*localsig*breitWigner(W,bwnorm); 
          dndWdY = Egamma*photonFlux(Egamma)*localsig*breitWigner(W,bwnorm); 
        }else if (A_2 ==1 && A_1 !=1){
          // localbmin = getbbs().beam1().nuclearRadius() + 0.7; 
          // localz = getbbs().beam1().Z(); 
	  //   dndWdY = Egamma*localz*localz*nepoint(Egamma,localbmin)*localsig*breitWigner(W,bwnorm); 
          dndWdY = Egamma*photonFlux(Egamma)*localsig*breitWigner(W,bwnorm); 
        }else{ 
          double csVN = sigma_N(Wgp);         
          double csVA = sigma_A(csVN); 
          double csgA= (csVA/csVN)*sigmagp(Wgp); 
          dndWdY = Egamma*photonFlux(Egamma)*csgA*breitWigner(W,bwnorm); 
        }
      }

      wylumfile << dndWdY << endl;
    }
  }

  wylumfile.close();
  
  wylumfile.open("slight.txt",ios::app);
  cout << "bwnorm: "<< bwnorm <<endl;
  wylumfile << bwnorm << endl;
  wylumfile.close();
}


//______________________________________________________________________________
double incoherentPhotonNucleusLuminosity::nofe(double Egamma, double bimp)
{
  //Function for the calculation of the "photon density".
  //nofe=numberofgammas/(energy*area)
  //Assume beta=1.0 and gamma>>1, i.e. neglect the (1/gamma^2)*K0(x) term
  
  double X=0.,nofex=0.,factor1=0.,factor2=0.,factor3=0.;
  
  X = (bimp*Egamma)/(_inputgammaa.beamLorentzGamma()*starlightConstants::hbarc);
  
  if( X <= 0.0) 
    cout<<"In nofe, X= "<<X<<endl;
  
  factor1 = (double(getbbs().beam1().Z()*getbbs().beam1().Z())
	     *starlightConstants::alpha)/(starlightConstants::pi*starlightConstants::pi);

  factor2 = 1./(Egamma*bimp*bimp);
  factor3 = X*X*(bessel::dbesk1(X))*(bessel::dbesk1(X));
  nofex    = factor1*factor2*factor3;
  return nofex;
}
