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

#include "inputparameters.h"
#include "starlightconstants.h"


using namespace std;


//______________________________________________________________________________
inputParameters::inputParameters() :
   _configFileName("slight.in")
   ,_Z1(0)
   ,_A1(0)
   ,_Z2(0)
   ,_A2(0)
   ,_gamma_em(0.)
   ,_Wmax(0.)
   ,_Wmin(0.)
   ,_numw(0)
   ,_Ymax(0.)
   ,_numy(0)
   ,_gg_or_gP(0)
   ,_ievents(0)
   ,_pid(0)
   ,_iseed(0)
   ,_iout(0)
   ,_ibreakup(0)
   ,_iinterfere(0)
   ,_xinterfere(0.)
   ,_in_or_co(0)
   ,_incoherentFactor(0.)
   ,_bford(0.)
   ,_ptMax(0.)
   ,_NPT(0)
   ,_dpt(0.)
   ,_protonEnergy(0.)
{ }


//______________________________________________________________________________
inputParameters::~inputParameters()
{ }


//______________________________________________________________________________
int inputParameters::init(const std::string& filename)
{
  _configFileName = filename;
  ifstream inputfile;
  inputfile.open(_configFileName.c_str());

  double Wmininput=0.;
  double Wmaxinput=0.;
  
  cout<<"CONSTRUCTING THE INPUT PARAMETERS"<<endl;
  
  inputfile >> _Z1;	//Number of protons in #1 nucleus
  inputfile >> _A1;	//Number of nucleons in #1 nucleus
  inputfile >> _Z2;	//Number of protons in #2 nucleus
  inputfile >> _A2;	//Number of nucleons in #2 nucleus
  inputfile >> _gamma_em;  //Lorentz gamma of the collider for one beam...
  inputfile >> Wmaxinput;  //Maximum value for center of mass energy--might change depending upon which channel
  inputfile >> Wmininput;	//Minimum center of mass energy, for default values use -1
  inputfile >> _numw;	//Number of center of mass bins to use
  inputfile >> _Ymax;	//Sets the maximum rapidity value...more like abs. value
  inputfile >> _numy;	//Number of rapidity bins
  inputfile >> _gg_or_gP;  //Is it a photonphoton, photonA-narrow, or photonA-wide:1,2,3
  inputfile >> _ievents;  //Number of events to generate
  inputfile >> _pid;	//Channel one wants to produce.
  inputfile >> _iseed;	//Number seed for the random number gen.
  inputfile >> _iout;	//Output format,1:txt,2:formatted text for gstar,3:ROOT?
  inputfile >> _ibreakup;	//Breakup modes, 1-...
  inputfile >> _iinterfere; //Interference 0-no,1-yes
  inputfile >> _xinterfere; //Percentage of interference
  inputfile >> _in_or_co;	//dAu,A-A, coherent or incoherent
  inputfile >> _incoherentFactor;//factor times the A*sigma(gp->vmp)
  inputfile >> _bford;	//
  inputfile >> _ptMax;	//
  inputfile >> _NPT;	//
  inputfile.close();
  _dpt = _ptMax/_NPT;

  _protonEnergy=_gamma_em*starlightConstants::mp;

  //Defining Interaction Type
  if(_gg_or_gP==1){
    _interactionTest=starlightConstants::PHOTONPHOTON;
  }
  if(_gg_or_gP==2){
    _interactionTest=starlightConstants::PHOTONPOMERONNARROW;
  }
  if(_gg_or_gP==3){
    _interactionTest=starlightConstants::PHOTONPOMERONWIDE;
  }
  
  //Trying to define the proper Wmins and Wmaxs. a TEMPORARY fix....Better solution=??
  
  double mass         = 0.0;
  double Wmin_default = 0.0; //This is the dault for _Wmin, unless defined later
  double width        = 0.0;
  if(_pid==11){
    //electron
    Wmin_default = 0.01; //_Wmin is settable in constants--default is 0.01; _Wmin up to 0.15 is safe for Summer 2000 triggering for e+e- pairs
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::ELECTRON;
    _decayTest=starlightConstants::LEPTONPAIR;
  }
  else if(_pid==13){
    //muon
    Wmin_default   =  2.*starlightConstants::mmu;
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::MUON;
    _decayTest=starlightConstants::LEPTONPAIR;
  }
  else if(_pid==15){
    //tauon
    Wmin_default   =  2*starlightConstants::mtau;
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::TAUON;
    _decayTest=starlightConstants::LEPTONPAIR;
  }
  else if(_pid==115){
    //a2(1320)
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::A2;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==221){
    //eta
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::ETA;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==225){
    //f2(1270)
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::F2;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==331){
    //eta'(958)
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::ETAPRIME;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==335){
    //f2'(1525)
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::F2PRIME;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==441){
    //eta-c
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::ETAC;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==9010221){
    //f0(980)--was orginally called 10221? updated to standard number
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::F0;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==33){
    //Z"/Z03
    _Wmax=Wmaxinput;
    _pidTest=starlightConstants::ZOVERZ03;
    _decayTest=starlightConstants::SINGLEMESON;
  }
  else if(_pid==113)   {
    mass        =  0.7685;                 //rho(770)0
    width       =  0.1507;
    Wmin_default=  2.*starlightConstants::mpi;
    _Wmax        =  mass+5.*width;
    _pidTest=starlightConstants::RHO;
    _decayTest=starlightConstants::WIDEVMDEFAULT;
  }
  else if(_pid==913)   {
    mass        =  0.7685;                      //rho(770)0 direct pi+pi- decay, interference given by ZEUS
    width       =  0.1507;
    Wmin_default=  2.*starlightConstants::mpi;
    _Wmax        =  mass+5.*width;    //Use the same (ZEUS) 1.5GeV max mass
    _pidTest=starlightConstants::RHOZEUS;
    _decayTest=starlightConstants::WIDEVMDEFAULT;
  }
  else if (_pid == 999) {  // pi+pi-pi+pi- decay
    mass         = 1.350;
    width        = 0.360;
    Wmin_default = 4 * starlightConstants::mpi;
    _Wmax         = 3;
    _pidTest      = starlightConstants::FOURPRONG;
    _decayTest    = starlightConstants::WIDEVMDEFAULT;
  }
  else if(_pid==223){
    mass          =  0.78194;          //omega(782)
    width         =  0.00843;
    _Wmax          =  mass+5.*width;
    Wmin_default  =  mass - 5.*width;
    _pidTest=starlightConstants::OMEGA;
    _decayTest=starlightConstants::NARROWVMDEFAULT;//will probably be shifted to 3body decay
  }
  else if(_pid==333) {
    mass          =  1.019413;     //phi(1020)
    width         =  0.00443;
    Wmin_default  =  2*starlightConstants::mK;
    _Wmax          =  mass+5.*width;
    _pidTest=starlightConstants::PHI;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==443){
    mass          =  3.09692;  // JN  3.09688;      //J/psi
    width         =  0.000091; //JN  0.000087;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::JPSI;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==444){
    mass          =  3.686093;      //J/psi
    width         =  0.000337;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::JPSI2S;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==444011){
    mass          =  3.686093;      //J/psi
    width         =  0.000337;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::JPSI2S_ee;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==444013){
    mass          =  3.686093;      //J/psi
    width         =  0.000337;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::JPSI2S_mumu;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==443011){
    mass          =  3.09692;  // JN  3.09688;      //J/psi
    width         =  0.000091; //JN  0.000087;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::JPSI_ee;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==443013){
    mass          =  3.09692;  // JN  3.09688;      //J/psi
    width         =  0.000091; //JN  0.000087;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::JPSI_mumu;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==553){
    mass          =  9.46030;      //Upsilon
    width         =  0.00005402;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==553011){
    mass          =  9.46030;      //Upsilon
    width         =  0.00005402;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON_ee;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==553013){
    mass          =  9.46030;      //Upsilon
    width         =  0.00005402;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON_mumu;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==554){
    mass          =  10.02326;      //Upsilon(2S)
    width         =  0.00003198;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON2S;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==554011){
    mass          =  10.02326;      //Upsilon(2S)
    width         =  0.00003198;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON2S_ee;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==554013){
    mass          =  10.02326;      //Upsilon(2S)
    width         =  0.00003198;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON2S_mumu;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==555){
    mass          =  10.3552;      //Upsilon(3S)
    width         =  0.00002032;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON3S;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==555011){
    mass          =  10.3552;      //Upsilon(3S)
    width         =  0.00002032;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON3S_ee;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else if(_pid==555013){
    mass          =  10.3552;      //Upsilon(3S)
    width         =  0.00002032;
    Wmin_default  =  mass - 5.*width;
    _Wmax          =  mass + 5.*width;
    _pidTest=starlightConstants::UPSILON3S_mumu;
    _decayTest=starlightConstants::NARROWVMDEFAULT;
  }
  else{
    _Wmax=Wmaxinput;
  }
  if (Wmininput==-1){ 
    _Wmin = Wmin_default;
  }
  else{  
    _Wmin = Wmininput;
  }

  return 0;
}


//______________________________________________________________________________
int inputParameters::getZ1()
{
  return _Z1;
}


//______________________________________________________________________________
int inputParameters::getA1()
{
  return _A1;
}


//______________________________________________________________________________
int inputParameters::getZ2()
{
  return _Z2;
}


//______________________________________________________________________________
int inputParameters::getA2()
{
  return _A2;
}


//______________________________________________________________________________
double inputParameters::getgamma_em()
{
  return _gamma_em;
}


//______________________________________________________________________________
double inputParameters::getWmax()
{
  return _Wmax;
}


//______________________________________________________________________________
double inputParameters::getWmin()
{
  return _Wmin;
}


//______________________________________________________________________________
int inputParameters::getnumw()
{
  return _numw;
}


//______________________________________________________________________________
double inputParameters::getYmax()
{
  return _Ymax;
}


//______________________________________________________________________________
int inputParameters::getnumy()
{
  return _numy;
}


//______________________________________________________________________________
int inputParameters::getgg_or_gP()
{
  return _gg_or_gP;
}


//______________________________________________________________________________
int inputParameters::getNumberOfEvents()
{
  return _ievents;
}


//______________________________________________________________________________
int inputParameters::getSeed()
{
  return _iseed;
}


//______________________________________________________________________________
int inputParameters::getParticleId()
{
  return _pid;
}


//______________________________________________________________________________
int inputParameters::getOutputMode()
{
  return _iout;
}


//______________________________________________________________________________
int inputParameters::getBreakupMode()
{
  return _ibreakup;
}


//______________________________________________________________________________
int inputParameters::getInterferenceMode()
{
  return _iinterfere;
}


//______________________________________________________________________________
double inputParameters::getInterferencePercent()
{
  return _xinterfere;
}


//______________________________________________________________________________
int inputParameters::getIncoherentOrCoherent()
{
  return _in_or_co;
}


//______________________________________________________________________________
double inputParameters::getIncoherentFactor()
{
  return _incoherentFactor;
}


//______________________________________________________________________________
double inputParameters::getbford()
{
  return _bford;
}


//______________________________________________________________________________
double inputParameters::getMaximumInterPt()
{
  return _ptMax;
}


//______________________________________________________________________________
int inputParameters::getNPT()
{
  return _NPT;
}


//______________________________________________________________________________
double inputParameters::getdpt()
{
  return _dpt;
}


//______________________________________________________________________________
starlightConstants::particle inputParameters::getPidTest()
{
  return _pidTest;
}


//______________________________________________________________________________
starlightConstants::decaytype inputParameters::getDecayTest()
{
  return _decayTest;
}


//______________________________________________________________________________
starlightConstants::interactiontype inputParameters::getInteractionTest()
{
  return _interactionTest;
}


//______________________________________________________________________________
double inputParameters::getProtonEnergy()
{
  return _protonEnergy;
}
