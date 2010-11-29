// inputparameters.cpp
/*
 * $Id: inputparameters.cpp,v 1.0 2010/07/04   $
 *
 * /author Joseph Butterwoth
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
 * $Log: $
 *
 */
#include <iostream>
#include <fstream>

using namespace std;

#include "inputparameters.h"
#include "starlightconstants.h"
//______________________________________________________________________________
Inputparameters::Inputparameters() :
   configFileName("slight.in")
   ,Z1(0)
   ,A1(0)
   ,Z2(0)
   ,A2(0)
   ,gamma_em(0.)
   ,Wmax(0.)
   ,Wmin(0.)
   ,numw(0)
   ,Ymax(0.)
   ,numy(0)
   ,gg_or_gP(0)
   ,ievents(0)
   ,pid(0)
   ,iseed(0)
   ,iout(0)
   ,ibreakup(0)
   ,iinterfere(0)
   ,xinterfere(0.)
   ,in_or_co(0)
   ,incoherentfactor(0.)
   ,bford(0.)
   ,ptmax(0.)
   ,NPT(0)
   ,dpt(0.)
   ,ProtonEnergy(0.)
{
}
//______________________________________________________________________________
int Inputparameters::Init(std::string filename) 
{


  configFileName = filename;
  ifstream inputfile;
  inputfile.open(configFileName.c_str());

  double Wmininput=0.;
  double Wmaxinput=0.;
  
  cout<<"CONSTRUCTING THE INPUT PARAMETERS"<<endl;
  
  inputfile >> Z1;	//Number of protons in #1 nucleus
  inputfile >> A1;	//Number of nucleons in #1 nucleus
  inputfile >> Z2;	//Number of protons in #2 nucleus
  inputfile >> A2;	//Number of nucleons in #2 nucleus
  inputfile >> gamma_em;  //Lorentz gamma of the collider for one beam...
  inputfile >> Wmaxinput;  //Maximum value for center of mass energy--might change depending upon which channel
  inputfile >> Wmininput;	//Minimum center of mass energy, for default values use -1
  inputfile >> numw;	//Number of center of mass bins to use
  inputfile >> Ymax;	//Sets the maximum rapidity value...more like abs. value
  inputfile >> numy;	//Number of rapidity bins
  inputfile >> gg_or_gP;  //Is it a photonphoton, photonA-narrow, or photonA-wide:1,2,3
  inputfile >> ievents;  //Number of events to generate
  inputfile >> pid;	//Channel one wants to produce.
  inputfile >> iseed;	//Number seed for the random number gen.
  inputfile >> iout;	//Output format,1:txt,2:formatted text for gstar,3:ROOT?
  inputfile >> ibreakup;	//Breakup modes, 1-...
  inputfile >> iinterfere; //Interference 0-no,1-yes
  inputfile >> xinterfere; //Percentage of interference
  inputfile >> in_or_co;	//dAu,A-A, coherent or incoherent
  inputfile >> incoherentfactor;//factor times the A*sigma(gp->vmp)
  inputfile >> bford;	//
  inputfile >> ptmax;	//
  inputfile >> NPT;	//
  inputfile.close();
  dpt = ptmax/NPT;

  ProtonEnergy=gamma_em*StarlightConstants::mp;


  //Defining Interaction Type
  if(gg_or_gP==1){
    interactiontest=StarlightConstants::PHOTONPHOTON;
  }
  if(gg_or_gP==2){
    interactiontest=StarlightConstants::PHOTONPOMERONNARROW;
  }
  if(gg_or_gP==3){
    interactiontest=StarlightConstants::PHOTONPOMERONWIDE;
  }
  
  //Trying to define the proper Wmins and Wmaxs. a TEMPORARY fix....Better solution=??
  
  double mass         = 0.0;
  double Wmin_default = 0.0; //This is the dault for Wmin, unless defined later
  double width        = 0.0;
  if(pid==11){
    //electron
    Wmin_default = 0.01; //Wmin is settable in constants--default is 0.01; Wmin up to 0.15 is safe for Summer 2000 triggering for e+e- pairs
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::ELECTRON;
    decaytest=StarlightConstants::LEPTONPAIR;
  }
  else if(pid==13){
    //muon
    Wmin_default   =  2.*StarlightConstants::mmu;
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::MUON;
    decaytest=StarlightConstants::LEPTONPAIR;
  }
  else if(pid==15){
    //tauon
    Wmin_default   =  2*StarlightConstants::mtau;
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::TAUON;
    decaytest=StarlightConstants::LEPTONPAIR;
  }
  else if(pid==115){
    //a2(1320)
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::A2;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==221){
    //eta
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::ETA;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==225){
    //f2(1270)
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::F2;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==331){
    //eta'(958)
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::ETAPRIME;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==335){
    //f2'(1525)
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::F2PRIME;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==441){
    //eta-c
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::ETAC;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==9010221){
    //f0(980)--was orginally called 10221? updated to standard number
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::F0;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==33){
    //Z"/Z03
    Wmax=Wmaxinput;
    pidtest=StarlightConstants::ZOVERZ03;
    decaytest=StarlightConstants::SINGLEMESON;
  }
  else if(pid==113)   {
    mass        =  0.7685;                 //rho(770)0
    width       =  0.1507;
    Wmin_default=  2.*StarlightConstants::mpi;
    Wmax        =  mass+5.*width;
    pidtest=StarlightConstants::RHO;
    decaytest=StarlightConstants::WIDEVMDEFAULT;
  }
  else if(pid==913)   {
    mass        =  0.7685;                      //rho(770)0 direct pi+pi- decay, interference given by ZEUS
    width       =  0.1507;
    Wmin_default=  2.*StarlightConstants::mpi;
    Wmax        =  mass+5.*width;    //Use the same (ZEUS) 1.5GeV max mass
    pidtest=StarlightConstants::RHOZEUS;
    decaytest=StarlightConstants::WIDEVMDEFAULT;
  }
  else if (pid == 999) {  // pi+pi-pi+pi- decay
    mass         = 1.350;
    width        = 0.360;
    Wmin_default = 4 * StarlightConstants::mpi;
    Wmax         = 3;
    pidtest      = StarlightConstants::FOURPRONG;
    decaytest    = StarlightConstants::WIDEVMDEFAULT;
  }
  else if(pid==223){
    mass          =  0.78194;          //omega(782)
    width         =  0.00843;
    Wmax          =  mass+5.*width;
    Wmin_default  =  mass - 5.*width;
    pidtest=StarlightConstants::OMEGA;
    decaytest=StarlightConstants::NARROWVMDEFAULT;//will probably be shifted to 3body decay
  }
  else if(pid==333) {
    mass          =  1.019413;     //phi(1020)
    width         =  0.00443;
    Wmin_default  =  2*StarlightConstants::mK;
    Wmax          =  mass+5.*width;
    pidtest=StarlightConstants::PHI;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==443){
    mass          =  3.09692;  // JN  3.09688;      //J/psi
    width         =  0.000091; //JN  0.000087;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::JPSI;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==444){
    mass          =  3.686093;      //J/psi
    width         =  0.000337;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::JPSI2S;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==444011){
    mass          =  3.686093;      //J/psi
    width         =  0.000337;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::JPSI2S_ee;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==444013){
    mass          =  3.686093;      //J/psi
    width         =  0.000337;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::JPSI2S_mumu;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==443011){
    mass          =  3.09692;  // JN  3.09688;      //J/psi
    width         =  0.000091; //JN  0.000087;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::JPSI_ee;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==443013){
    mass          =  3.09692;  // JN  3.09688;      //J/psi
    width         =  0.000091; //JN  0.000087;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::JPSI_mumu;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==553){
    mass          =  9.46030;      //Upsilon
    width         =  0.00005402;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==553011){
    mass          =  9.46030;      //Upsilon
    width         =  0.00005402;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON_ee;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==553013){
    mass          =  9.46030;      //Upsilon
    width         =  0.00005402;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON_mumu;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==554){
    mass          =  10.02326;      //Upsilon(2S)
    width         =  0.00003198;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON2S;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==554011){
    mass          =  10.02326;      //Upsilon(2S)
    width         =  0.00003198;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON2S_ee;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==554013){
    mass          =  10.02326;      //Upsilon(2S)
    width         =  0.00003198;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON2S_mumu;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==555){
    mass          =  10.3552;      //Upsilon(3S)
    width         =  0.00002032;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON3S;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==555011){
    mass          =  10.3552;      //Upsilon(3S)
    width         =  0.00002032;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON3S_ee;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else if(pid==555013){
    mass          =  10.3552;      //Upsilon(3S)
    width         =  0.00002032;
    Wmin_default  =  mass - 5.*width;
    Wmax          =  mass + 5.*width;
    pidtest=StarlightConstants::UPSILON3S_mumu;
    decaytest=StarlightConstants::NARROWVMDEFAULT;
  }
  else{
    Wmax=Wmaxinput;
  }
  if (Wmininput==-1){ 
    Wmin = Wmin_default;
  }
  else{  
    Wmin = Wmininput;
  }

  return 0;
}
//______________________________________________________________________________
int Inputparameters::getZ1()
{
  return Z1;
}
//______________________________________________________________________________
int Inputparameters::getA1()
{
  return A1;
}
//______________________________________________________________________________
int Inputparameters::getZ2()
{
  return Z2;
}
//______________________________________________________________________________
int Inputparameters::getA2()
{
  return A2;
}
//______________________________________________________________________________
double Inputparameters::getgamma_em()
{
  return gamma_em;
}
//______________________________________________________________________________
double Inputparameters::getWmax()
{
  return Wmax;
}
//______________________________________________________________________________
double Inputparameters::getWmin()
{
  return Wmin;
}
//______________________________________________________________________________
int Inputparameters::getnumw()
{
  return numw;
}
//______________________________________________________________________________
double Inputparameters::getYmax()
{
  return Ymax;
}
//______________________________________________________________________________
int Inputparameters::getnumy()
{
  return numy;
}
//______________________________________________________________________________
int Inputparameters::getgg_or_gP()
{
  return gg_or_gP;
}
//______________________________________________________________________________
int Inputparameters::getnumberofevents()
{
  return ievents;
}
//______________________________________________________________________________
int Inputparameters::getseed()
{
  return iseed;
}
//______________________________________________________________________________
int Inputparameters::getparticleid()
{
  return pid;
}
//______________________________________________________________________________
int Inputparameters::getoutputmode()
{
  return iout;
}
//______________________________________________________________________________
int Inputparameters::getbreakupmode()
{
  return ibreakup;
}
//______________________________________________________________________________
int Inputparameters::getinterferencemode()
{
  return iinterfere;
}
//______________________________________________________________________________
double Inputparameters::getinterferencepercent()
{
  return xinterfere;
}
//______________________________________________________________________________
int Inputparameters::getincoherentorcoherent()
{
  return in_or_co;
}
//______________________________________________________________________________
double Inputparameters::getincoherentfactor()
{
  return incoherentfactor;
}
//______________________________________________________________________________
double Inputparameters::getbford()
{
  return bford;
}
//______________________________________________________________________________
double Inputparameters::getmaximuminterpt()
{
  return ptmax;
}
//______________________________________________________________________________
int Inputparameters::getNPT()
{
  return NPT;
}
//______________________________________________________________________________
double Inputparameters::getdpt()
{
  return dpt;
}
//______________________________________________________________________________
StarlightConstants::particle Inputparameters::getpidtest()
{
  return pidtest;
}
//______________________________________________________________________________
StarlightConstants::decaytype Inputparameters::getdecaytest()
{
  return decaytest;
}
//______________________________________________________________________________
StarlightConstants::interactiontype Inputparameters::getinteractiontest()
{
  return interactiontest;
}
//______________________________________________________________________________
double Inputparameters::getProtonEnergy()
{
  return ProtonEnergy;
}
//______________________________________________________________________________
Inputparameters::~Inputparameters()
{

}
