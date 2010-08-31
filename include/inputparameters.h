#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H


#include "starlightconstants.h"
//This is where we read in our input values.

class Inputparameters
{

 public:
  Inputparameters();
  ~Inputparameters();

  int Init(std::string filename = "config/slight.in");

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
  int getnumberofevents();
  int getparticleid();
  int getseed();
  int getoutputmode();
  int getbreakupmode();
  int getinterferencemode();
  double getinterferencepercent();
  int getincoherentorcoherent();
  double getincoherentfactor();
  double getbford();
  double getmaximuminterpt();
  int getNPT();
  double getdpt();
  StarlightConstants::particle getpidtest();
  StarlightConstants::decaytype getdecaytest();
  StarlightConstants::interactiontype getinteractiontest();
  double getf2o4pi();
  double getbslope();
  double getProtonEnergy();
  
 private:

  /** Name of the configuration file (default: slight.in) */
  std::string configFileName; 

  int    Z1;          // Atomic Number
  int    A1;          //Atomic Mass
  int    Z2;
  int    A2;
  double gamma_em;   //gamma for the colliding ions
  double Wmax;       //gamma-gamma center of mass Energy, max
  double Wmin;       //min for Energy, -1 is default
  int    numw;          //number of bins for energy in the calculations
  double Ymax;       //max value for the rapidity
  int    numy;       //number of bins used
  int    gg_or_gP;   //1=2 photon channels,2=vector meson with narrow resonance,3=vector meson with wide resonance(B-W)
  int    ievents;          //number of events to generate
  int    pid;           //pdg ID --channel choice //ip
  int    iseed;        //random number seed
  int    iout;             //output format,1=text,2=GSTARtext,3=ntuple
  int    ibreakup;             //1=hard sphere nuclei (b>2),2=both nuclei break up (XnXn),3=a single neutron from each nucleus (1n1n),
  //4=require that neither nucleon break up (with b>2R),5=require that there be no hadronic break up (This is similar to option one, but with the actual hadronic interaction)
  int    iinterfere;             //0= no interference, 1= interference
  double xinterfere;            //When there is interference, this gives the %age of interference, 0=none to 1=full
  int    in_or_co;//0 for incoherent 1 for coherent
  double incoherentfactor; //Variable allows the user to scale the incoherent contribution in VM production
  double bford;
  double ptmax;          //When there is interference, this is the max pt considered
  int    NPT;           //When there is interference, this is the number of pt bins
  double dpt;
  StarlightConstants::particle pidtest;   //Testing the starlightconstants particle enumeration here...
  StarlightConstants::decaytype decaytest;
  StarlightConstants::interactiontype interactiontest;
  double ProtonEnergy;
};
#endif // INPUTPARAMETERS_H
