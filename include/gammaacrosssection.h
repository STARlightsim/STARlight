#ifndef GAMMAACROSSSECTION_H
#define GAMMAACROSSSECTION_H
                                                                                                                                                      
#include "starlightconstants.h"
#include "beambeamsystem.h"

class Gammaacrosssection
{
  
 public:
  Gammaacrosssection(Inputparameters& input,Beambeamsystem& bbsystem);
  ~Gammaacrosssection();
  
  double getcsgA(double Egamma,double W);
  double getbslope();
  //Made public since it is used in checking values for pp elsewhere

  void crosssectioncalculation(double bwnormsave);
  // Will think about it...For VMs we just calculate it
  // So just use the wide or narrow constructor to calculate it
  // wide/narrow will inherit this.
		
  double getchannelmass();
  double getBNORM();
  double getlum();
  double photonflux(double Egamma);
  double getMaxPhotonEnergy();
  double getdefaultC();
  double breitwigner(double W,double C);
  Beambeamsystem getbbs();	
  double sigmagp(double Wgp);
  double sigma_A(double sig_N);
  double getf2o4pi();

 private:
  Beambeamsystem bbs;
  
  double nepoint(double Egamma, double bmin);
  
  StarlightConstants::particle SigmaPID;
  double SigmaProtonEnergy;
  double SigmaGamma_em;
  int SigmaBreakup;
  double bslope;
  double f2o4pi;
  double ANORM;
  double BNORM;
  double defaultC;
  double channelmass;
  double lum;
  double EgMax;
  double width;
  int SigmaCoherence; //1=coherent, 0=incoherent
  double SigmaCoherenceFactor;
  int SigmaNucleus;
};

//Now let's define narrow and wide, which derive from here.
class Wideresonancesigma:public Gammaacrosssection
{

 public:
  Wideresonancesigma(Inputparameters& input,Beambeamsystem& bbsystem);
  ~Wideresonancesigma();
  void crosssectioncalculation(double bwnormsave);
 private:
  double Ep;//Proton Energy
  double WideWmax;
  double WideWmin;
  double WideYmax;
  double WideYmin;		
};

class Narrowresonancesigma:public Gammaacrosssection
{
 public:
  Narrowresonancesigma(Inputparameters& input,Beambeamsystem& bbsystem);
  ~Narrowresonancesigma();
  void crosssectioncalculation(double bwnormsave);
 private:
  double NarrowYmax;
  double NarrowYmin;
  int NarrowNumY;
  double Ep;
  
};

#endif //GAMMAACROSSSECTION_H

