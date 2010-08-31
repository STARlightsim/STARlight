#ifndef GAMMAGAMMALEPTONPAIR_H
#define GAMMAGAMMALEPTONPAIR_H

#include <vector>
#include "starlightconstants.h"
#include "readinluminosity.h"
#include "starlightlimits.h"
#include "eventchannel.h"

class Gammagammaleptonpair:  public Eventchannel
{
  
 public:
  Gammagammaleptonpair(Inputparameters& input,Beambeamsystem& bbsystem);
  ~Gammagammaleptonpair();
  
  void twoleptoncrosssection();
  void tablecalc();

  StarlightConstants::event produceevent(int &ievent);
  UPCEvent ProduceEvent();

 private:
  double sigmax[StarlightLimits::MAXWBINS][StarlightLimits::MAXYBINS];//=new double[500][500];   //decreased from 1000*1000; too big! causes fault!
  double sigmasum;
  double sigfint[StarlightLimits::MAXWBINS];
  double sigofw[StarlightLimits::MAXWBINS];
  double signormw;
  double wdelta;  //Added 7/26/07 for passing sigmadelta to pickw
  double remainwd;// "
  int ivalwd;     // "
  double dgammade[1000];
  double tautolangle[100];
  
  double twomuoncrosssection(double w);
  void pickw(double &w);
  void picky(double &y);
  
  void pairmomentum(double w,double y,double &E,double &px,double &py,double&pz);
  double pp(double E);
  void twodecay(StarlightConstants::particle &ipid,double E,double W,double px0,double py0,double pz0,double &px1,double &py1,double&pz1,double &px2,double &py2,/*double &py2,*/double &pz2,int &iFbadevent);
  double thetalep(double W,double theta);
  void taudecay(double &px1,double &py1,double &pz1,double &E1,double &px2,double &py2,double &pz2,double &E2);
  
  double getmass();
  double getwidth();
  double getspin();
  
  StarlightConstants::particle GGlepInputpidtest;
  int GGlepInputnumw;
  int GGlepInputnumy;
  double GGlepInputGamma_em;
};
#endif //GAMMAGAMMALEPTONPAIR_H
