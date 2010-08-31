class PythiaStarlight;


#ifndef GAMMAGAMMASINGLE_H
#define GAMMAGAMMASINGLE_H

#include <vector>
#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "eventchannel.h"
class Gammagammasingle:  public Eventchannel
{
  
 public:
  Gammagammasingle(Inputparameters& input,Beambeamsystem& bbsystem);
  ~Gammagammasingle();
  
  void singlecrosssection();
  StarlightConstants::event produceevent(int &ievent);

  UPCEvent ProduceEvent();

 private:
  double sigmax[StarlightLimits::MAXWBINS][StarlightLimits::MAXYBINS];//=new double[500][500];   //decreased from 1000*1000; too big! causes fault!
  double sigmasum;
  double wdelta;  //Added 7/26/07 for passing sigmadelta to pickw
  double remainwd;// "
  int ivalwd;     // "
  
  void pickw(double &w);
  void picky(double &y);
  
  void parentmomentum(double w,double y,double &E,double &px,double &py,double&pz);
  double pp(double E);
  void twodecay(StarlightConstants::particle &ipid,double E,double W,double px0,double py0,double pz0,double &px1,double &py1,double&pz1,double &px2,double &py2,/*double &py2,*/double &pz2,int &iFbadevent);
  // void transform(double betax,double betay,double betaz,double &E,double &px,double &py,double &pz,int &iFbadevent);
  void thephi(double W,double px,double py,double pz,double E,double &theta,double &phi);
  
  double getmass();
  double getwidth();
  double getspin();
  
  StarlightConstants::particle GGsingInputpidtest;
  int GGsingInputnumw;
  int GGsingInputnumy;
  double GGsingInputGamma_em;
  
  PythiaStarlight *pythia;

};
#endif //GAMMAGAMMASINGLE_H

