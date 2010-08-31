#ifndef GAMMAAVM_H
#define GAMMAAVM_H

#include <vector>
#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "eventchannel.h"
#include "upcevent.h"

class Gammaavectormeson:public Eventchannel//:public Readinluminosity
{
  
 public:
  Gammaavectormeson(Inputparameters& input,Beambeamsystem& bbsystem);
  virtual ~Gammaavectormeson();
  StarlightConstants::event produceevent(int &ievent);
  
   UPCEvent ProduceEvent();

  void pickwy(double &W, double &Y);
  void momenta(double W,double Y,double &E,double &px,double &py,double &pz,int &tcheck);
  void vmpt(double W,double Y,double &E,double &px,double &py, double &pz,int &tcheck);
  void twodecay(StarlightConstants::particle &ipid,double E,double W,double px0,double py0,double pz0,double &px1,double &py1,double&pz1,double &px2,double &py2,double &pz2,int &iFbadevent);
  double getmass();
  double getwidth();
  virtual double gettheta(StarlightConstants::particle ipid);
  double getspin();
  double VMbslope;
  virtual double getdaughtermass(StarlightConstants::particle &ipid);                
  
 private:
  StarlightConstants::particle VMpidtest;
  int VMnumw;
  int VMnumy;
  int VMinterferencemode;
  int VMCoherence;
  double VMCoherenceFactor;
  double VMgamma_em;
  double VMNPT;
  double VMWmax;
  double VMWmin;
  double VMYmax;
  double VMYmin;
  double mass;
  double width;
  double VMptmax;
  double VMdpt;
};
class Gammaanarrowvm:public Gammaavectormeson
{
 public:
  Gammaanarrowvm(Inputparameters& input,Beambeamsystem& bbsystem);
  virtual ~Gammaanarrowvm();
                                                                                                                                
};
class Gammaawidevm:public Gammaavectormeson
{  
 public:
  Gammaawidevm(Inputparameters& input,Beambeamsystem& bbsystem);
  virtual ~Gammaawidevm();
};
#endif //GAMMAAVM_H
