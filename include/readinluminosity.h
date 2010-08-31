#ifndef READINLUMINOSITY_H
#define READINLUMINOSITY_H

#include "inputparameters.h"
#include "starlightlimits.h"

class Readinluminosity
{

 public:
  Readinluminosity(Inputparameters& input);
  ~Readinluminosity();
  
  void read();
  double Warray[StarlightLimits::MAXWBINS];   //decreased from 1000; too big! causes fault!
  double Yarray[StarlightLimits::MAXYBINS];    
  double Farray[StarlightLimits::MAXWBINS][StarlightLimits::MAXYBINS];
  double f_max;
  double fptarray[500][500];
  //		Inputparameters inputread;
  double bwnormsave;

 protected:
  int ReadInputNPT;
  int ReadInputnumy;
  int ReadInputnumw;
  int ReadInputgg_or_gP;
  int ReadInputinterferencemode;
};
#endif //READINLUMINOSITY_H
