#ifndef GAMMAALUMINOSITY_H
#define GAMMAALUMINOSITY_H

#include "beambeamsystem.h"
#include "inputparameters.h"
#include "gammaacrosssection.h"                                                                                                                            


class Gammaaluminosity: Gammaacrosssection
{

                                                                                                                                           
 public:
  Gammaaluminosity(Inputparameters& input, Beambeamsystem& bbsystem);
  ~Gammaaluminosity();
  
 private:
  Inputparameters inputgammaa;
  void gammaadifferentialluminosity();
  double *vmsigmapt(double W,double Egamma,double *SIGMAPT);
  double nofe(double Egamma,double bimp);
  void pttablegen();
};

#endif //GAMMAALUMINOSITY_H

