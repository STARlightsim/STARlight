#ifndef TWOPHOTONLUMINOSITY_H
#define TWOPHOTONLUMINOSITY_H

#include "nucleus.h"
#include "beam.h"
#include "beambeamsystem.h"
#include "starlightlimits.h"
class Twophotonluminosity: Beambeamsystem
{
  
 public:
  Twophotonluminosity(Beam beam_1,Beam beam_2,int mode,double luminosity,Inputparameters& input);
  Twophotonluminosity(Beam beam_1, Beam beam_2,int mode,Inputparameters& input);
  ~Twophotonluminosity();
  
 private:
  void twophotondifferentialluminosity();
  double D2LDMDY(double M,double Y,double &Normalize);
  double Integral(double Normalize);
  double radmul(int N,double *Lower,double *Upper,int NIterMin,int NIterMax,double EPS,double *WK,int NIter,double &Result,double &ResErr,double &NFNEVL,double &Summary);
  double integrand(double N,double X[15]);
  double Nphoton(double W,double gamma,double Rho);
  
  double W1; //Energy of photon #1
  double W2; //Energy of photon #2
  double gamma; //Gamma of the system
  Inputparameters input2photon;
};
#endif //TWOPHOTONLUMINOSITY_H
