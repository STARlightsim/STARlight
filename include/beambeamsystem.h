#ifndef BEAMBEAMSYSTEM_H
#define BEAMBEAMSYSTEM_H

//This class covers a coliding beam system SK
#include "nucleus.h"
#include "beam.h"

class Beambeamsystem
{

 public:
	//Better way to do this? Memory issues creating all of theses Beams?
	Beambeamsystem(Beam& beam_1,Beam& beam_2,double luminosity,Inputparameters& input);
	Beambeamsystem(Beam& beam_1, Beam& beam_2,Inputparameters& input);
	Beambeamsystem(Inputparameters &input);
//	Beam Beam1;
//	Beam Beam2;
    Beam getBeam1();
    Beam getBeam2();
//	double getluminosity();
    double probabilityofbreakup(double D);
    /*double probabilityofhadronbreakup(double impactparameter);
    double probabilityofphotonbreakup(double impactparameter,int mode);
    */ //mode is for which type of breakup desired, 1=hardsphere.
    ~Beambeamsystem();
    //Inputparameters inputbbs;

 private:
//	int ibreakup;//temporary solution until read in parameters are done
    double probabilityofhadronbreakup(double impactparameter);
    double probabilityofphotonbreakup(double impactparameter,int mode);

    double phadronbreakup;
    double pphotonbreakup;
    //Inputparameters inputbbs;
    double BBSInputGamma_em;
    int BBSInputBreakupmode;
    Beam Beam1;
    Beam Beam2;
//		double luminosity;
};

#endif //BEAMBEAMSYSTEM_H
