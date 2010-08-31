#ifndef BEAM_H
#define BEAM_H
//This calls inclues a single beam of nucleons
#include "nucleus.h"
#include "inputparameters.h"

class Beam: public Nucleus
{
	public:
		Beam(int Zin, int Ain, double bdeuteron, int in_or_co, Inputparameters& input);
		double nofe(double impactparameter);//photon density
		double Egamma;
		~Beam();
	protected:
		//Inputparameters inputbeam;
		double BeamInputGamma_em;
};

#endif //BEAM_H
