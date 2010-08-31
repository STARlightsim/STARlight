#ifndef PSIFAMILY_H
#define PSIFAMILY_H

#include <vector>
#include "starlightconstants.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "gammaavm.h"

class Psifamily:public Gammaanarrowvm
{
	public:
		Psifamily(Inputparameters& input, Beambeamsystem& bbsystem);
		~Psifamily();
		double gettheta(StarlightConstants::particle ipid);
		double getdaughtermass(StarlightConstants::particle &ipid);
	private:
		double width;
		double mass;
};

#endif //PSIFAMILY_H
