#ifndef EVENTCHANNEL_H
#define EVENTCHANNEL_H

#include <vector>
#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "upcevent.h"



class Eventchannel:public Readinluminosity
{

	public:
		Eventchannel(Inputparameters& input, Beambeamsystem& bbsystem);
		virtual ~Eventchannel();
		virtual StarlightConstants::event produceevent(int &ievent)= 0;

		virtual UPCEvent ProduceEvent() = 0;
 
		void transform(double betax,double betay,double betaz,double &E,
                              double &px,double &py,double &pz,int &iFbadevent);
		Randomgenerator Randy;
		Beambeamsystem bbs;
};

#endif //EVENTCHANNEL_H
