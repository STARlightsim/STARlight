///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef EVENTCHANNEL_H
#define EVENTCHANNEL_H

#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "upcevent.h"



class eventChannel : public readLuminosity
{
	public:
		eventChannel(inputParameters& input, beamBeamSystem& bbsystem);
		virtual ~eventChannel();

		virtual starlightConstants::event produceEvent(int &ievent)= 0;

		virtual upcEvent produceEvent() = 0;
 
		void transform(double betax,double betay,double betaz,double &E,
                              double &px,double &py,double &pz,int &iFbadevent);
		randomGenerator _randy;
		beamBeamSystem _bbs;
		unsigned long nmbAttempts() const {return _nmbAttempts;}
		unsigned long nmbAccepted() const {return _nmbAccepted;}
		double pseudoRapidity(double px, double py, double pz);
 protected:
		unsigned long _nmbAttempts;
		unsigned long _nmbAccepted;
		bool _accCutPt; //Acceptance cut in pT
		bool _accCutEta; //Acceptance cut in pseudorapidity
		double _ptMin; //Min pT if cut
		double _ptMax; //Max pT if cut
		double _etaMin; //Min eta if cut
		double _etaMax; //Max eta if cut
		
};


#endif  // EVENTCHANNEL_H
