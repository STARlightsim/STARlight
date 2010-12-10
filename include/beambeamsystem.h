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


#ifndef BEAMBEAMSYSTEM_H
#define BEAMBEAMSYSTEM_H


// this class covers a coliding beam system SK
#include "nucleus.h"
#include "beam.h"


class beamBeamSystem
{
 public:
	//Better way to do this? Memory issues creating all of theses Beams?
	beamBeamSystem(beam& beam_1,beam& beam_2,double luminosity,inputParameters& input);
	beamBeamSystem(beam& beam_1, beam& beam_2,inputParameters& input);
	beamBeamSystem(inputParameters &input);
    ~beamBeamSystem();

    beam getBeam1();
    beam getBeam2();
//	double getluminosity();
    double probabilityOfBreakup(double D);
    /*double probabilityOfHadronBreakup(double impactparameter);
    double probabilityOfPhotonBreakup(double impactparameter,int mode);
    */ //mode is for which type of breakup desired, 1=hardsphere.

 private:
//	int _ibreakup;//temporary solution until read in parameters are done
    double probabilityOfHadronBreakup(double impactparameter);
    double probabilityOfPhotonBreakup(double impactparameter,int mode);

    double _pHadronBreakup;
    double _pPhotonBreakup;
    //inputParameters inputbbs;
    double _BBSInputGamma_em;
    int _BBSInputBreakupmode;
    beam _beam1;
    beam _beam2;
//		double luminosity;
};


#endif  // BEAMBEAMSYSTEM_H
