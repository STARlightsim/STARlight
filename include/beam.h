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


#ifndef BEAM_H
#define BEAM_H


//This calls inclues a single beam of nucleons
#include "nucleus.h"
#include "inputParameters.h"


class beam : public nucleus
{
	public:
		beam(int Zin, int Ain, double bdeuteron, int in_or_co, inputParameters& input);
		~beam();
		double nofe(double impactparameter);  //photon density
		double _photonEnergy;
	protected:
		//inputParameters inputbeam;
		double _beamInputGamma_em;
};


#endif  // BEAM_H
