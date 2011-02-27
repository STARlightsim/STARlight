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


#ifndef PHOTONNUCLEUSCROSSSECTION_H
#define PHOTONNUCLEUSCROSSSECTION_H


#include "starlightconstants.h"
#include "beambeamsystem.h"


class photonNucleusCrossSection {

public:

	photonNucleusCrossSection(const inputParameters& input,
	                          const beamBeamSystem&  bbsystem);
	~photonNucleusCrossSection();
  
	double                getbslope         () const { return _bSlope;      }
	double                getChannelMass    () const { return _channelMass; }
	double                getBNORM          () const { return _BNORM;       }
	double                getLum            () const { return _lum;         }
	beamBeamSystem        getbbs            () const { return _bbs;         }
	double                getf2o4pi         () const { return _f2o4pi;      }
	double                getDefaultC       () const { return _defaultC;    }
	double                getMaxPhotonEnergy() const { return _EgMax;       }

	void crossSectionCalculation(const double bwnormsave);
	// Will think about it...For VMs we just calculate it
	// So just use the wide or narrow constructor to calculate it
	// wide/narrow will inherit this.
	double getcsgA(const double Egamma,
	               const double W);
	double photonFlux(const double Egamma);
	double sigmagp(const double Wgp);
	double sigma_A(const double sig_N);
	double breitWigner(const double W,
	                   const double C);

private:

	beamBeamSystem _bbs;
  
	double nepoint(const double Egamma,
	               const double bmin);
  
	starlightConstants::particleTypeEnum _sigmaPID;
	double                               _sigmaProtonEnergy;
	double                               _sigmaGamma_em;
	int                                  _sigmaBreakup;
	double                               _bSlope;
	double                               _f2o4pi;
	double                               _ANORM;
	double                               _BNORM;
	double                               _defaultC;
	double                               _channelMass;
	double                               _lum;
	double                               _EgMax;
	double                               _width;
	int                                  _sigmaCoherence;  // 1=coherent, 0=incoherent
	double                               _sigmaCoherenceFactor;
	int                                  _sigmaNucleus;
};


#endif  // PHOTONNUCLEUSCROSSSECTION_H
