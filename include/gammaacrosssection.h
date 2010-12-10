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


#ifndef GAMMAACROSSSECTION_H
#define GAMMAACROSSSECTION_H


#include "starlightconstants.h"
#include "beambeamsystem.h"


class photonNucleusCrossSection
{
 public:

  photonNucleusCrossSection(inputParameters& input, beamBeamSystem& bbsystem);
  ~photonNucleusCrossSection();
  
  double getcsgA(double Egamma,double W);
  double getbslope();
  // Made public since it is used in checking values for pp elsewhere

  void crossSectionCalculation(double bwnormsave);
  // Will think about it...For VMs we just calculate it
  // So just use the wide or narrow constructor to calculate it
  // wide/narrow will inherit this.
		
  double getChannelMass();
  double getBNORM();
  double getLum();
  double photonFlux(double Egamma);
  double getMaxPhotonEnergy();
  double getDefaultC();
  double breitWigner(double W,double C);
  beamBeamSystem getbbs();	
  double sigmagp(double Wgp);
  double sigma_A(double sig_N);
  double getf2o4pi();

 private:
  beamBeamSystem _bbs;
  
  double nepoint(double Egamma, double bmin);
  
  starlightConstants::particle _sigmaPID;
  double _sigmaProtonEnergy;
  double _sigmaGamma_em;
  int _sigmaBreakup;
  double _bSlope;
  double _f2o4pi;
  double _ANORM;
  double _BNORM;
  double _defaultC;
  double _channelMass;
  double _lum;
  double _EgMax;
  double _width;
  int _sigmaCoherence;  // 1=coherent, 0=incoherent
  double _sigmaCoherenceFactor;
  int _sigmaNucleus;
};


// Now let's define narrow and wide, which derive from here.
class wideResonanceCrossSection : public photonNucleusCrossSection
{
 public:
  wideResonanceCrossSection(inputParameters& input, beamBeamSystem& bbsystem);
  ~wideResonanceCrossSection();
  void crossSectionCalculation(double bwnormsave);

 private:
  double _Ep;  // Proton Energy
  double _wideWmax;
  double _wideWmin;
  double _wideYmax;
  double _wideYmin;		
};


class narrowResonanceCrossSection : public photonNucleusCrossSection
{
 public:
  narrowResonanceCrossSection(inputParameters& input,beamBeamSystem& bbsystem);
  ~narrowResonanceCrossSection();
  void crossSectionCalculation(double bwnormsave);

 private:
  double _narrowYmax;
  double _narrowYmin;
  int _narrowNumY;
  double _Ep;
};


#endif  // GAMMAACROSSSECTION_H
