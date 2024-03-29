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
// $Rev:: 293                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//    Class needed for root output
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "eventchannel.h"


using namespace std;


//______________________________________________________________________________
eventChannel::eventChannel(const inputParameters& inputParametersInstance, randomGenerator* randy, beamBeamSystem& bbsystem)
	: readLuminosity(inputParametersInstance),
	  _bbs(bbsystem),
	  _randy(randy),
	  _nmbAttempts(0),
	  _nmbAccepted(0),
	  _totalChannelCrossSection(0)
{
  _ptCutEnabled  = inputParametersInstance.ptCutEnabled();
  _ptCutMin      = inputParametersInstance.ptCutMin();
  _ptCutMax      = inputParametersInstance.ptCutMax();
  _etaCutEnabled = inputParametersInstance.etaCutEnabled();
  _etaCutMin     = inputParametersInstance.etaCutMin();
  _etaCutMax     = inputParametersInstance.etaCutMax();
}


//______________________________________________________________________________
eventChannel::~eventChannel()
{ }


//______________________________________________________________________________
void
eventChannel::transform(const double  betax,
                        const double  betay,
                        const double  betaz,
                        double&       E,
                        double&       px,
                        double&       py,
                        double&       pz,
                        int&          iFbadevent)
{
  // carries out a lorentz transform of the frame.  (Not a boost!)???
  const double E0  = E;
  const double px0 = px;
  const double py0 = py;
  const double pz0 = pz;

  const double beta = sqrt(betax * betax + betay * betay + betaz * betaz);
  if (beta >= 1)
	  iFbadevent = 1;
  const double gamma = 1. / sqrt(1. - beta * beta);
  const double gob   = (gamma - 1) / (beta * beta);

  E   = gamma * (E0 - betax * px0 - betay * py0 - betaz*  pz0);
  px  = -gamma * betax * E0 + (1. + gob * betax * betax) * px0
	  + gob * betax * betay * py0 + gob * betax * betaz * pz0;
  py  = -gamma * betay * E0 + gob * betay * betax * px0
	  + (1. + gob * betay * betay) * py0 + gob * betay * betaz  *pz0;
  pz  = -gamma * betaz * E0 + gob * betaz * betax * px0
	  + gob * betaz * betay * py0 + (1. + gob * betaz * betaz) * pz0;
}


/**
 * @brief Determine the pseudorapidity of a particle
 * 
 * @param px The x-momentum of the particle
 * @param py The y-momentum
 * @param pz The z-momentum
 * @return [double]: The pseudorapidity of the particle.
 */
double
eventChannel::pseudoRapidity(const double px,
                             const double py,
                             const double pz)
{
  const double pT= sqrt(px * px + py * py);
  const double p = sqrt(pz * pz + pT * pT);
  double eta = 1.0 * pow(10.0,20);  // instead of special value, std::numeric_limits<double>::quiet_NaN() should be used
  if ((p - pz) != 0)//This is to avoid the division by zero error.
	  eta = 0.5 * log((p + pz)/(p - pz));
  else
  {
    if (pz<0)//This is to consider and distinguish between positive and negative infinity
      eta = eta*-1.0;
  }
  return eta;
}

/**
 * @brief Determines the pseudorapidity of a particle in Laboratory frame of reference.
 * 
 * @param px The x-momemtum of the particle in the CM frame
 * @param py The y-momentum of the particle in the CM frame
 * @param pz The z-momentum of the particle in the CM frame
 * @param E  The energy of the particle in the CM frame
 * @param beta The boost vector to transform the particle from CM frame to Lab Frame
 * @return double: The Lab frame pseudorapidity of the particle
 */
double eventChannel::pseudoRapidityLab(const double px,
								                       const double py,
								                       const double pz,
			                        				 const double E,
      								                 const vector3 beta)
{
  lorentzVector vecme(px,py,pz,E);// A lorentz vector is created from the particle
  vecme.Boost(beta);// The vector is boosted from CM frame to Lab frame
  return pseudoRapidity(vecme.GetPx(),vecme.GetPy(),vecme.GetPz());//The pseudorapidity in this new Frame (Lab) is computed and returned.
}