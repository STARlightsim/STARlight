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
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef GAMMAGAMMASINGLE_H
#define GAMMAGAMMASINGLE_H


#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "eventchannel.h"
#include "starlightconfig.h"
#include "inputParameters.h"
#ifdef ENABLE_PYTHIA
#include "pythiadecayer.h"
#endif

class Gammagammasingle : public eventChannel
{
 public:
  Gammagammasingle(const inputParameters& input, randomGenerator* randy, beamBeamSystem& bbsystem);
  ~Gammagammasingle();
  
  void singleCrossSection();
  starlightConstants::event produceEvent(int &ievent);

  //upcEvent produceEvent(vector3 beta);
  upcXEvent produceEvent(vector3 beta);

 private:
  double _sigmax[starlightLimits::MAXWBINS][starlightLimits::MAXYBINS];
  double _wdelta;  //Added 7/26/07 for passing sigmadelta to pickw
  double _remainwd;// "
  int _ivalwd;     // "
  
  void pickw(double &w);
  void picky(double &y);
  
  void parentMomentum(double w,double y,double &E,double &px,double &py,double&pz,
                    double &Eb1, double &pxb1, double &pyb1, double &pzb1,
								    double &Eb2, double &pxb2, double &pyb2, double &pzb2, double &t2,
								    double &Egam1, double&pxgam1, double &pygam1, double &pzgam1, double &Q2gam1,
                    double &Egam2, double&pxgam2, double &pygam2, double &pzgam2, double &Q2gam2);

  // Duplicated old function pp into pp1 and pp2 to handle the two nuclei separately, allowing for asymmetric species
  double pp1(double E);
  double pp2(double E);
  void twoBodyDecay(starlightConstants::particleTypeEnum &ipid,double W,double px0,double py0,double pz0,double &E1, double &px1,double &py1,double&pz1, double &E2, double &px2,double &py2,double &pz2,double &mass,int &iFbadevent);
  
  double getMass();
  double getWidth();
  double getSpin();
  
  starlightConstants::particleTypeEnum _GGsingInputpidtest;
  int _GGsingInputnumw;
  int _GGsingInputnumy;
  double _GGsingInputGamma_em;
  double _axionMass; // AXION HACK
//  double _axionWidth;// AXION HACK--unused
#ifdef ENABLE_PYTHIA 
  pythiaDecayer _pyDecayer;
#endif
  
};


#endif  // GAMMAGAMMASINGLE_H
