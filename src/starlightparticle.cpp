// starlightparticle.cpp
/*
 * $Id: starlightiparticle.cpp,v 1.0 2010/07/04  $
 *
 *
 * /author 
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * $Log: $
 */

#include "starlightparticle.h"

StarlightParticle::StarlightParticle() :
   LorentzVector()
   ,fPdgCode(0)
   ,fCharge(-999)
   ,fMass(-1)
{

}

StarlightParticle::StarlightParticle(double px, double py, double pz, double e, double mass, int pdgCode, short charge) :
LorentzVector(px, py, pz, e)
,fPdgCode(pdgCode)
,fCharge(charge)
,fMass(mass)
{
}

StarlightParticle::~StarlightParticle()
{
      
}

