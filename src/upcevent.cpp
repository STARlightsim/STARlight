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
// $Rev:: 263                         $: revision of last commit
// $Author:: butter                   $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#include "upcevent.h"


upcEvent::upcEvent() :
        _particles(0)
        ,_vertices(0)
{ }
/**
 * @brief Construct a new upc Event from a starlightConstants::event struct.
 * 
 * @param ev The starlightConstants::event struct
 */
upcEvent::upcEvent(starlightConstants::event &ev) :
        _particles(0)
        ,_vertices(0)
{
  for(int i = 0; i < ev._numberOfTracks; i++)
    {
      starlightParticle p(
			  ev.px[i], 
			  ev.py[i], 
			  ev.pz[i], 
			  ev.E[i], 
			  ev.mass[i], 
			  ev._fsParticle[i],
			  ev._charge[i]
			  );
      addParticle(p);
    }
}
/**
 * @brief Construct a new upc Event::upc Event object => A copy Constructor
 * @details Created to handle a warning with >= gcc9 compiler.
 * @param [event]: The upcEvent to be copied
 */
upcEvent::upcEvent(const upcEvent& event){
  _particles = event._particles;
  _vertices = event._vertices;
  _gammaEnergies = event._gammaEnergies;
}

upcEvent::~upcEvent()
{ }


upcEvent& upcEvent::operator=(const upcEvent& rhs)
{

  if(this != &rhs)
  {
    this->_particles = rhs._particles;
    this->_vertices = rhs._vertices;
    this->_gammaEnergies = rhs._gammaEnergies;
  }
  return *this;
}

upcEvent& upcEvent::operator+(const upcEvent& ev)
{
  for(unsigned int n = 0; n < ev._particles.size(); n++)
  {
    this->_particles.push_back(ev._particles.at(n));
  }
  for(unsigned int n = 0; n < ev._vertices.size(); n++)
  {
    this->_vertices.push_back(ev._vertices.at(n));
  }
 for(unsigned int n = 0; n < ev._gammaEnergies.size(); n++)
  {
    this->_gammaEnergies.push_back(ev._gammaEnergies.at(n));
  }
  return *this;
}

void upcEvent::boost(double rapidity)
{
    vector3 boostVector(0, 0, tanh(rapidity));
    std::vector<starlightParticle>::iterator part = _particles.begin();
      
    for (part = _particles.begin(); part != _particles.end(); part++)
    {
      (*part).Boost(boostVector);
    }
}
