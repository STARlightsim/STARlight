
#include "upcXevent.h"


upcXEvent::upcXEvent() :
        _particles(0)
        ,_vertices(0)
        ,_beams(0)
        ,_neutrons(0)
        ,_vertext(0)
        ,_gamma(0)
        ,_sourceBeamIndex(0)
        ,_meson(0)        
        ,_gammaMasses(0)
        ,_gammaEnergies(0)
{ }

upcXEvent::upcXEvent(const upcXEvent& event){
  
  _particles = event._particles;
  _vertices = event._vertices;
  _beams = event._beams;
  _sourceBeamIndex = event._sourceBeamIndex;
  _neutrons = event._neutrons;
  _gamma = event._gamma;
  _meson = event._meson;

  _vertext = event._vertext;
  _gammaMasses = event._gammaMasses;
  _gammaEnergies = event._gammaEnergies;
}

upcXEvent::upcXEvent(starlightConstants::event &ev) :
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

upcXEvent::~upcXEvent()
{ }


upcXEvent& upcXEvent::operator=(const upcXEvent& rhs)
{

  if(this != &rhs)
  {
    this->_particles = rhs._particles;
    this->_vertices = rhs._vertices;
    
    this->_beams = rhs._beams;
    this->_sourceBeamIndex =rhs._sourceBeamIndex;
    this->_neutrons = rhs._neutrons;
    this->_gamma = rhs._gamma;
    this->_meson = rhs._meson;

    this->_vertext = rhs._vertext;
    this->_gammaMasses = rhs._gammaMasses;
    this->_gammaEnergies = rhs._gammaEnergies;
  }
  return *this;
}

/*I must say that this operator is meaningless and useless at the moment, if you want to use, think critically and modify how you want to concatenate weird object like: _beams, _sourcebeamIndex, ... */
upcXEvent& upcXEvent::operator+(const upcXEvent& ev)
{
  for(unsigned int n = 0; n < ev._particles.size(); n++)
  {
    this->_particles.push_back(ev._particles.at(n));
  }
  for(unsigned int n = 0; n < ev._vertices.size(); n++)
  {
    this->_vertices.push_back(ev._vertices.at(n));
  }
  for(unsigned int n = 0; n < ev._beams.size(); n++)
  {
    this->_beams.push_back(ev._beams.at(n));
  }
  for(unsigned int n = 0; n < ev._sourceBeamIndex.size(); n++)
  {
    this->_sourceBeamIndex.push_back(ev._sourceBeamIndex.at(n));
  }
  for(unsigned int n = 0; n < ev._neutrons.size(); n++)
  {
    this->_neutrons.push_back(ev._neutrons.at(n));
  }
  for(unsigned int n = 0; n < ev._vertext.size(); n++)
  {
    this->_vertext.push_back(ev._vertext.at(n));
  }
  for(unsigned int n = 0; n < ev._gamma.size(); n++)
  {
    this->_gamma.push_back(ev._gamma.at(n));
  }

  for(unsigned int n = 0; n < ev._gammaMasses.size(); n++)
  {
    this->_gammaMasses.push_back(ev._gammaMasses.at(n));
  }

  for(unsigned int n = 0; n < ev._gammaEnergies.size(); n++)
  {
    this->_gammaEnergies.push_back(ev._gammaEnergies.at(n));
  }
  return *this;
  for(unsigned int n=0; n< ev._meson.size(); n++)
  {
    this->_meson.push_back(ev._meson.at(n));
  }
}

void upcXEvent::boost(double rapidity)
{
    vector3 boostVector(0, 0, tanh(rapidity));
    std::vector<starlightParticle>::iterator part = _particles.begin();
      
    for (part = _particles.begin(); part != _particles.end(); part++)
    {
      (*part).Boost(boostVector);
    }
    std::vector<lorentzVector>::iterator beam = _beams.begin();
    for( beam = _beams.begin(); beam != _beams.end(); ++beam){
      (*beam).Boost(boostVector);
    }
    
    std::vector<lorentzVector>::iterator gamma = _gamma.begin();
    for( gamma = _gamma.begin(); gamma != _gamma.end(); ++gamma){
      (*gamma).Boost(boostVector);
    }
    
    std::vector<lorentzVector>::iterator vm = _meson.begin();
    for(vm = _meson.begin(); vm != _meson.end(); ++vm){
      (*vm).Boost(boostVector);
    }

    std::vector<lorentzVector>::iterator neutron = _neutrons.begin();
    for( neutron = _neutrons.begin(); neutron != _neutrons.end(); ++neutron){
      (*neutron).Boost(boostVector);
    }
}
