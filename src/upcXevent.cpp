
#include "upcXevent.h"


upcXEvent::upcXEvent() :
        _particles(0)
        ,_vertices(0)
{ }

upcXEvent::upcXEvent(const upcXEvent& event){
  
  _particles = event._particles;
  _vertices = event._vertices;
  _beams = event._beams;
  _beamNo = event._beamNo;
  _beamIsTarget = event._beamIsTarget; 
  _neutrons = event._neutrons;
  _gamma = event._gamma;

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
    this->_beamIsTarget = rhs._beamIsTarget;
    this->_beamNo = rhs._beamNo;
    this->_neutrons = rhs._neutrons;
    this->_gamma = rhs._gamma;

    this->_vertext = rhs._vertext;
    this->_gammaMasses = rhs._gammaMasses;
    this->_gammaEnergies = rhs._gammaEnergies;
  }
  return *this;
}

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
  for(unsigned int n = 0; n < ev._beamNo.size(); n++)
  {
    this->_beamNo.push_back(ev._beamNo.at(n));
  }
  for(unsigned int n = 0; n < ev._beamIsTarget.size(); n++)
  {
    this->_beamIsTarget.push_back(ev._beamIsTarget.at(n));
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
    std::vector<lorentzVector>::iterator neutron = _neutrons.begin();
    for( neutron = _neutrons.begin(); neutron != _neutrons.end(); ++neutron){
      (*neutron).Boost(boostVector);
    }
}
