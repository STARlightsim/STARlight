#ifndef UPCXEVENT_H
#define UPCXEVENT_H


#include <vector>

#include "starlightconstants.h"
#include "starlightparticle.h"


class upcXEvent 
{
   public:

      upcXEvent();
      //upcXEvent(starlightConstants::event &ev);
      upcXEvent(const upcXEvent& event);
      upcXEvent(starlightConstants::event &ev);
      ~upcXEvent();

      void addParticle(starlightParticle &part) { _particles.push_back(part); }
      void addVertex(vector3 &vertex) { _vertices.push_back(vertex); }
      void addGamma(lorentzVector gamma, float egamma, float Q2) { _gamma.push_back(gamma); _gammaEnergies.push_back(egamma); _gammaMasses.push_back(Q2);}
      void addOutgoingBeam1(lorentzVector &beam1, bool istarget){_beams.push_back(beam1); _beamNo.push_back(1); _beamIsTarget.push_back(istarget);}
      void addOutgoingBeam1(lorentzVector &beam1,  int _targetBeamNo){_beams.push_back(beam1); _beamNo.push_back(1); _beamIsTarget.push_back( _targetBeamNo == 1 ? true : false );}
      void addOutgoingBeam2(lorentzVector &beam2,  bool istarget){_beams.push_back(beam2); _beamNo.push_back(2); _beamIsTarget.push_back(istarget);}
      void addOutgoingBeam2(lorentzVector &beam2,  int _targetBeamNo){_beams.push_back(beam2); _beamNo.push_back(2); _beamIsTarget.push_back( _targetBeamNo == 2 ? true : false );}
      void addNeutrons(lorentzVector &neutrons){_neutrons.push_back(neutrons);}
      void addVertext(double t){_vertext.push_back(t);}

      const std::vector<starlightParticle> * getParticles() const { return &_particles; }
      const std::vector<vector3> * getVertices() const { return &_vertices; }
      const std::vector<lorentzVector>* getBeams() const {return &_beams;}
      const std::vector<int>* getBeamNo() const {return &_beamNo;}
      const std::vector<bool>* getBeamIsTarget() const {return &_beamIsTarget;}
      //const std::vector<lorentzVector>* getBeam1() const{return &_beam1;}//later
      //const std::vector<lorentzVector>* getBeam2() const{return &_beam2;}//later
      const std::vector<lorentzVector>* getNeutrons() const{return &_neutrons;}
      const std::vector<lorentzVector>* getGamma() const {return &_gamma;}
      const std::vector<double> * getVertext() const {return &_vertext;}
      //const std::vector<int> * getTargetBeam() const {return ;} //later a subset of beams containing only target beams to match perfectly with vertext
      //const std::vector<int> * getSourceBeam() const {return ;} //later
      const std::vector<float>* getGammaMasses() const {return &_gammaMasses;}
      const std::vector<float> * getGammaEnergies() const { return &_gammaEnergies; }

      upcXEvent & operator=(const upcXEvent&);
      upcXEvent & operator+(const upcXEvent&);
      
      void boost(double rapidity);
   private:
      
      std::vector<starlightParticle> _particles;
      std::vector<vector3> _vertices;
      std::vector<lorentzVector> _beams;//contains both beam 1 and 2 lumped together
      //std::vector<int> beamA;//specify A and Z for outgoing beam1 and beam2 - important when neutrons are emmited
      //std::vector<int> beamZ;
      std::vector<int> _beamNo;//keep tracks of the beam number of each beam whether beam 1 or 2
      std::vector<bool> _beamIsTarget;//maitains whether the beam is a photon target or source. needed to associate _vertext
      std::vector<lorentzVector> _neutrons;
      std::vector<double> _vertext;//in estarlight this is associated with the target beams only... it will be treated same here

      //std:vector<int>_breakupGammaIndices;
      std::vector<lorentzVector> _gamma;
      std::vector<float> _gammaMasses;
      std::vector<float> _gammaEnergies;
};


#endif  // UPCXEVENT_H
