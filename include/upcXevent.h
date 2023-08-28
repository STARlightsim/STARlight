#ifndef UPCXEVENT_H
#define UPCXEVENT_H


#include <vector>
#include <cassert>

#include "starlightconstants.h"
#include "starlightparticle.h"

struct gammaPack{
   lorentzVector gamma;
   int sourceBeamIndex; // index of the source beam that emmitted the photon
   //int gammaIndex;//index of the photon's lorentzVector location in the _gamma array//No need to store this, as user cant access the array directly.
   double gammaMass;
   double gammaEnergy;
};

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
      void addGammaFromBeam1(lorentzVector gamma, float egamma, float Q2) { _gamma.push_back(gamma); _gammaEnergies.push_back(egamma); _gammaMasses.push_back(Q2); _sourceBeamIndex.push_back(0);}
      void addGammaFromBeam2(lorentzVector gamma, float egamma, float Q2) { _gamma.push_back(gamma); _gammaEnergies.push_back(egamma); _gammaMasses.push_back(Q2); _sourceBeamIndex.push_back(1);}
      
      
      /*STARlight Assumes that there are only two beams: Beam1 and Beam2. Beam1 travels in the +z direction.
      N.B: This function helps enforce the rule that Beam1 must be added to the array before Beam2*/
      void addOutgoingBeams(lorentzVector &beam1, lorentzVector &beam2){
         if(_beams.size()==0){
            _beams.push_back(beam1); 
            _beams.push_back(beam2);
         }
         else{
            //cout error "you cannot add beams to the array more than once. We also did not make provision to replace the beams already added. This is most likely a bug" 
            assert(false);
            }
      }
      void addNeutrons(lorentzVector &neutrons){_neutrons.push_back(neutrons);}
      /*Adds the Mandelstam variable "t" - transferred momentum squared*/
      void addVertext(double t){_vertext.push_back(t);}
      /*Use only for Photonuclear Interaction, i.e. Exclusive Vector Meson
      Only allows a single Meson*/
      void addVectorMeson(lorentzVector &VM){if(_meson.size() == 0) _meson.push_back(VM); else assert(false);}
      /*Use only for Two-Photon Interactions, allows multiple Mesons*/
      void addMeson(lorentzVector &meson){_meson.push_back(meson);}

      const lorentzVector getVectorMeson()const{if(_meson.size() == 1) return _meson[0]; else assert(false); return lorentzVector();}
      const std::vector<lorentzVector> * getMesons() const{return &_meson;}
      const std::vector<starlightParticle> * getParticles() const { return &_particles; }
      const std::vector<vector3> * getVertices() const { return &_vertices; }
      
      /*returns the list of Neutrons produced during nuclear breakup - NOT IMPLEMENTED*/
      const std::vector<lorentzVector>* getNeutrons() const{return &_neutrons;}
      const std::vector<double> * getVertext() const {return &_vertext;}

      /*Returns the Photon Q^2 for each photon. Index associated with those of @see getGammaEnergies*/
      const std::vector<float>* getGammaMasses() const {return &_gammaMasses;}
      /*Returns the Photon Target Energies. Index is associated with those of @see getGammaMasses*/
      const std::vector<float> * getGammaEnergies() const { return &_gammaEnergies; }
      const lorentzVector getBeam1() const{
         return _beams[0];
      }
      const lorentzVector getBeam2() const{
         return _beams[1];
      }

      /*should only be used in photonuclear interactions*/
      const gammaPack getGamma() const{
         gammaPack outt;
         if(_gamma.size() ==1){
            outt.gamma = _gamma[0];
            outt.gammaEnergy = _gammaEnergies[0];
            outt.gammaMass = _gammaMasses[0];
            outt.sourceBeamIndex = _sourceBeamIndex[0];
            return outt;
         }      
         else
            assert(false);
         return outt;
      }
      /*Should only be used in two photon interactions*/
      const gammaPack getGammaFromBeam1() const{
         gammaPack outt;
         if(_gamma.size() ==2){
            if(_sourceBeamIndex[0] == 0 && _sourceBeamIndex[1] == 1)
            {
               outt.sourceBeamIndex = 0;
               outt.gamma = _gamma[0];
               outt.gammaEnergy = _gammaEnergies[0];
               outt.gammaMass = _gammaMasses[0];
               return outt;
            }
            else if(_sourceBeamIndex[0] == 1 && _sourceBeamIndex[1] == 0)
            {
               outt.sourceBeamIndex = 0;
               outt.gamma = _gamma[1];
               outt.gammaEnergy = _gammaEnergies[1];
               outt.gammaMass = _gammaMasses[1];
               return outt;
            }
            else assert(false);
         }
         else assert(false);
         return outt;
      }
      /*Should only be used in two photon interactions*/
      const gammaPack getGammaFromBeam2() const{
         gammaPack outt;
         if(_gamma.size() ==2){
            if(_sourceBeamIndex[0] == 0 && _sourceBeamIndex[1] == 1)
            {
               outt.sourceBeamIndex = 1;
               outt.gamma = _gamma[1];
               outt.gammaEnergy = _gammaEnergies[1];
               outt.gammaMass = _gammaMasses[1];
               return outt;
            }
            else if(_sourceBeamIndex[0] == 1 && _sourceBeamIndex[1] == 0)
            {
               outt.sourceBeamIndex = 1;
               outt.gamma = _gamma[0];
               outt.gammaEnergy = _gammaEnergies[0];
               outt.gammaMass = _gammaMasses[0];
               return outt;
            }
            else assert(false);
         }
         else assert(false);
         return outt;
      }

      int targetBeamNo() const{
         if(_sourceBeamIndex.size() ==1){
            return 2 - _sourceBeamIndex[0];//if sourceIndex = 0, sourceNo =1, and targetNo = 2; else if sourceIndex =1, source no = 2, and targetNo = 1; 
         }
         else if(_sourceBeamIndex.size() ==2) return -1;
         else assert(false);
         return -1;
      }
      int gammaCount() const{
         return _gamma.size();
      }


      upcXEvent & operator=(const upcXEvent&);
      upcXEvent & operator+(const upcXEvent&);
      
      void boost(double rapidity);
   private:
      
      std::vector<starlightParticle> _particles;
      std::vector<vector3> _vertices;
      std::vector<lorentzVector> _beams;//contains both beam 1 and 2 lumped together
      //std::vector<int> beamA;//specify A and Z for outgoing beam1 and beam2 - important when neutrons are emmited
      //std::vector<int> beamZ;
      std::vector<lorentzVector> _neutrons;
      std::vector<double> _vertext; //we cannot associate this with target beams because in TwoPhoton case, there are no target beams but t exist.

      //std:vector<int>_breakupGammaIndices;
      std::vector<lorentzVector> _gamma;
      std::vector<int> _sourceBeamIndex;
      std::vector<lorentzVector> _meson; // Meson may/may not be vector meson.
      std::vector<float> _gammaMasses; //Q^2 for the photon
      std::vector<float> _gammaEnergies; //Energy in the target frame of reference for the photon.
};


#endif  // UPCXEVENT_H
