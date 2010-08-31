#ifndef STARLIGHTPARTICLE_H
#define STARLIGHTPARTICLE_H

#include "lorentzvector.h"
class StarlightParticle : public LorentzVector
{
   public:
      
      StarlightParticle();
      StarlightParticle ( double px, double py, double pz, double e, double mass, int pdgCode, short charge);
      virtual ~StarlightParticle();
   
      void SetPdgCode(int pdgCode) { fPdgCode = pdgCode; }
      int GetPdgCode() const { return fPdgCode; }
      
      short SetCharge(short charge) { fCharge = charge; }
      short GetCharge() const { return fCharge; }
      
   private:
    
    int fPdgCode;
    short fCharge;
    double fMass;
};

#endif // STARLIGHTPARTICLE_H
