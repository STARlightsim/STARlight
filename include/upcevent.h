
#ifndef UPCEVENT_H
#define UPCEVENT_H

#include <vector>
#include "starlightconstants.h"
#include "starlightparticle.h"

class UPCEvent 
{
   public:

      UPCEvent();
      ~UPCEvent();

//      UPCEvent & operator  = (const StarlightConstants::event &);

      void AddParticle(StarlightParticle &part) { fParticles.push_back(part); }
      void AddVertex(Vector3 &vertex) { fVertices.push_back(vertex); }
      
      const std::vector<StarlightParticle> * GetParticles() const { return &fParticles; }
      const std::vector<Vector3> * GetVertices() const { return &fVertices; }
      
   private:
      
      int fNTracks;
      std::vector<StarlightParticle> fParticles;
      std::vector<Vector3> fVertices;
};

#endif // UPCEVENT_H
