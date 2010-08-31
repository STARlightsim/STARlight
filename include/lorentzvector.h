#ifndef LORENTZVECTOR_H
#define LORENTZVECTOR_H

#include "vector3.h"
#include <vector>

class LorentzVector
{
   public:
      
      LorentzVector();
      virtual ~LorentzVector();
      
      LorentzVector(double x, double y, double z, double t);
      //LorentzVector(double px, double py, double pz, double e);
      
      void SetXYZT(double x, double y, double z, double t);
      void SetPxPyPzE(double px, double py, double pz, double e);
      
      double GetPx() const { return fSpaceVec.GetVector()[0]; }
      double GetPy() const { return fSpaceVec.GetVector()[1]; }
      double GetPz() const { return fSpaceVec.GetVector()[2]; }
      double GetE() const { return fTime; }
      
   private:
      
      Vector3 fSpaceVec;
      double fTime;
      
};

#endif // LORENTZVECTOR_H
