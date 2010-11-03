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
	    void SetPxPyPzE(double px, double py, double pz, double e) { SetXYZT(px, py, pz, e); };
      
      double GetPx() const { return fSpaceVec.GetVector()[0]; }
      double GetPy() const { return fSpaceVec.GetVector()[1]; }
      double GetPz() const { return fSpaceVec.GetVector()[2]; }
      double GetE() const { return fTime; }

	    LorentzVector& operator +=(const LorentzVector& vec)
	    {
		    fSpaceVec += vec.fSpaceVec;
		    fTime     += vec.fTime;
		    return *this;
	    }
	    LorentzVector& operator -=(const LorentzVector& vec)
	    {
		    fSpaceVec -= vec.fSpaceVec;
		    fTime     -= vec.fTime;
		    return *this;
	    }

	    double M2() const { return fTime * fTime - fSpaceVec.Mag2(); }
      double M () const
	    {
	      const double mag2 = M2();
	      return (mag2 < 0) ? -sqrt(-mag2) : sqrt(mag2);
      }

	    Vector3 BoostVector() const
	    { return Vector3(fSpaceVec.X() / fTime, fSpaceVec.Y() / fTime, fSpaceVec.Z() / fTime); }
	    void Boost(const Vector3& beta)
	    {
		    const double beta2        = beta.Mag2();
		    const double gamma        = 1 / sqrt(1 - beta2);
		    const double betaTimesMom = beta.X() * fSpaceVec.X() + beta.Y() * fSpaceVec.Y() + beta.Z() * fSpaceVec.Z();
		    const double gamma2       = (beta2 > 0) ? (gamma - 1) / beta2 : 0;
		    SetXYZT(fSpaceVec.X() + gamma2 * betaTimesMom * beta.X() + gamma * beta.X() * fTime,
		            fSpaceVec.Y() + gamma2 * betaTimesMom * beta.Y() + gamma * beta.Y() * fTime,
		            fSpaceVec.Z() + gamma2 * betaTimesMom * beta.Z() + gamma * beta.Z() * fTime,
		            gamma * (fTime + betaTimesMom));
	    }
      
	    friend std::ostream& operator << (std::ostream&        out,
	                                      const LorentzVector& vec)
	    {
		    out << "(" << vec.GetPx() << ", " << vec.GetPy() << ", " << vec.GetPz()
		        << "; " << vec.GetE() << ")";
		    return out;
	    }

   private:
      
      Vector3 fSpaceVec;
      double fTime;
      
};

#endif // LORENTZVECTOR_H
