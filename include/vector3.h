#ifndef VECTOR3_H
#define VECTOR3_H


#include <iostream>
#include <cmath>


class Vector3
{
   public:
      Vector3();
      Vector3(double *vec);
      Vector3(double x, double y, double z);
      virtual ~Vector3();
      
      const double* GetVector() const { return fVec; }
      
      void SetVector(double x, double y, double z);
      void SetVector(double *vec);

	    Vector3& operator =(const Vector3& vec)
	    {
		    if (this != &vec)
			    for (unsigned int i = 0; i < 3; ++i)
				    fVec[i] = vec.fVec[i];
		    return *this;
	    }

	    Vector3& operator +=(const Vector3& vec)
	    {
		    for (unsigned int i = 0; i < 3; ++i)
			    fVec[i] += vec.fVec[i];
		    return *this;
	    }
	    Vector3& operator -=(const Vector3& vec)
	    {
		    for (unsigned int i = 0; i < 3; ++i)
			    fVec[i] -= vec.fVec[i];
		    return *this;
	    }

	    double X() const { return fVec[0]; }
	    double Y() const { return fVec[1]; }
	    double Z() const { return fVec[2]; }

	    double Mag2() const { return fVec[0] * fVec[0] + fVec[1] * fVec[1] + fVec[2] * fVec[2]; }
	    double Mag () const { return sqrt(Mag2()); }
      
	    friend std::ostream& operator << (std::ostream&  out,
	                                      const Vector3& vec)
	    {
		    out << "(" << vec.X() << ", " << vec.Y() << ", " << vec.Z() << ")";
		    return out;
	    }
	
   private:
      
      double fVec[3];
   
};

#endif // VECTOR3_H
