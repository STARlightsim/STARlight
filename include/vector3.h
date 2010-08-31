#ifndef VECTOR3_H
#define VECTOR3_H

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
      
   private:
      
      double fVec[3];
   
};

#endif // VECTOR3_H
