#ifndef NUCLEUS_H
#define NUCLEUS_H

//This class holds the information for a target nucleus

class Nucleus
{

 public:
  Nucleus(int Zin, int Ain, double bdeuteron, int in_or_co);
  ~Nucleus();
 
  int getZin();
  int getAin();
  double getWoodSaxonradius();
  double getWoodSaxonskindepth();
  double RNuc(); //was used as wood-saxon anyway?
  double fritiofr0(); //Fritiof r0 (rws)/formfactor
  double rws(double r); //Woodsaxon nuclear density
  double formfactor(double t);  //Formfactor function
  double getQ0();
  double getrho0();
  double Thickness(double b);
 
 private:
  double r0;
  double rho0;
  double Q0;
  int Zint;
  int Aint;
  double bford;
  int NUCin_or_co;
};
#endif //NUCLEUS_H
