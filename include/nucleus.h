///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef NUCLEUS_H
#define NUCLEUS_H


//This class holds the information for a target nucleus
class nucleus
{

 public:
  nucleus(int Zin, int Ain, double bdeuteron, int in_or_co);
  ~nucleus();
 
  int getZin();
  int getAin();
  double getWoodSaxonRadius();
  double getWoodSaxonSkinDepth();
  double RNuc(); //was used as wood-saxon anyway?
  double fritiofR0(); //Fritiof _r0 (rws)/formfactor
  double rws(double r); //Woodsaxon nuclear density
  double formFactor(double t);  //Formfactor function
  double getQ0();
  double getRho0();
  double thickness(double b);
 
 private:
  double _r0;
  double _rho0;
  double _Q0;
  int _Zint;
  int _Aint;
  double _bford;
  int _NUCin_or_co;
};


#endif  // NUCLEUS_H
