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


#include <iostream>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "nucleus.h"


using namespace std;


//______________________________________________________________________________
nucleus::nucleus(int Zin, int Ain, double bdeuteron,int in_or_co)
{
  _Zint=Zin;//Used within the class
  _Aint=Ain;//"
  _bford=bdeuteron;
  _NUCin_or_co=in_or_co;//dAu incoherent or coherent
  //cout<<"bdeuteron, _bford: "<<bdeuteron<<" "<<_bford<<endl;
  if (Zin==79){
    _Q0   =   0.060;
    _rho0 =   0.159407;
  }
  else if (Zin==53){
    _Q0   =   0.069;
    _rho0 =   0.161626;
  }
  else if (Zin==49){
    _Q0   =   0.069;
    _rho0 =   0.161626;
  }
  else if (Zin==29){
    _Q0   =   0.087;
    _rho0 =   0.166878;
  }
  else if (Zin==14){
    _Q0   =   0.115;
    _rho0 =   0.177128;
  }
  else if (Zin==8){
    _Q0   =   0.138;
    _rho0 =   0.188459;
  }
  else if (Zin==82){
    _Q0   =   0.059;
    _rho0 =   0.159176;
  }
  else if (Zin==20){
    _Q0   =	 0.102;
    _rho0 =	 0.171907;
  }
  else if (Zin==1){
    //_Q0 and rho) are not relevant for protons.
    _Q0   =   -1.0;
    _rho0 =   -1.0;
  }
  else{
    cout << "Warning:  Density not defined for this projectile!";
    _rho0  = 0.16;  //'typical' density
    _Q0    = 0.0;   // not used.
  }
  _r0=1.16*(1.-1.16*pow(_Aint,(-2./3.)));//For FRITIOF and FormFactor.
}


//______________________________________________________________________________
nucleus::~nucleus()
{ }


//______________________________________________________________________________
double nucleus::getRho0()
{
  return _rho0;
}


//______________________________________________________________________________
double nucleus::getQ0()
{
  return _Q0;
}


//______________________________________________________________________________
int nucleus::getZin()
{
  return _Zint;
}


//______________________________________________________________________________
int nucleus::getAin()
{
  return _Aint;
}


//______________________________________________________________________________
double nucleus::getWoodSaxonRadius()
{
  double WSradius = 0.;
  WSradius = 1.2*pow(_Aint,(1.0/3.0));
  return WSradius;
}


//______________________________________________________________________________
double nucleus::getWoodSaxonSkinDepth()
{
  double WSskindepth =0.;
  WSskindepth=.53;  //.53fermis, also known as diffuseness constant--Au
  return WSskindepth;
}


//______________________________________________________________________________
double nucleus::RNuc()
{
  double nuclearradius; //we use this for specific nuclei, au,pb,protons...
  if(_Zint==79)  
    nuclearradius = 6.38;
  else 
    if(_Zint==82)  
      nuclearradius = 6.62;
    else 
      if(_Zint==1&&_Aint==1) 
	nuclearradius = 0.7;
      else 
	if(_Zint==1&&_Aint==2) 
	  nuclearradius = 1.9;
	else 
	  nuclearradius = getWoodSaxonRadius();
  return nuclearradius;
}


//______________________________________________________________________________
double nucleus::fritiofR0()
{
  double fritiof;
  //_r0=1.16*(1.-1.16*pow(_Aint,(-2./3.)));
  fritiof=_r0*pow(_Aint,(1./3.));
  return fritiof;
}


//______________________________________________________________________________
double nucleus::rws(double r)
{
  double rws_r=1.0/( 1. + exp( (r-fritiofR0())/getWoodSaxonSkinDepth() ) );
  return rws_r;
}


//______________________________________________________________________________
double nucleus::formFactor(double t)
{
  double q,sph,yuk,a0,formf_r,rec;//R is just FRITIOF
  double st,st4; //sergey's stuff, also replaced b with _bford and dropped t02 since it wasnt used

  q   = sqrt(t);

  if(_Zint==1 && _Aint==1){
    //This is for protons
    //em formfactor of a proton
    rec = 1./(1.+(q*q)/0.71E0);
    formf_r = rec*rec;
  }
  else if(_Zint==1&&_Aint==2){   //careful with this line on dAu
    //This is for dAu//Sergey
    //incoherent f-f F(t)=0.34e(141.5t)+0.58e(26.1t)+0.08e(15.5t)
    st=0.34*exp(-141.5*t)+0.58*exp(-26.1*t)+0.08*exp(-15.5*t);
    st4=0.34*exp(-141.5*t/4)+0.58*exp(-26.1*t/4)+0.08*exp(-15.5*t/4);
    //st Paramters from Franco and Varma for st eqn PRL33 ...
    //formf_R from nuclear physics B 104 Eisenberg
    if(_NUCin_or_co==0){
      formf_r=exp(-_bford*t)*0.5*(1+st)-0.068*exp(-_bford*t*3/4)-st4*st4*exp(-_bford*t)+0.068*st4*exp(-_bford*t*3/4);
    }
    if(_NUCin_or_co==1){
      formf_r=(st4*st4*exp(-_bford*t)-0.068*st4*exp(-_bford*t*3/4));
    }
  }
  else{
    //     >> Use Parameterization from FRITIOF
    //     _r0=1.16*(1.-1.16*pow(A,(-2./3.)));
    //  R=_r0*pow(A,(1./3.));
    sph = sin(q*fritiofR0()/starlightConstants::hbarc) - (q*fritiofR0()/starlightConstants::hbarc)*cos(q*fritiofR0()/starlightConstants::hbarc);
    sph = ((3.*starlightConstants::hbarc*starlightConstants::hbarc*starlightConstants::hbarc)/(q*q*q*_r0*_r0*_r0))*sph;
    sph = sph/double(_Aint);
    
    a0  = 0.70;//fm
    yuk = 1./( 1. + ((a0*a0*q*q)/(starlightConstants::hbarc*starlightConstants::hbarc)) );
    formf_r = sph*yuk;
  }
  return formf_r;
}


//______________________________________________________________________________
double nucleus::thickness(double b_thick)
{
  //    JS      This code calculates the nuclear thickness function as per Eq. 4 in
  //    Klein and Nystrand, PRC 60.
  //    former DOUBLE PRECISION FUNCTION T(b)
                                                                                                                                                        
  double zsp,zmin,zmax,r,sum;
  int NGAUSS;
  double T;
  
  zmin=0.0;
  zmax=15.0;
  NGAUSS=5;
  
  //     >> DATA FOR GAUSS INTEGRATION
  double xg[6]={0.,0.1488743390,0.4333953941,0.6794095683,0.8650633667,0.9739065285};
  double ag[6]={0.,0.2955242247,0.2692667193,0.2190863625,0.1494513492,0.0666713443};
  
  sum=0.;
  
  for(int I=1;I<=NGAUSS;I++){
    zsp = 0.5*(zmax-zmin)*xg[I] + 0.5*(zmax+zmin);
    r   = sqrt(b_thick*b_thick+zsp*zsp);
    sum = sum + ag[I]*rws(r);
    zsp = 0.5*(zmax-zmin)*(-xg[I]) + 0.5*(zmax+zmin);
    r   = sqrt(b_thick*b_thick+zsp*zsp);
    sum = sum + ag[I]*rws(r);
    }
  
  sum = 0.5*(zmax-zmin)*2.*sum;
  T=sum;
  
  return T;
}
