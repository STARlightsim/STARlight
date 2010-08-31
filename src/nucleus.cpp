// nucleus.cpp
/*
 * $Id: nucleus.cpp,v 1.0 2010/07/04  $
 *
 * /author Joseph Butterwoth
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * $Log: $
 *
 */
#include <iostream>
#include <fstream>

using namespace std;

#include <math.h>
#include "starlightconstants.h"
#include "nucleus.h"
//______________________________________________________________________________
Nucleus::Nucleus(int Zin, int Ain, double bdeuteron,int in_or_co)
{
  Zint=Zin;//Used within the class
  Aint=Ain;//"
  bford=bdeuteron;
  NUCin_or_co=in_or_co;//dAu incoherent or coherent
  //cout<<"bdeuteron, bford: "<<bdeuteron<<" "<<bford<<endl;
  if (Zin==79){
    Q0   =   0.060;
    rho0 =   0.159407;
  }
  else if (Zin==53){
    Q0   =   0.069;
    rho0 =   0.161626;
  }
  else if (Zin==49){
    Q0   =   0.069;
    rho0 =   0.161626;
  }
  else if (Zin==29){
    Q0   =   0.087;
    rho0 =   0.166878;
  }
  else if (Zin==14){
    Q0   =   0.115;
    rho0 =   0.177128;
  }
  else if (Zin==8){
    Q0   =   0.138;
    rho0 =   0.188459;
  }
  else if (Zin==82){
    Q0   =   0.059;
    rho0 =   0.159176;
  }
  else if (Zin==20){
    Q0   =	 0.102;
    rho0 =	 0.171907;
  }
  else if (Zin==1){
    //Q0 and rho) are not relevant for protons.
    Q0   =   -1.0;
    rho0 =   -1.0;
  }
  else{
    cout << "Warning:  Density not defined for this projectile!";
    rho0  = 0.16;  //'typical' density
    Q0    = 0.0;   // not used.
  }
  r0=1.16*(1.-1.16*pow(Aint,(-2./3.)));//For FRITIOF and FormFactor.
}
//______________________________________________________________________________
double Nucleus::getrho0()
{
  return rho0;
}
//______________________________________________________________________________
double Nucleus::getQ0()
{
  return Q0;
}
//______________________________________________________________________________
int Nucleus::getZin()
{
  return Zint;
}
//______________________________________________________________________________
int Nucleus::getAin()
{
  return Aint;
}
//______________________________________________________________________________
double Nucleus::getWoodSaxonradius()
{
  double WSradius = 0.;
  WSradius = 1.2*pow(Aint,(1.0/3.0));
  return WSradius;
}
//______________________________________________________________________________
double Nucleus::getWoodSaxonskindepth()
{
  double WSskindepth =0.;
  WSskindepth=.53;  //.53fermis, also known as diffuseness constant--Au
  return WSskindepth;
}
//______________________________________________________________________________
double Nucleus::RNuc()
{
  double nuclearradius; //we use this for specific nuclei, au,pb,protons...
  if(Zint==79)  
    nuclearradius = 6.38;
  else 
    if(Zint==82)  
      nuclearradius = 6.62;
    else 
      if(Zint==1&&Aint==1) 
	nuclearradius = 0.7;
      else 
	if(Zint==1&&Aint==2) 
	  nuclearradius = 1.9;
	else 
	  nuclearradius = getWoodSaxonradius();
  return nuclearradius;
}
//______________________________________________________________________________
double Nucleus::fritiofr0()
{
  double fritiof;
  //r0=1.16*(1.-1.16*pow(Aint,(-2./3.)));
  fritiof=r0*pow(Aint,(1./3.));
  return fritiof;

}
//______________________________________________________________________________
double Nucleus::rws(double r)
{
  double rws_r=1.0/( 1. + exp( (r-fritiofr0())/getWoodSaxonskindepth() ) );
  return rws_r;
}
//______________________________________________________________________________
double Nucleus::formfactor(double t)
{
  
  double q,sph,yuk,a0,formf_r,rec;//R is just FRITIOF
  double st,st4; //sergey's stuff, also replaced b with bford and dropped t02 since it wasnt used

  q   = sqrt(t);

  if(Zint==1 && Aint==1){
    //This is for protons
    //em formfactor of a proton
    rec = 1./(1.+(q*q)/0.71E0);
    formf_r = rec*rec;
  }
  else if(Zint==1&&Aint==2){   //careful with this line on dAu
    //This is for dAu//Sergey
    //incoherent f-f F(t)=0.34e(141.5t)+0.58e(26.1t)+0.08e(15.5t)
    st=0.34*exp(-141.5*t)+0.58*exp(-26.1*t)+0.08*exp(-15.5*t);
    st4=0.34*exp(-141.5*t/4)+0.58*exp(-26.1*t/4)+0.08*exp(-15.5*t/4);
    //st Paramters from Franco and Varma for st eqn PRL33 ...
    //formf_R from nuclear physics B 104 Eisenberg
    if(NUCin_or_co==0){
      formf_r=exp(-bford*t)*0.5*(1+st)-0.068*exp(-bford*t*3/4)-st4*st4*exp(-bford*t)+0.068*st4*exp(-bford*t*3/4);
    }
    if(NUCin_or_co==1){
      formf_r=(st4*st4*exp(-bford*t)-0.068*st4*exp(-bford*t*3/4));
    }
  }
  else{
    //     >> Use Parameterization from FRITIOF
    //     r0=1.16*(1.-1.16*pow(A,(-2./3.)));
    //  R=r0*pow(A,(1./3.));
    sph = sin(q*fritiofr0()/StarlightConstants::hbarc) - (q*fritiofr0()/StarlightConstants::hbarc)*cos(q*fritiofr0()/StarlightConstants::hbarc);
    sph = ((3.*StarlightConstants::hbarc*StarlightConstants::hbarc*StarlightConstants::hbarc)/(q*q*q*r0*r0*r0))*sph;
    sph = sph/double(Aint);
    
    a0  = 0.70;//fm
    yuk = 1./( 1. + ((a0*a0*q*q)/(StarlightConstants::hbarc*StarlightConstants::hbarc)) );
    formf_r = sph*yuk;
  }
  return formf_r;
}
//______________________________________________________________________________
double Nucleus::Thickness(double b_thick)
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
//______________________________________________________________________________
Nucleus::~Nucleus()
{
        // insert your code here
}

