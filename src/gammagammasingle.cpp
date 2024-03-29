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
// $Rev:: 293                         $: revision of last commit
// $Author:: butter                   $: author of last commit
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
#include <vector>


#include "starlightconstants.h"
#include "gammagammasingle.h"
#include "starlightconfig.h"

using namespace std;


//______________________________________________________________________________
Gammagammasingle::Gammagammasingle(const inputParameters& inputParametersInstance, randomGenerator* randy, beamBeamSystem& bbsystem)
: eventChannel(inputParametersInstance, randy, bbsystem)
#ifdef ENABLE_PYTHIA
,_pyDecayer()
#endif
{

#ifdef ENABLE_PYTHIA
    _pyDecayer.init();
#endif
  
  //Storing inputparameters into protected members for use
  _GGsingInputnumw=inputParametersInstance.nmbWBins();
  _GGsingInputnumy=inputParametersInstance.nmbRapidityBins();
  _GGsingInputpidtest=inputParametersInstance.prodParticleType();
  _GGsingInputGamma_em=inputParametersInstance.beamLorentzGamma();
  _axionMass=inputParametersInstance.axionMass(); // AXION HACK
  cout<<"SINGLE MESON pid test: "<<_GGsingInputpidtest<<endl;
  //reading in luminosity tables
  read();
  //Now calculating crosssection
  singleCrossSection();
}


//______________________________________________________________________________
Gammagammasingle::~Gammagammasingle()
{ }


//______________________________________________________________________________
void Gammagammasingle::singleCrossSection()
{
  //This function carries out a delta function cross-section calculation. For reference, see STAR Note 243, Eq. 8
  //Multiply all _Farray[] by _f_max
  double _sigmaSum=0.,remainw=0.;//_remainwd=0.;
  int ivalw =0;//_ivalwd;
  //calculate the differential cross section and place in the sigma table
  cout << "MASS  " << getMass() << "\n"; // AXION HACK, optional
  cout << "WIDTH  " << getWidth() << "\n";// AXION HACK, optional
  _wdelta=getMass();
  for(int i=0;i<_GGsingInputnumw;i++){
    for(int j=0;j<_GGsingInputnumy;j++){
      // Eq. 1 of starnote 347
      _sigmax[i][j]=(getSpin()*2.+1.)*4*starlightConstants::pi*starlightConstants::pi*getWidth()/
	(getMass()*getMass()*getMass())*_f_max*_Farray[i][j]*starlightConstants::hbarc*starlightConstants::hbarc/100.;
    }
  }
  //find the index, i,for the value of w just less than the mass because we want to use the value from the sigma table that has w=mass

  for(int i=0;i<_GGsingInputnumw;i++){
    if(getMass()>_Warray[i]) ivalw=i;
  }

  remainw = (getMass()-_Warray[ivalw])/(_Warray[ivalw+1]-_Warray[ivalw]);
  _ivalwd = ivalw;
  _remainwd = remainw;
  //if we are interested rho pairs at threshold, the just set sigma to 100nb
  switch(_GGsingInputpidtest){
  case starlightConstants::ZOVERZ03:
    _sigmaSum =0.;
    for(int j=0;j<_GGsingInputnumy-1;j++){
                        _sigmaSum = _sigmaSum +(_Yarray[j+1]-_Yarray[j])*
			  100.0E-9*(.1/getMass())*((1.-remainw)*_f_max*
						   (_Farray[ivalw][j]+_Farray[ivalw][j])/2.+remainw*_f_max*
						   (_Farray[ivalw+1][j]+_Farray[ivalw+1][j+1])/2.);
    }
    break;
  default:
    //Sum to find the total cross-section
    _sigmaSum =0.;
    for(int j =0;j<_GGsingInputnumy-1;j++){
                        _sigmaSum = _sigmaSum+
			  (_Yarray[j+1]-_Yarray[j])*((1.-remainw)*
						   (_sigmax[ivalw][j]+_sigmax[ivalw][j+1])/2.+remainw*
						   (_sigmax[ivalw+1][j]+_sigmax[ivalw+1][j+1])/2.);
    }
  }
  // if(_sigmaSum > 0.1) cout <<"The total cross-section is: "<<_sigmaSum<<" barns."<<endl;
  // else if(_sigmaSum > 0.0001)cout <<"The total cross-section is: "<<_sigmaSum*1000<<" mb."<<endl;
  // else cout <<"The total cross-section is: "<<_sigmaSum*1000000<<" ub."<<endl;
  cout<<endl;
  if (_sigmaSum > 1.){
     cout << "Total cross section: "<<_sigmaSum<<" barn."<<endl;  
  } else if (1000.*_sigmaSum > 1.){
     cout << "Total cross section: "<<1000.*_sigmaSum<<" mb."<<endl;  
  } else if (1000000.*_sigmaSum > 1.){
    cout << "Total cross section: "<<1000000.*_sigmaSum<<" microbarn."<<endl;  
  } else if (1.E9*_sigmaSum > 1.){
    cout << "Total cross section: "<<1.E9*_sigmaSum<<" nanobarn."<<endl;  
  } else if (1.E12*_sigmaSum > 1.){
    cout << "Total cross section: "<<1.E12*_sigmaSum<<" picobarn."<<endl;  
  } else {
    cout << "Total cross section: "<<1.E15*_sigmaSum<<" femtobarn."<<endl;  
  }
  cout<<endl; 
  setTotalChannelCrossSection(_sigmaSum);
     
  return;
}


//______________________________________________________________________________
void Gammagammasingle::pickw(double &w)
{
  //This function picks a w for the 2-photon calculation. 
  double sgf=0.,signorm=0.,x=0.,remainarea=0.,remainw=0.,a=0.,b=0.,c=0.;
  int ivalw=0;
  
  double * _sigofw;
  double * sgfint;
  _sigofw = new double[starlightLimits::MAXWBINS];
  sgfint = new double[starlightLimits::MAXYBINS];
 
  if(_wdelta != 0){
    w=_wdelta;
    ivalw=_ivalwd;
    remainw=_remainwd;
  }
  else{
    //deal with the case where sigma is an array
    //_sigofw is simga integrated over y using a linear interpolation
    //sigint is the integral of sgfint, normalized
    
    //integrate sigma down to a function of just w
    for(int i=0;i<_GGsingInputnumw;i++){
      _sigofw[i]=0.;
      for(int j=0;j<_GGsingInputnumy-1;j++){
	_sigofw[i] = _sigofw[i]+(_Yarray[j+1]-_Yarray[j])*(_sigmax[i][j+1]+_sigmax[i][j])/2.;
      }
    }
    //calculate the unnormalized sgfint array
    sgfint[0]=0.;
    for(int i=0;i<_GGsingInputnumw-1;i++){
      sgf=(_sigofw[i+1]+_sigofw[i])*(_Warray[i+1]-_Warray[i])/2.;
      sgfint[i+1]=sgfint[i]+sgf;
    }
    //normalize sgfint array
    signorm=sgfint[_GGsingInputnumw-1];
    
    for(int i=0;i<_GGsingInputnumw;i++){
      sgfint[i]=sgfint[i]/signorm;
    }
    //pick a random number
    x = _randy->Rndom();
    //compare x and sgfint to find the ivalue which is just less than the random number x
    for(int i=0;i<_GGsingInputnumw;i++){
      if(x > sgfint[i]) ivalw=i;
    }
    //remainder above ivalw
    remainarea = x - sgfint[ivalw];
    
    //figure out what point corresponds to the excess area in remainarea
    c = -remainarea*signorm/(_Warray[ivalw+1]-_Warray[ivalw]);
    b = _sigofw[ivalw];
    a = (_sigofw[ivalw+1]-_sigofw[ivalw])/2.;
    if(a==0.){
      remainw = -c/b;
    }
    else{
      remainw = (-b+sqrt(b*b-4.*a*c))/(2.*a);
    }
    _ivalwd = ivalw;
    _remainwd = remainw;
    //calculate the w value
    w = _Warray[ivalw]+(_Warray[ivalw+1]-_Warray[ivalw])*remainw;
    }

  delete[] _sigofw;
  delete[] sgfint;
}


//______________________________________________________________________________
void Gammagammasingle::picky(double &y)
{
  double * sigofy;
  double * sgfint;
  sigofy = new double[starlightLimits::MAXYBINS];
  sgfint = new double[starlightLimits::MAXYBINS];
  
  double remainw =0.,remainarea=0.,remainy=0.,a=0.,b=0.,c=0.,sgf=0.,signorm=0.,x=0.;
  int ivalw=0,ivaly=0;
  
  ivalw=_ivalwd;
  remainw=_remainwd;
  //average over two colums to get y array
  for(int j=0;j<_GGsingInputnumy;j++){
    sigofy[j]=_sigmax[ivalw][j]+(_sigmax[ivalw+1][j]-_sigmax[ivalw][j])*remainw;
  }
  //calculate the unnormalized sgfint
  
  sgfint[0]=0.;
  for(int j=0;j<_GGsingInputnumy-1;j++){
    sgf = (sigofy[j+1]+sigofy[j])/2.;
    sgfint[j+1]=sgfint[j]+sgf*(_Yarray[j+1]-_Yarray[j]);
  }
  
  //normalize the sgfint array
  signorm = sgfint[_GGsingInputnumy-1];
  
  for(int j=0;j<_GGsingInputnumy;j++){
    sgfint[j]=sgfint[j]/signorm;
  }
  //pick a random number
  x = _randy->Rndom();
  //compare x and sgfint to find the ivalue which is just less then the random number x
  for(int i=0;i<_GGsingInputnumy;i++){
    if(x > sgfint[i]) 
      ivaly = i;
  }
  //remainder above ivaly
  remainarea = x - sgfint[ivaly];
  //figure what point corresponds to the leftover area in remainarea
  c = -remainarea*signorm/(_Yarray[ivaly+1]-_Yarray[ivaly]);
  b = sigofy[ivaly];
  a = (sigofy[ivaly+1]-sigofy[ivaly])/2.;
  if(a==0.){
    remainy = -c/b;
  }
  else{
    remainy = (-b + sqrt(b*b-4.*a*c))/(2.*a);
  }
  //calculate the y value
  y = _Yarray[ivaly]+(_Yarray[ivaly+1]-_Yarray[ivaly])*remainy;
  delete[] sigofy;
  delete[] sgfint;
}


//______________________________________________________________________________
void Gammagammasingle::parentMomentum(double w,double y,double &E,double &px,double &py,double&pz,
                    double &Eb1, double &pxb1, double &pyb1, double &pzb1,
								    double &Eb2, double &pxb2, double &pyb2, double &pzb2, double &t2,
								    double &Egam1, double&pxgam1, double &pygam1, double &pzgam1, double &Q2gam1,
                    double &Egam2, double&pxgam2, double &pygam2, double &pzgam2, double &Q2gam2)
{
  //this function calculates px,py,pz,and E given w and y
  double anglepp1=0.,anglepp2=0.,ppp1=0.,ppp2=0.,E1=0.,E2=0.,signpx=0.,pt=0.;
  
  //E1 and E2 are for the 2 photons in the CM frame
  E1 = w*exp(y)/2.;//from beam1
  E2 = w*exp(-y)/2.;//from beam2
  //calculate px and py
  //to get x and y components-- phi is random between 0 and 2*pi
  anglepp1 = _randy->Rndom();
  anglepp2 = _randy->Rndom();
  
  ppp1 = pp1(E1);
  ppp2 = pp2(E2);
  px = ppp1*cos(2.*starlightConstants::pi*anglepp1)+ppp2*cos(2.*starlightConstants::pi*anglepp2);
  py = ppp1*sin(2.*starlightConstants::pi*anglepp1)+ppp2*sin(2.*starlightConstants::pi*anglepp2);
  //Compute vector sum Pt=Pt1+Pt2 to find pt for the produced particle
  pt = sqrt(px*px+py*py);
  //W is the mass of the produced particle (not necessarily on-mass-shell).Now compute its energy and pz
  E = sqrt(w*w+pt*pt)*cosh(y);
  pz= sqrt(w*w+pt*pt)*sinh(y);
  signpx = _randy->Rndom();
  //pick the z direction
  if(signpx > 0.5) 
    pz = -pz;

  double E0b1 = _ip->protonEnergy()*_ip->beam1A();
	double px0b1 = 0, px0b2 =0, py0b1=0, py0b2=0;
	double pz0b1 = sqrt(E0b1*E0b1 - _ip->protonMass()*_ip->protonMass()*_ip->beam1A()*_ip->beam1A());
	double E0b2 = _ip->protonEnergy()*_ip->beam2A();
	double pz0b2 = -sqrt(E0b2*E0b2 - _ip->protonMass()*_ip->protonMass()*_ip->beam2A()*_ip->beam2A());
  pxgam1 = ppp1*cos(2.*starlightConstants::pi*anglepp1);
  pygam1 = ppp1*sin(2.*starlightConstants::pi*anglepp1);
  pxgam2 = ppp2*cos(2.*starlightConstants::pi*anglepp2);
  pygam2 = ppp2*sin(2.*starlightConstants::pi*anglepp2);
  Egam1 = E1;
  Egam2 = E2;

  Eb1 = E0b1 - Egam1;
  pxb1 = px0b1 - pxgam1;
  pyb1 = py0b1 - pygam1;
  pzb1 = sqrt(Eb1*Eb1 - (pxb1*pxb1 + pyb1*pyb1 + _ip->protonMass()*_ip->beam1A()*_ip->protonMass()*_ip->beam1A()));
  //pzgam1 = pz0b1 - pzb1;
  //Q2gam1 = Egam1*Egam1 - (pxgam1*pxgam1 + pygam1*pygam1 + pzgam1*pzgam1);

  Eb2 = E0b2 - Egam2;
  pxb2 = px0b2 - pxgam2;
  pyb2 = py0b2 - pygam2;
  pzb2 = -sqrt(Eb2*Eb2 - (pxb2*pxb2 + pyb2*pyb2 + _ip->protonMass()*_ip->beam2A()*_ip->protonMass()*_ip->beam2A()));
  //pzgam2 = pz0b2 - pzb2;
  //Q2gam2 = Egam2*Egam2 - (pxgam2*pxgam2 + pygam2*pygam2 + pzgam2*pzgam2);
  t2= pt*pt;//not sure

  double pzgamA1, pzgamB1, pzgamA2, pzgamB2;

  pzgamA1 = pz0b1- pzb1;
	pzgamA2 = pz + pzb2 - pz0b2;
  pzgamB1 = pz0b2 - pzb2;
  pzgamB2 = pz - pzgamA1;

  pzgam1 = (2*pzgamA1 + pzgamA2)/3.0;
  pzgam2 = (2*pzgamB1 + pzgamB2)/3.0;
  Q2gam1 = Egam1*Egam1 - (pxgam1*pxgam1 + pygam1*pygam1 + pzgam1*pzgam1);
  Q2gam2 = Egam2*Egam2 - (pxgam2*pxgam2 + pygam2*pygam2 + pzgam2*pzgam2);
}


//______________________________________________________________________________
double Gammagammasingle::pp1(double E)
{
  // First 'copy' of pp, for nucleus 1 form factor.  The split was needed to handle asymmetric beams.  SRK 4/2015
  //  will probably have to pass in beambeamsys? that way we can do beam1.formFactor(t) or beam2..., careful with the way sergey did it for asymmetry
  //  returns on random draw from pp(E) distribution
      
  double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
  double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
  int satisfy =0;
        
  ereds = (E/_GGsingInputGamma_em)*(E/_GGsingInputGamma_em);
  Cm = sqrt(3.)*E/_GGsingInputGamma_em;
  //the amplitude of the p_t spectrum at the maximum
  singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
  //Doing this once and then storing it as a double, for the beam 1 form factor.
  Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
  //pick a test value pp, and find the amplitude there
  x = _randy->Rndom();
  pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); //Will use nucleus #1, there should be two for symmetry//nextline
  singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
  test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

  while(satisfy==0){
    u = _randy->Rndom();
    if(u*Coef <= test){
      satisfy =1;
    }
    else{
      x =_randy->Rndom();
      pp = 5*starlightConstants::hbarc/_bbs.beam1().nuclearRadius()*x;
      singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);//Symmetry
      test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
    }
  }
  return pp;
}

//______________________________________________________________________________
double Gammagammasingle::pp2(double E)
{
  // Second 'copy' of pp, for nucleus 1 form factor.  The split was needed to handle asymmetric beams.  SRK 4/2015
  //  will probably have to pass in beambeamsys? that way we can do beam1.formFactor(t) or beam2..., careful with the way sergey did it for asymmetry
  //  returns on random draw from pp(E) distribution
      
  double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
  double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
  int satisfy =0;
        
  ereds = (E/_GGsingInputGamma_em)*(E/_GGsingInputGamma_em);
  Cm = sqrt(3.)*E/_GGsingInputGamma_em;
  //the amplitude of the p_t spectrum at the maximum
  singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
  //Doing this once and then storing it as a double, which we square later...SYMMETRY?using beam1 for now.
  Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
  //pick a test value pp, and find the amplitude there
  x = _randy->Rndom();
  pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); //Will use nucleus #1, there should be two for symmetry//nextline
  singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
  test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

  while(satisfy==0){
    u = _randy->Rndom();
    if(u*Coef <= test){
      satisfy =1;
    }
    else{
      x =_randy->Rndom();
      pp = 5*starlightConstants::hbarc/_bbs.beam2().nuclearRadius()*x;
      singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);//Symmetry
      test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
    }
  }
  return pp;
}

/**
 * @brief determines the decay products for a decay of this channel.
 * 
 * @param ipid [output reference] - The correct ipid of the daughter particles is outputed
 * @param W [input] - CM energy of collision
 * @param px0 [input] - CM x-momentum
 * @param py0  [input] - CM y-momentum
 * @param pz0 [input] - CM z-momentum
 * @param E1 [output reference] 1st daughter CM energy
 * @param px1 [output reference] 1st daughter CM x-momentum
 * @param py1 [output reference] 1st daughter CM y-momentum
 * @param pz1 [output reference] 1st daughter CM z-momentum
 * @param E2 [output reference] 2nd daughter CM Energy
 * @param px2 [output reference] 2nd daughter CM x-momentum
 * @param py2 [output reference] 2nd daughter CM y-momentum
 * @param pz2 [output reference] 2nd daughter CM z-momentum
 * @param mass [output reference] returns the correct mass of the decay products.
 * @param iFbadevent [output reference] - Sets to 1 if the decay is unsuccessful ONLY. N.B. It only sets, It never resets (but remains in its initial value) if decay is successful.
 */
void Gammagammasingle::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,double W,double px0,double py0,double pz0,double &E1,double &px1,double &py1,double &pz1,double& E2,double &px2,double &py2,double &pz2,double &mass,int &iFbadevent)
{
  //     This routine decays a particle into two particles of mass mdec,
  //     taking spin into account
  
  double mdec=0.;
  double pmag,ytest=0.;
  double phi,theta,xtest,dndtheta,Ecm;
  double  betax,betay,betaz;
  
  //    set the mass of the daughter particles
  switch(_GGsingInputpidtest){ 
  case starlightConstants::ZOVERZ03:
  case starlightConstants::F2:	
    mdec = _ip->pionChargedMass();
    break;
  case starlightConstants::AXION:       // AXION HACK
    mdec = 0;//axion decays to two photons, set mass of decay products to zero
    break;
  case starlightConstants::F2PRIME:
    //  decays 50% to K+/K-, 50% to K_0's
    ytest = _randy->Rndom();
    if(ytest >= 0.5){
      mdec = _ip->kaonChargedMass();
    }
    else{
      mdec = _ip->kaonNeutralMass();
    }
    break;
  default :
    cout<<"No default mass selected for single photon-photon particle, expect errant results"<<endl;
  }
  mass = mdec; //returns the determined mass to the output reference parameter.

  //Calculating the momentum's magnitude
    if(W < 2*mdec){
      cout<<" ERROR: W="<<W<<endl;
      iFbadevent = 1;
      return;
    }
    pmag = sqrt(W*W/4. - mdec*mdec);
//  }
  //     pick an orientation, based on the spin
  //      phi has a flat distribution in 2*pi
  phi = _randy->Rndom()*2.*starlightConstants::pi;
  
  //     find theta, the angle between one of the outgoing particles and
  //    the beamline, in the frame of the two photons
  //this will depend on spin, F2,F2' and z/z03 all have spin 2, all other photonphoton-single mesons are handled by jetset/pythia
  //Applies to spin2 mesons.
 L300td:
  theta = starlightConstants::pi*_randy->Rndom();
  xtest = _randy->Rndom();
  dndtheta = sin(theta)*sin(theta)*sin(theta)*sin(theta)*sin(theta);
  if(xtest > dndtheta)
    goto L300td;

  //     compute unboosted momenta
  px1 = sin(theta)*cos(phi)*pmag;
  py1 = sin(theta)*sin(phi)*pmag;
  pz1 = cos(theta)*pmag;
  px2 = -px1;
  py2 = -py1;
  pz2 = -pz1;
  //        compute energies
  //Changed mass to W 11/9/2000 SRK
  Ecm = sqrt(W*W+px0*px0+py0*py0+pz0*pz0);
  E1 = sqrt(mdec*mdec+px1*px1+py1*py1+pz1*pz1);
  E2 = sqrt(mdec*mdec+px2*px2+py2*py2+pz2*pz2);

  //     Lorentz transform into the lab frame
  // betax,betay,betaz are the boost of the complete system
  betax = -(px0/Ecm);
  betay = -(py0/Ecm);
  betaz = -(pz0/Ecm);
  
  transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
  transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);
  
  
  if(iFbadevent == 1)
    return;
  
  //       change particle id from that of parent to that of daughters

  switch(_GGsingInputpidtest){
    //These decay into a pi+ pi- pair
  case starlightConstants::ZOVERZ03:
  case starlightConstants::F2:
    ipid=starlightConstants::PION;
    break;
  case starlightConstants::AXION:// AXION HACK
    ipid=starlightConstants::PHOTON;  // AXION HACK
    break;      // AXION HACK
  case starlightConstants::F2PRIME:
    if( ytest >= 0.5 )
      {
	//Decays 50/50 into k+ k- or k_s k_l
	ipid=starlightConstants::KAONCHARGE;	
      }
    else
      {
	ipid=starlightConstants::KAONSHORT;
      }	
    break;
  default:
    cout<<"Rethink the daughter particles"<<endl;
  }
}


//______________________________________________________________________________
starlightConstants::event Gammagammasingle::produceEvent(int &/*ievent*/)
{
  // Not in use anymore, default event struct returned
  return starlightConstants::event();
}


/**
 * @brief Produces an event of the specified Two Photon Channel.
 * 
 * @details Approach: properties of decay/daughter particles are temporarily stored in a starlightConstant::event struct and later converted to a upcEvent object.
 * 
 * @param beta The boost vector needed to transform particles from CM to Lab frame. Needed for calculating pseudorapidity in Lab Frame.
 * @return The newly produced event 
 */
//upcEvent Gammagammasingle::produceEvent(vector3 beta)
upcXEvent Gammagammasingle::produceEvent(vector3 beta)
{//this function decays particles and writes events to a file

  upcXEvent event;//keeps record of the produced event and its properties.
  double comenergy = 0.;
  double rapidity = 0.;
  double parentE = 0.;
  double parentmomx=0.,parentmomy=0.,parentmomz=0.;

  
  double Pgam1[4] = {0.0,0.0,0.0,0.0};//Photon from beam1 - Egam1,pxgam1,pygam1,pzgam1
  double Pgam2[4] = {0.0,0.0,0.0,0.0};//Photon from beam2 Egam2,pxgam2,pygam2,pzgam2
  double Pb1[4] = {0.0,0.0,0.0,0.0};//Outgoing beam1 Eb1,pxb1,pyb1,pzb1
  double Pb2[4] = {0.0,0.0,0.0,0.0};//Outgoing beam2 Eb2,pxb2,pyb2,pzb2
  double Q2gam1 =0.,Q2gam2, t=0.;
  
  
  if(_GGsingInputpidtest != starlightConstants::F2 && _GGsingInputpidtest != starlightConstants::F2PRIME && _GGsingInputpidtest != starlightConstants::AXION)
  {
    pickw(comenergy);
    picky(rapidity);
    parentMomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz,
                    Pb1[0],Pb1[1],Pb1[2],Pb1[3],
                    Pb2[0],Pb2[1],Pb2[2],Pb2[3],t,
                    Pgam1[0],Pgam1[1],Pgam1[2],Pgam1[3],Q2gam1,
                    Pgam2[0],Pgam2[1],Pgam2[2],Pgam2[3], Q2gam2);
#ifdef ENABLE_PYTHIA
    starlightParticle particle(parentmomx,parentmomy,parentmomz, parentE, getMass(),_GGsingInputpidtest , 0);
  
    _pyDecayer.addParticle(particle);
  
    return _pyDecayer.execute();
#endif
  }


  int ievent = 0;
  int iFbadevent=0;
  double mass = 0.;
  double mass2 =0; // This is used for events where we have more than two different decays rather than just one. It keeps track of the masses of the species in case they differ.

  

  double ptCutMin2 = _ptCutMin*_ptCutMin;//Used to make Pt Cuts comparison without using square roots.
  double ptCutMax2 = _ptCutMax*_ptCutMax;//Used to make Pt Cuts comparison without using square roots.
  starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
  double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.,E2=0.,E1=0.;
  double px3=0.,px4=0.,py3=0.,py4=0.,pz3=0.,pz4=0.,E3=0.,E4=0.;
  double xtest=0.,ztest=0.;
  bool accepted;
  switch(_GGsingInputpidtest){
  case starlightConstants::ZOVERZ03:

    do{
      iFbadevent = 0; //Ensures that this is reset after every loop cycle.
      pickw(comenergy);
      picky(rapidity);
      parentMomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz,
                    Pb1[0],Pb1[1],Pb1[2],Pb1[3],
                    Pb2[0],Pb2[1],Pb2[2],Pb2[3],t,
                    Pgam1[0],Pgam1[1],Pgam1[2],Pgam1[3],Q2gam1,
                    Pgam2[0],Pgam2[1],Pgam2[2],Pgam2[3], Q2gam2);

  
    //Decays into two pairs.
    parentmomx=parentmomx/2.;
    parentmomy=parentmomy/2.;
    parentmomz=parentmomz/2.;
    accepted = true;//re-initialize the acceptance flag in the loop 
    _nmbAttempts++;
    //Pair #1	
    twoBodyDecay(ipid,comenergy/2.,parentmomx,parentmomy,parentmomz,E1,px1,py1,pz1,E2,px2,py2,pz2,mass,iFbadevent);//Decaying/Producing the first two daughter
    //Pair #2
    
    twoBodyDecay(ipid,comenergy/2.,parentmomx,parentmomy,parentmomz,E3,px3,py3,pz3,E4,px4,py4,pz4,mass2,iFbadevent);//Decaying/producing the 3rd and 4th daughter.
  
    
    if(_ptCutEnabled){//carries out pt_Cut checks
      double pt1chk2 = px1*px1 + py1*py1;// determines the transverse momentum (squared) for daughter 1
      double pt2chk2 = px2*px2 + py2*py2;//similar as above
      double pt3chk2 = px3*px3 + py3*py3;//similar as above
      double pt4chk2 = px4*px4 + py4*py4;//similar as above
      if(pt1chk2 < ptCutMin2 || pt1chk2 > ptCutMax2 || pt2chk2 < ptCutMin2 || pt2chk2 > ptCutMax2 || pt3chk2 < ptCutMin2 || pt3chk2 > ptCutMax2 || pt4chk2 < ptCutMin2 || pt4chk2 > ptCutMax2){//if any of the 4 particles fall outside the range of ptCuts.
        accepted = false;
        continue;
      }
    }
    if(_etaCutEnabled){//Carries out _etaCut checks
      double eta1chk = pseudoRapidityLab(px1,py1,pz1,E1,beta);//computes the pseudorapidity in the Lab Frame.
      double eta2chk = pseudoRapidityLab(px2,py2,pz2,E2,beta);//similar as above
      double eta3chk = pseudoRapidityLab(px3,py3,pz3,E3,beta);//similar as above
      double eta4chk = pseudoRapidityLab(px4,py4,pz4,E4,beta);//similar as above
      if(eta1chk < _etaCutMin || eta1chk >_etaCutMax || eta2chk < _etaCutMin || eta2chk >_etaCutMax || eta3chk < _etaCutMin || eta3chk >_etaCutMax || eta4chk < _etaCutMin || eta4chk >_etaCutMax){//if any of the 4 particles fall outside the range of etaCuts
        accepted = false;
        continue;
      }
    }
    if(accepted and (iFbadevent==0)){//keeps count of the successfully accepted events.
      _nmbAccepted++;
    }
    }while(!accepted || iFbadevent !=0);//repeats the loop if event does not satisfy ptCut etaCut or does not decay successfully

    if (iFbadevent==0){
      
      starlightParticle particle1;
      starlightParticle particle2;
      starlightParticle particle3;
      starlightParticle particle4;

      xtest = _randy->Rndom();
      ztest = _randy->Rndom();
      //Assigning charges randomly.
      if (xtest<0.5){
        particle1.setCharge(1);//q1=1
        particle2.setCharge(-1);//q2=-1
      }
      else{
        particle1.setCharge(-1);//q1=-1
        particle2.setCharge(1);//q2=1
      }
      if (ztest<0.5){
        particle3.setCharge(1);//q3=1
        particle4.setCharge(-1);//q4=-1
      }
      else{
        particle3.setCharge(-1);//q3=-1
        particle4.setCharge(1);//q4=1
      }
      //storing the remaining event's track details
      //Track #1
      particle1.SetPxPyPzE(px1,py1,pz1,E1);
      particle1.setMass(mass);
      particle1.setPdgCode(ipid*particle1.getCharge());
      event.addParticle(particle1);

      //Track #2
      particle2.SetPxPyPzE(px2,py2,pz2,E2);
      particle2.setMass(mass);
      particle2.setPdgCode(ipid*particle2.getCharge());
      event.addParticle(particle2);

      //Track #3
      particle3.SetPxPyPzE(px3,py3,pz3,E3);
      particle3.setMass(mass2);
      particle3.setPdgCode(ipid*particle3.getCharge());
      event.addParticle(particle3);

      //Track #4
      particle4.SetPxPyPzE(px4,py4,pz4,E4);
      particle4.setMass(mass2);
      particle4.setPdgCode(ipid*particle4.getCharge());
      event.addParticle(particle4);

      //add information about the meadiating photons and beams when possible

      if(_ip->giveExtraBeamInfo()){//is it possible to obtain extra information
        lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
        lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
        double targetEgamma1, targetEgamma2, rap1cm = acosh(_ip->beamLorentzGamma()), cmsEgam1 = Pgam1[0];
        double cmsEgam2 = Pgam2[0], Pzgam1 = Pgam1[3], Pzgam2 = Pgam2[3];
        lorentzVector gamma1(Pgam1[1],Pgam1[2],Pzgam1,cmsEgam1);
        lorentzVector gamma2(Pgam2[1],Pgam2[2],Pzgam2,cmsEgam2);
        lorentzVector vmeson(parentmomx,parentmomy,parentmomz,parentE/2);
        lorentzVector vmeson2(parentmomx,parentmomy,parentmomz,parentE/2);

        targetEgamma2 = cmsEgam2*cosh(rap1cm) - Pzgam2*sinh(rap1cm);//beam 1 is target - hence for gamma2
        targetEgamma1 = cmsEgam1*cosh(rap1cm) + Pzgam1*sinh(rap1cm);//beam2 is target - hence for gamma1

        event.addMeson(vmeson);//rho
        event.addMeson(vmeson2);//rho
        event.addGammaFromBeam1(gamma1,targetEgamma1,Q2gam1);//emmitted by beam1.
        event.addGammaFromBeam2(gamma2,targetEgamma2,Q2gam2);//emmitted by beam2
        event.addOutgoingBeams(beam1,beam2);
        event.addVertext(t);
      }
                
      ievent=ievent+1;
      return event;
    }	
    
    break;
  case starlightConstants::F2:
  case starlightConstants::F2PRIME:
    do{
      iFbadevent = 0; //Resets this after every loop cycle-to avoid an infinite loop
      pickw(comenergy);
      picky(rapidity);
      parentMomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz,
                    Pb1[0],Pb1[1],Pb1[2],Pb1[3],
                    Pb2[0],Pb2[1],Pb2[2],Pb2[3],t,
                    Pgam1[0],Pgam1[1],Pgam1[2],Pgam1[3],Q2gam1,
                    Pgam2[0],Pgam2[1],Pgam2[2],Pgam2[3], Q2gam2);
      
      accepted = true;//Reinitialize the acceptance flag after every loop cycle
      _nmbAttempts++;
      twoBodyDecay(ipid,comenergy,parentmomx,parentmomy,parentmomz,E1,px1,py1,pz1,E2,px2,py2,pz2,mass,iFbadevent);//Decaying/producing daughter particles

      if(_ptCutEnabled){//Carries out Momentum cut Checks
        double pt1chk2 = px1*px1 + py1*py1;//computes transverse momentum (squared) of daughter 1
        double pt2chk2 = px2*px2 + py2*py2;// similar as above
        
        if(pt1chk2 < ptCutMin2 || pt1chk2 > ptCutMax2 || pt2chk2 < ptCutMin2 || pt2chk2 > ptCutMax2 ){//If any of both particles fall outside the ptCut range
          accepted = false;
          continue;
        }
      }
      if(_etaCutEnabled){//Carries out eta Cut Checks
        double eta1chk = pseudoRapidityLab(px1,py1,pz1,E1,beta);//Computes eta in Lab frame for daughter 1.
        double eta2chk = pseudoRapidityLab(px2,py2,pz2,E2,beta);//Similar as above
        
        if(eta1chk < _etaCutMin || eta1chk >_etaCutMax || eta2chk < _etaCutMin || eta2chk >_etaCutMax ){//if any of both particles fall outside the etaCut range
          accepted = false;
          continue;
        }
      }
      if(accepted and (iFbadevent == 0)){
        _nmbAccepted++;//maintain count of successfully accepted events
      }
    }while(!accepted || iFbadevent != 0);//repeats loop if ptCut, etaCut or succesful decay requirements are not satisfied.
    

    if (iFbadevent==0){
      starlightParticle particle1;
      starlightParticle particle2;
      
      xtest = _randy->Rndom();
      if (xtest<0.5){// randomly sets the charge of the particles
	      particle1.setCharge(1);//q1=1
        particle2.setCharge(-1);//q2=-1
      }
      else{
        particle1.setCharge(-1);//q1=-1
        particle2.setCharge(1);//q2=1
      }
      //storing the event details	
      //Track #1
      particle1.SetPxPyPzE(px1,py1,pz1,E1);
      particle1.setMass(mass);
      particle1.setPdgCode(ipid*particle1.getCharge());
      event.addParticle(particle1);

      //Track #2
      particle2.SetPxPyPzE(px2,py2,pz2,E2);
      particle2.setMass(mass);
      particle2.setPdgCode(ipid*particle2.getCharge());
      event.addParticle(particle2);

      //adds information about mediating photons and outgoing beams when possible

      if(_ip->giveExtraBeamInfo()){//is it possible to obtain extra information?
        lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
        lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
        double targetEgamma1, targetEgamma2, rap1cm = acosh(_ip->beamLorentzGamma()), cmsEgam1 = Pgam1[0];
        double cmsEgam2 = Pgam2[0], Pzgam1 = Pgam1[3], Pzgam2 = Pgam2[3];
        lorentzVector gamma1(Pgam1[1],Pgam1[2],Pzgam1,cmsEgam1);
        lorentzVector gamma2(Pgam2[1],Pgam2[2],Pzgam2,cmsEgam2);
        lorentzVector vmeson(parentmomx,parentmomy,parentmomz,parentE);

        targetEgamma2 = cmsEgam2*cosh(rap1cm) - Pzgam2*sinh(rap1cm);//beam 1 is target - hence for gamma2
        targetEgamma1 = cmsEgam1*cosh(rap1cm) + Pzgam1*sinh(rap1cm);//beam2 is target - hence for gamma1

        event.addMeson(vmeson);//F2 or F2_Prime
        event.addGammaFromBeam1(gamma1,targetEgamma1,Q2gam1);//emmitted by beam1.
        event.addGammaFromBeam2(gamma2,targetEgamma2,Q2gam2);//emmitted by beam2
        event.addOutgoingBeams(beam1,beam2);
        event.addVertext(t);
      }
      
      ievent=ievent+1;
    }
    break;

  case starlightConstants::AXION:      // AXION HACK, start
    do{
      iFbadevent =0;//resets variable after every loop cycle - to avoid infinite loops
      pickw(comenergy);
      picky(rapidity);
      parentMomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz,
                    Pb1[0],Pb1[1],Pb1[2],Pb1[3],
                    Pb2[0],Pb2[1],Pb2[2],Pb2[3],t,
                    Pgam1[0],Pgam1[1],Pgam1[2],Pgam1[3],Q2gam1,
                    Pgam2[0],Pgam2[1],Pgam2[2],Pgam2[3], Q2gam2);
      accepted = true;//reinitializes the acceptance flag after every loop cycle - avoid inifinte loop
      _nmbAttempts++;
      twoBodyDecay(ipid,comenergy,parentmomx,parentmomy,parentmomz,E1,px1,py1,pz1,E2,px2,py2,pz2,mass,iFbadevent);//decaying/producing daughter particles

      if(_ptCutEnabled){//apply ptCut checks
        double pt1chk2 = px1*px1 + py1*py1;//computing the transverse momentum (squared) for daughter 1.
        double pt2chk2 = px2*px2 + py2*py2;//similar as above
        
        if(pt1chk2 < ptCutMin2 || pt1chk2 > ptCutMax2 || pt2chk2 < ptCutMin2 || pt2chk2 > ptCutMax2 ){//if any of the 2 particles fall outside the range
          accepted = false;
          continue;
        }
      }
      if(_etaCutEnabled){//apply eta Cut checks
        double eta1chk = pseudoRapidityLab(px1,py1,pz1,E1,beta);//Computes eta in Lab frame for daughter 1
        double eta2chk = pseudoRapidityLab(px2,py2,pz2,E2,beta);//similar as above
        
        if(eta1chk < _etaCutMin || eta1chk >_etaCutMax || eta2chk < _etaCutMin || eta2chk >_etaCutMax ){//if any of the two particles fall outside the etaCut range
          accepted = false;
          continue;
        }
      }
      if(accepted and (iFbadevent == 0)){
        _nmbAccepted++;//Keeps count of successfully accepted events
      }
    }while(!accepted || iFbadevent != 0);//repeats loop if ptCut,EtaCut or successful decay requirements are not satisfied.
    

    if (iFbadevent==0){
      //storing the event's details.
      starlightParticle particle1;
      starlightParticle particle2;


      particle1.setCharge(0);//q1=0;
      particle2.setCharge(0);//q2=0;

      //Track #1
      particle1.SetPxPyPzE(px1,py1,pz1,E1);
      particle1.setMass(mass);
      particle1.setPdgCode(ipid);
      event.addParticle(particle1);

      //Track #2
      particle2.SetPxPyPzE(px2,py2,pz2,E2);
      particle2.setMass(mass);
      particle2.setPdgCode(ipid);
      event.addParticle(particle2);

      //adds information about the mediating photon and the outgoing beams

      if(_ip->giveExtraBeamInfo()){//is it possible to obtain extra information
        lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
        lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
        double targetEgamma1, targetEgamma2, rap1cm = acosh(_ip->beamLorentzGamma()), cmsEgam1 = Pgam1[0];
        double cmsEgam2 = Pgam2[0], Pzgam1 = Pgam1[3], Pzgam2 = Pgam2[3];
        lorentzVector gamma1(Pgam1[1],Pgam1[2],Pzgam1,cmsEgam1);
        lorentzVector gamma2(Pgam2[1],Pgam2[2],Pzgam2,cmsEgam2);
        lorentzVector vmeson(parentmomx,parentmomy,parentmomz,parentE);

        targetEgamma2 = cmsEgam2*cosh(rap1cm) - Pzgam2*sinh(rap1cm);//beam 1 is target - hence for gamma2
        targetEgamma1 = cmsEgam1*cosh(rap1cm) + Pzgam1*sinh(rap1cm);//beam2 is target - hence for gamma1

        event.addVectorMeson(vmeson);//Axion
        event.addGammaFromBeam1(gamma1,targetEgamma1,Q2gam1);//emmitted by beam1.
        event.addGammaFromBeam2(gamma2,targetEgamma2,Q2gam2);//emmitted by beam2
        event.addOutgoingBeams(beam1,beam2);
        event.addVertext(t);        
      }

      ievent=ievent+1;

    }
    break;  // AXION HACK, end


  default:
    break;
  }
  
  
  return event;
}


//______________________________________________________________________________
double Gammagammasingle::getMass()
{
  using namespace starlightConstants;
  double singlemass=0.;
  switch(_GGsingInputpidtest){
  case starlightConstants::ETA:
    singlemass = _ip->etaMass();
    break;
  case starlightConstants::ETAPRIME:
    singlemass = _ip->etaPrimeMass();
    break;
  case starlightConstants::ETAC:
    singlemass = _ip->etaCMass();
    break;
  case starlightConstants::F0:
    singlemass = _ip->f0Mass();
    break;
  case starlightConstants::F2:
    singlemass = _ip->f2Mass();
    break;
  case starlightConstants::A2:
    singlemass = _ip->a2Mass();
    break;
  case starlightConstants::F2PRIME:
    singlemass = _ip->f2PrimeMass();
    break;
  case starlightConstants::ZOVERZ03:
    singlemass = _ip->zoverz03Mass();
    break;
  case starlightConstants::AXION: // AXION HACK
    singlemass = _axionMass;      // AXION HACK
    break; // AXION HACK
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getmass(), mass = 0."<<endl;
  }
  return singlemass;
}


//______________________________________________________________________________
double Gammagammasingle::getWidth()
{

  /* Partial widths(GAMMA(gammgamma)) taken from PDG 2014- Chinese Physics C 38, no 9, Sept. 2014.*/
  double singlewidth=0.;
  switch(_GGsingInputpidtest){
  case starlightConstants::ETA:
    singlewidth = _ip->etaPartialggWidth();
    break;
  case starlightConstants::ETAPRIME:
    singlewidth = _ip->etaPrimePartialggWidth();
    break;
  case starlightConstants::ETAC:
    singlewidth = _ip->etaCPartialggWidth();
    break;
  case starlightConstants::F0:
    singlewidth = _ip->f0PartialggWidth();
    break;
  case starlightConstants::F2:
    singlewidth = _ip->f2PartialggWidth();
    break;
  case starlightConstants::A2:
    singlewidth = _ip->a2PartialggWidth();
    break;
  case starlightConstants::F2PRIME:
    singlewidth = _ip->f2PrimePartialggWidth();
    break;
  case starlightConstants::ZOVERZ03:
    singlewidth = _ip->zoverz03PartialggWidth();
    break;
  case starlightConstants::AXION: // AXION HACK
    singlewidth = 1/(64*starlightConstants::pi)*_axionMass*_axionMass*_axionMass/(1000*1000);//Fix Lambda=1000 GeV,rescaling is trivial.    // AXION HACK
    break;
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getwidth(), width = 0."<<endl;
  }
  return singlewidth; 
}


//______________________________________________________________________________
double Gammagammasingle::getSpin()
{
  double singlespin=0.5;
  switch(_GGsingInputpidtest){
  case starlightConstants::ETA:
    singlespin = _ip->etaSpin();
    break;
  case starlightConstants::ETAPRIME:
    singlespin = _ip->etaPrimeSpin();
    break;
  case starlightConstants::ETAC:
    singlespin = _ip->etaCSpin();
    break;
  case starlightConstants::F0:
    singlespin = _ip->f0Spin();
    break;
  case starlightConstants::F2:
    singlespin = _ip->f2Spin();
    break;
  case starlightConstants::A2:
    singlespin = _ip->a2Spin();
    break;
  case starlightConstants::F2PRIME:
    singlespin = _ip->f2PrimeSpin();
    break;
  case starlightConstants::ZOVERZ03:
    singlespin = _ip->zoverz03Spin();
    break;
  case starlightConstants::AXION:// AXION HACK
    singlespin = _ip->axionSpin();// AXION HACK
    break;// AXION HACK
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getspin(), spin = 0."<<endl;
  }
  return singlespin;
}



