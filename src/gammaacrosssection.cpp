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
//    Converted all nuclear formfactors to come from the scattering nucleus(nucleus2)
//    Added incoherent condition to the cross-section that follows a similar approach as pp
//    Could not figure out the scaling exactly for incoherent(possibly units) so divided by
//    1E-4 and there is a incoherence factor that can be selected in the input file,
//    starlight.in--JWB
//    Also, it has not been implemented for interference.
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "gammaacrosssection.h"
#include "bessel.h"


using namespace std;


//______________________________________________________________________________
photonNucleusCrossSection::photonNucleusCrossSection (inputParameters& input, beamBeamSystem& bbsystem)
: _bbs(bbsystem)
{
  _sigmaProtonEnergy=input.getProtonEnergy();
  _sigmaGamma_em=input.getgamma_em();
  _sigmaPID=input.getPidTest();
  _sigmaBreakup=input.getBreakupMode();
  _sigmaCoherence=input.getIncoherentOrCoherent();
  _sigmaCoherenceFactor=input.getIncoherentFactor();
  _sigmaNucleus=_bbs.getBeam2().getAin();

  switch(_bbs.getBeam1().getZin())
    {
    case 79://Au
      _lum=2.0;
      break;
    case 53://I
      _lum=27.;
      break;
    case 49://Indium,uses same as Iodine
      _lum=27.;
      break;
    case 29://Cu
      _lum=95.;
      break;
    case 14://Si
      _lum=440.;
      break;
    case 8://O
      _lum=980.;
      break;
    case 82://Pb
      _lum=1.;
      break;
    case 20://Ca
      _lum=2000.;
      break;
    case 1://proton
      _lum=1.E8;
      break;
    default:
      cout <<"Warning:Luminosity not defined.Gammaacrosssection::getlum"<<endl;
    }
  switch(_sigmaPID)
    {
    case starlightConstants::RHO:
      _bSlope= 11.0;
      _f2o4pi= 2.02;
      _ANORM= -2.75;
      _BNORM= 0.0;
      _defaultC=     1.0;
      _channelMass=0.7685;
      _width=0.1507;
      break;
    case starlightConstants::RHOZEUS:
      _bSlope=11.0;
      _f2o4pi=2.02;
      _ANORM=-2.75;
      _BNORM=1.84;
      _defaultC=1.0;
      _channelMass=0.7685;
      _width=0.1507;
      break;
    case starlightConstants::FOURPRONG:
      _bSlope      = 11.0;
      _f2o4pi      = 2.02;
      _ANORM       = -2.75;
      _BNORM       = 0;  // no coherent background component is implemented for four-prong
      _defaultC    = 11.0;
      _channelMass = 1.350;
      _width       = 0.360;
      break;
    case starlightConstants::OMEGA:
      _bSlope=10.0;
      _f2o4pi=23.13;
      _ANORM=-2.75;
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=0.78194;
      _width=0.00843;
      break;
    case starlightConstants::PHI:
      _bSlope=7.0;
      _f2o4pi=13.71;
      _ANORM=-2.75;
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=1.019413;
      _width=0.00443;
      break;
    case starlightConstants::JPSI:
    case starlightConstants::JPSI_ee:
    case starlightConstants::JPSI_mumu:
      _bSlope=4.0;
      _f2o4pi=10.45;
      _ANORM=-2.75;//Artificial Breit-Wigner parameters--no direct pions
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=3.09692;//JN 3.09688
      _width=0.000091;//JN 0.000087
      break;
    case starlightConstants::JPSI2S:
    case starlightConstants::JPSI2S_ee:
    case starlightConstants::JPSI2S_mumu:
      _bSlope=4.3;
      _f2o4pi=26.39;
      _ANORM=-2.75;//Artificial
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=3.686093;
      _width=0.000337;
      break;
    case starlightConstants::UPSILON:
    case starlightConstants::UPSILON_ee:
    case starlightConstants::UPSILON_mumu:
      _bSlope=4.0;
      _f2o4pi=125.37;
      _ANORM=-2.75;//Artificial
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=9.46030;
      _width=0.00005402;
      break;
    case starlightConstants::UPSILON2S:
    case starlightConstants::UPSILON2S_ee:
    case starlightConstants::UPSILON2S_mumu:
      _bSlope=4.0;
      _f2o4pi=290.84;
      _ANORM=-2.75;
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=10.02326;
      _width=0.00003198;
      break;
    case starlightConstants::UPSILON3S:
    case starlightConstants::UPSILON3S_ee:
    case starlightConstants::UPSILON3S_mumu:
      _bSlope=4.0;
      _f2o4pi=415.10;
      _ANORM=-2.75;
      _BNORM=0.0;
      _defaultC=1.0;
      _channelMass=10.3552;
      _width=0.00002032;
      break;
    default:
      cout <<"No sigma constants parameterized for pid: "<<_sigmaPID
	   <<" GammaAcrosssection"<<endl;
    }

  _EgMax= 4.*_sigmaGamma_em*starlightConstants::hbarc/_bbs.getBeam1().RNuc(); 
  //Max photon energy( for VM only, in GeV, lab frame, use beam energy
  //, nuclear size cutoff)

}

//______________________________________________________________________________
photonNucleusCrossSection::~photonNucleusCrossSection()
{ }


//______________________________________________________________________________
beamBeamSystem photonNucleusCrossSection::getbbs()
{
  return _bbs;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getBNORM()
{
  return _BNORM;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getLum()
{
  return _lum;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getf2o4pi()
{
  return _f2o4pi;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getChannelMass()
{
  return _channelMass;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getbslope()
{
  return _bSlope;
}


//______________________________________________________________________________
void photonNucleusCrossSection::crossSectionCalculation(double)
{
  cout << "Neither narrow/wide resonance cross-section calculation.--Derived" << endl;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getcsgA(double Egamma, double W)
{
  //This function returns the cross-section for photon-nucleus interaction 
  //producing vectormesons
  
  double Av,Wgp,cs,cvma;
  double t,tmin,tmax;
  double csgA,ax,bx;
  int NGAUSS;                                                                                                                                       
  
  //     DATA FOR GAUSS INTEGRATION
  double xg[6]={0,0.1488743390,0.4333953941,0.6794095683,0.8650633667,0.9739065285};
  double ag[6]={0,0.2955242247,0.2692667193,0.2190863625,0.1494513492,0.0666713443};
  NGAUSS = 6;
  
  //       Find gamma-proton CM energy
  Wgp=sqrt(2.*Egamma*(_sigmaProtonEnergy+sqrt(_sigmaProtonEnergy*_sigmaProtonEnergy-
					     starlightConstants::mp*starlightConstants::mp))+starlightConstants::mp*starlightConstants::mp);
	
  //Used for d-A and A-A
  tmin   = (W*W/(4.*Egamma*_sigmaGamma_em) )*(W*W/(4.*Egamma*_sigmaGamma_em) );
  
  if(_bbs.getBeam1().getAin()==1&&_bbs.getBeam2().getAin()==1)
    {  //Proton-proton, no scaling needed.
      csgA = sigmagp(Wgp);
    }
  else if(_bbs.getBeam2().getZin()==1&&_bbs.getBeam2().getAin()==2)
    {  

      // Deuteron-A interaction
      Av = _bSlope*sigmagp(Wgp);
      
      tmax   = tmin + 0.64;   //0.64
      ax     = 0.5*(tmax-tmin);
      bx     = 0.5*(tmax+tmin);
      csgA   = 0.;
      
      for( int k=1;k<NGAUSS;k++){ 
        t    = ax*xg[k]+bx;
        //We use beam2 here since the input stores the deuteron as nucleus 2
        //and nucleus 2 is the pomeron field source
        //Also this is the way sergey formatted the formfactor.
        csgA = csgA + ag[k]*_bbs.getBeam2().formFactor(t); 
        t    = ax*(-xg[k])+bx;
        csgA = csgA + ag[k]*_bbs.getBeam2().formFactor(t);
      }
      csgA = 0.5*(tmax-tmin)*csgA;
      csgA = Av*csgA;
    }
  else if(_sigmaCoherence==0&&(!(_bbs.getBeam2().getZin()==1&&_bbs.getBeam2().getAin()==2)))
    {

      // For incoherent AA interactions , since incoherent treating it as gamma-p
      // Calculate the differential V.M.+proton cross section
      csgA = 1.E-4*_sigmaCoherenceFactor*_sigmaNucleus*_bSlope*sigmagp(Wgp);//artifical 1E-3 to scale down sigma

      //Calculating int |F(t)| dt
      //Using proton formfactor for this case
      //Note the coherence scaling factor being intergrated with the F(t)
      //Should it just be F(t)^2?   
      //Pay attention to the way the formfactor is implemented in nucleus class
      //Also, notice the tmin value.  It starts higher for dAu, should we proceed
      //in a similar fashion
      //Why don't we use formfactors for pp? Is it because it is incorporated in 
      //the gamma-p fits done for dsigma/dt? Yes?

    }
  else 
    {	

      // For typical AA interactions.
      // Calculate V.M.+proton cross section
      cs=sqrt(16.*starlightConstants::pi*_f2o4pi*_bSlope
	      *starlightConstants::hbarc*starlightConstants::hbarc*sigmagp(Wgp)
	      /starlightConstants::alpha);
    
      //  Calculate V.M.+nucleus cross section
      cvma=sigma_A(cs); 
 
      // Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
      Av=(starlightConstants::alpha*cvma*cvma)/
        (16.*starlightConstants::pi*_f2o4pi*starlightConstants::hbarc*starlightConstants::hbarc);
    
      tmax   = tmin + 0.25;
      ax     = 0.5*(tmax-tmin);
      bx     = 0.5*(tmax+tmin);
      csgA   = 0.;
      for( int k=1;k<NGAUSS;k++){ 
        t    = ax*xg[k]+bx;
        csgA = csgA + ag[k]*_bbs.getBeam2().formFactor(t)*_bbs.getBeam2().formFactor(t);
        t    = ax*(-xg[k])+bx;
        csgA = csgA + ag[k]*_bbs.getBeam2().formFactor(t)*_bbs.getBeam2().formFactor(t);
      }
    
      csgA = 0.5*(tmax-tmin)*csgA;
      csgA = Av*csgA;

    }
  return csgA;	
}


//______________________________________________________________________________
double photonNucleusCrossSection::photonFlux(double Egamma)
{
  // This routine gives the photon flux as a function of energy Egamma
  // It works for arbitrary nuclei and gamma; the first time it is
  // called, it calculates a lookup table which is used on
  // subsequent calls.
  // It returns dN_gamma/dE (dimensions 1/E), not dI/dE
  // energies are in GeV, in the lab frame
  // rewritten 4/25/2001 by SRK
  
  double lEgamma,Emin,Emax;
  static double lnEmax, lnEmin, dlnE;
  double stepmult,energy,rZ,rA;
  int nbstep,nrstep,nphistep,nstep;
  double bmin,bmax,bmult,biter,bold,integratedflux;
  double fluxelement,deltar,riter;
  double deltaphi,phiiter,dist;
  static double dide[401];
  double lnElt;
  double rA2, rZ2; //Added sergey
  double flux_r; //Returns the flux.
  double Xvar;
  int Ilt;
  double RNuc=0.,RNuc2=0.;

  RNuc=_bbs.getBeam1().RNuc();
  RNuc2=_bbs.getBeam2().RNuc();
  // static ->>> dide,lnEMax,lnEmin,dlnE
  static int  Icheck = 0;
  
  //Check first to see if pp (JN0705)
  if( _bbs.getBeam1().getAin()==1 && _bbs.getBeam2().getAin()==1 ){
    int nbsteps = 200;
    double bmin = 0.5;
    double bmax = 5.0 + (5.0*_sigmaGamma_em*starlightConstants::hbarc/Egamma);
    double dlnb = (log(bmax)-log(bmin))/(1.*nbsteps);

    double local_sum=0.0;

    // Impact parameter loop 
    for (int i = 0; i<=nbsteps;i++){

      double bnn0 = bmin*exp(i*dlnb);
      double bnn1 = bmin*exp((i+1)*dlnb);
      double db   = bnn1-bnn0;
        
      //      double PofB0 = 1.0; 
      //      if( bnn0 > 1.4 )PofB0=0.0;
      //      double PofB1 = 1.0; 
      //      if( bnn1 > 1.4 )PofB1=0.0;
      
      double ppslope = 19.0; 
      double GammaProfile = exp(-bnn0*bnn0/(2.*starlightConstants::hbarc*starlightConstants::hbarc*ppslope));  
      double PofB0 = 1. - (1. - GammaProfile)*(1. - GammaProfile);   
      GammaProfile = exp(-bnn1*bnn1/(2.*starlightConstants::hbarc*starlightConstants::hbarc*ppslope));  
      double PofB1 = 1. - (1. - GammaProfile)*(1. - GammaProfile);   

      double Xarg = Egamma*bnn0/(starlightConstants::hbarc*_sigmaGamma_em);
      double loc_nofe0 = (_bbs.getBeam1().getZin()*_bbs.getBeam1().getZin()*starlightConstants::alpha)/
	(starlightConstants::pi*starlightConstants::pi); 
      loc_nofe0 *= (1./(Egamma*bnn0*bnn0)); 
      loc_nofe0 *= Xarg*Xarg*(bessel::dbesk1(Xarg))*(bessel::dbesk1(Xarg)); 

      Xarg = Egamma*bnn1/(starlightConstants::hbarc*_sigmaGamma_em);
      double loc_nofe1 = (_bbs.getBeam1().getZin()*_bbs.getBeam1().getZin()*starlightConstants::alpha)/
	(starlightConstants::pi*starlightConstants::pi); 
      loc_nofe1 *= (1./(Egamma*bnn1*bnn1)); 
      loc_nofe1 *= Xarg*Xarg*(bessel::dbesk1(Xarg))*(bessel::dbesk1(Xarg)); 

      local_sum += loc_nofe0*(1. - PofB0)*bnn0*db; 
      local_sum += loc_nofe1*(1. - PofB1)*bnn1*db; 

    }
    // End Impact parameter loop 

    // Note: 2*pi --> pi because of no factor 2 above 
    double flux_r=local_sum*starlightConstants::pi; 
    return flux_r;

    //    bmin = RNuc+RNuc;
    //    flux_r = nepoint(Egamma,bmin);
    //    return flux_r;
  }

  //   first call?  - initialize - calculate photon flux
  Icheck=Icheck+1;
  if(Icheck > 1) goto L1000f;
  
  rZ=double(_bbs.getBeam1().getZin());
  rA=double(_bbs.getBeam1().getAin());
  rZ2=double(_bbs.getBeam2().getZin());  //Sergey--dAu
  rA2=double(_bbs.getBeam2().getAin());  //Sergey
  
  //  Nuclear breakup is done by PofB
  //  collect number of integration steps here, in one place
  
  nbstep=400;
  nrstep=60;
  nphistep=40;
  
  //  this last one is the number of energy steps
  nstep=100;
  
  // following previous choices, take Emin=10 keV at LHC, Emin = 1 MeV at RHIC
  
  Emin=1.E-5;
  if (_sigmaGamma_em < 500) 
    Emin=1.E-3;
  
  //  maximum energy is 12 times the cutoff
  //  25 GeV for gold at RHIC, 650 GeV for lead at LHC
  
  Emax=12.*starlightConstants::hbarc*_sigmaGamma_em/RNuc;
  //Will this be diff for dAu?
  
  //     >> lnEmin <-> ln(Egamma) for the 0th bin
  //     >> lnEmax <-> ln(Egamma) for the last bin
  
  lnEmin=log(Emin);
  lnEmax=log(Emax);
  dlnE=(lnEmax-lnEmin)/nstep;                                                                                                                  

  cout<<" Calculating flux for photon energies from E= "<<Emin 
      <<" to  "<<Emax<<"  GeV (lab frame) "<<endl;


  stepmult= exp(log(Emax/Emin)/double(nstep));
  energy=Emin;
  
  for (int j = 1; j<=nstep;j++){
    energy=energy*stepmult;
    
    //  integrate flux over 2R_A < b < 2R_A+ 6* gamma hbar/energy
    //  use exponential steps
    
    bmin=RNuc+RNuc2; //2.*RNuc; Sergey
    bmax=bmin + 6.*starlightConstants::hbarc*_sigmaGamma_em/energy;
    
    bmult=exp(log(bmax/bmin)/double(nbstep));
    biter=bmin;
    integratedflux=0.;
    
    if (_bbs.getBeam2().getZin()==1&&_bbs.getBeam1().getAin()==2){
      //This is for deuteron-gold
      Xvar = (RNuc+RNuc2)*energy/(starlightConstants::hbarc*(_sigmaGamma_em));
      
      fluxelement = (2.0/starlightConstants::pi)*rZ*rZ*starlightConstants::alpha/
	energy*(Xvar*bessel::dbesk0(Xvar)*bessel::dbesk1(Xvar)-(1/2)*Xvar*Xvar*
		(bessel::dbesk1(Xvar)*bessel::dbesk1(Xvar)-bessel::dbesk0(Xvar)*bessel::dbesk0(Xvar)));
      
      integratedflux=integratedflux+fluxelement;
                
    }//if dAu
    else{ 
      for (int jb = 1; jb<=nbstep;jb++){
	bold=biter;
	biter=biter*bmult;
	// When we get to b>20R_A change methods - just take the photon flux
	//  at the center of the nucleus.
	if (biter > (10.*RNuc))
	  {
	    // if there is no nuclear breakup or only hadronic breakup, which only
	    // occurs at smaller b, we can analytically integrate the flux from b~20R_A
	    // to infinity, following Jackson (2nd edition), Eq. 15.54
	    Xvar=energy*biter/(starlightConstants::hbarc*_sigmaGamma_em);
	    // Here, there is nuclear breakup.  So, we can't use the integrated flux
	    //  However, we can do a single flux calculation, at the center of the
	    //  nucleus
	    
	    // Eq. 41 of Vidovic, Greiner and Soff, Phys.Rev.C47,2308(1993), among other places
	    //  this is the flux per unit area
	    fluxelement  = (rZ*rZ*starlightConstants::alpha*energy)*
	      (bessel::dbesk1(Xvar))*(bessel::dbesk1(Xvar))/
	      ((starlightConstants::pi*_sigmaGamma_em*starlightConstants::hbarc)*
	       (starlightConstants::pi*_sigmaGamma_em*starlightConstants::hbarc));
	    
	  }//if biter>10
	else{
	  // integrate over nuclear surface. n.b. this assumes total shadowing -
	  // treat photons hitting the nucleus the same no matter where they strike
	  fluxelement=0.;
	  deltar=RNuc/double(nrstep);
	  riter=-deltar/2.;
          
	  for (int jr =1; jr<=nrstep;jr++){
	    riter=riter+deltar;
	    // use symmetry;  only integrate from 0 to pi (half circle)
	    deltaphi=starlightConstants::pi/double(nphistep);
	    phiiter=0.;
            
	    for( int jphi=1;jphi<= nphistep;jphi++){
	      phiiter=(double(jphi)-0.5)*deltaphi;
	      //  dist is the distance from the center of the emitting nucleus to the point in question
	      dist=sqrt((biter+riter*cos(phiiter))*(biter+riter*
						    cos(phiiter))+(riter*sin(phiiter))*(riter*sin(phiiter)));
	      
	      Xvar=energy*dist/(starlightConstants::hbarc*_sigmaGamma_em);				
	      
	      flux_r = (rZ*rZ*starlightConstants::alpha*energy)*
		(bessel::dbesk1(Xvar)*bessel::dbesk1(Xvar))/
		((starlightConstants::pi*_sigmaGamma_em*starlightConstants::hbarc)*
		 (starlightConstants::pi*_sigmaGamma_em*starlightConstants::hbarc));
	      
	      //  The surface  element is 2.* delta phi* r * delta r
	      //  The '2' is because the phi integral only goes from 0 to pi
	      fluxelement=fluxelement+flux_r*2.*deltaphi*riter*deltar;
	      //  end phi and r integrations
	    }//for(jphi)
	  }//for(jr)
	  //  average fluxelement over the nuclear surface
	  fluxelement=fluxelement/(starlightConstants::pi*RNuc*RNuc);
	}//else
	//  multiply by volume element to get total flux in the volume element
	fluxelement=fluxelement*2.*starlightConstants::pi*biter*(biter-bold);
	//  modulate by the probability of nuclear breakup as f(biter)
	if (_sigmaBreakup > 1){
	  fluxelement=fluxelement*_bbs.probabilityOfBreakup(biter);
	}
	integratedflux=integratedflux+fluxelement;
	
      }//end of for
    }  //end of else
    // end energy integration
    // nobody going here any more
    
    //  In lookup table, store k*dN/dk because it changes less
    //  so the interpolation should be better
    
    dide[j]=integratedflux*energy;
                                     
  }//end of for.
       
  //  for 2nd and subsequent calls, use lookup table immediately
  
 L1000f:
  
  lEgamma=log(Egamma);
  if (lEgamma < (lnEmin+dlnE) ||  lEgamma  > lnEmax){
    flux_r=0.0;
    cout<<"  ERROR: Egamma outside defined range. Egamma= "<<Egamma
	<<"   "<<lnEmax<<" "<<(lnEmin+dlnE)<<endl;
  }
  else{
    //       >> Egamma between Ilt and Ilt+1
    Ilt = int((lEgamma-lnEmin)/dlnE);
    //       >> ln(Egamma) for first point 
    lnElt = lnEmin + Ilt*dlnE; 
    //       >> Interpolate
    flux_r = dide[Ilt] + ((lEgamma-lnElt)/dlnE)*(dide[Ilt+1]- dide[Ilt]);
    flux_r = flux_r/Egamma;
  }
  
  return flux_r;
}


//______________________________________________________________________________
double photonNucleusCrossSection::nepoint(double Egamma, double bmin)
{
  //     >> Function for the spectrum of virtual photons,
  //     >> dn/dEgamma, for a point charge q=Ze sweeping
  //     >> past the origin with velocity gamma
  //     >> (=1/SQRT(1-(V/c)**2)) integrated over impact
  //     >> parameter from bmin to infinity
  //     >> See Jackson eq15.54 Classical Electrodynamics
  //     >> Declare Local Variables
  double beta,X,C1,bracket,nepoint_r;
  
  beta = sqrt(1.-(1./(_sigmaGamma_em*_sigmaGamma_em)));
  X = (bmin*Egamma)/(beta*_sigmaGamma_em*starlightConstants::hbarc);
  
  bracket = -0.5*beta*beta*X*X*(bessel::dbesk1(X)*bessel::dbesk1(X)
				-bessel::dbesk0(X)*bessel::dbesk0(X));

  bracket = bracket+X*bessel::dbesk0(X)*bessel::dbesk1(X);
  
  C1=(2.*double((_bbs.getBeam1().getZin())*(_bbs.getBeam1().getZin()))*
      starlightConstants::alpha)/starlightConstants::pi;
  
  //Looks like this is only used in photon flux for the case of pp collisions..
  //might be able to remove the Zs.
  nepoint_r = C1*(1./beta)*(1./beta)*(1./Egamma)*bracket;
  
  return nepoint_r;
  
}


//______________________________________________________________________________
double photonNucleusCrossSection::sigmagp(double Wgp)
{
  //     >> Function for the gamma-proton --> VectorMeson
  //     >> cross section. Wgp is the gamma-proton CM energy.
  //     >> Unit for cross section: fm**2
  
  double sigmagp_r=0.;
  
  switch(_sigmaPID)
    { 
    case starlightConstants::RHO:
    case starlightConstants::RHOZEUS:
    case starlightConstants::FOURPRONG:
      sigmagp_r=1.E-4*(5.0*exp(0.22*log(Wgp))+26.0*exp(-1.23*log(Wgp)));
      break;
    case starlightConstants::OMEGA:
      sigmagp_r=1.E-4*(0.55*exp(0.22*log(Wgp))+18.0*exp(-1.92*log(Wgp)));
      break;                                                      
    case starlightConstants::PHI:
      sigmagp_r=1.E-4*0.34*exp(0.22*log(Wgp));
      break;
    case starlightConstants::JPSI:
    case starlightConstants::JPSI_ee:
    case starlightConstants::JPSI_mumu:
      sigmagp_r=(1.0-((_channelMass+starlightConstants::mp)*(_channelMass+starlightConstants::mp))/(Wgp*Wgp));
      sigmagp_r*=sigmagp_r;
      sigmagp_r*=1.E-4*0.00406*exp(0.65*log(Wgp));
      // sigmagp_r=1.E-4*0.0015*exp(0.80*log(Wgp));
      break;
    case starlightConstants::JPSI2S:
    case starlightConstants::JPSI2S_ee:
    case starlightConstants::JPSI2S_mumu:
      sigmagp_r=(1.0-((_channelMass+starlightConstants::mp)*(_channelMass+starlightConstants::mp))/(Wgp*Wgp));
      sigmagp_r*=sigmagp_r;
      sigmagp_r*=1.E-4*0.00406*exp(0.65*log(Wgp));
      sigmagp_r*=0.166;  
      //      sigmagp_r=0.166*(1.E-4*0.0015*exp(0.80*log(Wgp)));
      break;
    case starlightConstants::UPSILON:
    case starlightConstants::UPSILON_ee:
    case starlightConstants::UPSILON_mumu:
      //       >> This is W**1.7 dependence from QCD calculations
      sigmagp_r=1.E-10*(0.060)*exp(1.70*log(Wgp));
      break;
    case starlightConstants::UPSILON2S:
    case starlightConstants::UPSILON2S_ee:
    case starlightConstants::UPSILON2S_mumu:
      sigmagp_r=1.E-10*(0.0259)*exp(1.70*log(Wgp));
      break;
    case starlightConstants::UPSILON3S:
    case starlightConstants::UPSILON3S_ee:
    case starlightConstants::UPSILON3S_mumu:
      sigmagp_r=1.E-10*(0.0181)*exp(1.70*log(Wgp));
      break;
    default: cout<< "!!!  ERROR: Unidentified Vector Meson: "<< _sigmaPID <<endl;
    }                                                                  
  return sigmagp_r;
}


//______________________________________________________________________________
double photonNucleusCrossSection::sigma_A(double sig_N)
{                                                         
  // Nuclear Cross Section
  // sig_N,sigma_A in (fm**2) 

  double sum;
  double b,bmax,Pint,arg,sigma_A_r;
  
  int NGAUSS;
  
  double xg[17]=
    {.0,
     .0483076656877383162,.144471961582796493,
     .239287362252137075, .331868602282127650,
     .421351276130635345, .506899908932229390,
     .587715757240762329, .663044266930215201,
     .732182118740289680, .794483795967942407,
     .849367613732569970, .896321155766052124,
     .934906075937739689, .964762255587506430,
     .985611511545268335, .997263861849481564
    };
  
  double ag[17]=
    {.0,
     .0965400885147278006, .0956387200792748594,
     .0938443990808045654, .0911738786957638847,
     .0876520930044038111, .0833119242269467552,
     .0781938957870703065, .0723457941088485062,
     .0658222227763618468, .0586840934785355471,
     .0509980592623761762, .0428358980222266807,
     .0342738629130214331, .0253920653092620595,
     .0162743947309056706, .00701861000947009660
    };
  
  NGAUSS=16;
  
  // CALCULATE P(int) FOR b=0.0 - bmax (fm)
  bmax = 25.0;
  sum  = 0.;
  for(int IB=1;IB<=NGAUSS;IB++){
    
    b = 0.5*bmax*xg[IB]+0.5*bmax;
    
    arg=-sig_N*_bbs.getBeam1().getRho0()*_bbs.getBeam1().thickness(b);
    
    Pint=1.0-exp(arg);
    sum=sum+2.*starlightConstants::pi*b*Pint*ag[IB];
    
    b = 0.5*bmax*(-xg[IB])+0.5*bmax;
    arg=-sig_N*_bbs.getBeam1().getRho0()*_bbs.getBeam1().thickness(b);
    Pint=1.0-exp(arg);
    sum=sum+2.*starlightConstants::pi*b*Pint*ag[IB];

  }

  sum=0.5*bmax*sum;
  
  sigma_A_r=sum;
 
  return sigma_A_r;
}


//______________________________________________________________________________
double photonNucleusCrossSection::getDefaultC()
{
  return _defaultC;
}


//______________________________________________________________________________
double photonNucleusCrossSection::breitWigner(double W, double C)
{
	// use simple fixed-width s-wave Breit-Wigner without coherent backgorund for rho'
	// (PDG '08 eq. 38.56)
	if(_sigmaPID==starlightConstants::FOURPRONG) {
		if (W < 4.01 * starlightConstants::mpi)
			return 0;
		const double termA  = _channelMass * _width;
		const double termA2 = termA * termA;
		const double termB  = W * W - _channelMass * _channelMass;
		return C * _ANORM * _ANORM * termA2 / (termB * termB + termA2);
	}

  // Relativistic Breit-Wigner according to J.D. Jackson,
  // Nuovo Cimento 34, 6692 (1964), with nonresonant term. A is the strength
  // of the resonant term and b the strength of the non-resonant
  // term. C is an overall normalization.

  double ppi=0.,ppi0=0.,GammaPrim,rat;
  double aa,bb,cc;
  
  double nrbw_r;

  // width depends on energy - Jackson Eq. A.2
  // if below threshold, then return 0.  Added 5/3/2001 SRK
  // 0.5% extra added for safety margin
  if( _sigmaPID==starlightConstants::RHO ||_sigmaPID==starlightConstants::RHOZEUS){  
    if (W < 2.01*starlightConstants::mpi){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt( ((W/2.)*(W/2.)) - starlightConstants::mpi*starlightConstants::mpi );
    ppi0=0.358;
  }
  
  // handle phi-->K+K- properly
  if (_sigmaPID  ==  starlightConstants::PHI){
    if (W < 2.*starlightConstants::mK){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt( ((W/2.)*(W/2.))- starlightConstants::mK*starlightConstants::mK);
    ppi0=sqrt( ((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mK*starlightConstants::mK);
  }

  //handle J/Psi-->e+e- properly
  if (_sigmaPID==starlightConstants::JPSI || _sigmaPID==starlightConstants::JPSI2S){
    if(W<2.*starlightConstants::mel){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mel*starlightConstants::mel);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mel*starlightConstants::mel);
  }
  if (_sigmaPID==starlightConstants::JPSI_ee){
    if(W<2.*starlightConstants::mel){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mel*starlightConstants::mel);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mel*starlightConstants::mel);   
  }
  if (_sigmaPID==starlightConstants::JPSI_mumu){
    if(W<2.*starlightConstants::mmu){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mmu*starlightConstants::mmu);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mmu*starlightConstants::mmu); 
  }
  if (_sigmaPID==starlightConstants::JPSI2S_ee){
    if(W<2.*starlightConstants::mel){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mel*starlightConstants::mel);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mel*starlightConstants::mel);   
  }
  if (_sigmaPID==starlightConstants::JPSI2S_mumu){
    if(W<2.*starlightConstants::mmu){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mmu*starlightConstants::mmu);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mmu*starlightConstants::mmu); 
  }

  if(_sigmaPID==starlightConstants::UPSILON || _sigmaPID==starlightConstants::UPSILON2S ||_sigmaPID==starlightConstants::UPSILON3S ){ 
    if (W<2.*starlightConstants::mmu){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mmu*starlightConstants::mmu);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mmu*starlightConstants::mmu);
  }
  
  if(_sigmaPID==starlightConstants::UPSILON_mumu || _sigmaPID==starlightConstants::UPSILON2S_mumu ||_sigmaPID==starlightConstants::UPSILON3S_mumu ){ 
    if (W<2.*starlightConstants::mmu){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mmu*starlightConstants::mmu);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mmu*starlightConstants::mmu);
  }
  
  if(_sigmaPID==starlightConstants::UPSILON_ee || _sigmaPID==starlightConstants::UPSILON2S_ee ||_sigmaPID==starlightConstants::UPSILON3S_ee ){ 
    if (W<2.*starlightConstants::mel){
      nrbw_r=0.;
      return nrbw_r;
    }
    ppi=sqrt(((W/2.)*(W/2.))-starlightConstants::mel*starlightConstants::mel);
    ppi0=sqrt(((_channelMass/2.)*(_channelMass/2.))-starlightConstants::mel*starlightConstants::mel);
  }
  
  if(ppi==0.&&ppi0==0.) 
    cout<<"Improper Gammaacrosssection::breitwigner, ppi&ppi0=0."<<endl;
  
  rat=ppi/ppi0;
  GammaPrim=_width*(_channelMass/W)*rat*rat*rat;
  
  aa=_ANORM*sqrt(GammaPrim*_channelMass*W);
  bb=W*W-_channelMass*_channelMass;
  cc=_channelMass*GammaPrim;
  
  // First real part squared 
  nrbw_r = (( (aa*bb)/(bb*bb+cc*cc) + _BNORM)*( (aa*bb)/(bb*bb+cc*cc) + _BNORM));
  
  // Then imaginary part squared 
  nrbw_r = nrbw_r + (( (aa*cc)/(bb*bb+cc*cc) )*( (aa*cc)/(bb*bb+cc*cc) ));

  //  Alternative, a simple, no-background BW, following J. Breitweg et al.
  //  Eq. 15 of Eur. Phys. J. C2, 247 (1998).  SRK 11/10/2000
  //      nrbw_r = (_ANORM*_mass*GammaPrim/(bb*bb+cc*cc))**2
  
  nrbw_r = C*nrbw_r;
  
  return nrbw_r;
    
}


//______________________________________________________________________________
double photonNucleusCrossSection::getMaxPhotonEnergy()
{
  return _EgMax;
}


//_____________________/~~WIDE~~_________________________________________________________
wideResonanceCrossSection::wideResonanceCrossSection(inputParameters& input, beamBeamSystem& bbsystem)
  : photonNucleusCrossSection(input, bbsystem)//hrm
{
  _wideWmax=input.getWmax();
  _wideWmin=input.getWmin();
  _wideYmax=input.getYmax();
  _wideYmin=-1.0*_wideYmax;
  _Ep=input.getProtonEnergy();
}


//______________________________________________________________________________
wideResonanceCrossSection::~wideResonanceCrossSection()
{

}


//______________________________________________________________________________
void wideResonanceCrossSection::crossSectionCalculation(double bwnormsave)
{
  //     This subroutine calculates the cross-section assuming a wide
  //     (Breit-Wigner) resonance.

  // double Av,Wgp,cs,cvma;
  double W,dW,dY;
  double y1,y2,y12,ega1,ega2,ega12;
  // double t,tmin,tmax;
  double csgA1,csgA2,csgA12,int_r,dR,rate;
  double dsigdW,dsigdWalt,dndW,tmp;
  double dsigdW2;
  // double ax,bx;
  double Eth;
  int    I,J,NW,NY;
  // int    K,NGAUSS;
                                                                                                                                                      
  // ----------------- !!!!!!!!!!!!!!!!!!!! -----------------------------
                                                                                                                                                      
  double bwnorm =bwnormsave;//used to transfer the bwnorm from the luminosity tables

  // --------------------------------------------------------------------
  //gamma+nucleon threshold.

  Eth=0.5*(((_wideWmin+starlightConstants::mp)*(_wideWmin+starlightConstants::mp)
	    -starlightConstants::mp*starlightConstants::mp)/(_Ep+sqrt(_Ep*_Ep-starlightConstants::mp*starlightConstants::mp)));
                                                                                                                                                      
  NW   = 100;
  dW   = (_wideWmax-_wideWmin)/double(NW);
  
  NY   =  1200;
  dY   = (_wideYmax-_wideYmin)/double(NY);
  
  if (getBNORM()  ==  0.){
    cout<<" Using Breit-Wigner Resonance Profile."<<endl;
  }
  else{
    cout<<" Using Breit-Wigner plus direct pi+pi- profile."<<endl;
  }
  
  cout<<" Integrating over W from "<<_wideWmin<<" to "<<_wideWmax<<endl;
                                                                                                                                                      
  int_r=0.;
 
  for(I=0;I<=NW-1;I++){
    
    W = _wideWmin + double(I)*dW + 0.5*dW;
    
    tmp = 0.0;
    dsigdW=0.0;
    dsigdW2=0.0;
    dsigdWalt=0.0;
    dndW=0.0;
    
    for(J=0;J<=NY-1;J++){
      
      y1  = _wideYmin + double(J)*dY;
      y2  = _wideYmin + double(J+1)*dY;
      y12 = 0.5*(y1+y2);
      
      ega1  = 0.5*W*exp(y1);
      ega2  = 0.5*W*exp(y2);
      ega12 = 0.5*W*exp(y12);
      
      if(ega1 < Eth) continue;
      if(ega2 > getMaxPhotonEnergy()) continue;
      // check it !!
          
      if(J == 0){
	// >> 1st Point (Calculated only the first time)     =====>>>
	//ega1 used.                                                        
	csgA1=getcsgA(ega1,W);
      }
      else{
	csgA1 = csgA2;
      }
          
      //         >> Middle Point                      =====>>>
      csgA12=getcsgA(ega12,W);         

      //         >> Second Point                      =====>>>
      csgA2=getcsgA(ega2,W);
      
      //>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
      dR  = ega1*photonFlux(ega1)*csgA1;
      dR  = dR + 4.*ega12*photonFlux(ega12)*csgA12;
      dR  = dR + ega2*photonFlux(ega2)*csgA2;
      tmp = tmp+2.*dR*(dY/6.);
      dR  = dR*(dY/6.)*breitWigner(W,bwnorm)*dW;
      
      //For identical beams, we double.  Either may emit photon/scatter
      //For large differences in Z, we approx, that only beam1 emits photon
      //and beam2 scatters, eg d-Au where beam1=au and beam2=d
      if(getbbs().getBeam1().getAin()==getbbs().getBeam2().getAin()){
	dR  = 2.*dR;
      }
      int_r = int_r+dR;  
    }
  }
                                                                                                                                                      
  rate=getLum()*int_r;
  
  cout<<" Cross section (mb): "<<10.*int_r<<endl;
  cout<<" Production rate   : "<<rate<<" Hz"<<endl;
}


//______________________________________________________________________________
/////////////~~NARROW~~////////
narrowResonanceCrossSection::narrowResonanceCrossSection(inputParameters& input,beamBeamSystem& bbsystem)
  :photonNucleusCrossSection(input,bbsystem)
{
  _narrowYmax=input.getYmax();
  _narrowYmin=-1.0*_narrowYmax;
  _narrowNumY=input.getnumy();
  _Ep=input.getProtonEnergy();	
}


//______________________________________________________________________________
narrowResonanceCrossSection::~narrowResonanceCrossSection()
{ }


//______________________________________________________________________________
void narrowResonanceCrossSection::crossSectionCalculation(double)  // _bwnormsave (unused)
{
  // This subroutine calculates the vector meson cross section assuming
  // a narrow resonance.  For reference, see STAR Note 386.
  
  // double Av,Wgp,cs,cvma;
  double W,dY;
  double y1,y2,y12,ega1,ega2,ega12;
  // double t,tmin,tmax;
  double csgA1,csgA2,csgA12,int_r,dR,rate;
  double tmp;
  // double ax,bx;
  double Eth;
  int    J,NY;
  // int    K,NGAUSS;
  
  NY   =  _narrowNumY;
  dY   = (_narrowYmax-_narrowYmin)/double(NY);
  
  cout<<" Using Narrow Resonance ..."<<endl;
  
  W = getChannelMass();
  Eth=0.5*(((W+starlightConstants::mp)*(W+starlightConstants::mp)-
	    starlightConstants::mp*starlightConstants::mp)/(_Ep+sqrt(_Ep*_Ep-starlightConstants::mp*starlightConstants::mp)));
  
  cout<<" gamma+nucleon  Threshold: "<<Eth<<endl;
  int_r=0.;
  
  tmp = 0.0;
  
  for(J=0;J<=(NY-1);J++){
    
    y1  = _narrowYmin + double(J)*dY;
    y2  = _narrowYmin + double(J+1)*dY;
    y12 = 0.5*(y1+y2);
    
    ega1  = 0.5*W*exp(y1);
    ega2  = 0.5*W*exp(y2);
    ega12 = 0.5*W*exp(y12);
    
    if(ega1 < Eth)   
      continue;
    if(ega2 > getMaxPhotonEnergy()) 
      continue;

    csgA1=getcsgA(ega1,W);
    
    // Middle Point                      =====>>>
    csgA12=getcsgA(ega12,W); 

    // Second Point                      =====>>>
    csgA2=getcsgA(ega2,W);
    
    // Sum the contribution for this W,Y.
    dR  = ega1*photonFlux(ega1)*csgA1;
    dR  = dR + 4.*ega12*photonFlux(ega12)*csgA12;
    dR  = dR + ega2*photonFlux(ega2)*csgA2;
    tmp = tmp+2.*dR*(dY/6.);
    dR  = dR*(dY/6.);

    // The 2 accounts for the 2 beams
    //Not dAu
    //Double for identical beams.  
    //If the beams are different in Z by a large amount
    //eg dAu, then only beam1(Au) emits the photon and beam2(d) scatters it
    
    if(getbbs().getBeam1().getAin()==getbbs().getBeam2().getAin()){
      dR  = 2.*dR;
    }
    int_r = int_r+dR;
  }
  rate=getLum()*int_r;
  cout<<" Cross section (mb): " <<10.*int_r<<endl;
  cout<<" Production rate   : "<<rate<<" Hz"<<endl;
}
