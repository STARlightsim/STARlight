// gammaacrosssection.cpp
/*
 * $Id: gammaacrosssection.cpp,v 1.0 2010/07/04   $
 *
 * /author Joseph Butterwoth 
 *
 * $Log: $
 *
 *
 *Converted all nuclear formfactors to come from the scattering nucleus(nucleus2)
 *Added incoherent condition to the cross-section that follows a similar approach as pp
 *Could not figure out the scaling exactly for incoherent(possibly units) so divided by 
 *1E-4 and there is a incoherence factor that can be selected in the input file, 
 *starlight.in--JWB
 *Also, it has not been implemented for interference.
 *
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
 */
#include <iostream>
#include <fstream>

using namespace std;

#include <math.h>
#include "starlightconstants.h"
#include "gammaacrosssection.h"
#include "bessel.h"
//______________________________________________________________________________
Gammaacrosssection::Gammaacrosssection (Inputparameters& input,
                                       Beambeamsystem& bbsystem):bbs(bbsystem)
{
	SigmaProtonEnergy=input.getProtonEnergy();
	SigmaGamma_em=input.getgamma_em();
	SigmaPID=input.getpidtest();
	SigmaBreakup=input.getbreakupmode();
	SigmaCoherence=input.getincoherentorcoherent();
	SigmaCoherenceFactor=input.getincoherentfactor();
	SigmaNucleus=bbs.getBeam2().getAin();

	switch(bbs.getBeam1().getZin())
        {
                case 79://Au
                        lum=2.0;
                        break;
                case 53://I
                        lum=27.;
                        break;
                case 49://Indium,uses same as Iodine
                        lum-27.;
                        break;
                case 29://Cu
                        lum=95.;
                        break;
                case 14://Si
                        lum=440.;
                        break;
                case 8://O
                        lum=980.;
                        break;
                case 82://Pb
                        lum=1.;
                        break;
                case 20://Ca
                        lum=2000.;
                        break;
                case 1://proton
                        lum=1.E8;
                        break;
                default:
        cout <<"Warning:Luminosity not defined.Gammaacrosssection::getlum"<<endl;
	}
	switch(SigmaPID)
        {
                case StarlightConstants::RHO:
                        bslope= 11.0;
			f2o4pi= 2.02;
			ANORM= -2.75;
			BNORM= 0.0;
			defaultC=     1.0;
			channelmass=0.7685;
			width=0.1507;
                        break;
                case StarlightConstants::RHOZEUS:
                        bslope=11.0;
			f2o4pi=2.02;
			ANORM=-2.75;
			BNORM=1.84;
			defaultC=1.0;
			channelmass=0.7685;
			width=0.1507;
                        break;
                case StarlightConstants::OMEGA:
                        bslope=10.0;
			f2o4pi=23.13;
			ANORM=-2.75;
			BNORM=0.0;
			defaultC=1.0;
			channelmass=0.78194;
			width=0.00843;
                        break;
                case StarlightConstants::PHI:
                        bslope=7.0;
			f2o4pi=13.71;
			ANORM=-2.75;
			BNORM=0.0;
			defaultC=1.0;
			channelmass=1.019413;
			width=0.00443;
                        break;
                case StarlightConstants::JPSI:
                        bslope=4.0;
			f2o4pi=10.45;
			ANORM=-2.75;//Artificial Breit-Wigner parameters--no direct pions
			BNORM=0.0;
			defaultC=1.0;
			channelmass=3.09692;//JN 3.09688
			width=0.000091;//JN 0.000087
                        break;
                case StarlightConstants::JPSI2S:
                        bslope=4.3;
			f2o4pi=26.39;
			ANORM=-2.75;//Artificial
			BNORM=0.0;
			defaultC=1.0;
			channelmass=3.686093;
			width=0.000337;
                        break;
                case StarlightConstants::UPSILON:
                        bslope=4.0;
			f2o4pi=125.37;
			ANORM=-2.75;//Artificial
			BNORM=0.0;
			defaultC=1.0;
			channelmass=9.46030;
			width=0.00005402;
                        break;
                case StarlightConstants::UPSILON2S:
                        bslope=4.0;
			f2o4pi=290.84;
			ANORM=-2.75;
			BNORM=0.0;
			defaultC=1.0;
			channelmass=10.02326;
			width=0.00003198;
                        break;
                case StarlightConstants::UPSILON3S:
                        bslope=4.0;
			f2o4pi=415.10;
			ANORM=-2.75;
			BNORM=0.0;
			defaultC=1.0;
			channelmass=10.3552;
			width=0.00002032;
                        break;
                default:
   cout <<"No sigma constants parameterized for pid: "<<SigmaPID
   <<" GammaAcrosssection"<<endl;
        }

	EgMax= 4.*SigmaGamma_em*StarlightConstants::hbarc/bbs.getBeam1().RNuc(); 
	//Max photon energy( for VM only, in GeV, lab frame, use beam energy
    //, nuclear size cutoff)

}
//______________________________________________________________________________
Gammaacrosssection::~Gammaacrosssection()
{

}
//______________________________________________________________________________
Beambeamsystem Gammaacrosssection::getbbs()
{
return bbs;
}
//______________________________________________________________________________
double Gammaacrosssection::getBNORM()
{
return BNORM;
}
//______________________________________________________________________________
double Gammaacrosssection::getlum()
{
return lum;
}
//______________________________________________________________________________
double Gammaacrosssection::getf2o4pi()
{
return f2o4pi;
}
//______________________________________________________________________________
double Gammaacrosssection::getchannelmass()
{
return channelmass;
}
//______________________________________________________________________________
double Gammaacrosssection::getbslope()
{
return bslope;
}
//______________________________________________________________________________
void Gammaacrosssection::crosssectioncalculation(double bwnormsave)
{
cout <<"Neither narrow/wide resonance cross-section calculation.--Derived"<<endl;
}
//______________________________________________________________________________
double Gammaacrosssection::getcsgA(double Egamma,double W)
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
  Wgp=sqrt(2.*Egamma*(SigmaProtonEnergy+sqrt(SigmaProtonEnergy*SigmaProtonEnergy-
					     StarlightConstants::mp*StarlightConstants::mp))+StarlightConstants::mp*StarlightConstants::mp);
	
  //Used for d-A and A-A
  tmin   = (W*W/(4.*Egamma*SigmaGamma_em) )*(W*W/(4.*Egamma*SigmaGamma_em) );
  
  //how to make this more generic???
  if(bbs.getBeam1().getAin()==1&&bbs.getBeam2().getAin()==1)
    {  //Proton-proton, no scaling needed.
      csgA = sigmagp(Wgp);
    }
  else if(bbs.getBeam2().getZin()==1&&bbs.getBeam2().getAin()==2)
    {  //Deuteron-A interaction
      Av = bslope*sigmagp(Wgp);
      
      tmax   = tmin + 0.64;   //0.64
      ax     = 0.5*(tmax-tmin);
      bx     = 0.5*(tmax+tmin);
      csgA   = 0.;
      
      for( int k=1;k<NGAUSS;k++){                                                                                                                   
	t    = ax*xg[k]+bx;
	//We use beam2 here since the input stores the deuteron as nucleus 2
	//and nucleus 2 is the pomeron field source
	//Also this is the way sergey formatted the formfactor.
	csgA = csgA + ag[k]*bbs.getBeam2().formfactor(t); 
	t    = ax*(-xg[k])+bx;
	csgA = csgA + ag[k]*bbs.getBeam2().formfactor(t);
      }
      csgA = 0.5*(tmax-tmin)*csgA;
      csgA = Av*csgA;
    }
        else if(SigmaCoherence==0&&(!(bbs.getBeam2().getZin()==1&&bbs.getBeam2().getAin()==2)))
    {//For incoherent AA interactions , since incoherent treating it as gamma-p
    //      Calculate the differential V.M.+proton cross section
    csgA = 1.E-4*SigmaCoherenceFactor*SigmaNucleus*bslope*sigmagp(Wgp);//artifical 1E-3 to scale down sigma
    //cout<<"sigma for inco AA: "<<cs<<endl;
    //Calculating int |F(t)| dt
    //Using proton formfactor for this case
    //Note the coherence scaling factor being intergrated with the F(t)
    //Should it just be F(t)^2?   
    //Pay attention to the way the formfactor is implemented in Nucleus class
    //Also, notice the tmin value.  It starts higher for dAu, should we proceed
    //in a similar fashion
    //Why don't we use formfactors for pp? Is it because it is incorporated in 
    //the gamma-p fits done for dsigma/dt? Yes?


      /*  tmax   = tmin + 0.25;
        ax     = 0.5*(tmax-tmin);
        bx     = 0.5*(tmax+tmin);
        csgA   = 0.;

        for( int k=1;k<NGAUSS;k++)
        {                                                                                                                   
          t    = ax*xg[k]+bx;
              csgA = csgA + ag[k]*SigmaCoherenceFactor*SigmaNucleus*bbs.getBeam2().formfactor(t);
          t    = ax*(-xg[k])+bx;
              csgA = csgA + ag[k]*SigmaCoherenceFactor*SigmaNucleus*bbs.getBeam2().formfactor(t);
          }

                csgA = 0.5*(tmax-tmin)*csgA;
                csgA = Av*csgA;
        */

    }
  else {	//For typical AA interactions.
		//      Calculate V.M.+proton cross section
    cs=sqrt(16.*StarlightConstants::pi*f2o4pi*bslope
	    *StarlightConstants::hbarc*StarlightConstants::hbarc*sigmagp(Wgp)
	    /StarlightConstants::alpha);
    
    
    //  Calculate V.M.+Nucleus cross section
    cvma=sigma_A(cs);                                                                                                                                      
    
    // Calculate Av = dsigma/dt(t=0) Note Units: fm**s/Gev**2
    Av=(StarlightConstants::alpha*cvma*cvma)/
      (16.*StarlightConstants::pi*f2o4pi*StarlightConstants::hbarc*StarlightConstants::hbarc);
    
    tmax   = tmin + 0.25;
    ax     = 0.5*(tmax-tmin);
    bx     = 0.5*(tmax+tmin);
    csgA   = 0.;
    
    for( int k=1;k<NGAUSS;k++){                                                                                                                   
      t    = ax*xg[k]+bx;
      csgA = csgA + ag[k]*bbs.getBeam2().formfactor(t)*bbs.getBeam2().formfactor(t);
      t    = ax*(-xg[k])+bx;
      csgA = csgA + ag[k]*bbs.getBeam2().formfactor(t)*bbs.getBeam2().formfactor(t);
    }
    
    csgA = 0.5*(tmax-tmin)*csgA;
    csgA = Av*csgA;
    
  }
  return csgA;	
}
//______________________________________________________________________________
double Gammaacrosssection::photonflux(double Egamma)
{
  // DOUBLE PRECISION FUNCTION flux(Egamma)
  // This routine gives the photon flux as a function of energy Egamma
  //  It works for arbitrary nuclei and gamma; the first time it is
  //  called, it calculates a lookup table which is used on
  // subsequent calls
  //  x
  // it returns dN_gamma/dE (dimensions 1/E), not dI/dE
  // energies are in GeV, in the lab frame
  //  rewritten 4/25/2001 by SRK
  
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

  RNuc=bbs.getBeam1().RNuc();
  RNuc2=bbs.getBeam2().RNuc();
  // static ->>> dide,lnEMax,lnEmin,dlnE
  static int  Icheck = 0;
  
  //Check first to see if pp (JN0705)
  if( bbs.getBeam1().getZin()==1 && bbs.getBeam2().getAin()==1 ){
    bmin = RNuc+RNuc;
    flux_r = nepoint(Egamma,bmin);
    return flux_r;
  }
  //   first call?  - initialize - calculate photon flux
  
  Icheck=Icheck+1;
  if(Icheck > 1) goto L1000f;
  
  rZ=double(bbs.getBeam1().getZin());
  rA=double(bbs.getBeam1().getAin());
  rZ2=double(bbs.getBeam2().getZin());  //Sergey--dAu
  rA2=double(bbs.getBeam2().getAin());  //Sergey
  
  //  Nuclear breakup is done by PofB
  //  collect number of integration steps here, in one place
  
  nbstep=400;
  nrstep=60;
  nphistep=40;
  
  //  this last one is the number of energy steps
  nstep=100;
  
  // following previous choices, take Emin=10 keV at LHC, Emin = 1 MeV at RHIC
  
  Emin=1.E-5;
  if (SigmaGamma_em < 500) 
    Emin=1.E-3;
  
  //  maximum energy is 12 times the cutoff
  //  25 GeV for gold at RHIC, 650 GeV for lead at LHC
  
  Emax=12.*StarlightConstants::hbarc*SigmaGamma_em/RNuc;
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
    bmax=bmin + 6.*StarlightConstants::hbarc*SigmaGamma_em/energy;
    
    bmult=exp(log(bmax/bmin)/double(nbstep));
    biter=bmin;
    integratedflux=0.;
    
    if (bbs.getBeam2().getZin()==1&&bbs.getBeam1().getAin()==2){
      //This is for deuteron-gold
      Xvar = (RNuc+RNuc2)*energy/(StarlightConstants::hbarc*(SigmaGamma_em));
      
      fluxelement = (2.0/StarlightConstants::pi)*rZ*rZ*StarlightConstants::alpha/
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
	    Xvar=energy*biter/(StarlightConstants::hbarc*SigmaGamma_em);
	    // Here, there is nuclear breakup.  So, we can't use the integrated flux
	    //  However, we can do a single flux calculation, at the center of the
	    //  nucleus
	    
	    // Eq. 41 of Vidovic, Greiner and Soff, Phys.Rev.C47,2308(1993), among other places
	    //  this is the flux per unit area
	    fluxelement  = (rZ*rZ*StarlightConstants::alpha*energy)*
	      (bessel::dbesk1(Xvar))*(bessel::dbesk1(Xvar))/
	      ((StarlightConstants::pi*SigmaGamma_em*StarlightConstants::hbarc)*
	       (StarlightConstants::pi*SigmaGamma_em*StarlightConstants::hbarc));
	    
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
	    deltaphi=StarlightConstants::pi/double(nphistep);
	    phiiter=0.;
            
	    for( int jphi=1;jphi<= nphistep;jphi++){
	      phiiter=(double(jphi)-0.5)*deltaphi;
	      //  dist is the distance from the center of the emitting nucleus to the point in question
	      dist=sqrt((biter+riter*cos(phiiter))*(biter+riter*
						    cos(phiiter))+(riter*sin(phiiter))*(riter*sin(phiiter)));
	      
	      Xvar=energy*dist/(StarlightConstants::hbarc*SigmaGamma_em);				
	      
	      flux_r = (rZ*rZ*StarlightConstants::alpha*energy)*
		(bessel::dbesk1(Xvar)*bessel::dbesk1(Xvar))/
		((StarlightConstants::pi*SigmaGamma_em*StarlightConstants::hbarc)*
		 (StarlightConstants::pi*SigmaGamma_em*StarlightConstants::hbarc));
	      
	      //  The surface  element is 2.* delta phi* r * delta r
	      //  The '2' is because the phi integral only goes from 0 to pi
	      fluxelement=fluxelement+flux_r*2.*deltaphi*riter*deltar;
	      //  end phi and r integrations
	    }//for(jphi)
	  }//for(jr)
	  //  average fluxelement over the nuclear surface
	  fluxelement=fluxelement/(StarlightConstants::pi*RNuc*RNuc);
	}//else
	//  multiply by volume element to get total flux in the volume element
	fluxelement=fluxelement*2.*StarlightConstants::pi*biter*(biter-bold);
	//  modulate by the probability of nuclear breakup as f(biter)
	if (SigmaBreakup > 1){
	  fluxelement=fluxelement*bbs.probabilityofbreakup(biter);
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
double Gammaacrosssection::nepoint(double Egamma, double bmin)
{
  
  //     >> Function for the spectrum of virtual photons,
  //     >> dn/dEgamma, for a point charge q=Ze sweeping
  //     >> past the origin with velocity gamma
  //     >> (=1/SQRT(1-(V/c)**2)) integrated over impact
  //     >> parameter from bmin to infinity
  //     >> See Jackson eq15.54 Classical Electrodynamics
  //     >> Declare Local Variables
  double beta,X,C1,bracket,nepoint_r;
  
  beta = sqrt(1.-(1./(SigmaGamma_em*SigmaGamma_em)));
  X = (bmin*Egamma)/(beta*SigmaGamma_em*StarlightConstants::hbarc);
  
  bracket = -0.5*beta*beta*X*X*(bessel::dbesk1(X)*bessel::dbesk1(X)
				-bessel::dbesk0(X)*bessel::dbesk0(X));

  bracket = bracket+X*bessel::dbesk0(X)*bessel::dbesk1(X);
  
  C1=(2.*double((bbs.getBeam1().getZin())*(bbs.getBeam1().getZin()))*
      StarlightConstants::alpha)/StarlightConstants::pi;
  
  //Looks like this is only used in photonflux for the case of pp collisions..
  //might be able to remove the Zs.
  nepoint_r = C1*(1./beta)*(1./beta)*(1./Egamma)*bracket;
  
  return nepoint_r;
  
}
//______________________________________________________________________________
double Gammaacrosssection::sigmagp(double Wgp)
{
  //     >> Function for the gamma-proton --> VectorMeson
  //     >> cross section. Wgp is the gamma-proton CM energy.
  //     >> Unit for cross section: fm**2
  
  double sigmagp_r=0.;
  
  switch(SigmaPID)
    { 
    case StarlightConstants::RHO:
    case StarlightConstants::RHOZEUS:
      sigmagp_r=1.E-4*(5.0*exp(0.22*log(Wgp))+26.0*exp(-1.23*log(Wgp)));
      break;
    case StarlightConstants::OMEGA:
      sigmagp_r=1.E-4*(0.55*exp(0.22*log(Wgp))+18.0*exp(-1.92*log(Wgp)));
      break;                                                      
    case StarlightConstants::PHI:
      sigmagp_r=1.E-4*0.34*exp(0.22*log(Wgp));
      break;
    case StarlightConstants::JPSI:
      sigmagp_r=(1.0-((channelmass+StarlightConstants::mp)*(channelmass+StarlightConstants::mp))/(Wgp*Wgp));
      sigmagp_r*=sigmagp_r;
      sigmagp_r*=1.E-4*0.00406*exp(0.65*log(Wgp));
      // sigmagp_r=1.E-4*0.0015*exp(0.80*log(Wgp));
      break;
    case StarlightConstants::JPSI2S:
      sigmagp_r=(1.0-((channelmass+StarlightConstants::mp)*(channelmass+StarlightConstants::mp))/(Wgp*Wgp));
      sigmagp_r*=sigmagp_r;
      sigmagp_r*=1.E-4*0.00406*exp(0.65*log(Wgp));
      sigmagp_r*=0.166;  
      //      sigmagp_r=0.166*(1.E-4*0.0015*exp(0.80*log(Wgp)));
      break;
    case StarlightConstants::UPSILON:
      //       >> This is W**1.7 dependence from QCD calculations
      sigmagp_r=1.E-10*(0.060)*exp(1.70*log(Wgp));
      break;
    case StarlightConstants::UPSILON2S:
      sigmagp_r=1.E-10*(0.0259)*exp(1.70*log(Wgp));
      break;
    case StarlightConstants::UPSILON3S:
      sigmagp_r=1.E-10*(0.0181)*exp(1.70*log(Wgp));
      break;
    default: cout<< "!!!  ERROR: Unidentified Vector Meson: "<< SigmaPID <<endl;
    }                                                                  
  return sigmagp_r;

}
//______________________________________________________________________________
double Gammaacrosssection::sigma_A(double sig_N)
{                                                                                                                                                        
  //     >> Nuclear Cross Section
  //     >> sig_N,sigma_A in (fm**2)
                                                                                                                                                        
  
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
  
  bmax = 25.0;
  sum  = 0.;
  //     >> CALCULATE P(int) FOR b=0.0 - bmax (fm)
  
  
  for(int IB=1;IB<=NGAUSS;IB++){
    
    b = 0.5*bmax*xg[IB]+0.5*bmax;
    
    arg=-sig_N*bbs.getBeam1().getrho0()*bbs.getBeam1().Thickness(b);
    
    Pint=1.0-exp(arg);
    sum=sum+2.*StarlightConstants::pi*b*Pint*ag[IB];
    
    b = 0.5*bmax*(-xg[IB])+0.5*bmax;
    arg=-sig_N*bbs.getBeam1().getrho0()*bbs.getBeam1().Thickness(b);
    Pint=1.0-exp(arg);
    sum=sum+2.*StarlightConstants::pi*b*Pint*ag[IB];
  }
                                                                                                                                                        
  sum=0.5*bmax*sum;
  
  sigma_A_r=sum;
 
  return sigma_A_r;
      
}
//______________________________________________________________________________
double Gammaacrosssection::getdefaultC()
{
  return defaultC;
}
//______________________________________________________________________________
double Gammaacrosssection::breitwigner(double W,double C)
{
  
  
  //Relativistic Breit-Wigner according to J.D. Jackson,
  //    Nuovo Cimento 34, 6692 (1964), with nonresonant term. A is the strength
  //    of the resonant term and b the strength of the non-resonant
  //     term. C is an overall normalization.
  //I think the page is 1644, you can find a copy on Jackson's page on lbl.gov?
  
  
                                                                                                                                               
  double ppi=0.,ppi0=0.,GammaPrim,rat;
  double aa,bb,cc;
  
  // return value
  
  double nrbw_r;
  
                                                                                                                                               
  // width depends on energy - Jackson Eq. A.2
  // if below threshold, then return 0.  Added 5/3/2001 SRK
  // 0.5% extra added for safety margin
  if( SigmaPID==StarlightConstants::RHO ||SigmaPID==StarlightConstants::RHOZEUS){                                                                         
 
    if (W < 2.01*StarlightConstants::mpi){
      nrbw_r=0.;
      return nrbw_r;
    }
    
    ppi=sqrt( ((W/2.)*(W/2.)) - StarlightConstants::mpi*StarlightConstants::mpi );
    ppi0=0.358;
  }                                                                                                                                               
  
  // handle phi-->K+K- properly
  if (SigmaPID  ==  StarlightConstants::PHI){
      if (W < 2.*StarlightConstants::mK){
	nrbw_r=0.;
	return nrbw_r;
      }
      ppi=sqrt( ((W/2.)*(W/2.))- StarlightConstants::mK*StarlightConstants::mK);
      ppi0=sqrt( ((channelmass/2.)*(channelmass/2.))-StarlightConstants::mK*StarlightConstants::mK);
  }
  //handle J/Psi-->e+e- properly
  if (SigmaPID==StarlightConstants::JPSI || SigmaPID==StarlightConstants::JPSI2S){
    
    if(W<2.*StarlightConstants::mel){
      nrbw_r=0.;
      return nrbw_r;
    }
    
    ppi=sqrt(((W/2.)*(W/2.))-StarlightConstants::mel*StarlightConstants::mel);
    ppi0=sqrt(((channelmass/2.)*(channelmass/2.))-StarlightConstants::mel*StarlightConstants::mel);
    
  }
  if(SigmaPID==StarlightConstants::UPSILON || SigmaPID==StarlightConstants::UPSILON2S ||SigmaPID==StarlightConstants::UPSILON3S ){
    
    if (W<2.*StarlightConstants::mmu){
      nrbw_r=0.;
      return nrbw_r;
    }
	

    ppi=sqrt(((W/2.)*(W/2.))-StarlightConstants::mmu*StarlightConstants::mmu);
    ppi0=sqrt(((channelmass/2.)*(channelmass/2.))-StarlightConstants::mmu*StarlightConstants::mmu);
    
  }
  
  if(ppi==0.&&ppi0==0.) 
    cout<<"Improper Gammaacrosssection::breitwigner, ppi&ppi0=0."<<endl;
  
  rat=ppi/ppi0;
  GammaPrim=width*(channelmass/W)*rat*rat*rat;
  
  aa=ANORM*sqrt(GammaPrim*channelmass*W);
  bb=W*W-channelmass*channelmass;
  cc=channelmass*GammaPrim;
                                                                                                                                               
  //  real part^2
  
  nrbw_r = (( (aa*bb)/(bb*bb+cc*cc) + BNORM)*( (aa*bb)/(bb*bb+cc*cc) + BNORM));
  
  // imaginary part^2
                                                                                                                                               
  nrbw_r = nrbw_r + (( (aa*cc)/(bb*bb+cc*cc) )*( (aa*cc)/(bb*bb+cc*cc) ));
                                                                                                                                               
                                                                                                                                               
  //  Alternative, a simple, no-background BW, following J. Breitweg et al.
  //  Eq. 15 of Eur. Phys. J. C2, 247 (1998).  SRK 11/10/2000
  //      nrbw_r = (ANORM*mass*GammaPrim/(bb*bb+cc*cc))**2
  
  nrbw_r = C*nrbw_r;
                                                                                                                                               
  return nrbw_r;
    
}
//______________________________________________________________________________
double Gammaacrosssection::getMaxPhotonEnergy()
{
	return EgMax;
}
//_____________________/~~WIDE~~_________________________________________________________
Wideresonancesigma::Wideresonancesigma(Inputparameters& input,Beambeamsystem& bbsystem)
                                       :Gammaacrosssection(input,bbsystem)//hrm
{
	WideWmax=input.getWmax();
	WideWmin=input.getWmin();
	WideYmax=input.getYmax();
	WideYmin=-1.0*WideYmax;
	Ep=input.getProtonEnergy();

}
//______________________________________________________________________________
Wideresonancesigma::~Wideresonancesigma()
{

}
//______________________________________________________________________________
void Wideresonancesigma::crosssectioncalculation(double bwnormsave)
{
  //     This subroutine calculates the cross-section assuming a wide
  //     (Breit-Wigner) resonance.
                                                                                                                                                      

  double Av,Wgp,cs,cvma;
  double W,dW,dY;
  double y1,y2,y12,ega1,ega2,ega12;
  double t,tmin,tmax;
  double csgA1,csgA2,csgA12,int_r,dR,rate;
  double dsigdW,dsigdWalt,dndW,tmp;
  double dsigdW2;
  double ax,bx, Eth;
  int    I,J,K,NW,NY,NGAUSS;
                                                                                                                                                      
  // ----------------- !!!!!!!!!!!!!!!!!!!! -----------------------------
                                                                                                                                                      
  double bwnorm =bwnormsave;//used to transfer the bwnorm from the luminosity tables

  // --------------------------------------------------------------------
	//gamma+nucleon threshold.

  Eth=0.5*(((WideWmin+StarlightConstants::mp)*(WideWmin+StarlightConstants::mp)
	    -StarlightConstants::mp*StarlightConstants::mp)/(Ep+sqrt(Ep*Ep-StarlightConstants::mp*StarlightConstants::mp)));
                                                                                                                                                      
  NW   = 100;
  dW   = (WideWmax-WideWmin)/double(NW);
  
  NY   =  1200;
  dY   = (WideYmax-WideYmin)/double(NY);
  
  if (getBNORM()  ==  0.){
    cout<<" Using Breit-Wigner Resonance Profile."<<endl;
  }
  else{
    cout<<" Using Breit-Wigner plus direct pi+pi- profile."<<endl;
  }
  
  cout<<" Integrating over W from "<<WideWmin<<" to "<<WideWmax<<endl;
                                                                                                                                                      
  int_r=0.;
 
  for(I=0;I<=NW-1;I++){
    
    W = WideWmin + double(I)*dW + 0.5*dW;
    
    tmp = 0.0;
    dsigdW=0.0;
    dsigdW2=0.0;
    dsigdWalt=0.0;
    dndW=0.0;
    
    for(J=0;J<=NY-1;J++){
      
      y1  = WideYmin + double(J)*dY;
      y2  = WideYmin + double(J+1)*dY;
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
      dR  = ega1*photonflux(ega1)*csgA1;
      dR  = dR + 4.*ega12*photonflux(ega12)*csgA12;
      dR  = dR + ega2*photonflux(ega2)*csgA2;
      tmp = tmp+2.*dR*(dY/6.);
      dR  = dR*(dY/6.)*breitwigner(W,bwnorm)*dW;
      
      //For identical beams, we double.  Either may emit photon/scatter
      //For large differences in Z, we approx, that only beam1 emits photon
      //and beam2 scatters, eg d-Au where beam1=au and beam2=d
      if(getbbs().getBeam1().getAin()==getbbs().getBeam2().getAin()){
	dR  = 2.*dR;
      }
      int_r = int_r+dR;  
    }
  }
                                                                                                                                                      
  rate=getlum()*int_r;
  
  cout<<" Cross section (mb): "<<10.*int_r<<endl;
  cout<<" Production rate   : "<<rate<<" Hz"<<endl;
                                                                                                                                                     
}
//______________________________________________________________________________
/////////////~~NARROW~~////////
Narrowresonancesigma::Narrowresonancesigma(Inputparameters& input,Beambeamsystem& bbsystem)
  :Gammaacrosssection(input,bbsystem)
{
	
        NarrowYmax=input.getYmax();
        NarrowYmin=-1.0*NarrowYmax;
	NarrowNumY=input.getnumy();
        Ep=input.getProtonEnergy();	

}
//______________________________________________________________________________
Narrowresonancesigma::~Narrowresonancesigma()
{

}
//______________________________________________________________________________
void Narrowresonancesigma::crosssectioncalculation(double bwnormsave)
{
                                                                                                                                     
  
  //     This subroutine calculates the vector meson cross section assuming
  //     a narrow resonance.  For reference, see STAR Note 386.
  
  double Av,Wgp,cs,cvma;
  double W,dY;
  double y1,y2,y12,ega1,ega2,ega12;
  double t,tmin,tmax;
  double csgA1,csgA2,csgA12,int_r,dR,rate;
  double tmp;
  double ax,bx, Eth;
  int          J,K,NY,NGAUSS;
  
  NY   =  NarrowNumY;
  dY   = (NarrowYmax-NarrowYmin)/double(NY);
  
  cout<<" Using Narrow Resonance ..."<<endl;
  
  W = getchannelmass();
  Eth=0.5*(((W+StarlightConstants::mp)*(W+StarlightConstants::mp)-
	    StarlightConstants::mp*StarlightConstants::mp)/(Ep+sqrt(Ep*Ep-StarlightConstants::mp*StarlightConstants::mp)));
  
  cout<<" gamma+nucleon  Threshold: "<<Eth<<endl;
  int_r=0.;
  
  tmp = 0.0;
  
  for(J=0;J<=(NY-1);J++){
    
    y1  = NarrowYmin + double(J)*dY;
    y2  = NarrowYmin + double(J+1)*dY;
    y12 = 0.5*(y1+y2);
    
    ega1  = 0.5*W*exp(y1);
    ega2  = 0.5*W*exp(y2);
    ega12 = 0.5*W*exp(y12);
    
    if(ega1 < Eth)   
      continue;
    if(ega2 > getMaxPhotonEnergy()) 
      continue;
                                                                                                                                     
    csgA1=getcsgA(ega1,W);
    
    //       >> Middle Point                      =====>>>
    
    csgA12=getcsgA(ega12,W);                                                                                                                          
    //       >> Second Point                      =====>>>
    csgA2=getcsgA(ega2,W);
    
    //       >> Sum the contribution for this W,Y.
    dR  = ega1*photonflux(ega1)*csgA1;
    dR  = dR + 4.*ega12*photonflux(ega12)*csgA12;
    dR  = dR + ega2*photonflux(ega2)*csgA2;
    tmp = tmp+2.*dR*(dY/6.);
    dR  = dR*(dY/6.);
    //       >> The 2 accounts for the 2 beams
    //Not dAu
    //Double for identical beams.  
    //If the beams are different in Z by a large amount
    //eg dAu, then only beam1(Au) emits the photon and beam2(d) scatters it
    
    if(getbbs().getBeam1().getAin()==getbbs().getBeam2().getAin()){
      dR  = 2.*dR;
    }
    int_r = int_r+dR;
  }
  rate=getlum()*int_r;
  cout<<" Cross section (mb): " <<10.*int_r<<endl;
  cout<<" Production rate   : "<<rate<<" Hz"<<endl;
}


