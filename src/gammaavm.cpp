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
//    Added incoherent t2-> pt2 selection.  Following pp selection scheme
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "gammaavm.h"
#include "photonNucleusCrossSection.h"
#include "wideResonanceCrossSection.h"
#include "narrowResonanceCrossSection.h"
#include "incoherentVMCrossSection.h"

using namespace std;


//______________________________________________________________________________
Gammaavectormeson::Gammaavectormeson(beamBeamSystem& bbsystem):eventChannel(bbsystem), _phaseSpaceGen(0)  //:readLuminosity(input),_bbs(bbsystem)
{
	_VMNPT=inputParametersInstance.nmbPtBinsInterference();
	_VMWmax=inputParametersInstance.maxW();
	_VMWmin=inputParametersInstance.minW();
	_VMYmax=inputParametersInstance.maxRapidity();
	_VMYmin=-1.*_VMYmax;
	_VMnumw=inputParametersInstance.nmbWBins();
	_VMnumy=inputParametersInstance.nmbRapidityBins();
	_VMgamma_em=inputParametersInstance.beamLorentzGamma();
	_VMinterferencemode=inputParametersInstance.interferenceEnabled();
	_VMbslope=0.;//Will define in wide/narrow constructor
	_VMpidtest=inputParametersInstance.prodParticleType();
	_VMptmax=inputParametersInstance.maxPtInterference();
	_VMdpt=inputParametersInstance.ptBinWidthInterference();
	_VMCoherence=inputParametersInstance.coherentProduction();
	_VMCoherenceFactor=inputParametersInstance.coherentProduction();
        _ProductionMode=inputParametersInstance.productionMode();

        N1 = 0; N2 = 0; 

	switch(_VMpidtest){
	case starlightConstants::RHO:
	case starlightConstants::RHOZEUS:
		_width=0.1507;
		_mass=0.7685;
		break;
	case starlightConstants::FOURPRONG:
		// create n-body phase-space generator instance
		_phaseSpaceGen = new nBodyPhaseSpaceGen();
		_width = 0.360;
		_mass  = 1.350;
		break;
	case starlightConstants::OMEGA:
		_width=0.00843;
		_mass=0.78194;
		break;
	case starlightConstants::PHI:
		_width=0.00443;
		_mass=1.019413;
		break;
	case starlightConstants::JPSI:
	case starlightConstants::JPSI_ee:
	case starlightConstants::JPSI_mumu:
		_width=0.000091;
		_mass=3.09692;
		break;
	case starlightConstants::JPSI2S:
	case starlightConstants::JPSI2S_ee:
	case starlightConstants::JPSI2S_mumu:
		_width=0.000337;
		_mass=3.686093;
		break;
	case starlightConstants::UPSILON:
	case starlightConstants::UPSILON_ee:
	case starlightConstants::UPSILON_mumu:
		_width=0.00005402;
		_mass=9.46030;
		break;
	case starlightConstants::UPSILON2S:
	case starlightConstants::UPSILON2S_ee:
	case starlightConstants::UPSILON2S_mumu:
		_width=0.00003198;
		_mass=10.02326;
		break;
	case starlightConstants::UPSILON3S:
	case starlightConstants::UPSILON3S_ee:
	case starlightConstants::UPSILON3S_mumu:
		_width=0.00002032;
		_mass=10.3552;
		break;
	default: cout<<"No proper vector meson defined, gammaavectormeson::gammaavectormeson"<<endl;
	}  
}


//______________________________________________________________________________
Gammaavectormeson::~Gammaavectormeson()
{
	if (_phaseSpaceGen)
		delete _phaseSpaceGen;
}


//______________________________________________________________________________
void Gammaavectormeson::pickwy(double &W, double &Y)
{
        double dW, dY, xw,xy,xtest,btest;
	int  IW,IY;
  
	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	dY = (_VMYmax-_VMYmin)/double(_VMnumy);
  
 L201pwy:

	xw = randyInstance.Rndom();// random()/(RAND_MAX+1.0);
	W = _VMWmin + xw*(_VMWmax-_VMWmin);

	if (W < 2 * starlightConstants::pionChargedMass)
		goto L201pwy;
  
	IW = int((W-_VMWmin)/dW); //+ 1;
	xy = randyInstance.Rndom();//random()/(RAND_MAX+1.0);
	Y = _VMYmin + xy*(_VMYmax-_VMYmin);
	IY = int((Y-_VMYmin)/dY); //+ 1;
	xtest = randyInstance.Rndom();//random()/(RAND_MAX+1.0);

	if( xtest > _Farray[IW][IY] )
		goto L201pwy;

	// Determine the target nucleus 
	// For pA this is given, for all other cases use the relative probabilities in _Farray1 and _Farray2 
        if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
           if( _ProductionMode == 2 ){
	     _TargetBeam = 2;
	   } else {
             _TargetBeam = 1;
           }
        } else if(  _bbs.beam1().A() != 1 && _bbs.beam2().A()==1 ){
           if( _ProductionMode == 2 ){
	     _TargetBeam = 1;
	   } else {
             _TargetBeam = 2;
           }
        } else {
          btest = randyInstance.Rndom();
          if ( btest < _Farray1[IW][IY]/_Farray[IW][IY] ){
            _TargetBeam = 2;
            N2++;
          }  else {
            _TargetBeam = 1;
            N1++; 
          }
        }
	// cout<<"N1: "<<N1<<" N2; "<<N2<<endl;
}         


//______________________________________________________________________________                                               
void Gammaavectormeson::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,
                                     double,  // E (unused)
                                     double  W,
                                     double  px0, double  py0, double  pz0,
                                     double& px1, double& py1, double& pz1,
                                     double& px2, double& py2, double& pz2,
                                     int&    iFbadevent)
{
	// This routine decays a particle into two particles of mass mdec,
	// taking spin into account

	double pmag;
	// double anglelep[20001],xtest,ytest=0.,dndtheta;
	double phi,theta,Ecm;
	double betax,betay,betaz;
	double mdec=0.0;
	double E1=0.0,E2=0.0;

	//    set the mass of the daughter particles
	mdec=getDaughterMass(ipid);

	//     calculate the magnitude of the momenta
	if(W < 2*mdec){
		cout<<" ERROR: W="<<W<<endl;
		iFbadevent = 1;
		return;
	}
	pmag = sqrt(W*W/4. - mdec*mdec);
  
	//     pick an orientation, based on the spin
	//      phi has a flat distribution in 2*pi
	phi = randyInstance.Rndom()*2.*starlightConstants::pi;//(random()/(RAND_MAX+1.0))* 2.*pi;
                                                                                                                
	//     find theta, the angle between one of the outgoing particles and
	//    the beamline, in the frame of the two photons

	theta=getTheta(ipid);
 
	//     compute unboosted momenta
	px1 = sin(theta)*cos(phi)*pmag;
	py1 = sin(theta)*sin(phi)*pmag;
	pz1 = cos(theta)*pmag;
	px2 = -px1;
	py2 = -py1;
	pz2 = -pz1;

	Ecm = sqrt(W*W+px0*px0+py0*py0+pz0*pz0);
	E1 = sqrt(mdec*mdec+px1*px1+py1*py1+pz1*pz1);
	E2 = sqrt(mdec*mdec+px2*px2+py2*py2+pz2*pz2);

	betax = -(px0/Ecm);
	betay = -(py0/Ecm);
	betaz = -(pz0/Ecm);

	transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
	transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);

	if(iFbadevent == 1)
		return;

}


//______________________________________________________________________________                                               
// decays a particle into four particles with isotropic angular distribution
bool Gammaavectormeson::fourBodyDecay
(starlightConstants::particleTypeEnum& ipid,
 const double                  ,           // E (unused)
 const double                  W,          // mass of produced particle
 const double*                 p,          // momentum of produced particle; expected to have size 3
 lorentzVector*                decayVecs,  // array of Lorentz vectors of daughter particles; expected to have size 4
 int&                          iFbadevent)
{
	const double parentMass = W;

	// set the mass of the daughter particles
	const double daughterMass = getDaughterMass(ipid);
	if (parentMass < 4 * daughterMass){
		cout << " ERROR: W=" << parentMass << " GeV too small" << endl;
		iFbadevent = 1;
		return false;
	}

	// construct parent four-vector
	const double        parentEnergy = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
	                                        + parentMass * parentMass);
	const lorentzVector parentVec(p[0], p[1], p[2], parentEnergy);

	// setup n-body phase-space generator
	assert(_phaseSpaceGen);
	static bool firstCall = true;
	if (firstCall) {
		const double m[4] = {daughterMass, daughterMass, daughterMass, daughterMass};
		_phaseSpaceGen->setDecay(4, m);
		// estimate maximum phase-space weight
		_phaseSpaceGen->setMaxWeight(1.01 * _phaseSpaceGen->estimateMaxWeight(_VMWmax));
		firstCall = false;
	}

	// generate phase-space event
	if (!_phaseSpaceGen->generateDecayAccepted(parentVec))
		return false;

	// set Lorentzvectors of decay daughters
	for (unsigned int i = 0; i < 4; ++i)
		decayVecs[i] = _phaseSpaceGen->daughter(i);
	return true;
}


//______________________________________________________________________________
double Gammaavectormeson::getDaughterMass(starlightConstants::particleTypeEnum &ipid)
{
	//This will return the daughter particles mass, and the final particles outputed id...
	double ytest=0.,mdec=0.;
  
	switch(_VMpidtest){
	case starlightConstants::RHO:
	case starlightConstants::RHOZEUS:
	case starlightConstants::FOURPRONG:
	case starlightConstants::OMEGA:
		mdec = starlightConstants::pionChargedMass;
		ipid = starlightConstants::PION;
		break;
	case starlightConstants::PHI:
		mdec = starlightConstants::kaonChargedMass;
		ipid = starlightConstants::KAONCHARGE;
		break;
	case starlightConstants::JPSI:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break; 
	case starlightConstants::JPSI2S_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI2S_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break; 

	case starlightConstants::JPSI2S:
	case starlightConstants::UPSILON:
	case starlightConstants::UPSILON2S:
	case starlightConstants::UPSILON3S:
		//  decays 50% to e+/e-, 50% to mu+/mu-
		ytest = randyInstance.Rndom();//random()/(RAND_MAX+1.0);
    
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break;
	case starlightConstants::UPSILON_ee:
	case starlightConstants::UPSILON2S_ee:
	case starlightConstants::UPSILON3S_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break;
	case starlightConstants::UPSILON_mumu:
	case starlightConstants::UPSILON2S_mumu:
	case starlightConstants::UPSILON3S_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;   
		break;
	default: cout<<"No daughtermass defined, gammaavectormeson::getdaughtermass"<<endl;
	}
  
	return mdec;
}


//______________________________________________________________________________
double Gammaavectormeson::getTheta(starlightConstants::particleTypeEnum ipid)
{
	//This depends on the decay angular distribution
	//Valid for rho, phi, omega.
	double theta=0.;
	double xtest=0.;
	double dndtheta=0.;
 L200td:
                                                                                                                                                 
	theta = starlightConstants::pi*randyInstance.Rndom();//random()/(RAND_MAX+1.0);
	xtest = randyInstance.Rndom();//random()/(RAND_MAX+1.0);
	//  Follow distribution for helicity +/-1
	//  Eq. 19 of J. Breitweg et al., Eur. Phys. J. C2, 247 (1998)
	//  SRK 11/14/2000
  
	switch(ipid){
	  
	case starlightConstants::MUON:
	case starlightConstants::ELECTRON:
		//primarily for upsilon/j/psi.  VM->ee/mumu
		dndtheta = sin(theta)*(1.+((cos(theta))*(cos(theta))));
		break;
    
	case starlightConstants::PION:
	case starlightConstants::KAONCHARGE:
		//rhos etc
		dndtheta= sin(theta)*(1.-((cos(theta))*(cos(theta))));
		break;
    
	default: cout<<"No proper theta dependence defined, check gammaavectormeson::gettheta"<<endl;
	}//end of switch
  
	if(xtest > dndtheta)
		goto L200td;
  
	return theta;
  
}


//______________________________________________________________________________
double Gammaavectormeson::getWidth()
{
	return _width;
}


//______________________________________________________________________________
double Gammaavectormeson::getMass()
{
	return _mass;
}


//______________________________________________________________________________
double Gammaavectormeson::getSpin()
{
	return 1.0; //VM spins are the same
}


//______________________________________________________________________________
void Gammaavectormeson::momenta(double W,double Y,double &E,double &px,double &py,double &pz,int &tcheck)
{
	//     This subroutine calculates momentum and energy of vector meson
	//     given W and Y,   without interference.  Subroutine vmpt.f handles
	//     production with interference
 
	double dW,dY;
	double Egam,Epom,tmin,pt1,pt2,phi1,phi2;
	double px1,py1,px2,py2;
	double pt,xt,xtest,ytest;
	//	double photon_spectrum;
	double t2;

	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	dY  = (_VMYmax-_VMYmin)/double(_VMnumy);
  
	//Find Egam,Epom in CM frame
        if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
          if( _ProductionMode == 2 ){
    	    Egam = 0.5*W*exp(Y);
  	    Epom = 0.5*W*exp(-Y);
          }else{
    	    Egam = 0.5*W*exp(-Y);
  	    Epom = 0.5*W*exp(Y);
          }  
        } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
          if( _ProductionMode == 2 ){
  	    Egam = 0.5*W*exp(-Y);
  	    Epom = 0.5*W*exp(Y);
          }else{
    	    Egam = 0.5*W*exp(Y);
  	    Epom = 0.5*W*exp(-Y);
          }  
        } else if ( _TargetBeam == 1 ) {
	  Egam = 0.5*W*exp(-Y);
	  Epom = 0.5*W*exp(Y);
        } else { 
          Egam = 0.5*W*exp(Y);
	  Epom = 0.5*W*exp(-Y);
	}

        // cout<<" Y: "<<Y<<" W: "<<W<<" TargetBeam; "<<_TargetBeam<<" Egam: "<<Egam<<" Epom: "<<Epom<<endl; 
        pt1 = pTgamma(Egam);  
	phi1 = 2.*starlightConstants::pi*randyInstance.Rndom();

	if( (_bbs.beam1().A()==1 && _bbs.beam2().A()==1) || 
            (_ProductionMode == 4) ) {
	    if( (_VMpidtest == starlightConstants::RHO) || (_VMpidtest == starlightConstants::RHOZEUS) || (_VMpidtest == starlightConstants::OMEGA)){
	      // Use dipole form factor for light VM
	      L555vm:
	      xtest = 2.0*randyInstance.Rndom();
              double ttest = xtest*xtest; 
              ytest = randyInstance.Rndom();
              double t0 = 1./2.23; 
              double yprob = xtest*_bbs.beam1().dipoleFormFactor(ttest,t0)*_bbs.beam1().dipoleFormFactor(ttest,t0); 
              if( ytest > yprob ) goto L555vm; 
              t2 = ttest; 
              pt2 = xtest;              
	    }else{
		//Use dsig/dt= exp(-_VMbslope*t) for heavy VM
		xtest = randyInstance.Rndom(); 
		t2 = (-1./_VMbslope)*log(xtest);
		pt2 = sqrt(1.*t2);
	    }
	} else {
		//       >> Check tmin
		tmin = ((Epom/_VMgamma_em)*(Epom/_VMgamma_em));
	
		if(tmin > 0.5){
			cout<<" WARNING: tmin= "<<tmin<<endl;
                        cout<< " Y = "<<Y<<" W = "<<W<<" Epom = "<<Epom<<" gamma = "<<_VMgamma_em<<endl; 
			cout<<" Will pick a new W,Y "<<endl;
			tcheck = 1;
			return;
		}
 L203vm:
		xt = randyInstance.Rndom(); 
                if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
                  if( _ProductionMode == 2 ){
                    pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
                  }else{
                    pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
                  }   
                } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
                  if( _ProductionMode == 2 ){
                    pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
                  }else{
                    pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
                  }  
                } else if (_TargetBeam==1) {
                    pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
                } else {
                    pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
                }

		xtest = randyInstance.Rndom();
		t2 = tmin + pt2*pt2;
    
		if(_bbs.beam2().Z()==1&&_bbs.beam2().A()==2){
			if(1.0 < _bbs.beam2().formFactor(t2)*pt2)  cout <<"POMERON"<<endl;
			if( xtest > _bbs.beam2().formFactor(t2)*pt2) goto L203vm;
		}
		else{

		  double comp=0.0; 
                  if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
                    if( _ProductionMode == 2 ){
                      comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
                    }else{
                      comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
                    }   
                  } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
                    if( _ProductionMode == 2 ){
                      comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
                    }else{
                      comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
                    }  
                  } else if (_TargetBeam==1) {
                    comp = _bbs.beam2().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
                  } else {
                    comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2; 
                  }
	      
                  if( xtest > comp ) goto L203vm;
       		}

	}//else end from pp
	phi2 = 2.*starlightConstants::pi*randyInstance.Rndom();//random()/(RAND_MAX+1.0);

	px1 = pt1*cos(phi1);
	py1 = pt1*sin(phi1);
	px2 = pt2*cos(phi2);
	py2 = pt2*sin(phi2);
        
	// Compute vector sum Pt = Pt1 + Pt2 to find pt for the vector meson
	px = px1 + px2;
	py = py1 + py2;
	pt = sqrt( px*px + py*py );
       
	E  = sqrt(W*W+pt*pt)*cosh(Y);
	pz = sqrt(W*W+pt*pt)*sinh(Y);
 
        // cout<< " Y = "<<Y<<" W = "<<W<<" Egam = "<<Egam<<" gamma = "<<_VMgamma_em<<endl; 

	// Keep this special case for d+A 
	if(_bbs.beam2().Z()==1&&_bbs.beam2().A()==2){
		pz = -pz;
	}

        /*
        else if( (_bbs.beam1().A()==1 && _bbs.beam2().A() != 1) || (_bbs.beam2().A()==1 && _bbs.beam1().A() != 1) ){
	  // Don't switch      
        }
	else{
		if (randyInstance.Rndom() >= 0.5) pz = -pz;
	}
        */
}

//______________________________________________________________________________
double Gammaavectormeson::pTgamma(double E)
{
    // returns on random draw from pp(E) distribution
    double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
    double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
    int satisfy =0;
        
    ereds = (E/_VMgamma_em)*(E/_VMgamma_em);
    //sqrt(3)*E/gamma_em is p_t where the distribution is a maximum
    Cm = sqrt(3.)*E/_VMgamma_em;
    // If E is very small, the drawing of a pT below is extremely slow. 
    // ==> Set pT = sqrt(3.)*E/_VMgamma_em for very small E. 
    // Should have no observable consequences (JN, SRK 11-Sep-2014)
    if( E < 0.0005 )return Cm; 
 
    //the amplitude of the p_t spectrum at the maximum

    if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
      if( _ProductionMode == 2 ){
         singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
      }  
    } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
      if( _ProductionMode == 2 ){
         singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
      }  
    } else if (_TargetBeam == 1) {
      singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
    } else {
      singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
    }

    Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
    //pick a test value pp, and find the amplitude there
    x = randyInstance.Rndom();

    if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
      if( _ProductionMode == 2 ){
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
      }  
    } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
      if( _ProductionMode == 2 ){
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
      }  
    } else if (_TargetBeam == 1) {
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
    } else {
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
    }

    test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

    while(satisfy==0){
	u = randyInstance.Rndom();
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =randyInstance.Rndom();
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 ){
                pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);
              }  
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 ){
                pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
              }  
            } else if (_TargetBeam == 1) {
              pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
              singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
            } else {
              pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
              singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
            }
	    test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
	}
    }

    return pp;
}


//______________________________________________________________________________
void Gammaavectormeson::vmpt(double W,double Y,double &E,double &px,double &py, double &pz,
                             int&) // tcheck (unused)
{
	//    This function calculates momentum and energy of vector meson
	//    given W and Y, including interference.
	//    It gets the pt distribution from a lookup table.
	double dW=0.,dY=0.,yleft=0.,yfract=0.,xpt=0.,pt1=0.,ptfract=0.,pt=0.,pt2=0.,theta=0.;
	int IY=0,j=0;
  
	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	dY  = (_VMYmax-_VMYmin)/double(_VMnumy);
  
	//  Y is already fixed; choose a pt
	//  Follow the approavh in pickwy.f
	//  in  _fptarray(IY,pt) IY=1 corresponds to Y=0, IY=numy/2 corresponds to +y
  
	IY=int(fabs(Y)/dY);//+1;
	if (IY > (_VMnumy/2)-1){
		IY=(_VMnumy/2)-1;
	}
  
	yleft=fabs(Y)-(IY)*dY;
	yfract=yleft*dY;
  
	xpt=randyInstance.Rndom(); //random()/(RAND_MAX+1.0);
        
	for(j=0;j<_VMNPT+1;j++){
		if (xpt < _fptarray[IY][j]) goto L60;
	}
 L60:
  
	//  now do linear interpolation - start with extremes
  	if (j == 0){
		pt1=xpt/_fptarray[IY][j]*_VMdpt/2.;
		goto L80;
	}
	if (j == _VMNPT){
		pt1=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY][j])/(1.-_fptarray[IY][j]);
		goto L80;
	}
  
	//  we're in the middle
  	ptfract=(xpt-_fptarray[IY][j])/(_fptarray[IY][j+1]-_fptarray[IY][j]);
	pt1=(j+1)*_VMdpt+ptfract*_VMdpt;
  
	//  at an extreme in y?
	if (IY == (_VMnumy/2)-1){
		pt=pt1;
		goto L120;
	}
 L80:

	//  interpolate in y repeat for next fractional y bin      
	for(j=0;j<_VMNPT+1;j++){
		if (xpt < _fptarray[IY+1][j]) goto L90;
	}
 L90:
  
	//  now do linear interpolation - start with extremes
	if (j == 0){
		pt2=xpt/_fptarray[IY+1][j]*_VMdpt/2.;
		goto L100;
	}
	if (j == _VMNPT){
		pt2=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY+1][j])/(1.-_fptarray[IY+1][j]);
		goto L100;
	}
  
	//  we're in the middle
	ptfract=(xpt-_fptarray[IY+1][j])/(_fptarray[IY+1][j+1]-_fptarray[IY+1][j]);
	pt2=(j+1)*_VMdpt+ptfract*_VMdpt;
 L100:

	//  now interpolate in y  
	pt=yfract*pt2+(1-yfract)*pt1;

 L120:

	//  we have a pt 
	theta=2.*starlightConstants::pi*randyInstance.Rndom();//(random()/(RAND_MAX+1.0))*2.*pi;
	px=pt*cos(theta);
	py=pt*sin(theta);

	E  = sqrt(W*W+pt*pt)*cosh(Y);
	pz = sqrt(W*W+pt*pt)*sinh(Y);
	//      randomly choose to make pz negative 50% of the time
	if(randyInstance.Rndom()>=0.5) pz = -pz;
}


//______________________________________________________________________________
starlightConstants::event Gammaavectormeson::produceEvent(int&)
{
	// Not used; return default event
	return starlightConstants::event();
}


//______________________________________________________________________________
upcEvent Gammaavectormeson::produceEvent()
{
	// The new event type
	upcEvent event;

	int iFbadevent=0;
	int tcheck=0;
	starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
        starlightConstants::particleTypeEnum vmpid = starlightConstants::UNKNOWN; 

	if (_VMpidtest == starlightConstants::FOURPRONG) {
		double        comenergy = 0;
		double        mom[3]    = {0, 0, 0};
		double        E         = 0;
		lorentzVector decayVecs[4];
		do {
			double rapidity = 0;
			pickwy(comenergy, rapidity);
			if (_VMinterferencemode == 0)
				momenta(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
			else if (_VMinterferencemode==1)
				vmpt(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
		} while (!fourBodyDecay(ipid, E, comenergy, mom, decayVecs, iFbadevent));
		if ((iFbadevent == 0) and (tcheck == 0))
			for (unsigned int i = 0; i < 4; ++i) {
				starlightParticle daughter(decayVecs[i].GetPx(),
				                           decayVecs[i].GetPy(),
				                           decayVecs[i].GetPz(),
				                           starlightConstants::UNKNOWN,  // energy 
				                           starlightConstants::UNKNOWN,  // _mass
				                           ipid,
				                           (i < 2) ? -1 : +1);
				event.addParticle(daughter);
			}
	} else {
		double comenergy = 0.;
		double rapidity = 0.;
		double E = 0.;
		double momx=0.,momy=0.,momz=0.;

		double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
		bool accepted = false;
		//  if(_accCut){
		do{
			pickwy(comenergy,rapidity);

			if (_VMinterferencemode==0){
				momenta(comenergy,rapidity,E,momx,momy,momz,tcheck);
			} else if (_VMinterferencemode==1){
				vmpt(comenergy,rapidity,E,momx,momy,momz,tcheck);
			}
	   
			// cout << "_ptCutMin: " << _ptCutMin << " _ptCutMax: " << _ptCutMax << " _etaCutMin: " << _etaCutMin << " _etaCutMax: " << _etaCutMax << endl;
			_nmbAttempts++;
			//cout << "n tries: " << _nmbAttempts<< endl;
                        vmpid = ipid; 
			twoBodyDecay(ipid,E,comenergy,momx,momy,momz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
			double pt1chk = sqrt(px1*px1+py1*py1);
			double pt2chk = sqrt(px2*px2+py2*py2);
    
			//cout << "pt1: " << pt1chk  << " pt2: " << pt2chk << endl;
			double eta1 = pseudoRapidity(px1, py1, pz1);
			double eta2 = pseudoRapidity(px2, py2, pz2);
			//cout << "eta1: " << eta1 << " eta2: " << eta2 << endl;
			if(_ptCutEnabled && !_etaCutEnabled){
				if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
					accepted = true;
					_nmbAccepted++;
				}
			}
			else if(!_ptCutEnabled && _etaCutEnabled){
				if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
					accepted = true;
					_nmbAccepted++;
				}
			}
			else if(_ptCutEnabled && _etaCutEnabled){
				if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
					if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
						accepted = true;
						_nmbAccepted++;
					}
				}
			}
			else if(!_ptCutEnabled && !_etaCutEnabled)
				_nmbAccepted++;
	      
		}while((_ptCutEnabled || _etaCutEnabled) && !accepted);
		/*  }else{
		    twoBodyDecay(ipid,E,comenergy,momx,momy,momz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
		    }*/
		if (iFbadevent==0&&tcheck==0) {
			int q1=0,q2=0;
                        int ipid1,ipid2=0;

			double xtest = randyInstance.Rndom(); 
			if (xtest<0.5)
				{
					q1=1;
					q2=-1;
				}
			else {
				q1=-1;
				q2=1;
			}

                        if ( ipid == 11 || ipid == 13 ){
                          ipid1 = -q1*ipid;
                          ipid2 = -q2*ipid;
                        } else {
                          ipid1 = q1*ipid;
                          ipid2 = q2*ipid;
                        }
			//     The new stuff
			double md = getDaughterMass(vmpid); 
                        double Ed1 = sqrt(md*md+px1*px1+py1*py1+pz1*pz1); 
			starlightParticle particle1(px1, py1, pz1, Ed1, starlightConstants::UNKNOWN, ipid1, q1);
			event.addParticle(particle1);

                        double Ed2 = sqrt(md*md+px2*px2+py2*py2+pz2*pz2); 
			starlightParticle particle2(px2, py2, pz2, Ed2, starlightConstants::UNKNOWN, ipid2, q2);
			event.addParticle(particle2);
			//     End of the new stuff

		}
	}

	return event;

}
double Gammaavectormeson::pseudoRapidity(double px, double py, double pz)
{
	double pT = sqrt(px*px + py*py);
	double p = sqrt(pz*pz + pT*pT);
	double eta = -99.9; if((p-pz) != 0){eta = 0.5*log((p+pz)/(p-pz));}
	return eta;
}

//______________________________________________________________________________
Gammaanarrowvm::Gammaanarrowvm(beamBeamSystem& bbsystem):Gammaavectormeson(bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaanarrowvm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaanarrowvm()"<<endl;
	narrowResonanceCrossSection sigma(bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	_VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaanarrowvm::~Gammaanarrowvm()
{ }


//______________________________________________________________________________
Gammaaincoherentvm::Gammaaincoherentvm(beamBeamSystem& bbsystem):Gammaavectormeson(bbsystem)
{
        cout<<"Reading in luminosity tables. Gammaainkoherentvm()"<<endl;
        read();
        cout<<"Creating and calculating crosssection. Gammaainkoherentvm()"<<endl;
        incoherentVMCrossSection sigma(bbsystem); sigma.crossSectionCalculation(_bwnormsave);
        _VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaaincoherentvm::~Gammaaincoherentvm()
{ }


//______________________________________________________________________________
Gammaawidevm::Gammaawidevm(beamBeamSystem& bbsystem):Gammaavectormeson(bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaawidevm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaawidevm()"<<endl;
	wideResonanceCrossSection sigma(bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	_VMbslope=sigma.slopeParameter();
}


//______________________________________________________________________________
Gammaawidevm::~Gammaawidevm()
{ }


