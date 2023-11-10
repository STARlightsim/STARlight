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
// $Rev:: 313                         $: revision of last commit
// $Author:: aaronstanek              $: author of last commit
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
Gammaavectormeson::Gammaavectormeson(const inputParameters& inputParametersInstance, randomGenerator* randy, beamBeamSystem& bbsystem):eventChannel(inputParametersInstance, randy, bbsystem), _phaseSpaceGen(0)
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
        _bslopeDef=inputParametersInstance.bslopeDefinition();
	_bslopeVal=inputParametersInstance.bslopeValue();
	_pEnergy= inputParametersInstance.protonEnergy();
	_VMpidtest=inputParametersInstance.prodParticleType();
	_VMptmax=inputParametersInstance.maxPtInterference();
	_VMdpt=inputParametersInstance.ptBinWidthInterference();
        _ProductionMode=inputParametersInstance.productionMode();

        N0 = 0; N1 = 0; N2 = 0; 
	  if (_VMpidtest == starlightConstants::FOURPRONG || _VMpidtest == starlightConstants::OMEGA_pipipi){
		// create n-body phase-spage generator
		_phaseSpaceGen = new nBodyPhaseSpaceGen(_randy);
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

	xw = _randy->Rndom();
	W = _VMWmin + xw*(_VMWmax-_VMWmin);

	if (W < 2 * _ip->pionChargedMass())
		goto L201pwy;
  
	IW = int((W-_VMWmin)/dW);
	xy = _randy->Rndom();
	Y = _VMYmin + xy*(_VMYmax-_VMYmin);
	IY = int((Y-_VMYmin)/dY); 
	xtest = _randy->Rndom();

	if( xtest > _Farray[IW][IY] )
		goto L201pwy;

        N0++; 
	// Determine the target nucleus 
	// For pA this is given, for all other cases use the relative probabilities in _Farray1 and _Farray2 
        if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
           if( _ProductionMode == 2 || _ProductionMode ==3){
	     _TargetBeam = 2;
	   } else {
             _TargetBeam = 1;
           }
        } else if(  _bbs.beam1().A() != 1 && _bbs.beam2().A()==1 ){
           if( _ProductionMode == 2 || _ProductionMode ==3){
	     _TargetBeam = 1;
	   } else {
             _TargetBeam = 2;
           }
        } else {
          btest = _randy->Rndom();
	  if ( btest < _Farray1[IW][IY]/_Farray[IW][IY] ){
            _TargetBeam = 2;
            N2++;
          }  else {
            _TargetBeam = 1;
            N1++; 
          }
        }
}         


/**
 * @brief Creates the respective daughter particles for two body decay channels.
 * 
 * @param ipid [output reference] ipid of the daughter particle
 * @param W [input] Parent mass in CM Frame
 * @param px0 [input] Parent x-momentum in CM Frame
 * @param py0 [input] Parent y-momentum in CM Frame
 * @param pz0 [input] Parent z-momentum in CM Frame
 * @param E1 [output reference] Daughter 1's Energy in CM Frame
 * @param px1 [output reference] Daughter 1's x-momentum in CM Frame
 * @param py1 [output reference] Daughter 1's y-momentum in CM Frame
 * @param pz1 [output reference] Daughter 1's z-momentum in CM Frame
 * @param E2 [output reference] Daugter 2's Energy in CM Frame
 * @param px2 [output reference] Daughter 2's x-momentum in CM Frame
 * @param py2 [output reference] Daughter 2's y-momentum in CM Frame
 * @param pz2 [output reference] Daughter 2's z-momentum in CM Frame
 * @param iFbadevent 
 */
void Gammaavectormeson::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,
                                     double  W,
                                     double  px0, double  py0, double  pz0,
                                     double &E1, double& px1, double& py1, double& pz1,
                                     double &E2, double& px2, double& py2, double& pz2,
                                     int&    iFbadevent)
{
	// This routine decays a particle into two particles of mass mdec,
	// taking spin into account

	double pmag;
	double phi,theta,Ecm;
	double betax,betay,betaz;
	double mdec=0.0;
	//double E1=0.0,E2=0.0; not needed any more as references to these variables are now provided in parameter

	//    set the mass of the daughter particles
	mdec=getDaughterMass(ipid);

	//  calculate the magnitude of the momenta
	if(W < 2*mdec){
		cout<<" ERROR: W="<<W<<endl;
		iFbadevent = 1;
		return;
	}
	pmag = sqrt(W*W/4. - mdec*mdec);
  
	//  pick an orientation, based on the spin
	//  phi has a flat distribution in 2*pi
	phi = _randy->Rndom()*2.*starlightConstants::pi;
                                                                                                                
	//  find theta, the angle between one of the outgoing particles and
	//  the beamline, in the outgoing particles' rest frame

	theta=getTheta(ipid);
 
	//  compute unboosted momenta
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

/**
 * @brief decays a particle into three particles with isotropic angular distribution.
 * 
 * @param ipid [output reference] Holds the ipid of the charged particle
 * @param ipid2 [output reference] Holds the ipid of the uncharged particle
 * @param W [input] parent mass
 * @param p [input] parent momentum vector.
 * @param decayVecs [output reference] Array of daughter particles.
 * 						The first two particles are charged and the 3rd/last is uncharged.
 * @param iFbadevent [output reference] is set to 1 when decay is unsuccessful ONLY. N.B. Its value remains as initial when decay is successful - It does not reset to 0.
 * @return true if decay is successful
 * @return false if decay is unsuccessful
 */
bool Gammaavectormeson::omega3piDecay
(starlightConstants::particleTypeEnum& ipid,//This holds the ipid of the charged pion
starlightConstants::particleTypeEnum& ipid2,// This holds the ipid of the uncharged pion
 const double                  ,           // E (unused)
 const double                  W,          // mass of produced particle
 const double*                 p,          // momentum of produced particle; expected to have size 3
 lorentzVector*                decayVecs,  // array of Lorentz vectors of daughter particles; expected to have size 3
 int&                          iFbadevent)
{
	const double parentMass = W;

	// set the mass of the daughter particles
	const double daughterMass = getDaughterMass(ipid);//mass of the two charged particle
	const double _2ndDaughterMass = get2ndDaughterMass(ipid2);//mass of the 3rd uncharged particle
	if (parentMass < (2 * daughterMass + _2ndDaughterMass) ){
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
		const double m[3] = {_ip->pionChargedMass(), _ip->pionChargedMass(), _ip->pionNeutralMass()};
		_phaseSpaceGen->setDecay(3, m);
		// estimate maximum phase-space weight
		_phaseSpaceGen->setMaxWeight(1.01 * _phaseSpaceGen->estimateMaxWeight(_VMWmax));
		firstCall = false;
	}

	// generate phase-space event
	if (!_phaseSpaceGen->generateDecayAccepted(parentVec))
		return false;

	// set Lorentzvectors of decay daughters
	for (unsigned int i = 0; i < 3; ++i)
		decayVecs[i] = _phaseSpaceGen->daughter(i);
	return true;
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
/**
 * @brief This is used to determine the identity of the 2nd daughter particles for channels with more than one daughter type e.g. OMEGA_PIPIPI
 * 
 * @param ipid2 [output reference] holds the ipid of the second daughter specie.
 * @return [double]  The mass of the daughter specie.
 */
double Gammaavectormeson::get2ndDaughterMass(starlightConstants::particleTypeEnum &ipid2)
{
	// This is used only for channels producing more than one daughter type. Notably: Omega-3-pions
	//It returns both the mass of the 2nd daughter and the ipid of the second daughter.
	double mdec=0.;
  
	switch(_VMpidtest){
	case starlightConstants::OMEGA_pipipi:
		mdec = _ip->pionNeutralMass();
		ipid2 = starlightConstants::PIONNEUTRAL;
		break;
	default: cout<<"No 2nddaughtermass defined for this channel, try gammaavectormeson::getdaughtermass instead"<<endl;
	}
	return mdec;
}
//______________________________________________________________________________
double Gammaavectormeson::getDaughterMass(starlightConstants::particleTypeEnum &ipid)
{
	//This will return the daughter particles mass, and the final particles outputed id...
	double mdec=0.;
  
	switch(_VMpidtest){
	case starlightConstants::RHO:
	case starlightConstants::RHOZEUS:
	case starlightConstants::FOURPRONG:
	case starlightConstants::OMEGA:
	case starlightConstants::OMEGA_pipipi:
		mdec = _ip->pionChargedMass();
		ipid = starlightConstants::PION;
		break;
	case starlightConstants::PHI:
	case starlightConstants::PHIKK:
		mdec = _ip->kaonChargedMass();
		ipid = starlightConstants::KAONCHARGE;
		break;
	case starlightConstants::JPSI:
		mdec = _ip->mel();
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::RHO_ee:
	case starlightConstants::JPSI_ee:
	case starlightConstants::PHI_ee:
		mdec = _ip->mel();
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::RHO_mumu:
	case starlightConstants::JPSI_mumu:
		mdec = _ip->muonMass();
		ipid = starlightConstants::MUON;
		break; 
	case starlightConstants::JPSI_ppbar:
		mdec = _ip->protonMass();
		ipid = starlightConstants::PROTON;
		break; 
	case starlightConstants::JPSI2S_ee:
		mdec = _ip->mel();
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI2S_mumu:
		mdec = _ip->muonMass();
		ipid = starlightConstants::MUON;
		break; 

	case starlightConstants::JPSI2S:
	case starlightConstants::UPSILON:
	case starlightConstants::UPSILON2S:
	case starlightConstants::UPSILON3S:
		mdec = _ip->muonMass();
		ipid = starlightConstants::MUON;
		break;
	case starlightConstants::UPSILON_ee:
	case starlightConstants::UPSILON2S_ee:
	case starlightConstants::UPSILON3S_ee:
		mdec = _ip->mel();
		ipid = starlightConstants::ELECTRON;
		break;
	case starlightConstants::UPSILON_mumu:
	case starlightConstants::UPSILON2S_mumu:
	case starlightConstants::UPSILON3S_mumu:
		mdec = _ip->muonMass();
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
    
	theta = starlightConstants::pi*_randy->Rndom();
	xtest = _randy->Rndom();
	//  Follow distribution for helicity +/-1
	//  Eq. 19 of J. Breitweg et al., Eur. Phys. J. C2, 247 (1998)
	//  SRK 11/14/2000
  
	switch(ipid){
	  
	case starlightConstants::MUON:
	case starlightConstants::ELECTRON:
		//VM->ee/mumu
 	        dndtheta = sin(theta)*(1.+(cos(theta)*cos(theta)));
                break; 
		
	case starlightConstants::PROTON:
		//Pick this angular distribution for J/psi --> ppbar 
	        dndtheta = sin(theta)*(1.+(0.605*cos(theta)*cos(theta)));
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
void Gammaavectormeson::momenta(double W,double Y,
								double &E,double &px,double &py,double &pz,int &tcheck, //vector meson in cms frame
								double &Eb1, double &pxb1, double &pyb1, double &pzb1,//outgoing beam 1 in cms frame.
								double &Eb2, double &pxb2, double &pyb2, double &pzb2, double &t2, //outgoing beam 2 in cms frame.
								double &Egam, double&pxgam, double &pygam, double &pzgam, double &Q2gam) //photon in the cms frame.
{
	//     This subroutine calculates momentum and energy of vector meson
	//     given W and Y,   without interference.  Subroutine vmpt handles
	//     production with interference
 
	double Epom,tmin,pt1,pt2,phi1,phi2;
	double px1,py1,px2,py2;
	double pt,xt,xtest,ytest;
	//double t2, Egam;

  
	//Find Egam,Epom in CM frame
        if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
          // This is pA
          if( _ProductionMode == 2 || _ProductionMode ==3 ){
    	    Egam = 0.5*W*exp(Y);
  	    Epom = 0.5*W*exp(-Y);
          }else{
    	    Egam = 0.5*W*exp(-Y);
  	    Epom = 0.5*W*exp(Y);
          }  
        } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
          // This is Ap
          if( _ProductionMode == 2 || _ProductionMode == 3 ){
  	    Egam = 0.5*W*exp(-Y);
  	    Epom = 0.5*W*exp(Y);
          }else{
    	    Egam = 0.5*W*exp(Y);
  	    Epom = 0.5*W*exp (-Y);
          }
	} else {
          // This is pp or AA 
          if( _TargetBeam == 1 ){
            Egam = 0.5*W*exp(-Y);
	    Epom = 0.5*W*exp(Y);
	  }
          else {
            Egam = 0.5*W*exp(Y);
	    Epom = 0.5*W*exp(-Y);
	  }
	}

	//        } else if( _ProductionMode == 2 || _ProductionMode==3){
	//	  Egam = 0.5*W*exp(-Y);
	//	  Epom = 0.5*W*exp(Y);
	//        } else { 
	//          Egam = 0.5*W*exp(Y);
	//	  Epom = 0.5*W*exp(-Y);
	//	 }

        pt1 = pTgamma(Egam);  
	phi1 = 2.*starlightConstants::pi*_randy->Rndom();

	if( (_bbs.beam1().A()==1 && _bbs.beam2().A()==1) || 
            (_ProductionMode == 4) ) {
	    if( (_VMpidtest == starlightConstants::RHO) || (_VMpidtest == starlightConstants::RHOZEUS) || (_VMpidtest == starlightConstants::OMEGA)){
	      // Use dipole form factor for light VM
	      L555vm:
	      xtest = 2.0*_randy->Rndom();
              double ttest = xtest*xtest; 
              ytest = _randy->Rndom();
              double t0 = 1./2.23; 
              double yprob = xtest*_bbs.beam1().dipoleFormFactor(ttest,t0)*_bbs.beam1().dipoleFormFactor(ttest,t0); 
              if( ytest > yprob ) goto L555vm; 
              t2 = ttest; 
              pt2 = xtest;              
	    }else{
		//Use dsig/dt= exp(-_VMbslope*t) for heavy VM
                double bslope_tdist = _VMbslope; 
		double Wgammap = 0.0; 
                switch(_bslopeDef){
		  case 0:
		    //This is the default, as before
		    bslope_tdist = _VMbslope;
		    break;
		  case 1:
		    //User defined value of bslope. BSLOPE_VALUE default is 4.0 if not set. 
                    bslope_tdist = _bslopeVal;
		    if( N0 <= 1 )cout<<" ATTENTION: Using user defined value of bslope = "<<_bslopeVal<<endl;
                    break; 
		  case 2:
                    //This is Wgammap dependence of b from H1 (Eur. Phys. J. C 46 (2006) 585)
		    Wgammap = sqrt(4.*Egam*_pEnergy); 
		    bslope_tdist = 4.63 + 4.*0.164*log(Wgammap/90.0);
		    if( N0 <= 1 )cout<<" ATTENTION: Using energy dependent value of bslope!"<<endl; 
		    break;
		  default:
		    cout<<" Undefined setting for BSLOPE_DEFINITION "<<endl;
		}

	        xtest = _randy->Rndom(); 
		// t2 = (-1./_VMbslope)*log(xtest);
		t2 = (-1./bslope_tdist)*log(xtest);
		pt2 = sqrt(1.*t2);
	    }
	} else {
	    // >> Check tmin
	    tmin = ((Epom/_VMgamma_em)*(Epom/_VMgamma_em));
	
	    if(tmin > 0.5){
		cout<<" WARNING: tmin= "<<tmin<<endl;
                cout<< " Y = "<<Y<<" W = "<<W<<" Epom = "<<Epom<<" gamma = "<<_VMgamma_em<<endl; 
		cout<<" Will pick a new W,Y "<<endl;
		tcheck = 1;
		return;
	    }
 L203vm:
	    xt = _randy->Rndom(); 
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
		
 // Changed '8' to '32' 6 times below to extend the region of the p_T calculation up to 1 GeV.c  SRK May 28, 2019
//  'pt2' is the maximum vector meson momentum.  For heavy nuclei, the '32'coefficient corresonds to about 1 GeV/c
//  The downside of the larger coefficient is that the sampling runs more slowly.  This could be made into a parameter		
		

		pt2 = 32.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
              }else{
                pt2 = 32.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
              }   
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pt2 = 32.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
              }else{
                pt2 = 32.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
              }  
            } else if (_TargetBeam==1) {
                pt2 = 32.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
            } else {
                pt2 = 32.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
            }

	    xtest = _randy->Rndom();
	    t2 = tmin + pt2*pt2;

	    double comp=0.0; 
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
                comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
              }else{
                comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
              }   
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
                comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
              }else{
                comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
              }  
            } else if (_TargetBeam==1) {
              comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
            } else {
              comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2; 
            }	      
            if( xtest > comp ) goto L203vm;
       		
	}//else end from pp

	phi2 = 2.*starlightConstants::pi*_randy->Rndom();

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

	//pt1 is for gamma
	pxgam = px1;
	pygam = py1;//t2 and Egam naturally fixed._ip->protonMass()*_ip->beam2A()*_ip->protonMass()*_ip->beam2A()
	double E0b1 = _pEnergy*_ip->beam1A();
	double px0b1 = 0, px0b2 =0, py0b1=0, py0b2=0;
	double pz0b1 = sqrt(E0b1*E0b1 - _ip->protonMass()*_ip->protonMass()*_ip->beam1A()*_ip->beam1A());
	double E0b2 =_pEnergy*_ip->beam2A();
	double pz0b2 = -sqrt(E0b2*E0b2 - _ip->protonMass()*_ip->protonMass()*_ip->beam2A()*_ip->beam2A());

	double pzgam1, pzgam2, pzgamB, Q2gamB;
  	
	if(_TargetBeam == 2){
		Eb2 = Egam + E0b2 - E;//
		pxb2 = pxgam +px0b2 -px;//
		pyb2 = pygam + py0b2 - py;// correct nxt line
		pzb2 = -sqrt(Eb2*Eb2 - (pxb2*pxb2 + pyb2*pyb2 + _ip->protonMass()*_ip->beam2A()*_ip->protonMass()*_ip->beam2A()));
	
		Eb1 = E0b1 -Egam;
		pxb1 = px0b1 - pxgam;
		pyb1 = py0b1 - pygam;	
		pzb1 = sqrt(Eb1*Eb1 - (pxb1*pxb1 + pyb1*pyb1 + _ip->beam1A()*_ip->protonMass()*_ip->beam1A()*_ip->protonMass()));//correct
		
		pzgam1 = pz0b1- pzb1;
		pzgam2 = pz + pzb2 - pz0b2;
		
		//pzgamA = (pzgam1 + pzgam2)/2.0;//determining pzgam with the collapsed graph structure
		//Q2gamA = Egam*Egam - (pxgam*pxgam + pygam*pygam + pzgamA*pzgamA);//virtuality of photon in collapsed graph structure.

		pzgamB = (2*pzgam1 + pzgam2)/3.0;//z-momentum of photon using the  complete graph structure
		Q2gamB = Egam*Egam - (pxgam*pxgam + pygam*pygam + pzgamB*pzgamB);//virtuality of photon in complete graph structure.

		//pzgam3 = pz0b2 - pzb2;
		//pzPom = (2*pzgam3 + (pz-pzgam1))/3.0;//z-momentum of Pomeron using the  complete graph structure
	
		pzgam = pzgamB;
		Q2gam = Q2gamB;

	}
	else if(_TargetBeam == 1){
		Eb1 = Egam + E0b1 - E;//
		pxb1 = pxgam +px0b1 -px;//
		pyb1 = pygam + py0b1 - py;//
		pzb1 = sqrt(Eb1*Eb1 - (pxb1*pxb1 + pyb1*pyb1 + _ip->protonMass()*_ip->beam1A()*_ip->protonMass()*_ip->beam1A()));
		//pzgam = pzb1 + pz - pz0b1;
		//Q2gam = Egam*Egam - (pxgam*pxgam + pygam*pygam + pzgam*pzgam);
		Eb2 = E0b2 -Egam;
		pxb2 = px0b2 - pxgam;
		pyb2 = py0b2 - pygam;	
		//pzb2 = pz0b2 - pzgam;
		pzb2 = -sqrt(Eb2*Eb2 - (pxb2*pxb2 + pyb2*pyb2 + _ip->protonMass()*_ip->beam2A()*_ip->protonMass()*_ip->beam2A()));
		
		
		pzgam1 = pz0b2- pzb2;
		pzgam2 = pz + pzb1 - pz0b1;
		
		//pzgamA = (pzgam1 + pzgam2)/2.0;
		//Q2gamA = Egam*Egam - (pxgam*pxgam + pygam*pygam + pzgamA*pzgamA);//virtuality of photon in collapsed graph structure.

		pzgamB = (2.0*pzgam1 + pzgam2)/3.0;//z-momentum of photon using the  complete graph structure
		Q2gamB = Egam*Egam - (pxgam*pxgam + pygam*pygam + pzgamB*pzgamB);//virtuality of photon in complete graph structure.

		//pzgam3 = pz0b1 - pzb1;
		//pzPom = (2*pzgam3 + (pz-pzgam1))/3.0;//z-momentum of Pomeron using the  complete graph structure

		pzgam = pzgamB;
		Q2gam = Q2gamB;

	}
	else{
		cout << " ERROR: Target Beam Number: " << _TargetBeam << "is invalid" << endl;	
	}


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
    // If E is very small, the drawing of a pT below is extre_ip->mel()y slow. 
    // ==> Set pT = sqrt(3.)*E/_VMgamma_em for very small E. 
    // Should have no observable consequences (JN, SRK 11-Sep-2014)
    if( E < 0.0005 )return Cm; 
 
    //the amplitude of the p_t spectrum at the maximum

    if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
      if( _ProductionMode == 2 || _ProductionMode ==3 ){
         singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
      }  
    } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
      if( _ProductionMode == 2 || _ProductionMode ==3){
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
    x = _randy->Rndom();

    if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
      if( _ProductionMode == 2 || _ProductionMode ==3){
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
      }  
    } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
      if( _ProductionMode == 2 || _ProductionMode ==3){
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
	u = _randy->Rndom();
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =_randy->Rndom();
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);
              }  
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
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
	double dY=0.,yleft=0.,yfract=0.,xpt=0.,pt1=0.,ptfract=0.,pt=0.,pt2=0.,theta=0.;
	int IY=0,j=0;
  
	dY  = (_VMYmax-_VMYmin)/double(_VMnumy);
  
	//  Y is already fixed; choose a pt
	//  Follow the approach in pickwy
	//  in  _fptarray(IY,pt) IY=1 corresponds to Y=0, IY=numy/2 corresponds to +y
 	//  Changed,  now works -y to +y.
	IY=int((Y-_VMYmin)/dY);
	if (IY > (_VMnumy)-1){
        	IY=(_VMnumy)-1;
	}

	yleft=(Y-_VMYmin)-(IY)*dY;

	yfract=yleft*dY;
  
	xpt=_randy->Rndom();
	for(j=0;j<_VMNPT;j++){
		if (xpt < _fptarray[IY][j]) goto L60;
	}
	if(j == _VMNPT) j = _VMNPT-1;
 L60:
  
	//  now do linear interpolation - start with extremes
  	if (j == 0){
		pt1=xpt/_fptarray[IY][j]*_VMdpt/2.;
		goto L80;
	}
	if (j == _VMNPT-1){
		pt1=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY][j])/(1.-_fptarray[IY][j]);
		goto L80;
	}
  
	//  we're in the middle
  	ptfract=(xpt-_fptarray[IY][j])/(_fptarray[IY][j+1]-_fptarray[IY][j]);
	pt1=(j+1)*_VMdpt+ptfract*_VMdpt;
  
	//  at an extreme in y?
	if (IY == (_VMnumy)-1){
		pt=pt1;
		goto L120;
	}
 L80:

	//  interpolate in y repeat for next fractional y bin      
	for(j=0;j<_VMNPT;j++){
		if (xpt < _fptarray[IY+1][j]) goto L90;
	}
        if(j == _VMNPT) j = _VMNPT-1;
 L90:
  
	//  now do linear interpolation - start with extremes
	if (j == 0){
		pt2=xpt/_fptarray[IY+1][j]*_VMdpt/2.;
		goto L100;
	}
	if (j == _VMNPT-1){
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
	theta=2.*starlightConstants::pi*_randy->Rndom();
	px=pt*cos(theta);
	py=pt*sin(theta);

	E  = sqrt(W*W+pt*pt)*cosh(Y);
	pz = sqrt(W*W+pt*pt)*sinh(Y);
	//      randomly choose to make pz negative 50% of the time
	if(_randy->Rndom()>=0.5) pz = -pz;
}


//______________________________________________________________________________
starlightConstants::event Gammaavectormeson::produceEvent(int&)
{
	//Note used; return default event
	return starlightConstants::event();
}


/**
 * @brief 
 * 
 * @param beta 
 * @return upcEvent 
 */
//upcEvent Gammaavectormeson::produceEvent(vector3 beta)
upcXEvent Gammaavectormeson::produceEvent(vector3 beta)
{
	
	//upcEvent event;//former type output variable
	// The new event type
	upcXEvent event;//output variable

	int iFbadevent=0;//variable to track successful and unsuccessful vector Meson decay
	int tcheck=0;//variable to track successful and unsuccessful Vector Meson creation
	starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;//stores the ipid of the daughter particle ( or main daughter- if  there are more than one daughter)
	starlightConstants::particleTypeEnum ipid2 = starlightConstants::UNKNOWN;//used to store the ipid of 2nd daughters in channels that have more than one daughter particle.
        starlightConstants::particleTypeEnum vmpid = starlightConstants::UNKNOWN;//used for temporary local storage in twobodydecays(). 

	double ptCutMin2 = _ptCutMin*_ptCutMin;//used for ptCut comparison without using square roots - to reduce processing time
	double ptCutMax2 = _ptCutMax*_ptCutMax;//same as above
	double Pgam[4] = {0.0,0.0,0.0,0.0};//E,pxgam,pygam,pzgam
	double Pb1[4] = {0.0,0.0,0.0,0.0};//beam1 - Eb1, pxb1,pyb1,pzb1
	double Pb2[4] = {0.0,0.0,0.0,0.0};//beam2 - Eb2, pxb2,pyb2,pzb2
	double Q2gam =0., t=0.; //q2 for photon, t - transferred momenta squared
	if (_VMpidtest == starlightConstants::FOURPRONG) {
		double        comenergy = 0;
		double        mom[3]    = {0, 0, 0};
		double        E         = 0;
		lorentzVector decayVecs[4];
		bool accepted;
		double rapidity = 0;
		do {
			tcheck = 0;//reinitialized after every loop cycle - to avoid infinite loop
			iFbadevent = 0;//same as above.			
			pickwy(comenergy, rapidity);

			//Vector meson is created and its four momentum is determined below
			if (_VMinterferencemode == 0)
				momenta(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck,
						Pb1[0], Pb1[1],Pb1[2],Pb1[3],
						Pb2[0], Pb2[1],Pb2[2],Pb2[3],t,
						Pgam[0], Pgam[1],Pgam[2],Pgam[3],Q2gam);//without interference
			else if (_VMinterferencemode==1)
				vmpt(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);//with interference
			_nmbAttempts++;
			accepted = true;//re-initialized after every loop cycle -to avoid infinite loop

			
			if(tcheck != 0 || !fourBodyDecay(ipid, E, comenergy, mom, decayVecs, iFbadevent))
			{//if either vector meson creation, or further decay into four pions, is impossible
				accepted = false;
				continue;//this skips the etaCut and ptCut checks.
			}

			if (_ptCutEnabled) {
				for (int i = 0; i < 4; i++) {
					double pt_chk2 = 0;
					pt_chk2 += pow( decayVecs[i].GetPx() , 2);
					pt_chk2 += pow( decayVecs[i].GetPy() , 2);// compute transverse momentum (squared) for the particle, for ptCut checks.

					if (pt_chk2 < ptCutMin2 || pt_chk2 > ptCutMax2) {//if particle does not fall into ptCut range.
						accepted = false;
						break;//skips checking other daughter particles
					}
				}
			}
			if (_etaCutEnabled) {
				for (int i = 0; i < 4; i++) {
					double eta_chk = pseudoRapidityLab(
						decayVecs[i].GetPx(),
						decayVecs[i].GetPy(),
						decayVecs[i].GetPz(),
						decayVecs[i].GetE(),
						beta
					);//computes the pseudorapidity in the laboratory frame.
					if (eta_chk < _etaCutMin || eta_chk > _etaCutMax) {//if this particle does not fall into range
						accepted = false;
						break;//skips checking other daughter particles
					}
				}
			}
			if (accepted and (tcheck == 0)) {
				_nmbAccepted++;//maintain counts of accepted events.
			}

		} while (!accepted || tcheck != 0);//repeats loop if VM creation, decay, ptcut or etaCut criterias are not fulfilled. Important to avoid situations where events produced is less than requested.

		double md = getDaughterMass(ipid);
		if ((iFbadevent == 0) and (tcheck == 0)){
		//adds daughters as particles into the output event.
			for (unsigned int i = 0; i < 4; ++i) {
				starlightParticle daughter(decayVecs[i].GetPx(),
				                           decayVecs[i].GetPy(),
				                           decayVecs[i].GetPz(),
							   sqrt(decayVecs[i].GetPx()*decayVecs[i].GetPx()+decayVecs[i].GetPy()*decayVecs[i].GetPy()+decayVecs[i].GetPz()*decayVecs[i].GetPz()+md*md),//energy 
							   md,  // _mass
							   ipid*(2*(i/2)-1),   // make half of the particles pi^+, half pi^-
							   (2*(i/2)-1));//charge
				event.addParticle(daughter);
			}
			if (_ip->giveExtraBeamInfo()){
				lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
				lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
				double targetEgamma, rap1cm = acosh(_ip->beamLorentzGamma()),cmsEgam = Pgam[0], Pzgam = Pgam[3];
				lorentzVector gamma(Pgam[1],Pgam[2],Pzgam,cmsEgam);
				lorentzVector vmeson(mom[0],mom[1],mom[2],E);

				if(_TargetBeam == 1)
				targetEgamma = cmsEgam*cosh(rap1cm) - Pzgam*sinh(rap1cm);
				else
				targetEgamma = cmsEgam*cosh(rap1cm) + Pzgam*sinh(rap1cm);

				event.addVectorMeson(vmeson);
				if(_TargetBeam == 1)
					event.addGammaFromBeam2(gamma,targetEgamma,Q2gam);
				else if(_TargetBeam == 2)
					event.addGammaFromBeam1(gamma,targetEgamma,Q2gam);
				
				event.addOutgoingBeams(beam1,beam2);
				event.addVertext(t);
			}
		}
	} else if (_VMpidtest == starlightConstants::OMEGA_pipipi) {
		double        comenergy = 0;
		double        mom[3]    = {0, 0, 0};
		double        E         = 0;
		lorentzVector decayVecs[3];
		bool accepted;
		double mass, rapidity = 0;
		int charge;
		do {
			tcheck = 0;//reinitialized after every loop to avoid inifinite loop traps.
			iFbadevent = 0;//same as above.
			pickwy(comenergy, rapidity);
			//Creating the vector meson.
			if (_VMinterferencemode == 0)
				momenta(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck,
						Pb1[0], Pb1[1],Pb1[2],Pb1[3],
						Pb2[0], Pb2[1],Pb2[2],Pb2[3],t,
						Pgam[0], Pgam[1],Pgam[2],Pgam[3],Q2gam);
			else if (_VMinterferencemode==1)
				vmpt(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
			_nmbAttempts++;
			accepted = true;//re-initialized after every loop

			//VM decay excecuted and checked if it is successful. Otherwise restart the loop
			if(tcheck != 0 || !omega3piDecay(ipid, ipid2, E, comenergy, mom, decayVecs, iFbadevent)){
				accepted = false;
				continue;
			}
			//For etaCuts and Momentum cuts only charged particles are considered.
			//The first two particles are the charged particles while the last particle is uncharged.
			if (_ptCutEnabled) {
				for (int i = 0; i < 2/*//only the first two particles are charged particles*/; i++) {
					double pt_chk2 = 0;
					pt_chk2 += pow( decayVecs[i].GetPx() , 2);
					pt_chk2 += pow( decayVecs[i].GetPy() , 2);
					
					if (pt_chk2 < ptCutMin2 || pt_chk2 > ptCutMax2) {
						accepted = false;
						break;
					}
				}
			}
			if (_etaCutEnabled) {
				for (int i = 0; i < 2/*//only the first two particles are charged particles*/; i++) {
					
					double eta_chk = pseudoRapidityLab(
						decayVecs[i].GetPx(),
						decayVecs[i].GetPy(),
						decayVecs[i].GetPz(),
						decayVecs[i].GetE(),
						beta
					);
					if (eta_chk < _etaCutMin || eta_chk > _etaCutMax) {
						accepted = false;
						break;
					}
				}
			}
			if (accepted and (tcheck == 0)) {
				_nmbAccepted++;//maintain count of successfully accepted events.
			}
		} while (!accepted || tcheck != 0);//repeats loop if VM creation, decay,ptCut or etaCut requirements is not satisfied.

		if ((iFbadevent == 0) and (tcheck == 0)){
			for (unsigned int i = 0; i < 3; ++i) {
				if(i<2){
					mass = getDaughterMass(ipid);
					charge = 2*i - 1;//make the first negative and the second positive.
				}else{
					ipid = ipid2;//ensures that the right ipid is written for the neutral pion.
					charge = 0; // the last particle is neutral. A neutral pion
					mass = get2ndDaughterMass(ipid2);
				}
				starlightParticle daughter(decayVecs[i].GetPx(),
				                           decayVecs[i].GetPy(),
				                           decayVecs[i].GetPz(),
							   sqrt(decayVecs[i].GetPx()*decayVecs[i].GetPx()+decayVecs[i].GetPy()*decayVecs[i].GetPy()+decayVecs[i].GetPz()*decayVecs[i].GetPz()+mass*mass),//energy 
							   mass,  // _mass
							   ipid*(charge==0? 1: charge),   // make first two particles pi^-, pi^+ and last uncharged
							   charge);//charge
				event.addParticle(daughter);
			}
			if(_ip->giveExtraBeamInfo())
			{
				lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
				lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
				double targetEgamma, rap1cm = acosh(_ip->beamLorentzGamma()),cmsEgam = Pgam[0], Pzgam = Pgam[3];
				lorentzVector gamma(Pgam[1],Pgam[2],Pzgam,cmsEgam);
				lorentzVector vmeson(mom[0],mom[1],mom[2],E);

				if(_TargetBeam == 1)
				targetEgamma = cmsEgam*cosh(rap1cm) - Pzgam*sinh(rap1cm);
				else
				targetEgamma = cmsEgam*cosh(rap1cm) + Pzgam*sinh(rap1cm);

				event.addVectorMeson(vmeson);
				if(_TargetBeam == 1)
					event.addGammaFromBeam2(gamma,targetEgamma,Q2gam);
				else if(_TargetBeam == 2)
					event.addGammaFromBeam1(gamma,targetEgamma,Q2gam);
				
				event.addOutgoingBeams(beam1,beam2);
				event.addVertext(t);
			}
		}
	} else {
		//initializations
		double comenergy = 0.;
		double rapidity = 0.;
		double E = 0.;
		double momx=0.,momy=0.,momz=0.;

		double E2=0.,E1=0., px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
		bool accepted = false;
		do{
			tcheck = 0;//reset after every loop cycles - to avoid infinite loop traps
			iFbadevent = 0;//same as above

			pickwy(comenergy,rapidity);

			//Vector meson creation
			if (_VMinterferencemode==0){
				momenta(comenergy,rapidity,E,momx,momy,momz,tcheck,
						Pb1[0], Pb1[1],Pb1[2],Pb1[3],
						Pb2[0], Pb2[1],Pb2[2],Pb2[3],t,
						Pgam[0], Pgam[1],Pgam[2],Pgam[3],Q2gam);
			
			} else if (_VMinterferencemode==1){
				vmpt(comenergy,rapidity,E,momx,momy,momz,tcheck);
			}
	   
			_nmbAttempts++;

                        vmpid = ipid;
			twoBodyDecay(ipid,comenergy,momx,momy,momz,E1,px1,py1,pz1,E2,px2,py2,pz2,iFbadevent);//vector meson decay.
			if(tcheck !=  0 || iFbadevent != 0){
				continue;//jumps to loop repeat if either VM Creation or decay is unsuccessful.
			}
			double pt1chk2 = px1*px1+py1*py1;//computes transverse momentum (squared) of daughter 1 for ptCut checks
			double pt2chk2 = px2*px2+py2*py2;
			double eta1 = pseudoRapidityLab(px1, py1, pz1,E1,beta);//computes eta of daughter 1 in lab frame for etaCut checks
			double eta2 = pseudoRapidityLab(px2, py2, pz2,E2,beta);

			if(_ptCutEnabled && !_etaCutEnabled){//ptCut only enabled
				if(pt1chk2 > ptCutMin2 && pt1chk2 < ptCutMax2 &&  pt2chk2 > ptCutMin2 && pt2chk2 < ptCutMax2){//ptCut of both daughters within range
					accepted = true;
					_nmbAccepted++;
				}
			}
			else if(!_ptCutEnabled && _etaCutEnabled){//etaCut only enabled
				if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){//etaCut of both daughters within range
					accepted = true;
					_nmbAccepted++;
				}
			}
			else if(_ptCutEnabled && _etaCutEnabled){//both Cuts enabled
				if(pt1chk2 > ptCutMin2 && pt1chk2 < ptCutMax2 &&  pt2chk2 > ptCutMin2 && pt2chk2 < ptCutMax2){
					if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){//ptCut and etaCuts of both daughters within range
						accepted = true;
						_nmbAccepted++;
					}
				}
			}
			else if(!_ptCutEnabled && !_etaCutEnabled)//no Cuts enabled
				_nmbAccepted++;
		}while(tcheck !=  0 || iFbadevent != 0 || ((_ptCutEnabled || _etaCutEnabled) && !accepted));//repeats loop if either VM Creation, VM decay,ptCuts or etaCuts is unsuccessful.

		if (iFbadevent==0&&tcheck==0) {
			int q1=0,q2=0;
                        int ipid1,ipid2=0;

			//randomly assign charges to daughter particles
			double xtest = _randy->Rndom(); 
			if (xtest<0.5)
				{
					q1=1;
					q2=-1;
				}
			else {
				q1=-1;
				q2=1;
			}

                        //sets ipid of daughter particles based on their charges.
						if ( ipid == starlightConstants::ELECTRON || ipid == starlightConstants::MUON ){
                          ipid1 = -q1*ipid;
                          ipid2 = -q2*ipid;
                        } else {
                          ipid1 = q1*ipid;
                          ipid2 = q2*ipid;
                        }

			double md = getDaughterMass(vmpid);
                        double Ed1 = sqrt(md*md+px1*px1+py1*py1+pz1*pz1); 
			starlightParticle particle1(px1, py1, pz1, Ed1,md, ipid1, q1);
			event.addParticle(particle1);

                        double Ed2 = sqrt(md*md+px2*px2+py2*py2+pz2*pz2); 
			starlightParticle particle2(px2, py2, pz2, Ed2, md, ipid2, q2);
			event.addParticle(particle2);
			
			if(_ip->giveExtraBeamInfo() ){
				lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
				lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
				double targetEgamma, rap1cm = acosh(_ip->beamLorentzGamma()),cmsEgam = Pgam[0], Pzgam = Pgam[3];
				lorentzVector gamma(Pgam[1],Pgam[2],Pzgam,cmsEgam);
				lorentzVector vmeson(momx,momy,momz,E);


				if(_TargetBeam == 1)
				targetEgamma = cmsEgam*cosh(rap1cm) - Pzgam*sinh(rap1cm);
				else
				targetEgamma = cmsEgam*cosh(rap1cm) + Pzgam*sinh(rap1cm);

				event.addVectorMeson(vmeson);
				if(_TargetBeam == 1)
					event.addGammaFromBeam2(gamma,targetEgamma,Q2gam);
				else if(_TargetBeam == 2)
					event.addGammaFromBeam1(gamma,targetEgamma,Q2gam);
				
				event.addOutgoingBeams(beam1,beam2);
				event.addVertext(t);
			}

		}
	}

	return event;

}


//______________________________________________________________________________
Gammaanarrowvm::Gammaanarrowvm(const inputParameters& input, randomGenerator* randy, beamBeamSystem& bbsystem):Gammaavectormeson(input, randy, bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaanarrowvm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaanarrowvm()"<<endl;
	narrowResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaanarrowvm::~Gammaanarrowvm()
{ }


//______________________________________________________________________________
Gammaaincoherentvm::Gammaaincoherentvm(const inputParameters& input, randomGenerator* randy, beamBeamSystem& bbsystem):Gammaavectormeson(input, randy, bbsystem)
{
        cout<<"Reading in luminosity tables. Gammaainkoherentvm()"<<endl;
        read();
        cout<<"Creating and calculating crosssection. Gammaaincoherentvm()"<<endl;
        incoherentVMCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
        _VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaaincoherentvm::~Gammaaincoherentvm()
{ }


//______________________________________________________________________________
Gammaawidevm::Gammaawidevm(const inputParameters& input, randomGenerator* randy, beamBeamSystem& bbsystem):Gammaavectormeson(input, randy, bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaawidevm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaawidevm()"<<endl;
	wideResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter();
}


//______________________________________________________________________________
Gammaawidevm::~Gammaawidevm()
{ }
