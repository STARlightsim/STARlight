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
//    Nystrand 220710
//    Fixed bug which gave incorrect minv distribution in gammagammaleptonpair.
//    Moved loop over W and Y from pickw to twoLeptonCrossSection in
//    gammagammaleptonpair to speed up event generation.
//    Changed to Woods-Saxon radius in twophotonluminosity to be consistent
//    with old starligt.
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "starlightconstants.h"
#include "gammagammaleptonpair.h"


using namespace std;


//_____________________________________________________________________________
Gammagammaleptonpair::Gammagammaleptonpair(const inputParameters& inputParametersInstance, randomGenerator* randy, beamBeamSystem& bbsystem)
: eventChannel(inputParametersInstance, randy, bbsystem)
, _GGlepInputpidtest(inputParametersInstance.prodParticleType())
, _GGlepInputnumw(inputParametersInstance.nmbWBins())
, _GGlepInputnumy(inputParametersInstance.nmbRapidityBins())
, _GGlepInputGamma_em(inputParametersInstance.beamLorentzGamma())
{
    //Initialize randomgenerator with our seed.
    _randy->Rndom();
    //Let us read in the luminosity tables
    read();
    //Now we will calculate the crosssection
    twoLeptonCrossSection();
    //If it is a tauon, calculate its tables
    if(inputParametersInstance.prodParticleId()==starlightConstants::TAUONDECAY) calculateTable();
}


//______________________________________________________________________________
Gammagammaleptonpair::~Gammagammaleptonpair()
{ }


//______________________________________________________________________________
void Gammagammaleptonpair::twoLeptonCrossSection()
{
    //This function calculates section for 2-particle decay. For reference, see STAR Note 243, Eq. 9.
    //calculate the 2-lepton differential cross section
    //the 100 is to convert to barns
    //the 2 is to account for the fact that we integrate only over one half of the rapidity range
    //Multiply all _Farray[] by _f_max

    for(int i=0;i<_GGlepInputnumw;i++)
    {
	for(int j=0;j<_GGlepInputnumy;j++)
	{
            _sigmax[i][j]=twoMuonCrossSection(_Warray[i])*_f_max*_Farray[i][j]/(100.*_Warray[i]);
	}
    }
    //calculate the total two-lepton cross section
    double sigmasum =0.;
    for(int i=0;i<_GGlepInputnumw-1;i++)
    {
	for(int j=0;j<_GGlepInputnumy-1;j++)
	{
          sigmasum = sigmasum +(_sigmax[i][j]+_sigmax[i+1][j]+_sigmax[i][j+1]+_sigmax[i+1][j+1])/4.*(_Yarray[j+1]-_Yarray[j])*(_Warray[i+1]-_Warray[i]);
	}
    }
    //cout << "The total "<<_GGlepInputpidtest<<" cross-section is: "<<sigmasum<<" barns."<<endl;
    cout<<endl<<endl;
    if (sigmasum > 1.){
       cout << "Total cross section: "<<sigmasum<<" barn."<<endl;  
    } else if (1000.*sigmasum > 1.){
       cout << "Total cross section: "<<1000.*sigmasum<<" mb."<<endl;  
    } else if (1000000.*sigmasum > 1.){
      cout << "Total cross section: "<<1000000.*sigmasum<<" microbarn."<<endl;  
    } else if (1.E9*sigmasum > 1.){
       cout << "Total cross section: "<<1.E9*sigmasum<<" nanobarn."<<endl;  
    } else if (1.E12*sigmasum > 1.){
       cout << "Total cross section: "<<1.E12*sigmasum<<" picobarn."<<endl;  
    } else {
       cout << "Total cross section: "<<1.E15*sigmasum<<" femtobarn."<<endl;  
    }
    cout<<endl; 
    setTotalChannelCrossSection(sigmasum);
    
    // Do this integration here, once per run rather than once per event (JN 220710) 
    //integrate sigma down to a function of just w

    double sgf=0.;

    for(int i=0;i<_ReadInputnumw;i++)
    {
            _sigofw[i]=0.;
            for(int j=0;j<_ReadInputnumy-1;j++)
            {
                _sigofw[i] = _sigofw[i]+(_Yarray[j+1]-_Yarray[j])*(_sigmax[i][j+1]+_sigmax[i][j])/2.;
            }
    }

    //calculate the unnormalized sgfint array
    _sigfint[0]=0.;
    for(int i=0;i<_ReadInputnumw-1;i++)
    {
        sgf=(_sigofw[i+1]+_sigofw[i])*(_Warray[i+1]-_Warray[i])/2.;
        _sigfint[i+1]=_sigfint[i]+sgf;
     }

     //normalize sgfint array
     _signormw=_sigfint[_ReadInputnumw-1];
     for(int i=0;i<_ReadInputnumw;i++)
     {
          _sigfint[i]=_sigfint[i]/_signormw;
     }
    return;
}


//______________________________________________________________________________
double Gammagammaleptonpair::twoMuonCrossSection(double w)
{
    //This function gives the two muon cross section as a function of Y&W. 
    //Using the formula given in G.Soff et. al Nuclear Equation of State, part B, 579
    double s=0.,Etest=0.,deltat=0.,xL=0.,sigmuix=0.,alphasquared=0.,hbarcsquared=0.;
    s = w*w;
    Etest = 4.*getMass()*getMass()/s;  
    deltat = s * sqrt(1.-Etest);
    xL = 2.*log(sqrt(s)/(2.*getMass())+sqrt(1./Etest-1.));
    alphasquared = starlightConstants::alpha*starlightConstants::alpha;
    hbarcsquared = starlightConstants::hbarc*starlightConstants::hbarc;
    sigmuix = 4.*starlightConstants::pi*alphasquared/s*hbarcsquared*((1+Etest-0.5*Etest*Etest)*xL-(1./s+Etest/s)*deltat);
    if(Etest > 1.) 
	sigmuix = 0.;
    return sigmuix;
}


//______________________________________________________________________________
void Gammagammaleptonpair::pickw(double &w)
{
//  This function picks a w for the 2-photon calculation.

    double x=0.,remainarea=0.,remainw=0.,a=0.,b=0.,c=0.;
    int ivalw=0;

    if(_wdelta != 0)
    {
	w=_wdelta;
	ivalw=_ivalwd;
	remainw=_remainwd;
    }
    else{
	//deal with the case where sigma is an array
	//_sigofw is simga integrated over y using a linear interpolation
	//sigint is the integral of sgfint, normalized

	//pick a random number
	x = _randy->Rndom();
	//compare x and sgfint to find the ivalue which is just less than the random number x
	for(int i=0;i<_GGlepInputnumw;i++)
	{
	    if(x > _sigfint[i]) ivalw=i;
	}
	//remainder above ivalw
	remainarea = x - _sigfint[ivalw];

	//figure out what point corresponds to the excess area in remainarea
	c = -remainarea*_signormw/(_Warray[ivalw+1]-_Warray[ivalw]);
	b = _sigofw[ivalw];
	a = (_sigofw[ivalw+1]-_sigofw[ivalw])/2.;
	if(a==0.)
	{
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
}


//______________________________________________________________________________
void Gammagammaleptonpair::picky(double &y)
{
    // This function picks a y given a W 

    double * sigofy;
    double * sgfint;
    sigofy = new double[starlightLimits::MAXYBINS];
    sgfint = new double[starlightLimits::MAXWBINS];
	
    double remainw =0.,remainarea=0.,remainy=0.,a=0.,b=0.,c=0.,sgf=0.,signorm=0.,x=0.;
    int ivalw=0,ivaly=0;

    ivalw=_ivalwd;
    remainw=_remainwd;
    //average over two colums to get y array
    for(int j=0;j<_GGlepInputnumy;j++)
    {
	sigofy[j]=_sigmax[ivalw][j]+(_sigmax[ivalw+1][j]-_sigmax[ivalw][j])*remainw;
    }

    //calculate the unnormalized sgfint
    sgfint[0]=0.;
    for(int j=0;j<_GGlepInputnumy-1;j++)
    {
	sgf = (sigofy[j+1]+sigofy[j])/2.;
	sgfint[j+1]=sgfint[j]+sgf*(_Yarray[j+1]-_Yarray[j]);
    }

    //normalize the sgfint array
    signorm = sgfint[_GGlepInputnumy-1];

    for(int j=0;j<_GGlepInputnumy;j++)
    {
	sgfint[j]=sgfint[j]/signorm;
    }

    //pick a random number
    x = _randy->Rndom();
    //compare x and sgfint to find the ivalue which is just less then the random number x
    for(int i=0;i<_GGlepInputnumy;i++)
    {
	if(x > sgfint[i]) ivaly = i;
    }
    //remainder above ivaly
    remainarea = x - sgfint[ivaly];

    //figure what point corresponds to the leftover area in remainarea
    c = -remainarea*signorm/(_Yarray[ivaly+1]-_Yarray[ivaly]);
    b = sigofy[ivaly];
    a = (sigofy[ivaly+1]-sigofy[ivaly])/2.;
    if(a==0.)
    {
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
void Gammagammaleptonpair::pairMomentum(double w,double y,double &E,double &px,double &py,double &pz,
                                        double &Eb1, double &pxb1, double &pyb1, double &pzb1,
								        double &Eb2, double &pxb2, double &pyb2, double &pzb2, double &t2,
								        double &Egam1, double&pxgam1, double &pygam1, double &pzgam1, double &Q2gam1,
                                        double &Egam2, double&pxgam2, double &pygam2, double &pzgam2, double &Q2gam2)
{
    //this function calculates px,py,pz,and E given w and y

    double anglepp1=0.,anglepp2=0.,pp1=0.,pp2=0.,E1=0.,E2=0.;
    double pt=0.;

    //E1 and E2 are for the 2 photons in the CM frame
    E1 = w*exp(y)/2.;
    E2 = w*exp(-y)/2.;

    //calculate px and py
    //to get x and y components-- phi is random between 0 and 2*pi
    anglepp1 = _randy->Rndom();
    anglepp2 = _randy->Rndom();

    pp1 = pp_1(E1);//this is from beam1
    pp2 = pp_2(E2);//this is from beam2
    px = pp1*cos(2.*starlightConstants::pi*anglepp1)+pp2*cos(2.*starlightConstants::pi*anglepp2);
    py = pp1*sin(2.*starlightConstants::pi*anglepp1)+pp2*sin(2.*starlightConstants::pi*anglepp2);

    //Compute vector sum Pt=Pt1+Pt2 to find pt for the produced particle
    pt = sqrt(px*px+py*py);

    //W is the mass of the produced particle (not necessarily on-mass-shell).Now compute its energy and pz
    E = sqrt(w*w+pt*pt)*cosh(y);
    pz= sqrt(w*w+pt*pt)*sinh(y);

    double E0b1 = _ip->protonEnergy()*_ip->beam1A();
	double px0b1 = 0, px0b2 =0, py0b1=0, py0b2=0;
	double pz0b1 = sqrt(E0b1*E0b1 - _ip->protonMass()*_ip->protonMass()*_ip->beam1A()*_ip->beam1A());
	double E0b2 = _ip->protonEnergy()*_ip->beam2A();
	double pz0b2 = -sqrt(E0b2*E0b2 - _ip->protonMass()*_ip->protonMass()*_ip->beam2A()*_ip->beam2A());
    pxgam1 = pp1*cos(2.*starlightConstants::pi*anglepp1);
    pygam1 = pp1*sin(2.*starlightConstants::pi*anglepp1);
    pxgam2 = pp2*cos(2.*starlightConstants::pi*anglepp2);
    pygam2 = pp2*sin(2.*starlightConstants::pi*anglepp2);
    Egam1 = E1;
    Egam2 = E2;
    
    Eb1 = E0b1 - Egam1;
    pxb1 = px0b1 - pxgam1;
    pyb1 = py0b1 - pygam1;
    pzb1 = sqrt(Eb1*Eb1 - (pxb1*pxb1 + pyb1*pyb1 + _ip->protonMass()*_ip->beam1A()*_ip->protonMass()*_ip->beam1A()));
    pzgam1 = pz0b1 - pzb1;
    Q2gam1 = Egam1*Egam1 - (pxgam1*pxgam1 + pygam1*pygam1 + pzgam1*pzgam1);

    Eb2 = E0b2 - Egam2;
    pxb2 = px0b2 - pxgam2;
    pyb2 = py0b2 - pygam2;
    pzb2 = -sqrt(Eb2*Eb2 - (pxb2*pxb2 + pyb2*pyb2 + _ip->protonMass()*_ip->beam2A()*_ip->protonMass()*_ip->beam2A()));
    pzgam2 = pz0b2 - pzb2;
    Q2gam2 = Egam2*Egam2 - (pxgam2*pxgam2 + pygam2*pygam2 + pzgam2*pzgam2);
    t2= pt*pt;//not sure
}


//______________________________________________________________________________
double Gammagammaleptonpair::pp_1(double E)
{
    // This is for beam 1 
    // returns on random draw from pp(E) distribution
    double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
    double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
    int satisfy =0;
        
    ereds = (E/_GGlepInputGamma_em)*(E/_GGlepInputGamma_em);
    //sqrt(3)*E/gamma_em is p_t where the distribution is a maximum
    Cm = sqrt(3.)*E/_GGlepInputGamma_em;
    //the amplitude of the p_t spectrum at the maximum
    singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
    Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
    //pick a test value pp, and find the amplitude there
    x = _randy->Rndom();
    pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
    singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
    test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

    while(satisfy==0){
	u = _randy->Rndom();
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =_randy->Rndom();
	    pp = 5*starlightConstants::hbarc/_bbs.beam1().nuclearRadius()*x;
	    singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
	    test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
	}
    }

    return pp;
}

//______________________________________________________________________________
double Gammagammaleptonpair::pp_2(double E)
{

    // This is for beam 2 
    //returns on random draw from pp(E) distribution
    double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
    double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
    int satisfy =0;
        
    ereds = (E/_GGlepInputGamma_em)*(E/_GGlepInputGamma_em);
    //sqrt(3)*E/gamma_em is p_t where the distribution is a maximum
    Cm = sqrt(3.)*E/_GGlepInputGamma_em;
    //the amplitude of the p_t spectrum at the maximum
    singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
    Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
    //pick a test value pp, and find the amplitude there
    x = _randy->Rndom(); 
    pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); //Will use nucleus #1 
    singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
    test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

    while(satisfy==0){
	u = _randy->Rndom(); 
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =_randy->Rndom(); 
	    pp = 5*starlightConstants::hbarc/_bbs.beam2().nuclearRadius()*x;
	    singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds); 
	    test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
	}
    }

    return pp;
}



/**
 * @brief creates the daughter particles for the respective two photon decay channel in the CM Frame.
 * 
 * @param ipid [output reference] Sets to the ipid of the daughter particles.
 * @param W [input] Parent Energy in CM-Frame
 * @param px0 [input] Parent x-momentum in CM Frame
 * @param py0 [input] Parent y-momentum in CM Frame
 * @param pz0 [input] Parent z-momentum in CM Frame
 * @param E1 [output reference] Sets to daughter 1 Energy in CM Frame
 * @param px1 [output reference] Sets to daughter 1 x-momentum in CM Frame
 * @param py1 [output reference] Sets to daughter 1 y-momentum in CM Frame
 * @param pz1 [output reference] Sets to daughter 1 z-momentum in CM frame
 * @param E2 [output reference] Sets to daughter 2 Energy in CM Frame
 * @param px2 [output reference] Sets to daughter 2 x-momentum in CM Frame
 * @param py2 [output reference] Sets to daughter 2 y-momentum in CM Frame
 * @param pz2 [output reference] Sets to daughter 2 z-momentum in CM Frame
 * @param iFbadevent [output reference] It sets to 1 for unsuccessful decay ONLY. N.B. It does not reset
 * (but remains at initial value) for successful decays.
 */
void Gammagammaleptonpair::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,
                                    double  W,
                                    double  px0, double  py0, double  pz0,
                                    double& E1, double& px1, double& py1, double& pz1,
                                    double& E2, double& px2, double& py2, double& pz2,
                                    int&    iFbadevent)
{
    //     This routine decays a particle into two particles of mass mdec,
    //     taking spin into account

    double mdec=0.;
    double pmag, anglelep[20001];
    double phi,theta=0.,xtest,Ecm;
    double betax,betay,betaz;
    double hirestheta;
    double hiresw;  //added from JN->needed precision

    mdec = getMass();
    if(W < 2*mdec)
    {
	cout<<" ERROR: W="<<W<<endl;
	iFbadevent = 1;
	return;
    }
    pmag = sqrt(W*W/4. - mdec*mdec);

    //     pick an orientation, based on the spin
    //      phi has a flat distribution in 2*pi
    phi = _randy->Rndom()*2.*starlightConstants::pi;

    //     find theta, the angle between one of the outgoing particles and
    //    the beamline, in the frame of the two photons

    if(getSpin()==0.5){
	//  calculate a table of integrated angle values for leptons
	//  JN05: Go from 1000->20000bins, use double precision for anglelep and thetalep. needed when W >6Gev.
	hiresw = W;

	anglelep[0] = 0.;

	for(int i =1;i<=20000;i++)
	{
	    hirestheta = starlightConstants::pi * double(i) /20000.;

	    //  Added sin(theta) phase space factor (not in thetalep) and changed E to W in thetalep call
	    //  11/9/2000 SRK
	    //  Note that thetalep is form invariant, so it could be called for E, theta_lab just
	    //  as well as W,theta_cm.  Since there is a boost from cm to lab below, the former is fine.

	    anglelep[i] = anglelep[i-1] + thetalep(hiresw,hirestheta)*sin(hirestheta);
	}

	hirestheta = 0.;
	xtest = _randy->Rndom();
	for(int i =1;i<=20000;i++)
	{
	    if(xtest > (anglelep[i]/anglelep[20000]))
                hirestheta = starlightConstants::pi * double(i) / 20000.;
	}
	theta=hirestheta;

    }

    if(getSpin() != 0.5)
        cout<<" This model cannot yet handle this spin value for lepton pairs: "<<getSpin()<<endl; 


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
    //        decay tau to electrons
    //        note that after this routine px1, etc., refer to the electrons
    if(_GGlepInputpidtest == starlightConstants::TAUONDECAY)
        tauDecay(px1,py1,pz1,E1,px2,py2,pz2,E2);

    //     Lorentz transform into the lab frame
    // betax,betay,betaz are the boost of the complete system
    betax = -(px0/Ecm);
    betay = -(py0/Ecm);
    betaz = -(pz0/Ecm);

    transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
    transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);


    if(iFbadevent == 1)
        return;

    // change particle id from that of parent to that of daughters
    // change taoun id into electron id, already switched particles in tau decay
    if(_GGlepInputpidtest == starlightConstants::TAUONDECAY)
        ipid = starlightConstants::ELECTRON;
    //        electrons remain electrons; muons remain muons
    if ( (_GGlepInputpidtest == starlightConstants::ELECTRON) || (_GGlepInputpidtest == starlightConstants::MUON) || 
         (_GGlepInputpidtest == starlightConstants::TAUON) )
        ipid = _GGlepInputpidtest;
}


//______________________________________________________________________________
double Gammagammaleptonpair::thetalep(double W,double theta)
{
    //     This function calculates the cross-section as a function of
    //     angle for a given W and Y, for the production of two muons.
    //     (or tauons)
    //    expression is taken from Brodsky et al. PRD 1971, 1532
    //     equation 5.7
    //     factor that are not dependant on theta are scrapped, so the
    //     the absolute crosssections given by this function are inaccurate
    //     here we are working in the CM frame of the photons and the last
    //     term is 0

    //    real function thetalep (W,theta)

    double moverw=0., W1sq=0.;
    double thetalep_r=0.,denominator=0.;

    W1sq = (W / 2.)*(W/2.);
    moverw = getMass()*getMass() / W1sq;
    denominator = (1. - (1. - moverw)*(cos(theta)*cos(theta)));

    thetalep_r = 2. + 4.*(1.-moverw)*((1.-moverw)*sin(theta)*sin(theta)*cos(theta)*cos(theta) + moverw) / (denominator*denominator);

    return thetalep_r;

}


/**
 * @brief Creates an event with the decay particles inside.
 * 
 * @param ievent 
 * @return [starlightConstants::event] The produced event
 */
starlightConstants::event Gammagammaleptonpair::produceEvent(int &ievent)
{//returns the vector with the decay particles inside.
    starlightConstants::event leptonpair; //This object will store all the tracks for a single event
    double comenergy = 0.;
    double rapidity = 0.;
    double pairE = 0.;
    double pairmomx=0.,pairmomy=0.,pairmomz=0.;
    int iFbadevent=0;
    pairE = pairE + 1;
    starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
	
    double E1=0.,E2=0., px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
    //this function decays particles and writes events to a file
    //zero out the event structure
    leptonpair._numberOfTracks=0;
    for(int i=0;i<4;i++)
    {
	leptonpair.px[i]=0.;
	leptonpair.py[i]=0.;
	leptonpair.pz[i]=0.;
    leptonpair.E[i] =0.;
    leptonpair.mass[i] = 0.;
	leptonpair._fsParticle[i]=starlightConstants::UNKNOWN;
	leptonpair._charge[i]=0;
    }

    pickw(comenergy);

    picky(rapidity);

    //commented to avoid updating the pair momentum to accomodate beam and photon informations
    //hopes it has no unintended consequence.
    //pairMomentum(comenergy,rapidity,pairE,pairmomx,pairmomy,pairmomz);
    twoBodyDecay(ipid,comenergy,pairmomx,pairmomy,pairmomz,E1,px1,py1,pz1,E2,px2,py2,pz2,iFbadevent);//decaying/producing the daughter particles.
    if (iFbadevent==0){
	int q1=0,q2=0; 

	//charges of each daughter is randomly determined.
    double xtest = _randy->Rndom();
	if (xtest<0.5)
	{
	    q1=1;
	    q2=-1;
	}
	else{
	    q1=-1;
	    q2=1;
	}
    double mlepton = getMass();	
	leptonpair._numberOfTracks=2;//leptonpairs are two tracks...

    //storing the properties of each of the lepton pair
	leptonpair.px[0]=px1;
	leptonpair.py[0]=py1;
	leptonpair.pz[0]=pz1;
    leptonpair.E[0] = E1;
    leptonpair.mass[0] = mlepton;
	leptonpair._fsParticle[0]=ipid; 
	leptonpair._charge[0]=q1;

	leptonpair.px[1]=px2;
	leptonpair.py[1]=py2;
	leptonpair.pz[1]=pz2;
    leptonpair.E[1] =E2;
    leptonpair.mass[1] =mlepton;
	leptonpair._fsParticle[1]=ipid;
	leptonpair._charge[1]=q2;

	ievent=ievent+1;
    }

    return leptonpair;
}


/**
 * @brief Creates a upcEvent of the predefined channel containing the decay particles
 * 
 * @param beta The boost vector to transform from CM to Lab Frame. Important for pseudorapidity cuts.
 * @return [upcEvent] The produced event
 */
//upcEvent Gammagammaleptonpair::produceEvent(vector3 beta)
upcXEvent Gammagammaleptonpair::produceEvent(vector3 beta)
{
//returns the vector with the decay particles inside.

   //upcEvent event;
   upcXEvent event;
    //all important  variables are initialized.
   double comenergy = 0.;
   double rapidity = 0.;
   double pairE = 0.;
   double pairmomx=0.,pairmomy=0.,pairmomz=0.;
   int iFbadevent=0;//variable to track successful and unsuccessful decays.
   starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
   
   double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.,E2=0.,E1=0.;
   bool accepted = false;
   double ptCutMin2 = _ptCutMin*_ptCutMin;//used to make pt_Cut comparisons without using square_roots
   double ptCutMax2 = _ptCutMax*_ptCutMax;//used to make pt_Cut comparison without using square_roots

   double Pgam1[4] = {0.0,0.0,0.0,0.0};//Photon from beam1 - Egam1,pxgam1,pygam1,pzgam1
   double Pgam2[4] = {0.0,0.0,0.0,0.0};//Photon from beam2 Egam2,pxgam2,pygam2,pzgam2
   double Pb1[4] = {0.0,0.0,0.0,0.0};//Outgoing beam1 Eb1,pxb1,pyb1,pzb1
   double Pb2[4] = {0.0,0.0,0.0,0.0};//Outgoing beam2 Eb2,pxb2,pyb2,pzb2
   double Q2gam1 =0.,Q2gam2, t=0.;
   do{ 
     
     pickw(comenergy);
     
     picky(rapidity);
     
     pairMomentum(comenergy,rapidity,pairE,pairmomx,pairmomy,pairmomz,
                    Pb1[0],Pb1[1],Pb1[2],Pb1[3],
                    Pb2[0],Pb2[1],Pb2[2],Pb2[3],t,
                    Pgam1[0],Pgam1[1],Pgam1[2],Pgam1[3],Q2gam1,
                    Pgam2[0],Pgam2[1],Pgam2[2],Pgam2[3], Q2gam2);

   
  
     _nmbAttempts++;
     iFbadevent =0;
     twoBodyDecay(ipid,comenergy,pairmomx,pairmomy,pairmomz,E1,px1,py1,pz1,E2,px2,py2,pz2,iFbadevent);//Decaying/producing daughter particle
    
    if(iFbadevent != 0){//checks for successful decay - Important to avoid false accepted counts in the PtCut-EtaCut if-Block.
        continue;//skips all checks to repeat the do-while loop.
    } 
     double pt1chk2 = px1*px1+py1*py1;//used for ptCut comparison without using square roots.
     double pt2chk2 = px2*px2+py2*py2;//as above: computes transverse momentum (squared) for 2nd daughter.
     

     double eta1 = pseudoRapidityLab(px1,py1,pz1,E1,beta);//Determines pseudorapidity of  daughter particle 1 in the laboratory frame
     double eta2 = pseudoRapidityLab(px2,py2,pz2,E2,beta);//similar as above

    //ptCuts and Eta_Cuts are carried out below---->PtCut-EtaCut-if-Block

     if(_ptCutEnabled && !_etaCutEnabled){//Only ptCut is enabled
       if(pt1chk2 > ptCutMin2 && pt1chk2 < ptCutMax2 &&  pt2chk2 > ptCutMin2 && pt2chk2 < ptCutMax2){//if ALL daughter particles fall into ptCut range
	    accepted = true;
	    _nmbAccepted++;
       }
     }
     else if(!_ptCutEnabled && _etaCutEnabled){//only etaCut is enabled
       if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){//BOTH particles fall with EtaCut range.
	 accepted = true;
	 _nmbAccepted++;
       }
     }
     else if(_ptCutEnabled && _etaCutEnabled){//both cuts are enabled
       if(pt1chk2 > ptCutMin2 && pt1chk2 < ptCutMax2 &&  pt2chk2 > ptCutMin2 && pt2chk2 < ptCutMax2){
	 if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){//BOTH particles fall within BOTH ptCut and EtaCut range
	   accepted = true;
	    _nmbAccepted++;
	 }
       }
     }
     else if(!_ptCutEnabled && !_etaCutEnabled) //no cut is enabled
	_nmbAccepted++;
    
   }while((iFbadevent != 0) || ((_ptCutEnabled || _etaCutEnabled) && !accepted));//repeats loop if ptCuts, EtaCuts or successful decay requirements are not satisfied.
   
   if (iFbadevent==0){
     int q1=0,q2=0; 
     
     double xtest = _randy->Rndom();
     if (xtest<0.5)
       {
	 q1=1;
	 q2=-1;
       }
     else{
       q1=-1;
       q2=1;
     }

     // The new stuff
     double mlepton = getMass(); 
     E1 = sqrt( mlepton*mlepton + px1*px1 + py1*py1 + pz1*pz1 ); 
     E2 = sqrt( mlepton*mlepton + px2*px2 + py2*py2 + pz2*pz2 ); 

     starlightParticle particle1(px1, py1, pz1, E1, mlepton, -q1*ipid, q1);
     event.addParticle(particle1);
     
     starlightParticle particle2(px2, py2, pz2, E2, mlepton, -q2*ipid, q2);
     event.addParticle(particle2);

     if(_ip->giveExtraBeamInfo())
     {     
        lorentzVector beam1(Pb1[1],Pb1[2],Pb1[3],Pb1[0]);
        lorentzVector beam2(Pb2[1],Pb2[2],Pb2[3],Pb2[0]);
        double targetEgamma1, targetEgamma2, rap1cm = acosh(_ip->beamLorentzGamma()), cmsEgam1 = Pgam1[0];
        double cmsEgam2 = Pgam2[0], Pzgam1 = Pgam1[3], Pzgam2 = Pgam2[3];
        lorentzVector gamma1(Pgam1[1],Pgam1[2],Pzgam1,cmsEgam1);
        lorentzVector gamma2(Pgam2[1],Pgam2[2],Pzgam2,cmsEgam2);

        
        targetEgamma2 = cmsEgam2*cosh(rap1cm) - Pzgam2*sinh(rap1cm);//beam 1 is target - hence for gamma2
        targetEgamma1 = cmsEgam1*cosh(rap1cm) + Pzgam1*sinh(rap1cm);//beam2 is target - hence for gamma1

        event.addGamma(gamma1,targetEgamma1,Q2gam1);//emmitted by beam1. Order is important - write gamma1 b4 gamma2
        event.addGamma(gamma2,targetEgamma2,Q2gam2);//emmitted by beam2
        event.addOutgoingBeam1(beam1,false);//the order is important. Write beam1 before beam2 so that output can be consistent.
        event.addOutgoingBeam2(beam2,false);//and so that we can associate gamma1 to beam1 and gamma2 to beam2
        event.addVertext(t);
     }//end if

    }
   return event;
}


//______________________________________________________________________________
void Gammagammaleptonpair::calculateTable()
{
    //     this subroutine calculates the tables that are used
    //     elsewhere in the montecarlo
    //     the tauon decay is taken from V-A theory, 1 - 1/3 cos(theta)
    //     the energy of each of the two leptons in tau decay
    //     is calculated using formula 10.35 given
    //     in introduction to elementary particles by D. griffiths
    //     which assmes that the mass of the electron is 0.
    //     the highest energy attainable by the electron in such a system is
    //     .5 * mass of the tau

    //    subroutine calculateTable


    double E,theta;

    _tautolangle[0] = 0.;
    _dgammade[0] = 0.;


    for(int i =1;i<=100;i++)
    {
	//     calculate energy of tau decay
	E = double(i)/100. * .5 * _ip->tauMass();
	_dgammade[i] = _dgammade[i-1] + E*E * (1. - 4.*E/(3.*_ip->tauMass()));

	//     calculate angles for tau
	theta = starlightConstants::pi * double(i) / 100.;
	_tautolangle[i] = _tautolangle[i-1] + (1 + 0.3333333 * cos(theta));
    }


}


//______________________________________________________________________________
void Gammagammaleptonpair::tauDecay(double &px1,double &py1,double &pz1,double &E1,double &px2,double &py2,double &pz2,double &E2)
{
    //     this routine assumes that the tauons decay to electrons and
    //     calculates the directons of the decays

    double Ee1,Ee2,theta1,theta2,phi1,phi2, ran1, ran2 ;
    double pmag1,pex1,pey1,pez1,pex2,pey2,pez2,pmag2;
    double betax,betay,betaz,dir;

    int Tmp_Par=0; // temp variable for the transform function .. kind of starnge - being called with 7 parameter instead of 8

    //     the highest energy attainable by the electron in this system is
    //     .5 * mass of the tau

    //     get two random numbers to compare with


    ran1 = _randy->Rndom()*_dgammade[100];
    ran2 = _randy->Rndom()*_dgammade[100];

    //     compute the energies that correspond to those numbers
    Ee1 = 0.;
    Ee2 = 0.;

    for( int i =1;i<=100;i++)
    {
	if (ran1 > _dgammade[i])
	    Ee1 = double(i) /100. * .5 * getMass();
	if (ran2 > _dgammade[i])
	    Ee2 = double(i) /100. * .5 * getMass();
    }

    //     to find the direction of emmission, first
    //     we determine if the tauons have spin of +1 or -1 along the
    //     direction of the beam line
    dir = 1.;
    if ( _randy->Rndom() < 0.5 )
	dir = -1.;

    //     get two random numbers to compare with
    ran1 = _randy->Rndom()*_tautolangle[100];
    ran2 = _randy->Rndom()*_tautolangle[100];

    //     find the angles corrsponding to those numbers
    theta1 = 0.;
    theta2 = 0.;
    for( int i =1;i<=100;i++)
    {
	if (ran1 > _tautolangle[i]) theta1 = starlightConstants::pi * double(i) /100.;
	if (ran2 > _tautolangle[i]) theta2 = starlightConstants::pi * double(i) /100.;
    }

    //     grab another two random numbers to determine phi's
    phi1 = _randy->Rndom()*2.*starlightConstants::pi;
    phi2 = _randy->Rndom()*2.*starlightConstants::pi;
    //     figure out the momenta of the electron in the frames of the
    //     tauons from which they decayed, that is electron1 is in the
    //     rest frame of tauon1 and e2 is in the rest fram of tau2
    //     again the electrons are assumed to be massless
    pmag1 = Ee1;
    pex1 = cos(phi1)*sin(theta1)*pmag1;
    pey1 = sin(phi1)*sin(theta1)*pmag1;
    pez1 = cos(theta1)*pmag1*dir;
    pmag2 = Ee2;
    pex2 = cos(phi2)*sin(theta2)*pmag2;
    pey2 = sin(phi2)*sin(theta2)*pmag2;
    pez2 = cos(theta2)*pmag2*(-1.*dir);
    //     now Lorentz transform into the frame of each of the particles
    //     do particle one first
    betax = -(px1/E1);
    betay = -(py1/E1);
    betaz = -(pz1/E1);
    //cout<<"2decay betax,pex1= "<<betax<<" "<<pex1<<endl;
    transform (betax,betay,betaz,Ee1,pex1,pey1,pez1,Tmp_Par);
    //     then do particle two
    betax = -(px1/E2);
    betay = -(py1/E2);
    betaz = -(pz1/E2);

    transform (betax,betay,betaz,Ee2,pex2,pey2,pez2, Tmp_Par);
    //     finally dump the electron values into the approriate
    //     variables
    E1 = Ee1;
    E2 = Ee2;
    px1 = pex1;
    px2 = pex2;
    py1 = pey1;
    py2 = pey2;
    pz1 = pez1;
    pz2 = pez2;
}


//______________________________________________________________________________
double Gammagammaleptonpair::getMass()
{
    double leptonmass=0.;
    switch(_GGlepInputpidtest){
    case starlightConstants::ELECTRON:
	leptonmass=_ip->mel();
	break;
    case starlightConstants::MUON:
	leptonmass=_ip->muonMass();
	break;
    case starlightConstants::TAUON:
	leptonmass=_ip->tauMass();
	break;
    case starlightConstants::TAUONDECAY:
	leptonmass=_ip->tauMass();
	break;
    default:
	cout<<"Not a recognized lepton, Gammagammaleptonpair::getmass(), mass = 0."<<endl;
    }

    return leptonmass;
}


//______________________________________________________________________________
double Gammagammaleptonpair::getWidth()
{
    double leptonwidth=0.;
    return leptonwidth;

}


//______________________________________________________________________________
double Gammagammaleptonpair::getSpin()
{
    double leptonspin=0.5;
    return leptonspin;
}
