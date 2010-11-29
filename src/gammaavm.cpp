// gammaavm.cpp
/*
 * $Id: gammaavm.cpp,v 1.0 2010/07/04   $
 *
 * /author 
 *
 * $Log: $
 * Added incoherent t2-> pt2 selection.  Following pp selection scheme
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
#include <cassert>

using namespace std;
                                                                                
#include <math.h>
#include "gammaavm.h"
#include "gammaacrosssection.h"
//______________________________________________________________________________
Gammaavectormeson::Gammaavectormeson(Inputparameters& input,Beambeamsystem& bbsystem):Eventchannel(input,bbsystem), phaseSpaceGen(0)  //:Readinluminosity(input),bbs(bbsystem)
{
  VMNPT=input.getNPT();
  VMWmax=input.getWmax();
  VMWmin=input.getWmin();
  VMYmax=input.getYmax();
  VMYmin=-1.*VMYmax;
  VMnumw=input.getnumw();
  VMnumy=input.getnumy();
  VMgamma_em=input.getgamma_em();
  VMinterferencemode=input.getinterferencemode();
  VMbslope=0.;//Will define in wide/narrow constructor
  VMpidtest=input.getpidtest();
  VMptmax=input.getmaximuminterpt();
  VMdpt=input.getdpt();
  Randy.SetSeed(input.getseed());
  VMCoherence=input.getincoherentorcoherent();
  VMCoherenceFactor=input.getincoherentorcoherent();//probably not needed

  switch(VMpidtest){
  case StarlightConstants::RHO:
  case StarlightConstants::RHOZEUS:
    width=0.1507;
    mass=0.7685;
    break;
  case StarlightConstants::FOURPRONG:
	  // create n-body phase-space generator instance
	  phaseSpaceGen = new nBodyPhaseSpaceGen();
	  phaseSpaceGen->setSeed(input.getseed());
	  width = 0.360;
	  mass  = 1.350;
	  break;
  case StarlightConstants::OMEGA:
    width=0.00843;
    mass=0.78194;
    break;
  case StarlightConstants::PHI:
    width=0.00443;
    mass=1.019413;
    break;
  case StarlightConstants::JPSI:
  case StarlightConstants::JPSI_ee:
  case StarlightConstants::JPSI_mumu:
    width=0.000091;
    mass=3.09692;
    break;
  case StarlightConstants::JPSI2S:
  case StarlightConstants::JPSI2S_ee:
  case StarlightConstants::JPSI2S_mumu:
    width=0.000337;
    mass=3.686093;
    break;
  case StarlightConstants::UPSILON:
  case StarlightConstants::UPSILON_ee:
  case StarlightConstants::UPSILON_mumu:
    width=0.00005402;
    mass=9.46030;
    break;
  case StarlightConstants::UPSILON2S:
  case StarlightConstants::UPSILON2S_ee:
  case StarlightConstants::UPSILON2S_mumu:
    width=0.00003198;
    mass=10.02326;
    break;
  case StarlightConstants::UPSILON3S:
  case StarlightConstants::UPSILON3S_ee:
  case StarlightConstants::UPSILON3S_mumu:
    width=0.00002032;
    mass=10.3552;
    break;
  default: cout<<"No proper vector meson defined, gammaavectormeson::gammaavectormeson"<<endl;
  }  
}
//______________________________________________________________________________
Gammaavectormeson::~Gammaavectormeson()
{
	if (phaseSpaceGen)
		delete phaseSpaceGen;
}
//______________________________________________________________________________
void Gammaavectormeson::pickwy(double &W, double &Y)
{
  double dW, dY, xw,xy,xtest;
  int  IW,IY;
  
  dW = (VMWmax-VMWmin)/double(VMnumw);
  dY = (VMYmax-VMYmin)/double(VMnumy);
  
 L201pwy:

  xw = Randy.Rndom();// random()/(RAND_MAX+1.0);
  W = VMWmin + xw*(VMWmax-VMWmin);

  if (W < 2*StarlightConstants::mpi)
    goto L201pwy;
  
  IW = int((W-VMWmin)/dW); //+ 1;
  xy = Randy.Rndom();//random()/(RAND_MAX+1.0);
  Y = VMYmin + xy*(VMYmax-VMYmin);
  IY = int((Y-VMYmin)/dY); //+ 1;
  xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);

  if( xtest > Farray[IW][IY] )
    goto L201pwy;
  
}         
//______________________________________________________________________________                                               
void Gammaavectormeson::twodecay(StarlightConstants::particle &ipid,
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
  mdec=getdaughtermass(ipid);

  //     calculate the magnitude of the momenta
  if(W < 2*mdec){
      cout<<" ERROR: W="<<W<<endl;
      iFbadevent = 1;
      return;
  }
  pmag = sqrt(W*W/4. - mdec*mdec);
  
  //     pick an orientation, based on the spin
  //      phi has a flat distribution in 2*pi
  phi = Randy.Rndom()*2.*StarlightConstants::pi;//(random()/(RAND_MAX+1.0))* 2.*pi;
                                                                                                                
  //     find theta, the angle between one of the outgoing particles and
  //    the beamline, in the frame of the two photons

  theta=gettheta(ipid);
 
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
bool Gammaavectormeson::fourdecay
(StarlightConstants::particle& ipid,
 const double                  ,           // E (unused)
 const double                  W,          // mass of produced particle
 const double*                 p,          // momentum of produced particle; expected to have size 3
 LorentzVector*                decayVecs,  // array of Lorentz vectors of daughter particles; expected to have size 4
 int&                          iFbadevent)
{
	const double parentMass = W;

  // set the mass of the daughter particles
  const double daughterMass = getdaughtermass(ipid);
  if (parentMass < 4 * daughterMass){
	  cout << " ERROR: W=" << parentMass << " GeV too small" << endl;
	  iFbadevent = 1;
	  return false;
  }

  // construct parent four-vector
  const double        parentEnergy = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
                                          + parentMass * parentMass);
  const LorentzVector parentVec(p[0], p[1], p[2], parentEnergy);

  // setup n-body phase-space generator
  assert(phaseSpaceGen);
  static bool firstCall = true;
  if (firstCall) {
	  const double m[4] = {daughterMass, daughterMass, daughterMass, daughterMass};
	  phaseSpaceGen->setDecay(4, m);
	  // estimate maximum phase-space weight
	  phaseSpaceGen->setMaxWeight(1.01 * phaseSpaceGen->estimateMaxWeight(VMWmax));
	  firstCall = false;
  }

  // generate phase-space event
  if (!phaseSpaceGen->generateDecayAccepted(parentVec))
	  return false;

  // set Lorentzvectors of decay daughters
  for (unsigned int i = 0; i < 4; ++i)
	  decayVecs[i] = phaseSpaceGen->daughter(i);
  return true;
}
//______________________________________________________________________________
double Gammaavectormeson::getdaughtermass(StarlightConstants::particle &ipid)
{
  
  //This will return the daughter particles mass, and the final particles outputed id...
  double ytest=0.,mdec=0.;
  
  switch(VMpidtest){
  case StarlightConstants::RHO:
  case StarlightConstants::RHOZEUS:
  case StarlightConstants::FOURPRONG:
  case StarlightConstants::OMEGA:
    mdec = StarlightConstants::mpi;
    ipid = StarlightConstants::PION;
    break;
  case StarlightConstants::PHI:
    mdec = StarlightConstants::mK;
    ipid = StarlightConstants::KAONCHARGE;
    break;
  case StarlightConstants::JPSI:
    mdec = StarlightConstants::mel;
    ipid = StarlightConstants::ELECTRON;
    break; 
  case StarlightConstants::JPSI_ee:
    mdec = StarlightConstants::mel;
    ipid = StarlightConstants::ELECTRON;
    break; 
  case StarlightConstants::JPSI_mumu:
    mdec = StarlightConstants::mmu;
    ipid = StarlightConstants::MUON;
    break; 
  case StarlightConstants::JPSI2S_ee:
    mdec = StarlightConstants::mel;
    ipid = StarlightConstants::ELECTRON;
    break; 
  case StarlightConstants::JPSI2S_mumu:
    mdec = StarlightConstants::mmu;
    ipid = StarlightConstants::MUON;
    break; 

  case StarlightConstants::JPSI2S:
  case StarlightConstants::UPSILON:
  case StarlightConstants::UPSILON2S:
  case StarlightConstants::UPSILON3S:
    //  decays 50% to e+/e-, 50% to mu+/mu-
    ytest = Randy.Rndom();//random()/(RAND_MAX+1.0);
    
    mdec = StarlightConstants::mmu;
    ipid = StarlightConstants::MUON;
    break;
  case StarlightConstants::UPSILON_ee:
  case StarlightConstants::UPSILON2S_ee:
  case StarlightConstants::UPSILON3S_ee:
    mdec = StarlightConstants::mel;
    ipid = StarlightConstants::ELECTRON;
    break;
  case StarlightConstants::UPSILON_mumu:
  case StarlightConstants::UPSILON2S_mumu:
  case StarlightConstants::UPSILON3S_mumu:
    mdec = StarlightConstants::mmu;
    ipid = StarlightConstants::MUON;   
    break;
  default: cout<<"No daughtermass defined, gammaavectormeson::getdaughtermass"<<endl;
  }
  
  return mdec;
}
//______________________________________________________________________________
double Gammaavectormeson::gettheta(StarlightConstants::particle ipid)
{
  //This depends on the decay angular distribution
  //Valid for rho, phi, omega.
  double theta=0.;
  double xtest=0.;
  double dndtheta=0.;
 L200td:
                                                                                                                                                 
  theta = StarlightConstants::pi*Randy.Rndom();//random()/(RAND_MAX+1.0);
  xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
  //  Follow distribution for helicity +/-1
  //  Eq. 19 of J. Breitweg et al., Eur. Phys. J. C2, 247 (1998)
  //  SRK 11/14/2000
  
  switch(ipid){
	  
  case StarlightConstants::MUON:
  case StarlightConstants::ELECTRON:
    //primarily for upsilon/j/psi.  VM->ee/mumu
    dndtheta = sin(theta)*(1.+((cos(theta))*(cos(theta))));
    break;
    
  case StarlightConstants::PION:
  case StarlightConstants::KAONCHARGE:
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
double Gammaavectormeson::getwidth()
{
  return width;
}
//______________________________________________________________________________
double Gammaavectormeson::getmass()
{
  return mass;
}
//______________________________________________________________________________
double Gammaavectormeson::getspin()
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
  double pt,xt,xtest;
  double photon_spectrum;
  double t1,t2;

  dW = (VMWmax-VMWmin)/double(VMnumw);
  dY  = (VMYmax-VMYmin)/double(VMnumy);
  
  //Find Egam,Epom in CM frame
  Egam = 0.5*W*exp(Y);
  Epom = 0.5*W*exp(-Y);
  
 L202vm:
  xt = Randy.Rndom();//random()/(RAND_MAX+1.0);
  pt1  = 0.5*xt;

  tmin = ((Egam/VMgamma_em)*(Egam/VMgamma_em));
  if(tmin > 0.5)
    {
      cout<< " WARNING: tmin= "<<tmin<<endl;
      cout<< " Will pick a new W,Y "<<endl;
      tcheck = 1;
      return;
    }
 
  xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
  t1 = tmin + pt1*pt1;
  photon_spectrum = (bbs.getBeam1().formfactor(t1)*bbs.getBeam1().formfactor(t1)*pt1*pt1*pt1)/(t1*t1);
  
  photon_spectrum = 16.*sqrt(tmin)*photon_spectrum/(3.*sqrt(3.));
                                                                                                                                  
  if( photon_spectrum >  1.0 )
    {
      cout<< "WARNING: photon pt spectrum error "<<"  photon_spectrum="<<photon_spectrum<<endl;
    }
  if( xtest > photon_spectrum )
    goto L202vm;
  phi1 = 2.*StarlightConstants::pi*Randy.Rndom();//random()/(RAND_MAX+1.0);

  if( bbs.getBeam1().getZin()==1 && bbs.getBeam1().getAin()==1) {
    //dsig/dt= exp(-VMbslope*t)
    xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
    t2 = (-1./VMbslope)*log(xtest);
    pt2 = sqrt(1.*t2);
  }
  else{
  L203vm:
    xt = Randy.Rndom(); //random()/(RAND_MAX+1.0);
    //dAu--Sergey
    if(bbs.getBeam2().getZin()==1&&bbs.getBeam2().getAin()==2){
      pt2  = 0.8*xt; //it was 0.5,0.8, 1.0  (Sergey)
    }
    else{
      if(VMCoherence==1) pt2  = 0.5*xt;
    }
    //       >> Check tmin
    tmin = ((Epom/VMgamma_em)*(Epom/VMgamma_em));
	
    if(tmin > 0.5){
      cout<<" WARNING: tmin= "<<tmin<<endl;
      cout<<" Will pick a new W,Y "<<endl;
      tcheck = 1;
      return;
    }
    
    xtest = Randy.Rndom();
    t2 = tmin + pt2*pt2;
    
    if(bbs.getBeam2().getZin()==1&&bbs.getBeam2().getAin()==2){
      if(1.0 < bbs.getBeam2().formfactor(t2)*pt2)  cout <<"POMERON"<<endl;
      if( xtest > bbs.getBeam2().formfactor(t2)*pt2) goto L203vm;
    }
    else{
	if(VMCoherence==1){
      		if(1.0 < bbs.getBeam2().formfactor(t2)*bbs.getBeam2().formfactor(t2)*pt2) cout <<"POMERON:Sergey"<<endl;
      		if( xtest > bbs.getBeam2().formfactor(t2)*bbs.getBeam2().formfactor(t2)*pt2 )
		goto L203vm;
	}
    }//dAu else end

        if(VMCoherence==0 && (!(bbs.getBeam2().getZin()==1&&bbs.getBeam2().getAin()==2))){
		  //Incoherent pt2 selection
                  //dsig/dt= exp(-VMbslope*t)
                  xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
                  t2 = (-1./VMbslope)*log(xtest);//-1./(VMbslope*VMNucleus)?
                  pt2 = sqrt(1.*t2);
                  }

  }//else end from pp
  phi2 = 2.*StarlightConstants::pi*Randy.Rndom();//random()/(RAND_MAX+1.0);
  
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

  // Randomly choose to make pz negative 50% of the time
  if(bbs.getBeam2().getZin()==1&&bbs.getBeam2().getAin()==2){
    pz = -pz;
  }
  else{
    if (Randy.Rndom() >= 0.5) pz = -pz;
  }

}
//______________________________________________________________________________
void Gammaavectormeson::vmpt(double W,double Y,double &E,double &px,double &py, double &pz,
                             int&) // tcheck (unused)
{
  //    This function calculates momentum and energy of vector meson
  //     given W and Y, including interference.
  //     It gets the pt distribution from a lookup table.
  double dW=0.,dY=0.,yleft=0.,yfract=0.,xpt=0.,pt1=0.,ptfract=0.,pt=0.,pt2=0.,theta=0.;
  int IY=0,j=0;
  
  dW = (VMWmax-VMWmin)/double(VMnumw);
  dY  = (VMYmax-VMYmin)/double(VMnumy);
  
  //  Y is already fixed; choose a pt
  //  Follow the approavh in pickwy.f
  // in   fptarray(IY,pt) IY=1 corresponds to Y=0, IY=numy/2 corresponds to +y
  
  IY=int(fabs(Y)/dY);//+1;
  if (IY > (VMnumy/2)-1){
    IY=(VMnumy/2)-1;
  }
  
  yleft=fabs(Y)-(IY)*dY;
  yfract=yleft*dY;
                                                                                                                                  
  xpt=Randy.Rndom(); //random()/(RAND_MAX+1.0);
                                                                                                                                  
  for(j=0;j<VMNPT+1;j++){
    if (xpt < fptarray[IY][j]) goto L60;
  }
 L60:
  
  //  now do linear interpolation - start with extremes
  
  if (j == 0){
    pt1=xpt/fptarray[IY][j]*VMdpt/2.;
    goto L80;
  }
  if (j == VMNPT){
    pt1=(VMptmax-VMdpt/2.) + VMdpt/2.*(xpt-fptarray[IY][j])/(1.-fptarray[IY][j]);
    goto L80;
  }
  
  //  we're in the middle
  
  ptfract=(xpt-fptarray[IY][j])/(fptarray[IY][j+1]-fptarray[IY][j]);
  pt1=(j+1)*VMdpt+ptfract*VMdpt;
  
  //  at an extreme in y?
  if (IY == (VMnumy/2)-1){
    pt=pt1;
    goto L120;
  }
 L80:
  //  interpolate in y repeat for next fractional y bin
                                                                                                                                  
  for(j=0;j<VMNPT+1;j++){
    if (xpt < fptarray[IY+1][j]) goto L90;
  }
 L90:
  
  //  now do linear interpolation - start with extremes
                                                                                                                                  
  if (j == 0){
    pt2=xpt/fptarray[IY+1][j]*VMdpt/2.;
    goto L100;
  }
  if (j == VMNPT){
    pt2=(VMptmax-VMdpt/2.) + VMdpt/2.*(xpt-fptarray[IY+1][j])/(1.-fptarray[IY+1][j]);
    goto L100;
  }
  
  //  we're in the middle
                                                                                                                                  
  ptfract=(xpt-fptarray[IY+1][j])/(fptarray[IY+1][j+1]-fptarray[IY+1][j]);
  pt2=(j+1)*VMdpt+ptfract*VMdpt;
                                                                                                                                  
 L100:
                                                                                                                                  
  //  now interpolate in y
                                                                                                                                  
  pt=yfract*pt2+(1-yfract)*pt1;
                                                                                                                                  
 L120:
                                                                                                                                  
  //  we have a pt
                                                                                                                                  
  theta=2.*StarlightConstants::pi*Randy.Rndom();//(random()/(RAND_MAX+1.0))*2.*pi;
  px=pt*cos(theta);
  py=pt*sin(theta);
                                                                                                                                  
  //      I guess W is the mass of the vector meson (not necessarily
  //      on-mass-shell), and E is the energy
                                                                                                                                  
  E  = sqrt(W*W+pt*pt)*cosh(Y);
  pz = sqrt(W*W+pt*pt)*sinh(Y);
  //      randomly choose to make pz negative 50% of the time
  if(Randy.Rndom()>=0.5) pz = -pz;
}
//______________________________________________________________________________
StarlightConstants::event Gammaavectormeson::produceevent(int&)
{
// Not used; return default event
	return StarlightConstants::event();
}
//------------------------------------------------------------------------------
UPCEvent Gammaavectormeson::ProduceEvent()
{
    // The new event type
    UPCEvent event;

    int iFbadevent=0;
    int tcheck=0;
    StarlightConstants::particle ipid = StarlightConstants::UNKNOWN;

    if (VMpidtest == StarlightConstants::FOURPRONG) {
	    double        comenergy = 0;
	    double        mom[3]    = {0, 0, 0};
	    double        E         = 0;
	    LorentzVector decayVecs[4];
	    do {
		    double rapidity = 0;
		    pickwy(comenergy, rapidity);
		    if (VMinterferencemode == 0)
			    momenta(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
		    else if (VMinterferencemode==1)
			    vmpt(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
	    } while (!fourdecay(ipid, E, comenergy, mom, decayVecs, iFbadevent));
	    if ((iFbadevent == 0) and (tcheck == 0))
		    for (unsigned int i = 0; i < 4; ++i) {
			    StarlightParticle daughter(decayVecs[i].GetPx(),
			                               decayVecs[i].GetPy(),
			                               decayVecs[i].GetPz(),
			                               StarlightConstants::UNKNOWN,  // energy 
			                               StarlightConstants::UNKNOWN,  // mass
			                               ipid,
			                               (i < 2) ? -1 : +1);
			    event.AddParticle(daughter);
		    }
    } else {
	    double comenergy = 0.;
	    double rapidity = 0.;
	    double E = 0.;
	    double momx=0.,momy=0.,momz=0.;

	    double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;

	    pickwy(comenergy,rapidity);

	    if (VMinterferencemode==0){
		    momenta(comenergy,rapidity,E,momx,momy,momz,tcheck);
	    } else if (VMinterferencemode==1){
		    vmpt(comenergy,rapidity,E,momx,momy,momz,tcheck);
	    }

	    twodecay(ipid,E,comenergy,momx,momy,momz,px1,py1,pz1,px2,py2,pz2,iFbadevent);

	    if (iFbadevent==0&&tcheck==0) {
		    int q1=0,q2=0;

		    double xtest = Randy.Rndom(); 
		    if (xtest<0.5)
			    {
				    q1=1;
				    q2=-1;
			    }
		    else {
			    q1=-1;
			    q2=1;
		    }

		    //     The new stuff
		    StarlightParticle particle1(px1, py1, pz1, StarlightConstants::UNKNOWN, StarlightConstants::UNKNOWN, ipid, q1);
		    event.AddParticle(particle1);

		    StarlightParticle particle2(px2, py2, pz2, StarlightConstants::UNKNOWN, StarlightConstants::UNKNOWN, ipid, q2);
		    event.AddParticle(particle2);
		    //     End of the new stuff

	    }
    }

    return event;

}
//______________________________________________________________________________
Gammaanarrowvm::Gammaanarrowvm(Inputparameters& input,Beambeamsystem& bbsystem):Gammaavectormeson(input,bbsystem)
{
  //Need to make sigma object/run it and read in luminosity tables.
  //will just do that outside...of it?
  cout<<"Reading in luminosity tables. Gammaanarrowvm()"<<endl;
  read();
  cout<<"Creating and calculating crosssection. Gammaanarrowvm()"<<endl;
  Narrowresonancesigma sigma(input,bbsystem);
  sigma.crosssectioncalculation(bwnormsave);
  VMbslope=sigma.getbslope(); 
}
//______________________________________________________________________________
Gammaanarrowvm::~Gammaanarrowvm()
{
  
}
//______________________________________________________________________________
Gammaawidevm::Gammaawidevm(Inputparameters& input,Beambeamsystem& bbsystem):Gammaavectormeson(input,bbsystem)
{
  cout<<"Reading in luminosity tables. Gammaawidevm()"<<endl;
  read();
  cout<<"Creating and calculating crosssection. Gammaawidevm()"<<endl;
  Wideresonancesigma sigma(input,bbsystem);
  sigma.crosssectioncalculation(bwnormsave);
  VMbslope=sigma.getbslope();
}
//______________________________________________________________________________
Gammaawidevm::~Gammaawidevm()
{

}


