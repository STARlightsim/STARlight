// gammagammaleptonpair.cpp
/*
 * $Id: gammagammaleptonpair.cpp,v 1.0 2010/07/04   $
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
 * Nystrand 220710
 * Fixed bug which gave incorrect minv distribution in gammagammaleptonpair.
 * Moved loop over W and Y from pickw to twoleptoncrosssection in
 * gammagammaleptonpair to speed up event generation.
 * Changed to Woods-Saxon radius in twophotonluminosity to be consistent
 * with old starligt.
 */
#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <vector>
#include "starlightconstants.h"
#include "gammagammaleptonpair.h"
//_____________________________________________________________________________
Gammagammaleptonpair::Gammagammaleptonpair(Inputparameters& input,Beambeamsystem& bbsystem):Eventchannel(input,bbsystem)
{

    //Initialize randomgenerator with our seed.
    Randy.SetSeed(input.getseed());
    cout<<"Randy in leptonpair construction: "<<Randy.Rndom()<<endl;
    //Storing inputparameters into protected members for use
    GGlepInputnumw=input.getnumw();
    GGlepInputnumy=input.getnumy();
    GGlepInputpidtest=input.getpidtest();
    GGlepInputGamma_em=input.getgamma_em();
    //Let us read in the luminosity tables
    read();
    //Now we will calculate the crosssection
    twoleptoncrosssection();
    //If it is a tauon, calculate its tables
    if(input.getparticleid()==StarlightConstants::TAUON) tablecalc();
}
//______________________________________________________________________________
void Gammagammaleptonpair::twoleptoncrosssection()
{

    //This function calculates section for 2-particle decay. For reference, see STAR Note 243, Eq. 9.
    //calculate the 2-lepton differential cross section
    //the 100 is to convert to barns
    //the 2 is to account for the fact that we integrate only over one half of the rapidity range
    //Multiply all Farray[] by f_max

    for(int i=0;i<GGlepInputnumw;i++)
    {
	for(int j=0;j<GGlepInputnumy;j++)
	{
	    // sigmax[i][j]=2.*Gammagammaleptonpair::twomuoncrosssection(Warray[i])*f_max*Farray[i][j]/100.;
            sigmax[i][j]=2.*twomuoncrosssection(Warray[i])*f_max*Farray[i][j]/(100.*Warray[i]);
	}
    }
    //calculate the total two-lepton cross section
    double sigmasum =0.;
    for(int i=0;i<GGlepInputnumw-1;i++)
    {
	for(int j=0;j<GGlepInputnumy-1;j++)
	{
	    // sigmasum = sigmasum +2.*((sigmax[i][j]+sigmax[i+1][j]+sigmax[i][j+1]+sigmax[i+1][j+1])/4.*(Yarray[j+1]-Yarray[j])*(Warray[i+1]-Warray[i])/((Warray[i+1]+Warray[i])/2.));
	    // sigmasum = sigmasum +((sigmax[i][j]+sigmax[i+1][j]+sigmax[i][j+1]+sigmax[i+1][j+1])/4.*(Yarray[j+1]-Yarray[j])*(Warray[i+1]-Warray[i])/((Warray[i+1]+Warray[i])/2.));
          sigmasum = sigmasum +(sigmax[i][j]+sigmax[i+1][j]+sigmax[i][j+1]+sigmax[i+1][j+1])/4.*(Yarray[j+1]-Yarray[j])*(Warray[i+1]-Warray[i]);
	}
    }
    cout << "The total "<<GGlepInputpidtest<<" cross-section is: "<<sigmasum<<" barns."<<endl;

    // Do this integration here, once per run rather than once per event (JN 220710) 
    //integrate sigma down to a function of just w

    double sgf=0.;

    for(int i=0;i<ReadInputnumw;i++)
    {
            sigofw[i]=0.;
            for(int j=0;j<ReadInputnumy-1;j++)
            {
                sigofw[i] = sigofw[i]+(Yarray[j+1]-Yarray[j])*(sigmax[i][j+1]+sigmax[i][j])/2.;
            }
    }

    //calculate the unnormalized sgfint array
    sigfint[0]=0.;
    for(int i=0;i<ReadInputnumw-1;i++)
    {
        sgf=(sigofw[i+1]+sigofw[i])*(Warray[i+1]-Warray[i])/2.;
        sigfint[i+1]=sigfint[i]+sgf;
     }

     //normalize sgfint array
     signormw=sigfint[ReadInputnumw-1];
     for(int i=0;i<ReadInputnumw;i++)
     {
          sigfint[i]=sigfint[i]/signormw;
     }



    return;


}
//______________________________________________________________________________
double Gammagammaleptonpair::twomuoncrosssection(double w)
{
    //This function gives the two muon cross section as a function of Y&W. Using the formula given in G.Soff et. al Nuclear Equation of State, part B, 579
    double s=0.,Etest=0.,deltat=0.,xL=0.,sigmuix=0.,alphasquared=0.,hbarcsquared=0.;
    s = w*w;
    Etest = 4.*getmass()*getmass()/s;  
    deltat = s * sqrt(1.-Etest);
    xL = 2.*log(sqrt(s)/(2.*getmass())+sqrt(1./Etest-1.));
    alphasquared = StarlightConstants::alpha*StarlightConstants::alpha;
    hbarcsquared = StarlightConstants::hbarc*StarlightConstants::hbarc;
    sigmuix = 4.*StarlightConstants::pi*alphasquared/s*hbarcsquared*((1+Etest-0.5*Etest*Etest)*xL-(1./s+Etest/s)*deltat);
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

    if(wdelta != 0)
    {
	w=wdelta;
	ivalw=ivalwd;
	remainw=remainwd;
    }
    else{
	//deal with the case where sigma is an array
	//sigofw is simga integrated over y using a linear interpolation
	//sigint is the integral of sgfint, normalized

	//pick a random number
	x = Randy.Rndom();//random()/(RAND_MAX+1.0);
	//compare x and sgfint to find the ivalue which is just less than the random number x
	for(int i=0;i<GGlepInputnumw;i++)
	{
	    if(x > sigfint[i]) ivalw=i;
	}
	//remainder above ivalw
	remainarea = x - sigfint[ivalw];

	//figure out what point corresponds to the excess area in remainarea
	c = -remainarea*signormw/(Warray[ivalw+1]-Warray[ivalw]);
	b = sigofw[ivalw];
	a = (sigofw[ivalw+1]-sigofw[ivalw])/2.;
	if(a==0.)
	{
	    remainw = -c/b;
	}
	else{
	    remainw = (-b+sqrt(b*b-4.*a*c))/(2.*a);
	}
	ivalwd = ivalw;
	remainwd = remainw;
	//calculate the w value
	w = Warray[ivalw]+(Warray[ivalw+1]-Warray[ivalw])*remainw;

    }

}
//______________________________________________________________________________
void Gammagammaleptonpair::picky(double &y)
{

    // This function picks a y given a W 

    double * sigofy;
    double * sgfint;
    sigofy = new double[StarlightLimits::MAXYBINS];
    sgfint = new double[StarlightLimits::MAXWBINS];
	
    double remainw =0.,remainarea=0.,remainy=0.,a=0.,b=0.,c=0.,sgf=0.,signorm=0.,x=0.;
    int ivalw=0,ivaly=0;

    ivalw=ivalwd;
    remainw=remainwd;
    //average over two colums to get y array
    for(int j=0;j<GGlepInputnumy;j++)
    {
	sigofy[j]=sigmax[ivalw][j]+(sigmax[ivalw+1][j]-sigmax[ivalw][j])*remainw;
    }

    //calculate the unnormalized sgfint
    sgfint[0]=0.;
    for(int j=0;j<GGlepInputnumy-1;j++)
    {
	sgf = (sigofy[j+1]+sigofy[j])/2.;
	sgfint[j+1]=sgfint[j]+sgf*(Yarray[j+1]-Yarray[j]);
    }

    //normalize the sgfint array
    signorm = sgfint[GGlepInputnumy-1];

    for(int j=0;j<GGlepInputnumy;j++)
    {
	sgfint[j]=sgfint[j]/signorm;
    }

    //pick a random number
    x = Randy.Rndom();//random()/(RAND_MAX+1.0);
    //compare x and sgfint to find the ivalue which is just less then the random number x
    for(int i=0;i<GGlepInputnumy;i++)
    {
	if(x > sgfint[i]) ivaly = i;
    }
    //remainder above ivaly
    remainarea = x - sgfint[ivaly];

    //figure what point corresponds to the leftover area in remainarea
    c = -remainarea*signorm/(Yarray[ivaly+1]-Yarray[ivaly]);
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
    y = Yarray[ivaly]+(Yarray[ivaly+1]-Yarray[ivaly])*remainy;
    delete[] sigofy;
    delete[] sgfint;

}
//______________________________________________________________________________
void Gammagammaleptonpair::pairmomentum(double w,double y,double &E,double &px,double &py,double &pz)
{

    //this function calculates px,py,pz,and E given w and y

    double anglepp1=0.,anglepp2=0.,pp1=0.,pp2=0.,E1=0.,E2=0.,signpx=0.,pt=0.;

    //E1 and E2 are for the 2 photons in the CM frame
    E1 = w*exp(y)/2.;
    E2 = w*exp(-y)/2.;

    //calculate px and py
    //to get x and y components-- phi is random between 0 and 2*pi
    anglepp1 = Randy.Rndom();//random()/(RAND_MAX+1.0);
    anglepp2 = Randy.Rndom();//random()/(RAND_MAX+1.0);

    pp1 = pp(E1);
    pp2 = pp(E2);
    px = pp1*cos(2.*StarlightConstants::pi*anglepp1)+pp2*cos(2.*StarlightConstants::pi*anglepp2);
    py = pp1*sin(2.*StarlightConstants::pi*anglepp1)+pp2*sin(2.*StarlightConstants::pi*anglepp2);

    //Compute vector sum Pt=Pt1+Pt2 to find pt for the produced particle
    pt = sqrt(px*px+py*py);

    //W is the mass of the produced particle (not necessarily on-mass-shell).Now compute its energy and pz
    E = sqrt(w*w+pt*pt)*cosh(y);
    pz= sqrt(w*w+pt*pt)*sinh(y);
    signpx = Randy.Rndom();//random()/(RAND_MAX+1.0);

    //pick the z direction
    if(signpx > 0.5) pz = -pz;

}
//______________________________________________________________________________
double Gammagammaleptonpair::pp(double E)
{
//will probably have to pass in beambeamsys? that way we can do beam1.formfactor(t) or beam2..., careful with the way sergey did it for asymmetry

//returns on random draw from pp(E) distribution
    double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
    double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
    int satisfy =0;
        
    ereds = (E/GGlepInputGamma_em)*(E/GGlepInputGamma_em);
    //sqrt(3)*E/gamma_em is p_t where the distribution is a maximum
    Cm = sqrt(3.)*E/GGlepInputGamma_em;
    //the amplitude of the p_t spectrum at the maximum
    singleformfactorCm=bbs.getBeam1().formfactor(Cm*Cm+ereds);//Doing this once and then storing it as a double, which we square later...SYMMETRY?using beam1 for now.
    Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(StarlightConstants::pi)*(ereds+Cm*Cm))*(2.*(StarlightConstants::pi)*(ereds+Cm*Cm)));
        
    //pick a test value pp, and find the amplitude there
    x = Randy.Rndom();//random()/(RAND_MAX+1.0);
    pp = x*5.*StarlightConstants::hbarc/bbs.getBeam1().RNuc(); //Will use nucleus #1, there should be two for symmetry//nextline
    singleformfactorpp1=bbs.getBeam1().formfactor(pp*pp+ereds);
    test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*StarlightConstants::pi*(ereds+pp*pp))*(2.*StarlightConstants::pi*(ereds+pp*pp)));

    while(satisfy==0){
	u = Randy.Rndom();//random()/(RAND_MAX+1.0);
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =Randy.Rndom();//random()/(RAND_MAX+1.0);
	    pp = 5*StarlightConstants::hbarc/bbs.getBeam1().RNuc()*x;
	    singleformfactorpp2=bbs.getBeam1().formfactor(pp*pp+ereds);//Symmetry
	    test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*StarlightConstants::pi*(ereds+pp*pp)*2.*StarlightConstants::pi*(ereds+pp*pp));
	}
    }

    return pp;

}
//______________________________________________________________________________
void Gammagammaleptonpair::twodecay(StarlightConstants::particle &ipid,
                                    double  ,  // E (unused)
                                    double  W,
                                    double  px0, double  py0, double  pz0,
                                    double& px1, double& py1, double& pz1,
                                    double& px2, double& py2, double& pz2,
                                    int&    iFbadevent)
{


    //     This routine decays a particle into two particles of mass mdec,
    //     taking spin into account

    double mdec=0.,E1=0.,E2=0.;
    double pmag, anglelep[20001];
    // double ytest=0.,dndtheta;
    double phi,theta,xtest,Ecm;
    double betax,betay,betaz;
    double hirestheta,hirestest,hiresw;  //added from JN->needed precision

    //    set the mass of the daughter particles

    mdec = getmass();
    if(W < 2*mdec)
    {
	cout<<" ERROR: W="<<W<<endl;
	iFbadevent = 1;
	return;
    }
    pmag = sqrt(W*W/4. - mdec*mdec);

    //     pick an orientation, based on the spin
    //      phi has a flat distribution in 2*pi
    phi = Randy.Rndom()*2.*StarlightConstants::pi; //(random()/(RAND_MAX+1.0))* 2.*StarlightConstants::pi;

    //     find theta, the angle between one of the outgoing particles and
    //    the beamline, in the frame of the two photons

    if(getspin()==0.5){
	//  calculate a table of integrated angle values for leptons
	//  JN05: Go from 1000->20000bins, use double precision for anglelep and thetalep. needed when W >6Gev.
	hiresw = W;

	anglelep[0] = 0.;

	for(int i =1;i<=20000;i++)
	{
	    hirestheta = StarlightConstants::pi * double(i) /20000.;

	    //  Added sin(theta) phase space factor (not in thetalep) and changed E to W in thetalep call
	    //  11/9/2000 SRK
	    //  Note that thetalep is form invariant, so it could be called for E, theta_lab just
	    //  as well as W,theta_cm.  Since there is a boost from cm to lab below, the former is fine.

	    anglelep[i] = anglelep[i-1] + thetalep(hiresw,hirestheta)*sin(hirestheta);
	}

	hirestheta = 0.;
	xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
	hirestest = xtest;
	for(int i =1;i<=20000;i++)
	{
	    if(xtest > (anglelep[i]/anglelep[20000]))
                hirestheta = StarlightConstants::pi * double(i) / 20000.;
	}
	theta=hirestheta;

    }

    if(getspin() != 0.5)
        cout<<" This model cannot yet handle this spin value for lepton pairs: "<<getspin()<<endl; 


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
    if(GGlepInputpidtest == StarlightConstants::TAUON)
        taudecay(px1,py1,pz1,E1,px2,py2,pz2,E2);

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
    // change taoun id into electron id, already switched particles in taudecay
    if(GGlepInputpidtest == StarlightConstants::TAUON)
        ipid = StarlightConstants::ELECTRON;
    //        electrons remain electrons; muons remain muons
    if ((GGlepInputpidtest == StarlightConstants::ELECTRON) || (GGlepInputpidtest == StarlightConstants::MUON))
        ipid = GGlepInputpidtest;

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
    moverw = getmass()*getmass() / W1sq;
    denominator = (1. - (1. - moverw)*(cos(theta)*cos(theta)));

    thetalep_r = 2. + 4.*(1.-moverw)*((1.-moverw)*sin(theta)*sin(theta)*cos(theta)*cos(theta) + moverw) / (denominator*denominator);

    return thetalep_r;

}
//______________________________________________________________________________
StarlightConstants::event Gammagammaleptonpair::produceevent(int &ievent)
{//returns the vector with the decay particles inside.
    StarlightConstants::event leptonpair; //This object will store all the tracks for a single event
    double comenergy = 0.;
    double rapidity = 0.;
    double pairE = 0.;
    double pairmomx=0.,pairmomy=0.,pairmomz=0.;
    int iFbadevent=0;
    StarlightConstants::particle ipid = StarlightConstants::UNKNOWN;
	
    double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
//this function decays particles and writes events to a file
    //zero out the event structure
    leptonpair.numberoftracks=0;
    for(int i=0;i<4;i++)
    {
	leptonpair.px[i]=0.;
	leptonpair.py[i]=0.;
	leptonpair.pz[i]=0.;
	leptonpair.fsparticle[i]=StarlightConstants::UNKNOWN;
	leptonpair.charge[i]=0;
    }

    pickw(comenergy);

    picky(rapidity);

    pairmomentum(comenergy,rapidity,pairE,pairmomx,pairmomy,pairmomz);

    twodecay(ipid,pairE,comenergy,pairmomx,pairmomy,pairmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);

    if (iFbadevent==0){
	int q1=0,q2=0; 

	double xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
	if (xtest<0.5)
	{
	    q1=1;
	    q2=-1;
	}
	else{
	    q1=-1;
	    q2=1;
	}	
	leptonpair.numberoftracks=2;//leptonpairs are two tracks...
	leptonpair.px[0]=px1;
	leptonpair.py[0]=py1;
	leptonpair.pz[0]=pz1;
	leptonpair.fsparticle[0]=ipid; 
	leptonpair.charge[0]=q1;

	leptonpair.px[1]=px2;
	leptonpair.py[1]=py2;
	leptonpair.pz[1]=pz2;
	leptonpair.fsparticle[1]=ipid;
	leptonpair.charge[1]=q2;

	ievent=ievent+1;
    }

    return leptonpair;
}
//______________________________________________________________________________
UPCEvent Gammagammaleptonpair::ProduceEvent()
{
//returns the vector with the decay particles inside.

   UPCEvent event;

   double comenergy = 0.;
   double rapidity = 0.;
   double pairE = 0.;
   double pairmomx=0.,pairmomy=0.,pairmomz=0.;
   int iFbadevent=0;
   StarlightConstants::particle ipid = StarlightConstants::UNKNOWN;
   
   double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
   
   
   //this function decays particles and writes events to a file
   //zero out the event structure
   pickw(comenergy);

   picky(rapidity);
   
   pairmomentum(comenergy,rapidity,pairE,pairmomx,pairmomy,pairmomz);
   
   twodecay(ipid,pairE,comenergy,pairmomx,pairmomy,pairmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
   
   if (iFbadevent==0){
     int q1=0,q2=0; 
     
     double xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
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
     StarlightParticle particle1(px1, py1, pz1, StarlightConstants::UNKNOWN, StarlightConstants::UNKNOWN, ipid, q1);
     event.AddParticle(particle1);
     
     StarlightParticle particle2(px2, py2, pz2, StarlightConstants::UNKNOWN, StarlightConstants::UNKNOWN, ipid, q2);
     event.AddParticle(particle2);
     
    }
   return event;
}

//______________________________________________________________________________

void Gammagammaleptonpair::tablecalc()
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

    //    subroutine tablecalc


    double E,theta;

    tautolangle[0] = 0.;
    dgammade[0] = 0.;


    for(int i =1;i<=100;i++)
    {
	//     calculate energy of tau decay
	E = double(i)/100. * .5 * StarlightConstants::mtau;
	dgammade[i] = dgammade[i-1] + E*E * (1. - 4.*E/(3.*StarlightConstants::mtau));

	//     calculate angles for tau
	theta = StarlightConstants::pi * double(i) / 100.;
	tautolangle[i] = tautolangle[i-1] + (1 + 0.3333333 * cos(theta));
    }


}
//______________________________________________________________________________

void Gammagammaleptonpair::taudecay(double &px1,double &py1,double &pz1,double &E1,double &px2,double &py2,double &pz2,double &E2)
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


    ran1 = Randy.Rndom()*dgammade[100];//(random()/(RAND_MAX+1.0)) * dgammade[100];
    ran2 = Randy.Rndom()*dgammade[100];//(random()/(RAND_MAX+1.0)) * dgammade[100];

    //     compute the energies that correspond to those numbers
    Ee1 = 0.;
    Ee2 = 0.;

    for( int i =1;i<=100;i++)
    {
	if (ran1 > dgammade[i])
	    Ee1 = double(i) /100. * .5 * getmass();
	if (ran2 > dgammade[i])
	    Ee2 = double(i) /100. * .5 * getmass();
    }

    //     to find the direction of emmission, first
    //     we determine if the tauons have spin of +1 or -1 along the
    //     direction of the beam line
    dir = 1.;
    if ( Randy.Rndom() < 0.5 )//(random()/(RAND_MAX+1.0)) < 0.5)
	dir = -1.;

    //     get two random numbers to compare with
    ran1 = Randy.Rndom()*tautolangle[100];//(random()/(RAND_MAX+1.0))  * tautolangle[100];
    ran2 = Randy.Rndom()*tautolangle[100];//(random()/(RAND_MAX+1.0))  * tautolangle[100];

    //     find the angles corrsponding to those numbers
    theta1 = 0.;
    theta2 = 0.;
    for( int i =1;i<=100;i++)
    {
	if (ran1 > tautolangle[i]) theta1 = StarlightConstants::pi * double(i) /100.;
	if (ran2 > tautolangle[i]) theta2 = StarlightConstants::pi * double(i) /100.;
    }

    //     grab another two random numbers to determine phi's
    phi1 = Randy.Rndom()*2.*StarlightConstants::pi;// (random()/(RAND_MAX+1.0))* 2. * StarlightConstants::pi;
    phi2 = Randy.Rndom()*2.*StarlightConstants::pi;// (random()/(RAND_MAX+1.0))* 2. * StarlightConstants::pi;
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
double Gammagammaleptonpair::getmass()
{

    double leptonmass=0.;
    switch(GGlepInputpidtest){
    case StarlightConstants::ELECTRON:
	leptonmass=StarlightConstants::mel;
	break;
    case StarlightConstants::MUON:
	leptonmass=StarlightConstants::mmu;
	break;
    case StarlightConstants::TAUON:
	leptonmass=StarlightConstants::mtau;
	break;
    default:
	cout<<"Not a recognized lepton, Gammagammaleptonpair::getmass(), mass = 0."<<endl;
    }

    return leptonmass;
}
//______________________________________________________________________________
double Gammagammaleptonpair::getwidth()
{

    double leptonwidth=0.;
    return leptonwidth;

}
//______________________________________________________________________________
double Gammagammaleptonpair::getspin()
{
    double leptonspin=0.5;

    return leptonspin;
}
//______________________________________________________________________________
Gammagammaleptonpair::~Gammagammaleptonpair()
{

}
