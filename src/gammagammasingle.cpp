// gammagammasingle.cpp
#ifdef ENABLE_PYTHIA

/*
 * $Id: gammagammasingle.cpp,v 1.0 2010/07/04   $
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
#include <vector>
#include "starlightconstants.h"
#include "gammagammasingle.h"
#include "PythiaStarlight.h"


//______________________________________________________________________________
Gammagammasingle::Gammagammasingle(Inputparameters& input,Beambeamsystem& bbsystem):Eventchannel(input,bbsystem)
{
  pythia = new PythiaStarlight();
  std::cout << "Initialising pythia" << std::endl;
  //TODO: remove the env var stupidity...	
  const char* pythiaPath = getenv("PYTHIADIR");
  if(pythiaPath != 0)
    {
      char path[128];
      sprintf(path, "%s/xmldoc\0", pythiaPath);
      pythia->Init(std::string(path));
    }
  else
    {
      std::cerr << "ERROR: Trying to initialise pythia but cannot find the PYTHIADIR environment variable, please set it accordingly" << std::endl;
    }
  

  //Initialize randomgenerator with our seed.
  Randy.SetSeed(input.getseed());
  cout<<"Randy in Single Meson construction: "<<Randy.Rndom()<<endl;
  //Storing inputparameters into protected members for use
  GGsingInputnumw=input.getnumw();
  GGsingInputnumy=input.getnumy();
  GGsingInputpidtest=input.getpidtest();
  GGsingInputGamma_em=input.getgamma_em();
  cout<<"SINGLE MESON pid test: "<<GGsingInputpidtest<<endl;
  //reading in luminosity tables
  read();
  //Now calculating crosssection
  singlecrosssection();
  
}
//______________________________________________________________________________
void Gammagammasingle::singlecrosssection()
{

  //This function carries out a delta function cross-section calculation. For reference, see STAR Note 243, Eq. 8
  //Multiply all Farray[] by f_max
  double sigmasum=0.,remainw=0.;//remainwd=0.;
  int ivalw =0;//ivalwd;
  //calculate the differential cross section and place in the sigma table
  wdelta=getmass();
  for(int i=0;i<GGsingInputnumw;i++){
    for(int j=0;j<GGsingInputnumy;j++){
      //Change 8 to 4 in the following equation, to fix integration of a delta function.
      //Now, it matches the standard literature(cf. Eq. 67 of G.Baur et al., Phys. Rep. 364, 359(2002).
      //STAR Note 243 gives the incorrect 4.
      sigmax[i][j]=(getspin()*2.+1.)*4*StarlightConstants::pi*StarlightConstants::pi*getwidth()/
	(getmass()*getmass()*getmass())*f_max*Farray[i][j]*StarlightConstants::hbarc*StarlightConstants::hbarc/100.;
    }
  }
  //find the index, i,for the value of w just less than the mass because we want to use the value from the sigma table that has w=mass

  for(int i=0;i<GGsingInputnumw;i++){
    if(getmass()>Warray[i]) ivalw=i;
  }

  remainw = (getmass()-Warray[ivalw])/(Warray[ivalw+1]-Warray[ivalw+1]-Warray[ivalw]);
  ivalwd = ivalw;
  remainwd = remainw;
  //if we are interested rho pairs at threshold, the just set sigma to 100nb
  switch(GGsingInputpidtest){
  case StarlightConstants::ZOVERZ03:
    sigmasum =0.;
    for(int j=0;j<GGsingInputnumy-1;j++){
                        sigmasum = sigmasum +2.0*(Yarray[j+1]-Yarray[j])*
			  100.0E-9*(.1/getmass())*((1.-remainw)*f_max*
						   (Farray[ivalw][j]+Farray[ivalw][j])/2.+remainw*f_max*
						   (Farray[ivalw+1][j]+Farray[ivalw+1][j+1])/2.);
    }
    break;
  default:
    //Sum to find the total cross-section
    //The two iss to account for the fact that we integrate over just 1/2 of the rapidity range
    sigmasum =0.;
    for(int j =0;j<GGsingInputnumy-1;j++){
                        sigmasum = sigmasum+2.*
			  (Yarray[j+1]-Yarray[j])*((1.-remainw)*
						   (sigmax[ivalw][j]+sigmax[ivalw][j+1])/2.+remainw*
						   (sigmax[ivalw+1][j]+sigmax[ivalw+1][j+1])/2.);
    }
  }
  cout <<"The total cross-section is: "<<sigmasum<<" barns."<<endl;
  return;
}
//______________________________________________________________________________
void Gammagammasingle::pickw(double &w)
{
  //This function picks a w for the 2-photon calculation. 
  double sgf=0.,signorm=0.,x=0.,remainarea=0.,remainw=0.,a=0.,b=0.,c=0.;
  int ivalw=0;
  
  double * sigofw;
  double * sgfint;
  sigofw = new double[StarlightLimits::MAXWBINS];
  sgfint = new double[StarlightLimits::MAXYBINS];
 
  if(wdelta != 0){
    w=wdelta;
    ivalw=ivalwd;
    remainw=remainwd;
  }
  else{
    //deal with the case where sigma is an array
    //sigofw is simga integrated over y using a linear interpolation
    //sigint is the integral of sgfint, normalized
    
    //integrate sigma down to a function of just w
    for(int i=0;i<GGsingInputnumw;i++){
      sigofw[i]=0.;
      for(int j=0;j<GGsingInputnumy-1;j++){
	sigofw[i] = sigofw[i]+(Yarray[j+1]-Yarray[j])*(sigmax[i][j+1]+sigmax[i][j])/2.;
      }
    }
    //calculate the unnormalized sgfint array
    sgfint[0]=0.;
    for(int i=0;i<GGsingInputnumw-1;i++){
      sgf=(sigofw[i+1]+sigofw[i])*(Warray[i+1]-Warray[i])/2.;
      sgfint[i+1]=sgfint[i]+sgf;
    }
    //normalize sgfint array
    signorm=sgfint[GGsingInputnumw-1];
    
    for(int i=0;i<GGsingInputnumw;i++){
      sgfint[i]=sgfint[i]/signorm;
    }
    //pick a random number
    x = Randy.Rndom();//random()/(RAND_MAX+1.0);
    //compare x and sgfint to find the ivalue which is just less than the random number x
    for(int i=0;i<GGsingInputnumw;i++){
      if(x > sgfint[i]) ivalw=i;
    }
    //remainder above ivalw
    remainarea = x - sgfint[ivalw];
    
    //figure out what point corresponds to the excess area in remainarea
    c = -remainarea*signorm/(Warray[ivalw+1]-Warray[ivalw]);
    b = sigofw[ivalw];
    a = (sigofw[ivalw+1]-sigofw[ivalw])/2.;
    if(a==0.){
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

  delete[] sigofw;
  delete[] sgfint;
}
//______________________________________________________________________________
void Gammagammasingle::picky(double &y)
{
  
  double * sigofy;
  double * sgfint;
  sigofy = new double[StarlightLimits::MAXYBINS];
  sgfint = new double[StarlightLimits::MAXYBINS];
  
  double remainw =0.,remainarea=0.,remainy=0.,a=0.,b=0.,c=0.,sgf=0.,signorm=0.,x=0.;
  int ivalw=0,ivaly=0;
  
  ivalw=ivalwd;
  remainw=remainwd;
  //average over two colums to get y array
  for(int j=0;j<GGsingInputnumy;j++){
    sigofy[j]=sigmax[ivalw][j]+(sigmax[ivalw+1][j]-sigmax[ivalw][j])*remainw;
  }
  //calculate the unnormalized sgfint
  
  sgfint[0]=0.;
  for(int j=0;j<GGsingInputnumy-1;j++){
    sgf = (sigofy[j+1]+sigofy[j])/2.;
    sgfint[j+1]=sgfint[j]+sgf*(Yarray[j+1]-Yarray[j]);
  }
  
  //normalize the sgfint array
  signorm = sgfint[GGsingInputnumy-1];
  
  for(int j=0;j<GGsingInputnumy;j++){
    sgfint[j]=sgfint[j]/signorm;
  }
  //pick a random number
  x = Randy.Rndom();//random()/(RAND_MAX+1.0);
  //compare x and sgfint to find the ivalue which is just less then the random number x
  for(int i=0;i<GGsingInputnumy;i++){
    if(x > sgfint[i]) 
      ivaly = i;
  }
  //remainder above ivaly
  remainarea = x - sgfint[ivaly];
  //figure what point corresponds to the leftover area in remainarea
  c = -remainarea*signorm/(Yarray[ivaly+1]-Yarray[ivaly]);
  b = sigofy[ivaly];
  a = (sigofy[ivaly+1]-sigofy[ivaly])/2.;
  if(a==0.){
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
void Gammagammasingle::parentmomentum(double w,double y,double &E,double &px,double &py,double &pz)
{
  
  //this function calculates px,py,pz,and E given w and y
  double anglepp1=0.,anglepp2=0.,pp1=0.,pp2=0.,E1=0.,E2=0.,signpx=0.,pt=0.;
  
  //E1 and E2 are for the 2 photons in the CM frame
  E1 = w*exp(y)/2.;
  E2 = w*exp(-y)/2.;
  //pz = E1-E2;
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
  if(signpx > 0.5) 
    pz = -pz;	
}
//______________________________________________________________________________
double Gammagammasingle::pp(double E)
{
  //  will probably have to pass in beambeamsys? that way we can do beam1.formfactor(t) or beam2..., careful with the way sergey did it for asymmetry
  //  returns on random draw from pp(E) distribution
      
  double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
  double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
  int satisfy =0;
        
  ereds = (E/GGsingInputGamma_em)*(E/GGsingInputGamma_em);
  Cm = sqrt(3.)*E/GGsingInputGamma_em;
  //the amplitude of the p_t spectrum at the maximum
  singleformfactorCm=bbs.getBeam1().formfactor(Cm*Cm+ereds);
  //Doing this once and then storing it as a double, which we square later...SYMMETRY?using beam1 for now.
  Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(StarlightConstants::pi)*(ereds+Cm*Cm))*(2.*(StarlightConstants::pi)*(ereds+Cm*Cm)));
        
  //pick a test value pp, and find the amplitude there
  x = Randy.Rndom();//random()/(RAND_MAX+1.0);
  pp = x*5.*StarlightConstants::hbarc/bbs.getBeam1().RNuc(); //Will use nucleus #1, there should be two for symmetry//nextline
  singleformfactorpp1=bbs.getBeam1().formfactor(pp*pp+ereds);
  test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*StarlightConstants::pi*(ereds+pp*pp))*(2.*StarlightConstants::pi*(ereds+pp*pp)));

  while(satisfy==0){
    u = Randy.Rndom();//random()/(RAND_MAX+1.0);
    if(u*Coef <= test){
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
void Gammagammasingle::twodecay(StarlightConstants::particle &ipid,double E,double W,double px0,double py0,double pz0,double &px1,double &py1,double &pz1,double &px2,double &py2,double &pz2,int &iFbadevent)
{
  //     This routine decays a particle into two particles of mass mdec,
  //     taking spin into account
  
  double mdec=0.,E1=0.,E2=0.;
  double pmag,ytest=0.;
  double phi,theta,xtest,dndtheta,Ecm;
  double  betax,betay,betaz;
  
  //    set the mass of the daughter particles
  switch(GGsingInputpidtest){ 
  case StarlightConstants::ZOVERZ03:
  case StarlightConstants::F2:	
    mdec = StarlightConstants::mpi;
    break;
  case StarlightConstants::F2PRIME:
    //  decays 50% to K+/K-, 50% to K_0's
    ytest = Randy.Rndom();
    if(ytest >= 0.5){
      mdec = StarlightConstants::mK;
    }
    else{
      mdec = 0.493677;
    }
    break;
  default :
    cout<<"No default mass selected for single photon-photon particle, expect errant results"<<endl;
  }
  
  //Calculating the momentum's magnitude
  //add switch for rho pairs at threshold and everything else.
  switch(GGsingInputpidtest){
  case StarlightConstants::ZOVERZ03:	//the rho pairs produced at threshold
    pmag = sqrt(getmass()*getmass()/4. - mdec*mdec);
    break;
  default :
    if(W < 2*mdec){
      cout<<" ERROR: W="<<W<<endl;
      iFbadevent = 1;
      return;
    }
    pmag = sqrt(W*W/4. - mdec*mdec);
  }
  //     pick an orientation, based on the spin
  //      phi has a flat distribution in 2*pi
  phi = Randy.Rndom()*2.*StarlightConstants::pi; //(random()/(RAND_MAX+1.0))* 2.*StarlightConstants::pi;
  
  //     find theta, the angle between one of the outgoing particles and
  //    the beamline, in the frame of the two photons
  //this will depend on spin, F2,F2' and z/z03 all have spin 2, all other photonphoton-single mesons are handled by jetset
  //Applies to spin2 mesons.
 L300td:
  theta = StarlightConstants::pi*Randy.Rndom();
  xtest = Randy.Rndom();
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

  switch(GGsingInputpidtest){
    //These decay into a pi+ pi- pair
  case StarlightConstants::ZOVERZ03:
  case StarlightConstants::F2:
    ipid=StarlightConstants::PION;
    break;
  case StarlightConstants::F2PRIME:
    if( ytest >= 0.5 )
      {
	//Decays 50/50 into k+ k- or k_s k_l
	ipid=StarlightConstants::KAONCHARGE;	
      }
    else
      {
	ipid=StarlightConstants::KAONNEUTRAL;
      }	
    break;
  default:
    cout<<"Rethink the daughter particles"<<endl;
  }
}
//______________________________________________________________________________
StarlightConstants::event Gammagammasingle::produceevent(int &ievent)
{//returns the vector with the decay particles inside.
//	onedecayparticle single;
	StarlightConstants::event single;
	double comenergy = 0.;
	double rapidity = 0.;
	double parentE = 0.;
	double parentmomx=0.,parentmomy=0.,parentmomz=0.;
        int iFbadevent=0;
	StarlightConstants::particle ipid = StarlightConstants::UNKNOWN;
	double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
	double px3=0.,px4=0.,py3=0.,py4=0.,pz3=0.,pz4=0.;
	double theta=0.,phi=0.;//angles from jetset
//this function decays particles and writes events to a file
	//zeroing out the event structure
                single.numberoftracks=0;
        for(int i=0;i<4;i++)
        {
                single.px[i]=0.;
                single.py[i]=0.;
                single.pz[i]=0.;
                single.fsparticle[i]=StarlightConstants::UNKNOWN;
                single.charge[i]=0;
        }


		double xtest=0.,ztest=0.;


		pickw(comenergy);

                picky(rapidity);

                parentmomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz);

	switch(GGsingInputpidtest){
        case StarlightConstants::ZOVERZ03:
		//Decays into two pairs.
		parentmomx=parentmomx/2.;
		parentmomy=parentmomy/2.;
		parentmomz=parentmomz/2.;
		//Pair #1	
		twodecay(ipid,parentE,comenergy,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
		//Pair #2
		twodecay(ipid,parentE,comenergy,parentmomx,parentmomy,parentmomz,px3,py3,pz3,px4,py4,pz4,iFbadevent);
		//Now add them to vectors to be written out later.
		
		single.numberoftracks=4;//number of tracks per event
		if (iFbadevent==0){
                                                                                                                                                        
                        xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
			ztest = Randy.Rndom();
			//Assigning charges randomly.
                        if (xtest<0.5)
                        {
                                single.charge[0]=1;//q1=1;
                                single.charge[1]=-1;//q2=-1;
                        }
                        else{
                                single.charge[0]=1;//q1=-1;
                                single.charge[1]=-1;//q2=1;
                        }
			if (ztest<0.5)
			{
				single.charge[2]=1;//q3=1;
				single.charge[3]=-1;//q4=-1;
			}
			else{
				single.charge[2]=-1;//q3=-1;
				single.charge[3]=1;//q4=1;
			}
                        //Track #1
                        single.px[0]=px1;
                        single.py[0]=py1;
                        single.pz[0]=pz1;
                        single.fsparticle[0]=ipid;
                        //Track #2                                                                                                                      
                        single.px[1]=px2;
                        single.py[1]=py2;
                        single.pz[1]=pz2;
                        single.fsparticle[1]=ipid;
			//Track #3
			single.px[2]=px3;
                        single.py[2]=py3;
                        single.pz[2]=pz3;
                        single.fsparticle[2]=ipid;
			//Track #4
			single.px[3]=px4;
                        single.py[3]=py4;
                        single.pz[3]=pz4;
                        single.fsparticle[3]=ipid;

                        ievent=ievent+1;
		}	

		break;
	case StarlightConstants::F2:
	case StarlightConstants::F2PRIME:
	        twodecay(ipid,parentE,comenergy,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
		
		single.numberoftracks=2;
		if (iFbadevent==0){

			xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
			if (xtest<0.5)
			{
				single.charge[0]=1;//q1=1;
				single.charge[1]=-1;//q2=-1;
			}
			else{
				single.charge[0]=-1;//q1=-1;
				single.charge[1]=1;//q2=1;
			}	
			//Track #1
			single.px[0]=px1;
			single.py[0]=py1;
			single.pz[0]=pz1;
			single.fsparticle[0]=ipid; 
			//Track #2
                	single.px[1]=px2;
	                single.py[1]=py2;
        	        single.pz[1]=pz2;
                	single.fsparticle[1]=ipid;
			ievent=ievent+1;
		}
		break;
	default://This will be jetset stuff...just need to add link to jetset program
		
		//cout << "Submitting "<<GGsingInputpidtest<<" to jetset to decay."<<endl;
		//calculate theta and phi of the particle being submitted to jetset.
		//We do not need to do this anymore?  Pythia takes pid,px,py,pz,e, and m
		//thephi(comenergy,parentmomx,parentmomy,parentmomz,parentE,theta,phi);
		///lu1ent(0,pid,E,theta,phi);
		//Lets begin to link Pythia8 into Starlight.  Probably best to do this in main.cpp
		//However, we only do pythia here for the time being.  Lets be self contained in case
		//We want to remove it later or whatever other reason.
		/*#include "Pythia.h"
		using namespace Pythia8;
		// Generator; shorthand for event.
  		Pythia pythia;
  		Event& event = pythia.event;
		// Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  		pythia.readString("ProcessLevel:all = off");
  		//When limittau0 is on, only particles with tau0<tau0max are decayed.
  		pythia.readString("ParticleDecays:limitTau0 = on");
  		//Default tau0 max is 10mm/c, we are setting it to 1mm/c
  		pythia.readString("ParticleDecays:tau0Max = 1");
		// Provide printout of initial information.
 		pythia.settings.listChanged();//Useful for earlier stages.
  		// Initialize.
  		pythia.init();*/
		
		//Reseting the event.
		Pythia8::Event &event = pythia->GetPythia()->event;
		event.reset();
		//Adding an event to pythia
		//Adding Particles information (id,status,color,anticolor,px,py,pz,energy,restmass)
		double restmass=0.;
		restmass=getmass();
		double tempx=0.,tempy=0.,tempz=0.;
		event.append(GGsingInputpidtest,1,0,0,parentmomx,parentmomy,parentmomz,parentE,restmass);

		// Generate events. Quit if failure.
    		if (!pythia->GetPythia()->next()) {
      		cout << " Event generation aborted prematurely, owing to error! Gammagammasingle::produceevent\n";
      		break;
    		}
		single.vertx[0]=0.;
		single.verty[0]=0.;
		single.vertz[0]=0.;
		single.numberofvertices=1;
	
		//cout<<"Event statistics that I can output: "<<endl;
        	//cout<<"Event size: "<<event.size()<<endl;
		single.numberoftracks=(event.size()-1);
        	for(int j=1;j<event.size();j++){//skipping event[0], this just outputs the mother as if the beam/system
			int b=j-1;//start counting at 0 for arrays in single
			single.charge[b]=1;//Charge should be returned on the id
			single.px[b]=event[j].px();
			single.py[b]=event[j].py();
			single.pz[b]=event[j].pz();
			single.fsparticle[b]=event[j].id();
			single.mother1[b]=event[j].mother1();
			single.mother2[b]=event[j].mother2();
			single.daughter1[b]=event[j].daughter1();
			single.daughter2[b]=event[j].daughter2();
		//might have to change the type for fsparticle since pythia could return something we do not have defined.
			if(!(event[j].xProd()==tempx&&event[j].yProd()==tempy&&event[j].zProd()==tempz)){
				single.numberofvertices++;
				single.vertx[single.numberofvertices-1]=event[j].xProd();
				single.verty[single.numberofvertices-1]=event[j].yProd();
				single.vertz[single.numberofvertices-1]=event[j].zProd();
				tempx=event[j].xProd();
                                tempy=event[j].yProd();
                                tempz=event[j].zProd();
			}
/*	
        	cout<<"Event["<<j<<"] pid: "<<event[j].id()<<" p: "<<event[j].px()<<" "<<event[j].py()<<" "<<event[j].pz()<<endl;
        	cout<<" E: "<<event[j].e()<<" m: "<<event[j].m()<<endl;
        	cout<<"Vertex(x,y,z) : "<<setprecision(10)<<event[j].xProd()<<" "<<event[j].yProd()<<" "<<event[j].zProd()<<endl;
		cout<<"Mother1,2: "<<event[j].mother1()<<" "<<event[j].mother2()<<" Daughters1,2: "<<event[j].daughter1()<<" "<<event[j].daughter2()<<endl;
		}//forloop
		for(int c=0;c<single.numberofvertices;c++){
			cout<<"vertex: "<<single.vertx[c]<<" "<<single.verty[c]<<" "<<single.vertz[c]<<endl;*/


		}
		ievent=ievent+1;
		



	}
	return single;
}
//______________________________________________________________________________
// fix it ... lost functionality 
//StarlightConstants::event Gammagammasingle::produceevent(int &ievent)
UPCEvent Gammagammasingle::ProduceEvent()
{


	 cout << "NOT IMPLEMENTED!" << endl;
	 
	 return UPCEvent();
	 int ievent;

  //    returns the vector with the decay particles inside.
  //	onedecayparticle single;
  StarlightConstants::event single;
  double comenergy = 0.;
  double rapidity = 0.;
  double parentE = 0.;
  double parentmomx=0.,parentmomy=0.,parentmomz=0.;
  int iFbadevent=0;
  StarlightConstants::particle ipid = StarlightConstants::UNKNOWN;
  double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
  double px3=0.,px4=0.,py3=0.,py4=0.,pz3=0.,pz4=0.;
  double theta=0.,phi=0.;//angles from jetset
  double xtest=0.,ztest=0.;
 

  //this function decays particles and writes events to a file
  //zeroing out the event structure
  single.numberoftracks=0;
  for(int i=0;i<4;i++){
    single.px[i]=0.;
    single.py[i]=0.;
    single.pz[i]=0.;
    single.fsparticle[i]=StarlightConstants::UNKNOWN;
    single.charge[i]=0;
  }
  
  pickw(comenergy);
  picky(rapidity);
  parentmomentum(comenergy,rapidity,parentE,parentmomx,parentmomy,parentmomz);
  
  switch(GGsingInputpidtest){
  case StarlightConstants::ZOVERZ03:
    //Decays into two pairs.
    parentmomx=parentmomx/2.;
    parentmomy=parentmomy/2.;
    parentmomz=parentmomz/2.;
    //Pair #1	
    twodecay(ipid,parentE,comenergy,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
    //Pair #2
    twodecay(ipid,parentE,comenergy,parentmomx,parentmomy,parentmomz,px3,py3,pz3,px4,py4,pz4,iFbadevent);
    //Now add them to vectors to be written out later.
		
    single.numberoftracks=4;//number of tracks per event
    if (iFbadevent==0){
      xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
      ztest = Randy.Rndom();
      //Assigning charges randomly.
      if (xtest<0.5){
	single.charge[0]=1;//q1=1;
	single.charge[1]=-1;//q2=-1;
      }
      else{
	single.charge[0]=1;//q1=-1;
	single.charge[1]=-1;//q2=1;
      }
      if (ztest<0.5){
	single.charge[2]=1;//q3=1;
	single.charge[3]=-1;//q4=-1;
      }
      else{
	single.charge[2]=-1;//q3=-1;
	single.charge[3]=1;//q4=1;
      }
      //Track #1
      single.px[0]=px1;
      single.py[0]=py1;
      single.pz[0]=pz1;
      single.fsparticle[0]=ipid;
      //Track #2                                                                                                                      
      single.px[1]=px2;
      single.py[1]=py2;
      single.pz[1]=pz2;
      single.fsparticle[1]=ipid;
      //Track #3
      single.px[2]=px3;
      single.py[2]=py3;
      single.pz[2]=pz3;
      single.fsparticle[2]=ipid;
      //Track #4
      single.px[3]=px4;
      single.py[3]=py4;
      single.pz[3]=pz4;
      single.fsparticle[3]=ipid;
      
      ievent=ievent+1;
    }	
    
    break;
  case StarlightConstants::F2:
  case StarlightConstants::F2PRIME:
    twodecay(ipid,parentE,comenergy,parentmomx,parentmomy,parentmomz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
    
    single.numberoftracks=2;
    if (iFbadevent==0){
      xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
      if (xtest<0.5){
	single.charge[0]=1;//q1=1;
	single.charge[1]=-1;//q2=-1;
      }
      else{
	single.charge[0]=-1;//q1=-1;
	single.charge[1]=1;//q2=1;
      }	
      //Track #1
      single.px[0]=px1;
      single.py[0]=py1;
      single.pz[0]=pz1;
      single.fsparticle[0]=ipid; 
      //Track #2
      single.px[1]=px2;
      single.py[1]=py2;
      single.pz[1]=pz2;
      single.fsparticle[1]=ipid;
      ievent=ievent+1;
    }
    break;
  default:
    //This will be jetset stuff...just need to add link to jetset program
    
    //cout << "Submitting "<<GGsingInputpidtest<<" to jetset to decay."<<endl;
    //calculate theta and phi of the particle being submitted to jetset.
    //We do not need to do this anymore?  Pythia takes pid,px,py,pz,e, and m
    //thephi(comenergy,parentmomx,parentmomy,parentmomz,parentE,theta,phi);
    ///lu1ent(0,pid,E,theta,phi);
    //Lets begin to link Pythia8 into Starlight.  Probably best to do this in main.cpp
    //However, we only do pythia here for the time being.  Lets be self contained in case
    //We want to remove it later or whatever other reason.
    /*#include "Pythia.h"
      using namespace Pythia8;
      // Generator; shorthand for event.
      Pythia pythia;
      Event& event = pythia.event;
      // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
      pythia.readString("ProcessLevel:all = off");
      //When limittau0 is on, only particles with tau0<tau0max are decayed.
      pythia.readString("ParticleDecays:limitTau0 = on");
      //Default tau0 max is 10mm/c, we are setting it to 1mm/c
      pythia.readString("ParticleDecays:tau0Max = 1");
      // Provide printout of initial information.
      pythia.settings.listChanged();//Useful for earlier stages.
      // Initialize.
      pythia.init();*/
    
    //Reseting the event.
    
    Pythia8::Event &event = pythia->GetPythia()->event;
    event.reset();
    //Adding an event to pythia
    //Adding Particles information (id,status,color,anticolor,px,py,pz,energy,restmass)
    double restmass=0.;
    restmass=getmass();
    double tempx=0.,tempy=0.,tempz=0.;
    event.append(GGsingInputpidtest,1,0,0,parentmomx,parentmomy,parentmomz,parentE,restmass);
    
    // Generate events. Quit if failure.
    
    if (!pythia->GetPythia()->next()) {
      cout << " Event generation aborted prematurely, owing to error! Gammagammasingle::produceevent\n";
      break;
    }
    single.vertx[0]=0.;
    single.verty[0]=0.;
    single.vertz[0]=0.;
    single.numberofvertices=1;
    
    single.numberoftracks=(event.size()-1);
    for(int j=1;j<event.size();j++){//skipping event[0], this just outputs the mother as if the beam/system
      int b=j-1;//start counting at 0 for arrays in single
      single.charge[b]=1;//Charge should be returned on the id
      single.px[b]=event[j].px();
      single.py[b]=event[j].py();
      single.pz[b]=event[j].pz();
      single.fsparticle[b]=event[j].id();
      single.mother1[b]=event[j].mother1();
      single.mother2[b]=event[j].mother2();
      single.daughter1[b]=event[j].daughter1();
      single.daughter2[b]=event[j].daughter2();
      //might have to change the type for fsparticle since pythia could return something we do not have defined.
      if(!(event[j].xProd()==tempx&&event[j].yProd()==tempy&&event[j].zProd()==tempz)){
	single.numberofvertices++;
	single.vertx[single.numberofvertices-1]=event[j].xProd();
	single.verty[single.numberofvertices-1]=event[j].yProd();
	single.vertz[single.numberofvertices-1]=event[j].zProd();
	tempx=event[j].xProd();
	tempy=event[j].yProd();
	tempz=event[j].zProd();
      }
    }
    ievent=ievent+1;
  }
  // return single;
}
//______________________________________________________________________________
void Gammagammasingle::thephi(double W,double px,double py,double pz,double E,double &theta,double &phi)
{

  //     This subroutine calculates angles for channels decayed by jetset.
  //    subroutine thephi(W,px,py,pz,E,theta,phi)
  E = sqrt (W*W+px*px+py*py+pz*pz);

  theta = acos(pz/sqrt(px*px+py*py+pz*pz));
  phi = acos(px/sqrt(px*px+py*py));
  
  if ((px == 0)  && (py == 0))
    phi = 0.;
  if (py < 0)
    phi = 2*StarlightConstants::pi - phi;
}
//______________________________________________________________________________
double Gammagammasingle::getmass()
{
  
  double singlemass=0.;
  switch(GGsingInputpidtest){
  case StarlightConstants::ETA:
    singlemass=0.54730;
    break;
  case StarlightConstants::ETAPRIME:
    singlemass=0.95777;
    break;
  case StarlightConstants::ETAC:
    singlemass=2.980;
    break;
  case StarlightConstants::F0:
    singlemass=0.980;
    break;
  case StarlightConstants::F2:
    singlemass=1.2754;
    break;
  case StarlightConstants::A2:
    singlemass=1.318;
    break;
  case StarlightConstants::F2PRIME:
    singlemass=1.525;
    break;
  case StarlightConstants::ZOVERZ03:
    singlemass=1.540;
    break;
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getmass(), mass = 0."<<endl;
  }
  return singlemass;
}
//______________________________________________________________________________
double Gammagammasingle::getwidth()
{
  
  double singlewidth=0.;
  switch(GGsingInputpidtest){
  case StarlightConstants::ETA:
    singlewidth=1.E-6;
    break;
  case StarlightConstants::ETAPRIME:
    singlewidth=5.E-6;
    break;
  case StarlightConstants::ETAC:
    singlewidth=6.3E-6;
    break;
  case StarlightConstants::F0:
    singlewidth=0.56E-6;
    break;
  case StarlightConstants::F2:
    singlewidth=2.6E-6;
    break;
  case StarlightConstants::A2:
    singlewidth=1.04E-6;
    break;
  case StarlightConstants::F2PRIME:
    singlewidth=0.1E-6;
    break;
  case StarlightConstants::ZOVERZ03:
    singlewidth=0.1E-6;
    break;
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getwidth(), width = 0."<<endl;
  }
  return singlewidth; 
}
//______________________________________________________________________________
double Gammagammasingle::getspin()
{
  double singlespin=0.5;
  switch(GGsingInputpidtest){
  case StarlightConstants::ETA:
    singlespin=0.0;
    break;
  case StarlightConstants::ETAPRIME:
    singlespin=0.0;
    break;
  case StarlightConstants::ETAC:
    singlespin=0.0;
    break;
  case StarlightConstants::F0:
    singlespin=0.0;
    break;
  case StarlightConstants::F2:
    singlespin=2.0;
    break;
  case StarlightConstants::A2:
    singlespin=2.0;
    break;
  case StarlightConstants::F2PRIME:
    singlespin=2.0;
    break;
  case StarlightConstants::ZOVERZ03:
    singlespin=2.0;
    break;
  default:
    cout<<"Not a recognized single particle, Gammagammasingle::getspin(), spin = 0."<<endl;
  }
  return singlespin;
}
//______________________________________________________________________________
Gammagammasingle::~Gammagammasingle()
{

}
#endif
