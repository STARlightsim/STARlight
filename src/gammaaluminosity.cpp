// gammaaluminosity.cpp
/*
 * $Id: gammaaluminosity.cpp,v 1.0 2010/07/04   $
 *
 * /author Joseph Butterwoth
 *
 * $Log: $
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
#include "inputparameters.h"
#include "beambeamsystem.h"
#include "beam.h"
#include "starlightconstants.h"
#include "nucleus.h"
#include "bessel.h"
#include "gammaaluminosity.h"
//______________________________________________________________________________
Gammaaluminosity::Gammaaluminosity(Inputparameters& input, Beambeamsystem& bbsystem)
  :Gammaacrosssection(input,bbsystem),inputgammaa(input)
{
  cout <<"Creating Luminosity Tables."<<endl;
  gammaadifferentialluminosity();
  cout <<"Luminosity Tables created."<<endl;
}
//______________________________________________________________________________
void Gammaaluminosity::gammaadifferentialluminosity()
{

	// double Av,Wgp,cs,cvma;
  double W,dW,dY;
  double Egamma,Y;
  // double t,tmin,tmax;
  double testint,dndWdY;
  double csgA;
  // double ax,bx;
  double C;  

  ofstream wylumfile;
  wylumfile.precision(15);
  
  double  bwnorm,Eth;

  dW = (inputgammaa.getWmax()-inputgammaa.getWmin())/inputgammaa.getnumw();
  dY  = (inputgammaa.getYmax()-(-1.0)*inputgammaa.getYmax())/inputgammaa.getnumy();
    
  // Write the values of W used in the calculation to slight.txt.  
  wylumfile.open("slight.txt");
  wylumfile << getbbs().getBeam1().getZin() <<endl;
  wylumfile << getbbs().getBeam1().getAin() <<endl;
  wylumfile << getbbs().getBeam2().getZin() <<endl;
  wylumfile << getbbs().getBeam2().getAin() <<endl;
  wylumfile << inputgammaa.getgamma_em() <<endl;
  wylumfile << inputgammaa.getWmax() <<endl;
  wylumfile << inputgammaa.getWmin() <<endl;
  wylumfile << inputgammaa.getnumw() <<endl;
  wylumfile << inputgammaa.getYmax() <<endl;
  wylumfile << inputgammaa.getnumy() <<endl;
  wylumfile << inputgammaa.getgg_or_gP() <<endl;
  wylumfile << inputgammaa.getbreakupmode() <<endl;
  wylumfile << inputgammaa.getinterferencemode() <<endl;
  wylumfile << inputgammaa.getinterferencepercent() <<endl;
  wylumfile << inputgammaa.getincoherentorcoherent() <<endl;
  wylumfile << inputgammaa.getincoherentfactor() <<endl;
  wylumfile << inputgammaa.getbford() <<endl;
  wylumfile << inputgammaa.getmaximuminterpt() <<endl;
  wylumfile << inputgammaa.getNPT() <<endl;
  
  //     Normalize the Breit-Wigner Distribution and write values of W to slight.txt
  testint=0.0;
  //Grabbing default value for C in the breit-wigner calculation
  C=getdefaultC();
  for(int i=0;i<=inputgammaa.getnumw()-1;i++){    
    W = inputgammaa.getWmin() + double(i)*dW + 0.5*dW;
    testint = testint + breitwigner(W,C)*dW;
    wylumfile << W << endl;
  }
  bwnorm = 1./testint;
  
  //     Write the values of Y used in the calculation to slight.txt.
  for(int i=0;i<=inputgammaa.getnumy()-1;i++){
    Y = -1.0*inputgammaa.getYmax() + double(i)*dY + 0.5*dY;
    wylumfile << Y << endl;
  }
    
  Eth=0.5*(((inputgammaa.getWmin()+StarlightConstants::mp)*(inputgammaa.getWmin()
							    +StarlightConstants::mp)-StarlightConstants::mp*StarlightConstants::mp)/
	   (inputgammaa.getProtonEnergy()+sqrt(inputgammaa.getProtonEnergy()*
					       inputgammaa.getProtonEnergy()-StarlightConstants::mp*StarlightConstants::mp)));
  
  for(int i=0;i<=inputgammaa.getnumw()-1;i++){

    W = inputgammaa.getWmin() + double(i)*dW + 0.5*dW;
    
    for(int j=0;j<=inputgammaa.getnumy()-1;j++){

      Y = -1.0*inputgammaa.getYmax() + double(j)*dY + 0.5*dY;
      Egamma = 0.5*W*exp(Y);
      
      dndWdY = 0.; 

      if(Egamma > Eth){
        if(Egamma > getMaxPhotonEnergy())Egamma = getMaxPhotonEnergy();
        csgA=getcsgA(Egamma,W);
        dndWdY = Egamma*photonflux(Egamma)*csgA*breitwigner(W,bwnorm);
      }

      wylumfile << dndWdY << endl;
      
    }
  }

  wylumfile.close();
  
  if(inputgammaa.getinterferencemode()==1) 
    pttablegen();
 
  wylumfile.open("slight.txt",ios::app);
  cout << "bwnorm: "<< bwnorm <<endl;
  wylumfile << bwnorm << endl;
  wylumfile.close();
	
}
//______________________________________________________________________________
void Gammaaluminosity::pttablegen()
{
  //  Calculates the pt spectra for VM production with interference
  //  Follows S. Klein and J. Nystrand, Phys. Rev Lett. 84, 2330 (2000).
  //  Written by S. Klein, 8/2002
  
  //  fill in table pttable in one call
  //  Integrate over all y (using the same y values as in table yarray
  //  note that the cross section goes from ymin (<0) to ymax (>0), in numy points
  //  here,  we go from 0 to ymax in (numy/2)+1 points
  //  numy must be even.
  
  //  At each y, calculate the photon energies Egamma1 and Egamma2
  //  and the two photon-A cross sections
  
  //  loop over each p_t entry in the table.
  
  //  Then, loop over b and phi (the angle between the VM \vec{p_t} and \vec{b}
  //  and calculate the cross section at each step.
  //  Put the results in pttable
  
  ofstream wylumfile;
  wylumfile.precision(15);
  wylumfile.open("slight.txt",ios::app);
  
  
  double param1pt[500],param2pt[500];
  double  *ptparam1=param1pt;
  double  *ptparam2=param2pt;
  double dY=0.,Yp=0.,Egamma1=0.,Egamma2=0.,Wgp=0.,cs=0.,cvma=0.,Av=0.,tmin=0.,tmax=0.,ax=0.,bx=0.;
  double csgA=0.,t=0.,sig_ga_1=0.,sig_ga_2=0.,bmax=0.,bmin=0.,db=0.,pt=0.,sum1=0.,b=0.,A1=0.,A2=0.;
  double sumg=0.,theta=0.,amp_i_2=0.,sumint=0.;
  int NGAUSS=0,NBIN=0,NThetaBIN=0;
  
  double xg[16]={.0483076656877383162E0,.144471961582796493E0,
		 .239287362252137075E0, .331868602282127650E0,
		 .421351276130635345E0, .506899908932229390E0,
		 .587715757240762329E0, .663044266930215201E0,
		 .732182118740289680E0, .794483795967942407E0,
		 .849367613732569970E0, .896321155766052124E0,
		 .934906075937739689E0, .964762255587506430E0,
		 .985611511545268335E0, .997263861849481564E0};
  double ag[16]={.0965400885147278006E0, .0956387200792748594E0,
		 .0938443990808045654E0, .0911738786957638847E0,
		 .0876520930044038111E0, .0833119242269467552E0,
		 .0781938957870703065E0, .0723457941088485062E0,
		 .0658222227763618468E0, .0586840934785355471E0,
		 .0509980592623761762E0, .0428358980222266807E0,
		 .0342738629130214331E0, .0253920653092620595E0,
		 .0162743947309056706E0, .00701861000947009660E0};

  NGAUSS=16;

  //Setting input calls to variables/less calls this way.
  double Ymax=inputgammaa.getYmax();
  int numy = inputgammaa.getnumy();
  double Ep = inputgammaa.getProtonEnergy();
  int ibreakup = inputgammaa.getbreakupmode();
  double NPT = inputgammaa.getNPT();
  double gamma_em=inputgammaa.getgamma_em();
  double mass= getchannelmass();
  
  //  loop over y from 0 (not -ymax) to yma
  
  dY=(2.*Ymax)/numy;
  for(int jy=1;jy<=numy/2;jy++){
    Yp=(double(jy)-0.5)*dY;
    
    //  Find the photon energies.  Yp >= 0, so Egamma2 is smaller
    //  Use the vector meson mass for W here - neglect the width
    
    Egamma1 = 0.5*mass*exp(Yp);
    Egamma2 = 0.5*mass*exp(-Yp);
    
    //  Find the sigma(gammaA) for the two directions
    //  Photonuclear Cross Section 1
    //  Gamma-proton CM energy
    
    Wgp=sqrt(2.*Egamma1*(Ep+sqrt(Ep*Ep-StarlightConstants::mp*
				 StarlightConstants::mp))+StarlightConstants::mp*StarlightConstants::mp);
    
    // Calculate V.M.+proton cross section
    
    cs=sqrt(16.*StarlightConstants::pi*getf2o4pi()*getbslope()*
	    StarlightConstants::hbarc*StarlightConstants::hbarc*sigmagp(Wgp)
	    /StarlightConstants::alpha);
    
    // Calculate V.M.+Nucleus cross section
    
    cvma=sigma_A(cs);
    
    // Calculate Av = dsigma/dt(t=0) Note Units: fm**2/Gev**2
    
    Av=(StarlightConstants::alpha*cvma*cvma)/(16.*StarlightConstants::pi
					      *getf2o4pi()*StarlightConstants::hbarc*StarlightConstants::hbarc);
    
    tmin  = ((mass*mass)/(4.*Egamma1*gamma_em)*(mass*mass)/(4.*Egamma1*gamma_em));
    tmax  = tmin + 0.25;
    ax    = 0.5*(tmax-tmin);
    bx    = 0.5*(tmax+tmin);
    csgA  = 0.;
    
    for(int k=0;k<NGAUSS;k++){
      t     = sqrt(ax*xg[k]+bx);
      csgA  = csgA + ag[k]*getbbs().getBeam2().formfactor(t)*getbbs().getBeam2().formfactor(t);
      t     = sqrt(ax*(-xg[k])+bx);
      csgA  = csgA + ag[k]*getbbs().getBeam2().formfactor(t)*getbbs().getBeam2().formfactor(t);
    }
    
    csgA = 0.5*(tmax-tmin)*csgA;
    csgA = Av*csgA;
    sig_ga_1 = csgA;
	   
    // Photonuclear Cross Section 2
    
    Wgp=sqrt(2.*Egamma2*(Ep+sqrt(Ep*Ep-StarlightConstants::mp*
				 StarlightConstants::mp))+StarlightConstants::mp*StarlightConstants::mp);
    
    cs=sqrt(16.*StarlightConstants::pi*getf2o4pi()*getbslope()*
	    StarlightConstants::hbarc*StarlightConstants::hbarc*sigmagp(Wgp)/StarlightConstants::alpha);
    
    cvma=sigma_A(cs);
    
    Av=(StarlightConstants::alpha*cvma*cvma)/(16.*StarlightConstants::pi
					      *getf2o4pi()*StarlightConstants::hbarc*StarlightConstants::hbarc);
    
    tmin  = (((mass*mass)/(4.*Egamma2*gamma_em))*((mass*mass)/(4.*Egamma2*gamma_em)));
    tmax  = tmin + 0.25;
    ax    = 0.5*(tmax-tmin);
    bx    = 0.5*(tmax+tmin);
    csgA  = 0.;
    
    for(int k=0;k<NGAUSS;k++){
      t     = sqrt(ax*xg[k]+bx);
      csgA  = csgA + ag[k]*getbbs().getBeam2().formfactor(t)*getbbs().getBeam2().formfactor(t);
      t     = sqrt(ax*(-xg[k])+bx);
      csgA  = csgA + ag[k]*getbbs().getBeam2().formfactor(t)*getbbs().getBeam2().formfactor(t);
    }
	   
    csgA = 0.5*(tmax-tmin)*csgA;
    csgA = Av*csgA;
    sig_ga_2 = csgA;
    
    //  Set up pttables - they find the reduction in sigma(pt)
    //  due to the nuclear form factors.
    //  Use the vector meson mass for W here - neglect width in
    //  interference calculation
    
    ptparam1=vmsigmapt(mass,Egamma1,ptparam1);
    ptparam2=vmsigmapt(mass,Egamma2,ptparam2);
    
    //  set  bmax according to the smaller photon energy, following flux.f
    
    bmax=bmin+6.*StarlightConstants::hbarc*gamma_em/Egamma2;
    bmin = getbbs().getBeam1().RNuc()+getbbs().getBeam2().RNuc();
    //  if we allow for nuclear breakup, use a slightly smaller bmin
    
    if (ibreakup != 1) 
      bmin=0.95*bmin;
    //  set number of bins to a reasonable number to start
    NBIN = 2000;
    NThetaBIN = 1000;
    db   = (bmax-bmin)/float(NBIN);
    // loop over pt
    for(int i=1;i<=NPT;i++){
      
      pt = (float(i)-0.5)*inputgammaa.getdpt();
      sum1=0.0;
      // loop over b
      for(int j=1;j<=NBIN;j++){
	
	b = bmin + (float(j)-0.5)*db;
	//  nofe is the photon flux function
	A1 = Egamma1*nofe(Egamma1,b)*sig_ga_1*ptparam1[i];
	A2 = Egamma2*nofe(Egamma2,b)*sig_ga_2*ptparam2[i];
	sumg=0.0;
	//  do this as a Gaussian integral, from 0 to pi
	for(int k=0;k<NGAUSS;k++){
	  
	  theta=xg[k]*StarlightConstants::pi;
	  //  allow for a linear sum of interfering and non-interfering amplitudes
	  amp_i_2 = A1 + A2 - 2.*inputgammaa.getinterferencepercent()
	    *sqrt(A1*A2)*cos(pt*b*cos(theta)/StarlightConstants::hbarc);
	  sumg  = sumg+ag[k]*amp_i_2;
	}
	//  this is dn/dpt^2
	//  The factor of 2 is because the theta integral is only from 0 to pi
	sumint=2.*sumg*b*db;
	if (ibreakup > 1)
	  sumint=sumint*getbbs().probabilityofbreakup(b);
	sum1 = sum1 + sumint;
      }
      //  normalization is done in readDiffLum.f
      //  This is d^2sigma/dpt^2; convert to dsigma/dpt
      //CHECK THIS OUT---->      write(20,*)sum1*pt*dpt
      wylumfile << sum1*pt*inputgammaa.getdpt() <<endl;
      //  end of pt loop
    }
    //  end of y loop
  }
  wylumfile.close();
}
//______________________________________________________________________________
double *Gammaaluminosity::vmsigmapt(double W,double Egamma,double *SIGMAPT)
{
  //
  //  This subroutine calculates the effect of the nuclear form factor
  // on the pt spectrum, for use in interference calculations
  // For an interaction with mass W and photon energy Egamma,
  // it calculates the cross section suppression SIGMAPT(PT)
  // as a function of pt.
  // The input pt values come from pttable.inc

  
  double pxmax=0.,pymax=0.,dx=0.,Epom=0.,pt=0.,px0=0.,py0=0.,sum=0.,sumy=0.;
  double py=0.,px=0.,pt1=0.,pt2=0.,f1=0.,f2=0.,q1=0.,q2=0.,norm=0.;
  int NGAUSS =0,Nxbin=0;
  double xg[16]={.0483076656877383162e0,.144471961582796493e0,
		 .239287362252137075e0, .331868602282127650e0,
		 .421351276130635345e0, .506899908932229390e0,
		 .587715757240762329e0, .663044266930215201e0,
		 .732182118740289680e0, .794483795967942407e0,
		 .849367613732569970e0, .896321155766052124e0,
		 .934906075937739689e0, .964762255587506430e0,
		 .985611511545268335e0, .997263861849481564e0};
  double ag[16]={.0965400885147278006e0, .0956387200792748594e0,
		 .0938443990808045654e0, .0911738786957638847e0,
		 .0876520930044038111e0, .0833119242269467552e0,
		 .0781938957870703065e0, .0723457941088485062e0,
		 .0658222227763618468e0, .0586840934785355471e0,
		 .0509980592623761762e0, .0428358980222266807e0,
		 .0342738629130214331e0, .0253920653092620595e0,
		 .0162743947309056706e0, .00701861000947009660e0};
  NGAUSS=16;

  //     >> Initialize
  pxmax = 10.*(StarlightConstants::hbarc/getbbs().getBeam1().RNuc());
  pymax = 10.*(StarlightConstants::hbarc/getbbs().getBeam1().RNuc());
  
  Nxbin = 500;
  
  dx = 2.*pxmax/double(Nxbin);
  Epom   = W*W/(4.*Egamma);
  
  //     >> Loop over total Pt to find distribution
      
      for(int k=1;k<=inputgammaa.getNPT();k++){
	
              pt=inputgammaa.getdpt()*(double(k)-0.5);
              
              px0 = pt;
              py0 = 0.0;
              
              //  For each total Pt, integrate over Pt1, , the photon pt
              //  The pt of the Pomeron  is the difference
              //  pt1 is
              sum=0.;
              for(int i=1;i<=Nxbin;i++){
		
		px = -pxmax + (double(i)-0.5)*dx;
		sumy=0.0;
		for(int j=0;j<NGAUSS;j++){

		  py = 0.5*pymax*xg[j]+0.5*pymax;
		  //  photon pt
		  pt1 = sqrt( px*px + py*py );
		  //  pomeron pt
		  pt2 = sqrt( (px-px0)*(px-px0) + (py-py0)*(py-py0) );
		  q1  = sqrt( ((Egamma/inputgammaa.getgamma_em())*(Egamma/inputgammaa.getgamma_em())) + pt1*pt1 );        
		  q2  = sqrt( ((Epom/inputgammaa.getgamma_em())*(Epom/inputgammaa.getgamma_em()))   + pt2*pt2 );
		  
		  //  photon form factor
		  // add in phase space factor?
		  f1  = (getbbs().getBeam1().formfactor(q1*q1)*getbbs().getBeam1().formfactor(q1*q1)*pt1*pt1)/(q1*q1*q1*q1);
		  
		  //  Pomeron form factor
		  f2  = getbbs().getBeam1().formfactor(q2*q2)*getbbs().getBeam1().formfactor(q2*q2);
		  sumy= sumy + ag[j]*f1*f2;
		  
		  //  now consider other half of py phase space - why is this split?
		  py = 0.5*pymax*(-xg[j])+0.5*pymax;
		  pt1 = sqrt( px*px + py*py );
		  pt2 = sqrt( (px-px0)*(px-px0) + (py-py0)*(py-py0) );
		  q1  = sqrt( ((Egamma/inputgammaa.getgamma_em())*Egamma/inputgammaa.getgamma_em()) + pt1*pt1 );
		  q2  = sqrt( ((Epom/inputgammaa.getgamma_em())*(Epom/inputgammaa.getgamma_em()))   + pt2*pt2 );
		  //  add in phase space factor?
		  f1  = (getbbs().getBeam1().formfactor(q1*q1)*getbbs().getBeam1().formfactor(q1*q1)*pt1*pt1)/(q1*q1*q1*q1);
		  f2  = getbbs().getBeam1().formfactor(q2*q2)*getbbs().getBeam1().formfactor(q2*q2);
		  sumy= sumy + ag[j]*f1*f2;
      
		}
		//         >> This is to normalize the gaussian integration
		sumy = 0.5*pymax*sumy;
		//         >> The 2 is to account for py: 0 -- pymax
		sum  = sum + 2.*sumy*dx;
	      }
	      
	      if(k==1) norm=1./sum;
	      SIGMAPT[k]=sum*norm;
      }
      return (SIGMAPT);
}
//______________________________________________________________________________
double Gammaaluminosity::nofe(double Egamma,double bimp)
{
  //Function for the calculation of the "photon density".
  //nofe=numberofgammas/(energy*area)
  //Assume beta=1.0 and gamma>>1, i.e. neglect the (1/gamma^2)*K0(x) term
  
  double X=0.,nofex=0.,factor1=0.,factor2=0.,factor3=0.;
  
  X = (bimp*Egamma)/(inputgammaa.getgamma_em()*StarlightConstants::hbarc);
  
  if( X <= 0.0) 
    cout<<"In nofe, X= "<<X<<endl;
  
  factor1 = (double(getbbs().getBeam1().getZin()*getbbs().getBeam1().getZin())
	     *StarlightConstants::alpha)/(StarlightConstants::pi*StarlightConstants::pi);

  factor2 = 1./(Egamma*bimp*bimp);
  factor3 = X*X*(bessel::dbesk1(X))*(bessel::dbesk1(X));
  nofex    = factor1*factor2*factor3;
  return nofex;
}
//______________________________________________________________________________
Gammaaluminosity::~Gammaaluminosity()
{

}
