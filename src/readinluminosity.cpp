// readinluminosity.cpp
/*
 * $Id: readinluminosity.cpp,v 1.0 2010/07/04  $
 *
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
 * 
 * $Log: $
 *  Added 18->19 for reading in the luminosity table
 *  Incoherent factor added to table --Joey
 */
#include <iostream>
#include <fstream>

using namespace std;

#include "readinluminosity.h"
#include "starlightconstants.h"
#include "inputparameters.h"
//______________________________________________________________________________
Readinluminosity::Readinluminosity(Inputparameters& input)//:inputread(input)
{

  //storing inputparameters into protected variables for the object to use them
  ReadInputNPT=input.getNPT();
  ReadInputnumy=input.getnumy();
  ReadInputnumw=input.getnumw();
  ReadInputgg_or_gP=input.getgg_or_gP();
  ReadInputinterferencemode=input.getinterferencemode();
  
}
//______________________________________________________________________________
Readinluminosity::~Readinluminosity()
{

}
//______________________________________________________________________________
void Readinluminosity::read()
{
  double dummy[19]; //14//18
  double (*finterm)[StarlightLimits::MAXWBINS]=new double[StarlightLimits::MAXWBINS][StarlightLimits::MAXYBINS];  
  //decreased from 1000*1000; too big! causes fault!
  double fpart =0.;
  double fptsum=0.;
  ifstream wylumfile;

  f_max=0.0;


  wylumfile.open("slight.txt");
  for(int i=0;i < 19;i++){ // was 14; this is to account for sergei's additional parameters ie d-Au//was19
    wylumfile >> dummy[i];
  }
  for(int i=0;i<ReadInputnumw;i++){
    wylumfile >> Warray[i];
  }
  for(int i=0;i<ReadInputnumy;i++){
    wylumfile >> Yarray[i];
  }
  for(int i=0;i<ReadInputnumw;i++){
    for(int j=0;j<ReadInputnumy;j++){
      wylumfile >> Farray[i][j];
      if( Farray[i][j] > f_max ) f_max=Farray[i][j];
    }
  }
  //Normalize farray (JN 010402)
  for(int i=0;i<ReadInputnumw;i++){
    for(int j=0;j<ReadInputnumy;j++){
      Farray[i][j]=Farray[i][j]/f_max;
    }
  }

  if(ReadInputgg_or_gP == 1) goto L1000;
  if(ReadInputinterferencemode == 0) goto L1000;
  // only numy/2 y bins here, from 0 (not -ymax) to ymax
 
  for (int i=0;i<ReadInputnumy/2;i++){
    //fmax=0;
    //we want to convert fptarray to an integral array where fpt(i,j) is near 0, and fpt(j,NPT) ~1. This will facilitate a simple table loookup
    fptsum=0.;
    for(int j=0;j<ReadInputNPT;j++){
      wylumfile >> fpart;
      finterm[i][j] = fpart;
      fptarray[i][j]=0.;
      fptsum=fptsum+fpart;
    }
    //convert array to integral
    fptarray[i][0]=finterm[i][0]/fptsum;
    for(int j=1;j<ReadInputNPT;j++){
      for(int k=0;k<=j;k++){
	fptarray[i][j]=fptarray[i][j]+finterm[i][k];
      }
      fptarray[i][j]=fptarray[i][j]/fptsum;
    }
  }

 L1000:

  wylumfile >> bwnormsave;
  wylumfile.close();
  delete[] finterm;	
  return;
}

