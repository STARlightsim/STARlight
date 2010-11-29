// psifamily.cpp
/*
 * $Id: psifamily.cpp,v 1.0 2010/07/04  $
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
#include "psifamily.h"

Psifamily::Psifamily(Inputparameters& input,Beambeamsystem& bbsystem):Gammaanarrowvm(input,bbsystem)
{
//Defining width and mass...
 

	// switch(input.getpidtest()){
	// 	case StarlightConstants::JPSI:
	// 	  cout << "JPSI goddamnit!" << endl;
        //                 width=0.000091;
        //                 mass=3.09692;
        //         break;
        //         case StarlightConstants::JPSI2S:
        //                 width=0.000337;
        //                 mass=3.686093;
        //         break;
	// 	default: cout<<"This PSI Family Member Has Not Been Defined, Psifamily::Psifamily()"<<endl;
	// }

}

Psifamily::~Psifamily()
{
}

double Psifamily::gettheta(StarlightConstants::particle)
{
//should probably merge the psi fmaily back to the vm stuff.

//This depends on the decay angular distribution
//Valid for J/Psi, Psi(2s)?
//cout<<"Psi family theta"<<endl;
double theta=0.;
double xtest=0.;
double dndtheta=0.;
        L200td:
          theta = StarlightConstants::pi*Randy.Rndom();//random()/(RAND_MAX+1.0);
          xtest = Randy.Rndom();//random()/(RAND_MAX+1.0);
          //  Follow distribution for helicity +/-1
          //  Eq. 19 of J. Breitweg et al., Eur. Phys. J. C2, 247 (1998)//Does Not Apply for J/psi?
          //  SRK 11/14/2000

          dndtheta = sin(theta)*(1.+((cos(theta))*(cos(theta))));
          if(xtest > dndtheta)
            goto L200td;
        return theta;
}

double Psifamily::getdaughtermass(StarlightConstants::particle &ipid)
{
	double ytest=0.,mdec=0.;
	//  decays 50% to e+/e-, 50% to mu+/mu-
        ytest = Randy.Rndom();//random()/(RAND_MAX+1.0);
        if(ytest >= 0.5)
        {
	        mdec = StarlightConstants::mel;
        	ipid = StarlightConstants::ELECTRON;
        }
        else
        {
	        mdec = StarlightConstants::mmu;
        	ipid = StarlightConstants::MUON;
        }
	return mdec;


}
