// eventfilewriter.cpp
/*
 * $Id: eventfilewriter.cpp,v 1.0 2010/07/04   $
 *
 * /author 
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
#include "eventfilewriter.h"
#include "starlightparticlecodes.h"
//______________________________________________________________________________
EventFileWriter::EventFileWriter() : FileWriter()
{
}
//______________________________________________________________________________
EventFileWriter::EventFileWriter(std::string filename) : FileWriter(filename)
{

}
//______________________________________________________________________________
int EventFileWriter::WriteEvent(UPCEvent &event, int eventnumber)
{

   //TODO: Write code also for the pythia part!
   
    int numberoftracks = event.GetParticles()->size();
    int numberofvertices = event.GetVertices()->size();

    // sometimes we don't have tracks due to wrongly picked W , check it
    if(numberoftracks){
      eventnumber++;
      
      fFileStream << "EVENT: " << eventnumber << " " << numberoftracks << " " << 1 << std::endl;
      fFileStream <<"VERTEX: "<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<" "<<1<<" "<<0<<" "<<0<<" "<<numberoftracks<<std::endl;
      
      int ipart = 0;
      std::vector<StarlightParticle>::const_iterator part = (event.GetParticles())->begin();
      
      for (part = event.GetParticles()->begin(); part != event.GetParticles()->end(); part++, ipart++)
	{
	  // TODO:  convert the pdg codes to geant  codes
	  fFileStream << "TRACK: " << " " << StarlightParticleCodes::Jtog((*part).GetCharge() * (*part).GetPdgCode()) <<" "<< (*part).GetPx() << " " << (*part).GetPy()
		      << " "<< (*part).GetPz() << " " <<eventnumber << " " << ipart << " " << 0 << " "
		      << (*part).GetCharge() * (*part).GetPdgCode() <<std::endl;
	}
    }

    return 0;
}
//______________________________________________________________________________
// Output from the pythia based generation...
// if (PythiaOutput==true) {
//
//             for (int t=0;t<(*VD).numberofvertices;t++) {
//                 outputfile <<"VERTEX: "<<(*VD).vertx[t]/10.<<" "<<(*VD).verty[t]/10.<<" "<<(*VD).vertz[t]/10.<<" "<<0.<<" "<<1<<" "<<0<<" "<<0<<" "<<numberoftracks<<endl; //convert from mm to cm for Geant
//             }
//
//             for (int i=0;i<numberoftracks;i++) {
//                 outputfile << "TRACK: " <<" "<<jtog((*VD).charge[i]*(*VD).fsparticle[i])<<" "<<(*VD).px[i]<<" "<<(*VD).py[i]<<" "<<(*VD).pz[i]<<" "<<eventnumber<<" "<<i<<" "<<(*VD).mother1[i]<<" "<<(*VD).mother2[i]<<" "<<(*VD).daughter1[i]<<" "<<(*VD).daughter2[i]<<" "<<(*VD).charge[i]*(*VD).fsparticle[i]<<endl;
//             }//end of track for loop
//         }//if pythiatrue loop
