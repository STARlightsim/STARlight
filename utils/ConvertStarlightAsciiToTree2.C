// This macro reads a starlight output file (default name slight.out) and creates a root file 
// with  TLorentzVectors for the parent and a TClonesArray of TLorentzVectors for the daughter 
// particles.  The output is stored in a root file (default name starlight.root) with one branch 
// labeled "parent" and the other labeled "daughters". Any number of daughter tracks can be 
// accomodated.  Daughter species currently accomodated are:  electrons, muons, charged or neutral 
// pions, charged or neutral kaons, and protons.  
//
// To use this macro, open root and then 
// type .x convertStarlightAsciiToTree.C("inputfilename", "outputfilename")
// 
// The root file produced can be examined in a root TBrowser.
//
// A macro to read this root file and make some standard plots is also provided.  This macro is 
// called AnalyzeTree.cxx; it can be compiled and run with the anaTree.C macro by opening root 
// and typing .x anaTree.C()

#include <iostream>
#include <fstream>
#include <sstream>

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TFile.h"


using namespace std;
double IDtoMass(int particleCode);


void ConvertStarlightAsciiToTree(const char* inFileName  = "slight.out",
                        const char* outFileName = "starlight.root")
{

	// create output tree
	TFile* outFile = new TFile(outFileName, "RECREATE");
	if (!outFile) {
		cerr << "    error: could not create output file '" << outFileName << "'" << endl;
		return;
	}
	TTree*          outTree           = new TTree("starlightTree", "starlightTree");
	TLorentzVector* parentParticle    = new TLorentzVector();
  	TClonesArray*   daughterParticles = new TClonesArray("TLorentzVector");
	TLorentzVector* beam1             = new TLorentzVector();
	TLorentzVector* beam2             = new TLorentzVector();
	double			t  				  = 0;
	double			q2_gamma1		  = 0;
	double			q2_gamma2		  = 0;
	double			targetEgamma1	  = 0;
	double 			targetEgamma2	  = 0;
	TLorentzVector* target			  = new TLorentzVector();
	TClonesArray*   sources			  = new TClonesArray("TLorentzVector");
	
	outTree->Branch("parent",    "TLorentzVector", &parentParticle,    32000, -1);
	outTree->Branch("beam1", 	 "TLorentzVector", &beam1, 			   32000, -1);
	outTree->Branch("beam2", 	 "TLorentzVector", &beam1, 			   32000, -1);
	outTree->Branch("daughters", "TClonesArray",   &daughterParticles, 32000, -1);
	outTree->Branch("t",		 	 &t, 		     "value/D");
	outTree->Branch("q2_gamma1",	 &q2_gamma1, 	 "value/D");
	outTree->Branch("q2_gamma2",	 &q2_gamma2, 	 "value/D");
	outTree->Branch("targetEgamma1", &targetEgamma1, "value/D");
	outTree->Branch("targetEgamma2", &targetEgamma2, "value/D");
	outTree->Branch("target", "TLorentzVector",   &target, 32000, -1);
	outTree->Branch("sources", "TClonesArray",   &sources, 32000, -1);

	ifstream inFile;
	inFile.open(inFileName);
	unsigned int countLines = 0;
	while (inFile.good()) {
		string       line;
		stringstream lineStream;
		
		// read EVENT
		string label;
		int    eventNmb, nmbTracks;
		if (!getline(inFile, line))
			break;
		++countLines;
		lineStream.str(line);
		assert(lineStream >> label >> eventNmb >> nmbTracks);
		if (!(label == "EVENT:"))
			continue;
		
		// read vertex
		if (!getline(inFile, line))
			break;
		++countLines;
		lineStream.str(line);
		assert(lineStream >> label);
		assert(label == "VERTEX:");
			
		*parentParticle = TLorentzVector(0, 0, 0, 0);
		int itrack = 0;
		int igam = 0;
		int itarget =0;
		int isource =0;
		while (itrack < nmbTracks) {
			
			if (!getline(inFile, line))
				break;
			++countLines;
			lineStream.str(line);
			assert(lineStream >> label);
			
			if(label == "GAMMA:"){
				double in_targetEgamma, in_Q2;				
				assert(lineStream>>in_targetEgamma>>in_Q2);
				lineStream.clear();
				if (igam==0) {
					q2_gamma1 =in_Q2; 
					targetEgamma1 =in_targetEgamma;
				}else if(igam==1){
					q2_gamma2 = in_Q2;
					targetEgamma2 = in_targetEgamma;
				}else assert(false);
				igam++;
			}
			else if(label == "t:")
			{
				double t_origin;
				assert(lineStream >> t_origin);
				lineStream.clear();
				t = t_origin;
			}
			else if(label == "TARGET:" || label =="SOURCE:")
			{
				double momentum[4];
				string label2;
				assert(lineStream >> label2 >>momentum[0]>>momentum[1]>>momentum[2]>>momentum[3]);
				lineStream.clear();
				if(label2 == "BEAM1:")
				{
					*beam1 = TLorentzVector(momentum[0], momentum[1], momentum[2], momentum[3]);

				}else if(label2 == "BEAM2:"){
					*beam2 = TLorentzVector(momentum[0], momentum[1], momentum[2], momentum[3]);
				}
				if(label == "SOURCE:"){
					isource++;
					new ( (*sources)[isource] ) TLorentzVector(momentum[0], momentum[1], momentum[2], momentum[3]);			
				}
				else if(label == "TARGET:" && itarget == 0){
					*target = TLorentzVector(momentum[0], momentum[1], momentum[2], momentum[3]);
					itarget++;
				}else assert(false);//You can't have more than one targets
			}			
			else if(label == "TRACK:")// read tracks
			{
				int    particleCode;
				double momentum[3];
				
				assert(lineStream>> particleCode >> momentum[0] >> momentum[1] >> momentum[2]);
				lineStream.clear();
				itrack++;
				Double_t daughterMass = IDtoMass(particleCode);
				if (daughterMass < 0) {break;}
				const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
									+ momentum[2] * momentum[2] + daughterMass * daughterMass);
				new ( (*daughterParticles)[itrack] ) TLorentzVector(momentum[0], momentum[1], momentum[2], E);
				*parentParticle += *(static_cast<TLorentzVector*>(daughterParticles->At(itrack)));

			}
			 
		}
		sources->Compress();
		daughterParticles->Compress();
		outTree->Fill();
	}

	outTree->Write();
	if (outFile) {
		outFile->Close();
		delete outFile;
	}
}

	double IDtoMass(int particleCode){
		double mass;
		if (particleCode == 2 || particleCode==3) {mass = 0.00051099907;} // electron
		else if (particleCode == 5 || particleCode==6) {mass = 0.105658389;} // muon
		else if (particleCode == 8 || particleCode==9)  {mass = 0.13956995;} // charged pion
		else if (particleCode == 7) {mass = 0.1345766;} // neutral pion
		else if (particleCode == 11|| particleCode==12) {mass = 0.493677;} // charged kaon
		else if (particleCode == 10 || particleCode == 16)  {mass = 0.497614;} // neutral kaon
		else if (particleCode == 14)	{mass = 0.93827231;} // proton
		else {
			cout << "unknown daughter particle (ID = " << particleCode << "), please modify code to accomodate" << endl;
			mass = -1.0;
//			exit(0); 
		     } 

		return mass;
	}
