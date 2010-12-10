#include <iostream>
#include <fstream>
#include <sstream>

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TFile.h"


using namespace std;


void readStarlightAscii(const char* inFileName  = "slight.out",
                        const char* outFileName = "starlight.root")
{
	const double daughterMass = 0.13957018; // charged pion [GeV]

	// create output tree
	TFile* outFile = new TFile(outFileName, "RECREATE");
	if (!outFile) {
		cerr << "    error: could not create output file '" << outFileName << "'" << endl;
		return;
	}
	TTree*          outTree           = new TTree("starlightTree", "starlightTree");
	TLorentzVector* parentParticle    = new TLorentzVector();
	TClonesArray*   daughterParticles = new TClonesArray("TLorentzVector");
	outTree->Branch("parent",    "TLorentzVector", &parentParticle,    32000, -1);
	outTree->Branch("daughters", "TClonesArray",   &daughterParticles, 32000, -1);

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
		for (int i = 0; i < nmbTracks; ++i) {
			// read tracks
			int    particleCode;
			double momentum[3];
			if (!getline(inFile, line))
				break;
			++countLines;
			lineStream.str(line);
			assert(lineStream >> label >> particleCode >> momentum[0] >> momentum[1] >> momentum[2]);
			assert(label == "TRACK:");
			const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
			                      + momentum[2] * momentum[2] + daughterMass * daughterMass);
			new ( (*daughterParticles)[i] ) TLorentzVector(momentum[0], momentum[1], momentum[2], E);
			*parentParticle += *(static_cast<TLorentzVector*>(daughterParticles->At(i)));
		}

		daughterParticles->Compress();
		outTree->Fill();
	}

	outTree->Write();
	if (outFile) {
		outFile->Close();
		delete outFile;
	}
}
