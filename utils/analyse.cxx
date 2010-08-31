#include "analyse.h"

Analyse::Analyse() :
  fInfile("slight.root")
{
  //Constructor
  
  //Creating histograms
  fPtEl = new TH1F("PtEl", "Transverse momentum e+/e-", 100, 0, 2.);
  fPtMu = new TH1F("PtMu", "Transverse momentum mu+/mu-", 100, 0, 2.);
  fPtPi = new TH1F("PtPi", "Transverse momentum pi+/pi-", 100, 0, 2.);
  
  fRapEl = new TH1F("RapEl", "Rapidity e+/e-", 200, -10, 10);
  fRapMu = new TH1F("RapMu", "Rapidity mu+/mu-", 200, -10, 10);
  fRapPi = new TH1F("RapPi", "Rapidity pi+/pi-", 200, -10, 10);

  fInvMassEl = new TH1F("InvMassEl", "Invariant mass", 100, 0, 5);
  fInvMassMu = new TH1F("InvMassMu", "Invariant mass", 100, 0, 5);
  fInvMassPi = new TH1F("InvMassPi", "Invariant mass", 100, 0, 5);
}

Analyse::Analyse(TString infile) :
   fInfile(infile)
{
  //Special constructor

  //Creating histograms
  fPtEl = new TH1F("PtEl", "Transverse momentum e+/e-", 100, 0, 2.);
  fPtMu = new TH1F("PtMu", "Transverse momentum mu+/mu-", 100, 0, 2.);
  fPtPi = new TH1F("PtPi", "Transverse momentum pi+/pi-", 100, 0, 2.);
  
  fRapEl = new TH1F("RapEl", "Rapidity e+/e-", 200, -10, 10);
  fRapMu = new TH1F("RapMu", "Rapidity mu+/mu-", 200, -10, 10);
  fRapPi = new TH1F("RapPi", "Rapidity pi+/pi-", 200, -10, 10);

  fInvMassEl = new TH1F("InvMassEl", "Invariant mass", 100, 0, 5);
  fInvMassMu = new TH1F("InvMassMu", "Invariant mass", 100, 0, 5);
  fInvMassPi = new TH1F("InvMassPi", "Invariant mass", 100, 0, 5);
 
}

Analyse::~Analyse()
{
  //Destructor
  delete fPtEl;
  delete fPtMu;
  delete fPtPi;
 
  delete fRapEl;
  delete fRapMu;
  delete fRapPi;
  
  delete fInvMassEl;
  delete fInvMassMu;
  delete fInvMassPi;
 
}

void Analyse::doAnalysis()
{

  //Doing the analysis

  //Opening the file with the tree
  TFile *f = new TFile(fInfile);
 
  //Getting the tree from the file
  TTree *T = (TTree*)f->Get("outdata");
  
  //Creating a TClonesArray of TParticles
  TClonesArray *arr = new TClonesArray("TParticle");
  //Getting the branch
  T->GetBranch("branch")->SetAutoDelete(kFALSE);
  //Setting the branch address
  T->SetBranchAddress("branch", &arr);
  Int_t nentries = (Int_t)(T->GetEntries());

  //Looping over all events
  for(Int_t ev = 0; ev < nentries; ev++){
    	
    arr->Clear();
    T->GetEntry(ev);
    Int_t ntracks = arr->GetEntriesFast();
    //Array of TLorentzVectors. One vector for each tracks
    TLorentzVector* vecArr[ntracks];
    
    //Looping over the tracks of the event
    for(Int_t tr = 0; tr < ntracks; tr++){
      //Getting a TParticle from the TClonesArray
      TParticle *part = (TParticle*)arr->At(tr);
      //Creating a new TLorentzVector and setting px, py, pz and E.
      vecArr[tr] = new TLorentzVector;
      vecArr[tr]->SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
    }
    
    //Creating a new TLorentzVector, which is the sum of the elements in vecArr
    TLorentzVector sum;
    for(Int_t i = 0; i < ntracks; i++){
      sum += *vecArr[i];
    }
    //Filling the histograms depending on particle type
    TParticle *p = (TParticle*)arr->At(0);
    if(p->GetPdgCode() == 11 || p->GetPdgCode() == -11){
      fPtEl->Fill(sum.Pt());
      fRapEl->Fill(sum.Rapidity());
      fInvMassEl->Fill(sum.M());
    }
    else if(p->GetPdgCode() == 13 || p->GetPdgCode() == -13){
      fPtMu->Fill(sum.Pt());
      fRapMu->Fill(sum.Rapidity());
      fInvMassMu->Fill(sum.M());
    }
    else if(p->GetPdgCode() == 211 || p->GetPdgCode() == -211){
      fPtPi->Fill(sum.Pt());
      fRapPi->Fill(sum.Rapidity());
      fInvMassPi->Fill(sum.M());
    }
  
  }
  //Writing the histograms to file
  TFile file("histograms.root", "RECREATE");
  fPtEl->Write();
  fRapEl->Write();
  fInvMassEl->Write();
  fPtMu->Write();
  fRapMu->Write();
  fInvMassMu->Write();
  fPtPi->Write();
  fRapPi->Write();
  fInvMassPi->Write();
  
  
}
