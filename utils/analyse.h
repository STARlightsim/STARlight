#ifndef ANALYSE_H
#define ANALYSE_H

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH1F.h"


class Analyse 
{
 public:
  Analyse(); //Constructor
  Analyse(TString infile); //Special constructor
  ~Analyse(); //Destructor
  
  void doAnalysis(); //Function doing the analysis
 private:
 
  TH1F *fPtEl;       //Histogram of pt of electrons
  TH1F *fRapEl;      //Histogram of rapidity of electrons
  TH1F *fInvMassEl;  //Histogram of ivariant mass of electrons

  TH1F *fPtMu;       //Histogram of pt of muons
  TH1F *fRapMu;      //Histogram of rapidity of muons
  TH1F *fInvMassMu;  //Histogram of ivariant mass of muons

  TH1F *fPtPi;       //Histogram of pt of pions
  TH1F *fRapPi;      //Histogram of rapidity of pions
  TH1F *fInvMassPi;  //Histogram of ivariant mass of pions
  
  TString fInfile;

};

#endif
