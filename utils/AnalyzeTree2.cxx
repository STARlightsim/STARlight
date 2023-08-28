// This macro reads the starlight.root file produced by convertStarlightAsciiToTree.C, 
// which contains TLorentzVectors for the parents and a TClonesArray of TLorentzVectors 
// for the daughters. 
//
// It creates histograms of the p_T and rapidity of the daughters, as well as the p_T, 
// rapidity and mass of the parent.  While the parents may have been created as the 
// vector sum of any number of daughter particles, this macro currently produces 
// histograms for only the first two daughter particles.  The daughter histograms are
// called D1Pt, D2Pt, D1Rapidity, and D1Rapidity.  Parent histograms are 
// named ParentPt, ParentRapidity, and ParentMass.  The histograms are stored in 
// starlight_histos.root.  
//
// To use this macro, you must first run convertStarlightAsciiToTree.C to produce the 
// starlight.root file.  If needed, modify the file AnalyzeTree.h to call your input file
// (as downloaded, it calls starlight.root).  Then open root and type .x anaTree.C .

#define AnalyzeTree_cxx
#include "AnalyzeTree2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>

void AnalyzeTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MyClass.C
//      Root > MyClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

//  define histos
    TH1D* D1Rapidity     = new TH1D("D1Rapidity","Rapidity of Daughter 1",200, -10, 10);
    TH1D* D1Pt     = new TH1D("D1Pt","Transverse Momentum of Daughter 1",100, 0, 2.);
    TH1D* D2Rapidity     = new TH1D("D2Rapidity","Rapidity of Daughter 2",200, -10, 10);
    TH1D* D2Pt     = new TH1D("D2Pt","Transverse Momentum of Daughter 2",100, 0, 2.);
    TH1D* ParentRapidity     = new TH1D("ParentRapidity","Rapidity of Parent",200, -10, 10);
    TH1D* ParentPt     = new TH1D("ParentPt","Transverse Momentum of Parent",100, 0, 2.);
    TH1D* ParentMass     = new TH1D("ParentMass","Invariant Mass of Parent",100, 0, 5.);
    TH1D* B1Pt = new TH1D("B1Pt","Transverse Momentum of Beam 1",100, 0,2.);
    TH1D* B2Pt = new TH1D("B2Pt","Transverse Momentum of Beam 2",100, 0,2.);
    TH1D* H_t    = new TH1D("t","Squared Transfered Momentum",100,-0.002, 0.1);
    TH1D* Q2Gamma1 = new TH1D("Q2Gam1", "Q^2 for Photon 1",100, -20.,5.);
    TH1D* Q2Gamma2 = new TH1D("Q2Gam2", "Q^2 for Photon 2",100, -20.,5.);
    TH1D* TargetEgam1 = new TH1D("TarEGam1", "Photon Energy for Photon 1 in target frame",100,0,1000000);
    TH1D* TargetEgam2 = new TH1D("TarEGam2", "Photon Engergy for Photon 2 in target frame",100,0,1000000);
    TH1D* targetPt  = new TH1D("TPt", "Transverse momentum of target Beam", 100, 0, 2.);
    TH1D* source1Pt = new TH1D("S1Pt", "Transverse momentum of Source Beam", 100, 0, 2.);
    TH1D* source2Pt = new TH1D("S2Pt", "Transverse Momentum for 2nd Source Beam",100,0,2.);

// Fill histograms
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   TLorentzVector *source1;
   TLorentzVector *source2;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      daughters->Clear();
      sources->Clear();
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      TLorentzVector *D1 = (TLorentzVector*)daughters->At(0);
      TLorentzVector *D2 = (TLorentzVector*)daughters->At(1);

         
      
      source1 = (TLorentzVector*)sources->At(0);

      
      if(sources->GetEntriesFast() > 1){
         //source2 = (TLorentzVector*)sources->At(1);//Useful for two photon cases.
      }else source2 = nullptr;
      
      
      // if desired, acceptance or analysis cuts can be applied here, before filling histograms
         // if (Cut(ientry) < 0) continue;
      D1Rapidity->Fill(D1->Rapidity());
      D2Rapidity->Fill(D2->Rapidity());
      D1Pt->Fill(D1->Pt());
      D2Pt->Fill(D2->Pt());
      ParentRapidity->Fill(parent->Rapidity());
      ParentPt->Fill(parent->Pt());
      ParentMass->Fill(parent->M());
      B1Pt->Fill(beam1->Pt());
      B2Pt->Fill(beam2->Pt());
      H_t->Fill(t);
      Q2Gamma1->Fill(q2_gamma1);
      Q2Gamma2-> Fill(q2_gamma2);
      TargetEgam1-> Fill(targetEgamma1);
      TargetEgam2-> Fill(targetEgamma2);
      if(target) targetPt->Fill(target->Pt());
      if(source1) source1Pt->Fill(source1->Pt());
      if(source2) source2Pt->Fill(source2->Pt());
   }// jentry


// write out the histos
     TFile *histofile = new TFile("starlight_histos.root","recreate");
//     f->cd;
    
    D1Rapidity->Write();
    D2Rapidity->Write();
    D1Pt->Write();
    D2Pt->Write();
    ParentRapidity->Write();
    ParentPt->Write();
    ParentMass->Write();
    B1Pt->Write();
    B2Pt->Write();
    H_t->Write();
    Q2Gamma1->Write();
    Q2Gamma2-> Write();
    TargetEgam1-> Write();
    TargetEgam2-> Write();
    targetPt->Write();
    source1Pt->Write();
    source2Pt->Write();

     histofile->Save();
     histofile->Close();

}
