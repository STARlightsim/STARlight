//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 17 14:39:53 2014 by ROOT version 5.34/14
// from TTree starlightTree/starlightTree
// found on file: starlight.root
//////////////////////////////////////////////////////////

#ifndef AnalyzeTree_h
#define AnalyzeTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalyzeTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *parent;
   TLorentzVector  *beam1;
   TLorentzVector  *beam2;
   double           t;
   double           q2_gamma1;
   double           q2_gamma2;
   double           targetEgamma1;
   double           targetEgamma2;
   TClonesArray    *daughters;
   TClonesArray    *sources;
   TLorentzVector  *target;
   

   // List of branches
   TBranch        *b_parent;   //!
   TBranch        *b_beam1;
   TBranch        *b_beam2;
   TBranch        *b_t;
   TBranch        *b_q2_gamma1;
   TBranch        *b_q2_gamma2;
   TBranch        *b_targetEgamma1;
   TBranch        *b_targetEgamma2;
   TBranch        *b_daughters;   //!
   TBranch        *b_sources;
   TBranch        *b_target;

   AnalyzeTree(TTree *tree=0);
   virtual ~AnalyzeTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalyzeTree_cxx
AnalyzeTree::AnalyzeTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("starlight.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("starlight.root");
      }
      f->GetObject("starlightTree",tree);
   }
   Init(tree);
}

AnalyzeTree::~AnalyzeTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalyzeTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalyzeTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalyzeTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   parent = 0;
   beam1 = 0;
   beam2 = 0;
   daughters = 0;
   q2_gamma1 =0;
   q2_gamma2 =0;
   targetEgamma1 = 0;
   targetEgamma2 = 0;
   target =0;
   sources = 0;
   t =0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("parent", &parent, &b_parent);
   fChain->SetBranchAddress("beam1", &beam1, &b_beam1);
   fChain->SetBranchAddress("beam2", &beam2, &b_beam2);
   fChain->SetBranchAddress("t",     &t,     &b_t);
   fChain->SetBranchAddress("daughters", &daughters, &b_daughters);
   fChain->SetBranchAddress("q2_gamma1", &q2_gamma1, &b_q2_gamma1);
   fChain->SetBranchAddress("q2_gamma2", &q2_gamma2, &b_q2_gamma2);
   fChain->SetBranchAddress("targetEgamma1", &targetEgamma1, &b_targetEgamma1);
   fChain->SetBranchAddress("targetEgamma2", &targetEgamma2, &b_targetEgamma2);
   fChain->SetBranchAddress("target", &target, &b_target);
   fChain->SetBranchAddress("sources", &sources, &b_sources);
   Notify();
}

Bool_t AnalyzeTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalyzeTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalyzeTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalyzeTree_cxx
