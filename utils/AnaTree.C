// this macro compiles and runs AnalyzeTree.cxx, which takes as input the 
// starlight.root file produced by convertStarlightAsciiToTree.cxx
// output histograms are stored in starlight_histos.root 
//
void AnaTree(){
gROOT->ProcessLine(".L AnalyzeTree.cxx");
gROOT->ProcessLine("AnalyzeTree* l = new AnalyzeTree()");
gROOT->ProcessLine("l->Loop()");
}
