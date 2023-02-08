// this macro compiles and runs AnalyzeTree.cxx, which takes as input the 
// starlight.root file produced by convertStarlightAsciiToTree.cxx
// output histograms are stored in starlight_histos.root 
//
void AnaTree0(){

gROOT->ProcessLine("TFile f(\"starlightb4_histos.root\")");
gROOT->ProcessLine("TFile g(\"starlightnew_histos.root\")");
gROOT->ProcessLine("TCanvas *PR = new TCanvas(\"PR\",\"ParentRapidity\")");
gROOT->ProcessLine("TH1D* PR0 = (TH1D*) (f.FindObjectAny(\"ParentRapidity;1\"))");
gROOT->ProcessLine("PR0->Draw()");
gROOT->ProcessLine("gStyle->SetHistLineColor(kRed)");
gROOT->ProcessLine("ParentRapidity->UseCurrentStyle()");
gROOT->ProcessLine("ParentRapidity->Draw(\"SAME\")");
gROOT->ProcessLine("PR -> SaveAs(\"starlight_ParentRapidity.png\")");

gROOT->ProcessLine("TCanvas *PPT = new TCanvas(\"PPT\",\"ParentPT\")");
gROOT->ProcessLine("PPT->cd(1)");
gROOT->ProcessLine("TH1D* PPT0 = (TH1D*) (f.FindObjectAny(\"ParentPt;1\"))");
gROOT->ProcessLine("PPT0->Draw()");
gROOT->ProcessLine("gStyle->SetHistLineColor(kRed)");
gROOT->ProcessLine("ParentPt->UseCurrentStyle()");
gROOT->ProcessLine("ParentPt->Draw(\"SAME\")");
gROOT->ProcessLine("PPT -> SaveAs(\"starlight_ParentPt.png\")");

//gROOT->ProcessLine("ParentPt->Draw()");
gROOT->ProcessLine("TCanvas *PMass = new TCanvas(\"PMass\",\"ParentMass\")");
gROOT->ProcessLine("PMass->cd(1)");
gROOT->ProcessLine("TH1D* PM0 = (TH1D*) (f.FindObjectAny(\"ParentMass;1\"))");
gROOT->ProcessLine("PM0->Draw()");
gROOT->ProcessLine("gStyle->SetHistLineColor(kRed)");
gROOT->ProcessLine("ParentMass->UseCurrentStyle()");
gROOT->ProcessLine("ParentMass->Draw(\"SAME\")");
gROOT->ProcessLine("PMass -> SaveAs(\"starlight_ParentMass.png\")");

//gROOT->ProcessLine("ParentMass->Draw()");
}
