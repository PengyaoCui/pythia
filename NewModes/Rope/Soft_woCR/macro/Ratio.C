#include "./SourceFun.h"

void Ratio(){

  const TString sRatio = "p-Pb/pp";
  TString sPathNu[] = {"/home/cuipengyao/study/pythia/pyStrangeJets/pp13TeV/hard"};
  TString sFileNu[] = {"AnalysisResults_hard.root"};
  TString sListNu[] = {"list_results"};
  TString sHistNu[] = {"hLambda_Jet10_C04"};

  TString sPathDe = "/home/cuipengyao/study/pythia/pyStrangeJets/pp13TeV/hard";
  TString sFileDe = "AnalysisResults_hard.root";
  TString sListDe = "list_results";
  TString sHistDe = "hLambda_Jet10_C04";

  TString sLeg[] = {"p-Pb", "pp"};

  auto sLatex(Form(""));
  Double_t XMin = 0.6;
  Double_t XMax = 12.;
  //-----------------------------------
  TCanvas *can = nullptr;
  TPad *padT = nullptr;
  TPad *padB = nullptr;
  TLatex* tex = nullptr;
  TLegend *leg = nullptr;

  can = MakeCanvas("can");
  leg = new TLegend(0.7,0.9,1.0,0.6); SetLegend(leg);
  //-----------------------------------
  auto nHist = sizeof(sHistNu)/sizeof(TString);
  const auto nhist = nHist;

  TH1D* h2 = nullptr;
  TH1D* h1 = nullptr;
  TH1D* hRatio = nullptr;

  padB = MakePadB("padB");
  can->cd();
  padT = MakePadT("padT");
  padT->SetLogy(); 
  padT->SetLogy();

  h2 = (TH1D*)GetTH1D(sPathDe.Data(), sFileDe.Data(), sListDe.Data(), sHistDe.Data());
  h2->SetName("h2_De");
  h2->SetTitle("");
  h2->Rebin(10);
  padT->cd();
  //h2->GetYaxis()->SetRangeUser(1e-7, 1.);
  h2->SetMarkerSize(1);
  DrawHisto(h2, cLine[nHist], sMark[nHist], "same");
  SetFrame(h2, "#it{p}_{T}", "d#rho/d#it{p}_{T}");
  leg->AddEntry(h2, sLeg[nHist],"lp");


  for (Int_t i = 0; i< nHist; ++i){
    TString sMyPath = sPathNu[i]; 
    TString sMyFile = sFileNu[i]; 
    TString sMyList = sListNu[i]; 
    TString sMyHist = sHistNu[i]; 
    h1 = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHist.Data());
    h1->SetName(Form("hist_%d", i)); 
    h1->Rebin(10);
    //DeNormHistBinWidth(h1); 
    //h1=RebinTH1D(h1, h2);
    //NormHistBinWidth(h1); 
    hRatio = (TH1D*) h2->Clone(Form("hRatio_%d", i)); hRatio->Reset();
    hRatio->Divide(h1, h2);
    hRatio->SetName(Form("hratio_%d", i));
    //hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    padT->cd();
    DrawHisto(h1, cLine[i], sMark[i], "same");
    leg->AddEntry(h1, sLeg[i], "lp");
    padB->cd();
    padB->SetGridy();
    DrawHisto(hRatio, cLine[i], sMark[i], "same");
    SetFrame(hRatio, "#it{p}_{T} (GeV/#it{c})", sRatio, 0.08, 0.08, 0.1, 0.1, 1.10, 0.6);
  }


  padB->cd();
  //TLine* l = new TLine(XMin, 1., XMax, 1.);
  //l->SetLineColor(kRed);
  //l->SetLineWidth(2);
  //l->Draw("same");

  padT->cd();
  tex =  new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.2, 0.9, sLatex);
  leg->Draw();
  gStyle->SetOptStat("");
  //can->SaveAs(Form("./figure/LambdaNegPostoLambda_%s%s_Corr.eps", CentMin.Data(), CentMax.Data()));
  CanvasEnd(can);
  return;
}
