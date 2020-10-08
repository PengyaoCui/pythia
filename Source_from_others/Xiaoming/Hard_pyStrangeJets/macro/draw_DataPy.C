#include "./SourceFun.h"

void draw_DataPy(){
  auto Path = "";
  auto File = "";
  TString sDaPath[] = {"./Data", "./Data", "./Data"};
  TString sDaFile[] = {"ParToKRatio_DATA_pp.root", "ParToKRatio_DATA_pp.root", "ParToKRatio_DATA_pp.root"};
  TString sDaList[] = {"", "", ""};
  TString sDaHist[] = {"LKRatio_Incl", "LKRatio_JE", "LKRatio_PC"};
  
  TString sPyPath[] = {"./Data", "./Data", "./Data"};
  TString sPyFile[] = {"ParToKRatio_PYTHIA_pp.root", "ParToKRatio_PYTHIA_pp.root", "ParToKRatio_PYTHIA_pp.root"};
  TString sPyList[] = {"", "", ""};
  TString sPyHist[] = {"LKRatio_Incl", "LKRatio_JC", "LKRatio_PC"};

  TString sLeg[] = {"Inclsive", "JE", "PC"}; 

  auto sLatex("#Lambda to K^{0}_{S} ratio");
  //-----------------------------------
  TCanvas *can = nullptr;

  can = MakeCanvas("can");
  //can->SetLogy();
  //can->SetGridx(); can->SetGridy();
  
  //Double_t dBin[] = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15 };
  Double_t dBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5, 8, 12.};
  
  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  //-----------------------------------

  auto leg = new TLegend(0.66,0.9, 0.96,0.75);    
  //-----------------------------------
  auto nHist = sizeof(sDaHist)/sizeof(TString);
  for (Int_t i = 0; i< nHist; i++){
    TString sMyPath = sDaPath[i]; 
    TString sMyFile = sDaFile[i]; 
    TString sMyList = sDaList[i]; 
    TString sMyHist = sDaHist[i]; 
    
    TH1D* h1 = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHist.Data());
    h1->SetName(Form("hist_%d", i));
    auto dMini = h1->GetMinimum(); 
    auto dMaxi = h1->GetMaximum(); 
    //h1->GetXaxis()->SetRangeUser(0.6, 12.);
   // h1->GetYaxis()->SetRangeUser(0.1*dMini, 3.*dMaxi);
    h1->GetYaxis()->SetRangeUser(0., 1.8*dMaxi);
    h1->SetTitle(""); 
    SetFrame(h1, "#it{p}_{T}(GeV/#it{c})", "#Lambda/K^{0}_{S}");
    DrawHisto(h1, cLine[i], sMark[i], "same");
    leg->AddEntry(h1, sLeg[i],"lp");
  }
  SetLegend(leg);
  leg->Draw();

  leg = new TLegend(0.66, 0.7, 0.96, 0.6);

  nHist = sizeof(sPyHist)/sizeof(TString); 
  TGraph *g1 = nullptr;
  for (Int_t i = 0; i< nHist; i++){
    TString sMyPath = sPyPath[i]; 
    TString sMyFile = sPyFile[i]; 
    TString sMyList = sPyList[i]; 
    TString sMyHist = sPyHist[i]; 
    
    TH1D* h1 = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHist.Data());
    h1->SetName(Form("HIST_%d", i));
    h1->SetTitle(""); 
    g1 = new TGraph(h1); g1->SetLineStyle(9);
    g1->SetName(Form("graph_%d", i));
    DrawGraph(g1, cLine[i], "C");
  }
  TGraph *gTemp = new TGraph();
  gTemp->SetLineStyle(9);
  gTemp->SetLineWidth(2);

  leg->AddEntry(gTemp, "PYTHIA8", "LP")->SetTextSizePixels(24);
  SetLegend(leg);
  leg->Draw();


  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.91, sLatex);
  tex->DrawLatex(0.16, 0.8, "|#eta| < 0.75");
  tex->DrawLatex(0.16, 0.7, "pp #sqrt{#it{s}} = 13 TeV");
  gStyle->SetOptStat("");

  //can->SaveAs(Form("./figure/LKRatio_AllCone_DataPy.eps"));
  CanvasEnd(can);
  return;
}
