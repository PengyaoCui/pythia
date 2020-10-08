#include "./SourceFun.h"

void Test(){


  TString Path = "/home/cuipengyao/study/pythia/pyStrangeJets/pp13TeV";
  TString sPath[] = {Path};
  TString sFile[] = {"AnalysisResults_hard.root"};
  TString sList[] = {"list_results"};
  
  TString sK[] = {"hKshort"};
  TString sP[] = {"hOmegaNeg"};
  TString sA[] = {"hOmegaPos"};

  TString sLeg[] = {"Incl", "JC", "PC", "OC"};

  auto sLatex("#Omega to K^{0}_{S} ratio in pp at #sqrt{s} = 13 TeV");
  Double_t XMin = 0.6;
  Double_t XMax = 12;
  //-----------------------------------
 
  auto nHist = sizeof(sK)/sizeof(TString);
  
  TH1D* hP = nullptr;
  TH1D* hA = nullptr;
  TH1D* hK = nullptr;
  TH1D* hRatio = nullptr;

  //Double_t dBin[] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5.0, 6.0, 8.0, 12.};
  Double_t dBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5, 8, 12.};
  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  auto *h = new TH1D("h", "", nBin, dBin);

  for (Int_t i = 0; i< nHist; ++i){
    TString sMyPath = sPath[i]; 
    TString sMyFile = sFile[i]; 
    TString sMyList = sList[i]; 
    TString sMyHistP = sP[i]; 
    TString sMyHistA = sA[i]; 
    TString sMyHistK = sK[i]; 

    hK = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistK.Data());
    hK->SetName(Form("Kshort_%d", i)); 
    DeNormHistBinWidth(hK); 
    hK=RebinTH1D(hK, h);
    NormHistBinWidth(hK);
    hK->Scale(2.); 

    hP = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistP.Data());
    hP->SetName(Form("OmegaNeg_%d", i)); 
    DeNormHistBinWidth(hP); 
    hP=RebinTH1D(hP, hK);
    NormHistBinWidth(hP);
    
    hA = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistA.Data());
    hA->SetName(Form("OmegaPos_%d", i)); 
    DeNormHistBinWidth(hA); 
    hA=RebinTH1D(hA, hK);
    NormHistBinWidth(hA);
    hP->Add(hA);
   
    hRatio = (TH1D*) hK->Clone(Form("hRatio_%d", i)); hRatio->Reset();
    hRatio->Divide(hP, hK);
    hRatio->SetName(Form("hratio_%d", i));
    
    //auto dMini = hRatio->GetMinimum();
    //auto dMaxi = hRatio->GetMaximum();
    //hRatio->GetYaxis()->SetRangeUser(0.5*dMini, 1.5*dMaxi);
    //hRatio->GetXaxis()->SetRangeUser(XMin, XMax);
    //hRatio->Draw("LCP");
    //hRatio->SetOption("L");
    //DrawHisto(hRatio, cLine[0], sMark[0], "same");
    TGraph *gPySoft = new TGraph(hRatio); gPySoft->SetLineStyle(9);
    DrawGraph(gPySoft, cLine[0], "C");
    //gPySoft->Draw("C");
  }
  return;
}
