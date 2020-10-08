#include "./SourceFun.h"

void PartoKRatio(){


  TString sPath[] = {"/home/cuipengyao/study/pythia/Rope_pyStrangeJets/Soft/pp13TeV"};
  TString sFile[] = {"AnalysisResults_merged.root"};
  TString sList[] = {"list_results"};
  
  TString sK[] = {"hKshort", "hKshort_Jet10_C04", "hKshort_PCJet10_C04", "hKshort_OCJet10_C04"};
  TString sP[] = {"hLambda", "hLambda_Jet10_C04", "hLambda_PCJet10_C04", "hLambda_OCJet10_C04"};
  TString sA[] = {"hAntiLa", "hAntiLa_Jet10_C04", "hAntiLa_PCJet10_C04", "hAntiLa_OCJet10_C04"};

  TString sLeg[] = {"Incl", "JE", "PC", "OC"};

  auto sLatex("#Lambda to K^{0}_{S} ratio in pPb");
  Double_t XMin = 0.6;
  Double_t XMax = 12;
  //-----------------------------------
 
  TCanvas *can = nullptr;
  TLatex* tex = nullptr;
  TLegend *leg = nullptr;

  TFile *f = TFile::Open("./result/ParToKRatio_PYTHIA_pp.root", "UPDATE");
  can = MakeCanvas("can");

  leg = new TLegend(0.7,0.9,1.0,0.6); SetLegend(leg);
  //-----------------------------------
  auto nHist = sizeof(sK)/sizeof(TString);
  
  TH1D* hP = nullptr;
  TH1D* hA = nullptr;
  TH1D* hK = nullptr;
  TH1D* hRatio = nullptr;

  //Double_t dBin[] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5.0, 6.0, 8.0, 12.}; // LK
  //Double_t dBin[] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5.0, 6.0, 8.0}; // XK
  //Double_t dBin[] = {0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5.0}; // OK
  Double_t dBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5, 8, 12.};
  //Double_t dBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5, 8.};
  //Double_t dBin[] = { 0.6, 2.2, 2.8, 3.7, 5.};

  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  auto *h = new TH1D("h", "", nBin, dBin);

  for (Int_t i = 0; i< nHist; ++i){
    TString sMyPath = sPath[0]; 
    TString sMyFile = sFile[0]; 
    TString sMyList = sList[0]; 
    TString sMyHistP = sP[i]; 
    TString sMyHistA = sA[i]; 
    TString sMyHistK = sK[i]; 

    hK = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistK.Data());
    hK->SetName(Form("Kshort_%d", i)); 
    //DeNormHistBinWidth(hK); 
    hK=RebinTH1D(hK, h);
    //NormHistBinWidth(hK);
    hK->Scale(2.); 

    hP = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistP.Data());
    hP->SetName(Form("Lambda_%d", i)); 
    //DeNormHistBinWidth(hP); 
    hP=RebinTH1D(hP, hK);
    //NormHistBinWidth(hP);
    
    hA = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistA.Data());
    hA->SetName(Form("AntiLa_%d", i)); 
    //DeNormHistBinWidth(hA); 
    hA=RebinTH1D(hA, hK);
    //NormHistBinWidth(hA);
    hP->Add(hA);
    
    hRatio = (TH1D*) hK->Clone(Form("hRatio_%d", i)); hRatio->Reset();
    hRatio->Divide(hP, hK);
    hRatio->SetName(Form("hratio_%d", i));
    f->cd(); hRatio->Write(Form("LKRatio_%s", sLeg[i].Data()), TObject::kSingleKey);
    
    auto dMini = hRatio->GetMinimum();
    auto dMaxi = hRatio->GetMaximum();
    //hRatio->GetYaxis()->SetRangeUser(0.8*dMini, 1.2*dMaxi);
    hRatio->GetXaxis()->SetRangeUser(XMin, XMax);
    hRatio->GetYaxis()->SetRangeUser(0., 1.5);
    DrawHisto(hRatio, cLine[i], sMark[i], "same");
    leg->AddEntry(hRatio, sLeg[i], "lp");
    //SetFrame(hRatio, "#it{p}_{T} (GeV/#it{c})", "#frac{#Lambda + #bar{#Lambda}}{2K^{0}_{S}}");
    SetFrame(hRatio, "#it{p}_{T} (GeV/#it{c})", "Particles to {2K^{0}_{S}} Ratio");
  }
  f->Close();
  tex =  new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  //tex->DrawLatex(0.2, 0.9, sLatex);
  leg->Draw();
  gStyle->SetOptStat("");
  CanvasEnd(can);
  //can->SaveAs(Form("./figure/OKRatio_JE_pPb.eps"));
  return;
}
