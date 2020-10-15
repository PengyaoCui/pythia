#include "./SourceFun.h"

void PartoKRatio_JE(){


  TString sPath = "../";
  TString sFile = "AnalysisResults_merged.root";
  TString sList = "list_results";
  
  TString sKJC = "hKshort_Jet10_C04";
  TString sPJC[] = {"hLambda_Jet10_C04", "hXiNeg_Jet10_C04", "hOmegaNeg_Jet10_C04"};
  TString sAJC[] = {"hAntiLa_Jet10_C04", "hXiPos_Jet10_C04", "hOmegaPos_Jet10_C04"};
  
  TString sKPC = "hKshort_PCJet10_C04";
  TString sPPC[] = {"hLambda_PCJet10_C04", "hXiNeg_PCJet10_C04", "hOmegaNeg_PCJet10_C04"};
  TString sAPC[] = {"hAntiLa_PCJet10_C04", "hXiPos_PCJet10_C04", "hOmegaPos_PCJet10_C04"};

  TString sLeg[] = {"LKRatio", "XKRatio", "OKRatio"};

  auto sLatex("#Lambda to K^{0}_{S} ratio in pp at #sqrt{s} = 13 TeV");
  Double_t XMin = 0.6;
  Double_t XMax = 12;
  //-----------------------------------
 
  TCanvas *can = nullptr;
  TLatex* tex = nullptr;
  TLegend *leg = nullptr;

  TFile *f = TFile::Open("./result/ParToKRatio_PYTHIA_pp.root", "UPDATE");
  can = MakeCanvas("can");

  leg = new TLegend(0.7,0.9,1.0,0.7); SetLegend(leg);
  //-----------------------------------

  Int_t nHist = sizeof(sPPC)/sizeof(TString);

  TH1D* hPJC = nullptr;
  TH1D* hAJC = nullptr;
  TH1D* hKJC = nullptr;
  
  TH1D* hPPC = nullptr;
  TH1D* hAPC = nullptr;
  TH1D* hKPC = nullptr;
  
  TH1D* hRatio = nullptr;

  //Double_t dBin[] = {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2, 3.7, 4.2, 5.0, 6.0, 8.0, 12.};
  //Double_t dBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5, 8, 12.};
  Double_t dBin[] = { 0.6, 2.2, 3.7, 5, 8, 12., 16.};
  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  auto *h = new TH1D("h", "", nBin, dBin);

  hKJC = (TH1D*)GetTH1D(sPath.Data(), sFile.Data(), sList.Data(), sKJC.Data());
  hKPC = (TH1D*)GetTH1D(sPath.Data(), sFile.Data(), sList.Data(), sKPC.Data());
  hKJC->SetName(Form("Kshort_JC")); 
  hKPC->SetName(Form("Kshort_PC")); 
  DeNormHistBinWidth(hKJC); 
  DeNormHistBinWidth(hKPC); 
  hKJC->Add(hKPC, -1.*0.25);
  hKJC=RebinTH1D(hKJC, h);
  hKPC=RebinTH1D(hKPC, h);
  NormHistBinWidth(hKJC);
  NormHistBinWidth(hKPC);
  hKJC->Scale(2.); 

  for(Int_t i = 0; i< nHist; i++){
    hPJC = (TH1D*)GetTH1D(sPath.Data(), sFile.Data(), sList.Data(), sPJC[i].Data());
    hPJC->SetName(Form("Particle_JC_%d", i)); 
    hPPC = (TH1D*)GetTH1D(sPath.Data(), sFile.Data(), sList.Data(), sPPC[i].Data());
    hPPC->SetName(Form("Particle_PC_%d", i)); 
    
    DeNormHistBinWidth(hPJC); 
    DeNormHistBinWidth(hPPC); 
    hPJC=RebinTH1D(hPJC, hKJC);
    hPPC=RebinTH1D(hPPC, hKPC);
    NormHistBinWidth(hPJC);
    NormHistBinWidth(hPPC);
    hPJC->Add(hPPC, -1.*0.25);

    hAJC = (TH1D*)GetTH1D(sPath.Data(), sFile.Data(), sList.Data(), sAJC[i].Data());
    hAJC->SetName(Form("AntiPa_JC_%d", i)); 
    hAPC = (TH1D*)GetTH1D(sPath.Data(), sFile.Data(), sList.Data(), sAPC[i].Data());
    hAPC->SetName(Form("AntiPa_PC_%d", i)); 
    DeNormHistBinWidth(hAJC); 
    DeNormHistBinWidth(hAPC); 
    hAJC=RebinTH1D(hAJC, hKJC);
    hAPC=RebinTH1D(hAPC, hKPC);
    NormHistBinWidth(hAJC);
    NormHistBinWidth(hAPC);
    hAJC->Add(hAPC, -1.*0.25);
    
    hPJC->Add(hAJC);
    
    hRatio = (TH1D*) hKJC->Clone(Form("hRatio_%d", i)); hRatio->Reset();
    hRatio->Divide(hPJC, hKJC);
    hRatio->SetName(Form("hratio_%d", i));
    f->cd(); hRatio->Write(Form("%s_JE", sLeg[i].Data()), TObject::kSingleKey);
    auto dMini = hRatio->GetMinimum();
    auto dMaxi = hRatio->GetMaximum();
    hRatio->GetYaxis()->SetRangeUser(0.5*dMini, 1.5*dMaxi);
    hRatio->GetXaxis()->SetRangeUser(XMin, XMax);
    //hRatio->GetYaxis()->SetRangeUser(0., 0.2);
    DrawHisto(hRatio, cLine[i], sMark[i], "same");
    hRatio->SetOption("L");
    hRatio->Draw("same");
    leg->AddEntry(hRatio, sLeg[i], "lp");
  } 
  SetFrame(hRatio, "#it{p}_{T} (GeV/#it{c})", "#frac{#Lambda + #bar{#Lambda}}{2K^{0}_{S}}");
  //SetFrame(hRatio, "#it{p}_{T} (GeV/#it{c})", "Particles to K^{0}_{S} Ratio");
  f->Close();
  tex =  new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.15, 0.9, sLatex);
  leg->Draw();
  gStyle->SetOptStat("");
  CanvasEnd(can);
  //can->SaveAs(Form("./figure/LKRatio_HadUE_pp.eps"));
  return;
}
