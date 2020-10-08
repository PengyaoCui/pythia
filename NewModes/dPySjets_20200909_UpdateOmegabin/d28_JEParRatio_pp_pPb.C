#include "inc/PyJetUtils.h"

void PKRatio_pPb(TString sType = "Lambda_sum");
void PLRatio_pPb(TString sType = "Xi");
void PXRatio_pPb(TString sType = "Omega");

void d28_JEParRatio_pp_pPb(){
  PKRatio_pPb("Lambda_sum");
  PKRatio_pPb("Xi");
  PKRatio_pPb("Omega");
  PLRatio_pPb("Xi");
  PLRatio_pPb("Omega");
  PXRatio_pPb("Omega");
}
void PKRatio_pPb(TString sType = "Lambda_sum"){
//=============================================================================
  TList* l[5]; TFile *f[2];
  f[0] = TFile::Open("./data/Data/pp.root", "read");
  f[1] = TFile::Open("./data/Data/pPb.root", "read");
  l[4] = (TList*)f[0]->Get(Form("%s_toKRatio", sType.Data()));
  
  l[0] = (TList*)f[1]->Get(Form("%s_toKRatio_0100", sType.Data()));
  l[1] = (TList*)f[1]->Get(Form("%s_toKRatio_010", sType.Data()));
  l[2] = (TList*)f[1]->Get(Form("%s_toKRatio_1040", sType.Data()));
  l[3] = (TList*)f[1]->Get(Form("%s_toKRatio_40100", sType.Data()));
  f[0]->Close(); 
  f[1]->Close();
  
  TH1D *h[5]; TGraph *gE[5];
  for(Int_t i = 0; i<5; i++){
    h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));
    gE[i] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }
//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(1.0);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");
  if(sType == "Xi"){
    dfux = 8.;
    dfuy = 0.1;
    stny = "Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / 2K_{S}^{0}";
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 0.05;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / 2K_{S}^{0}";
  }
 
  SetStyle(kTRUE);

  auto can(MakeCanvas(Form("%s_KRatio", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);
  if(sType != "Omega"){
    DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
    DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
    DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
    DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");
    DrawHisto(h[4], wcl[4], wmk[4], "same"); DrawGraph(gE[4], wcl[4], "E2");
  }else{
    DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
    DrawHisto(h[4], wcl[4], wmk[4], "same"); DrawGraph(gE[4], wcl[4], "E2");
  }
  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  if(sType != "Omega"){
    leg->AddEntry(h[0], "0-100%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[1], "0-10%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[2], "10-40%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[3], "40-100%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[4], "pp", "LP")->SetTextSizePixels(24);
  }else{
    leg->AddEntry(h[0], "p-Pb", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[4], "pp", "LP")->SetTextSizePixels(24);
  }
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, "pp at #sqrt{#it{s}} = 13 TeV TeV");
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.72, "#frac{#Lambda + #bar{#Lambda}}{2K^{0}_{S}}  in jet");
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.72, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{2K^{0}_{S}} in jet", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toKRatioData_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toKRatioData_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toKRatioData_Cone3.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void PLRatio_pPb(TString sType = "Xi"){
//=============================================================================
  TList* l[5]; TFile *f[2];
  f[0] = TFile::Open("./data/Data/pp.root", "read");
  f[1] = TFile::Open("./data/Data/pPb.root", "read");
  l[4] = (TList*)f[0]->Get(Form("%s_toLRatio", sType.Data()));

  l[0] = (TList*)f[1]->Get(Form("%s_toLRatio_0100", sType.Data()));
  l[1] = (TList*)f[1]->Get(Form("%s_toLRatio_010", sType.Data()));
  l[2] = (TList*)f[1]->Get(Form("%s_toLRatio_1040", sType.Data()));
  l[3] = (TList*)f[1]->Get(Form("%s_toLRatio_40100", sType.Data()));
  f[0]->Close();
  f[1]->Close();

  TH1D *h[5]; TGraph *gE[5];
  for(Int_t i = 0; i<5; i++){
    h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));
    gE[i] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }


//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(1.0);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Lambda + #bar{#Lambda}");
  if(sType == "Xi"){
    dfux = 8.;
    dfuy = 0.5;
    stny = "Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / (#Lambda + #bar{#Lambda})";
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 0.1;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / (#Lambda + #bar{#Lambda})";
  }
 
  SetStyle(kTRUE);
//=============================================================================

  auto can(MakeCanvas(Form("%s_LRatio", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  if(sType != "Omega"){
    DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
    DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
    DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
    DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");
    DrawHisto(h[4], wcl[4], wmk[4], "same"); DrawGraph(gE[4], wcl[4], "E2");
  }else{
    DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
    DrawHisto(h[4], wcl[4], wmk[4], "same"); DrawGraph(gE[4], wcl[4], "E2");
  }
  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  if(sType != "Omega"){
    leg->AddEntry(h[0], "0-100%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[1], "0-10%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[2], "10-40%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[3], "40-100%", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[4], "pp", "LP")->SetTextSizePixels(24);
  }else{
    leg->AddEntry(h[0], "p-Pb", "LP")->SetTextSizePixels(24);
    leg->AddEntry(h[4], "pp", "LP")->SetTextSizePixels(24);
  }
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();


  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, "pp at #sqrt{#it{s}} = 13 TeV");
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.72, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Lambda + #bar{#Lambda}} in jet", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toLRatioData_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toLRatioData_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toLRatioData_Cone3.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void PXRatio_pPb(TString sType = "Omega"){
//=============================================================================
  TList* l[2]; TFile *f[2];
  f[0] = TFile::Open("./data/Data/pp.root", "read");
  f[1] = TFile::Open("./data/Data/pPb.root", "read");
  l[1] = (TList*)f[0]->Get(Form("%s_toXRatio", sType.Data()));

  l[0] = (TList*)f[1]->Get(Form("%s_toXRatio_0100", sType.Data()));
  f[0]->Close();
  f[1]->Close();

  TH1D *h[2]; TGraph *gE[2];
  for(Int_t i = 0; i<2; i++){
    h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));
    gE[i] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }

//=============================================================================
  auto dflx(0.), dfux(5.);
  auto dfly(0.), dfuy(1.);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Xi^{-} + #bar{#Xi}^{+}");
  
  SetStyle(kTRUE);
  
  auto can(MakeCanvas(Form("%s_toXRatio", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[4], wmk[4], "same"); DrawGraph(gE[1], wcl[4], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "p-Pb", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "pp", "LP")->SetTextSizePixels(24);
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.72, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Xi^{-} + #bar{#Xi}^{+}} in jet", sType.Data(), sType.Data()));
 
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toXRatioData_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toXRatioData_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/28_pp_pPb_%s_toXRatioData_Cone3.png", sType.Data()));
 
  CanvasEnd(can);

  return;
}

