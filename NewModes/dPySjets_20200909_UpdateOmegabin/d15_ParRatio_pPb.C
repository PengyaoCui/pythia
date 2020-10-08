#include "inc/PyJetUtils.h"

void PKRatio_pPb(TString sType = "Lambda_sum");
void PLRatio_pPb(TString sType = "Xi");
void PXRatio_pPb(TString sType = "Omega");

void d15_ParRatio_pPb(){
  PKRatio_pPb("Lambda_sum");
  PKRatio_pPb("Xi");
  PKRatio_pPb("Omega");
  PLRatio_pPb("Xi");
  PLRatio_pPb("Omega");
  PXRatio_pPb("Omega");
}
void PKRatio_pPb(TString sType = "Lambda_sum"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio_0100", sType.Data()));
  f->Close();
  TH1D *h[4]; TGraph *gE[4];
  h[0] = (TH1D*)l->FindObject(Form("hInR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJCR"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JCRerr"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UERerr"));
  h[3] = (TH1D*)l->FindObject(Form("hJER"));   gE[3] = (TGraphErrors*)l->FindObject(Form("JERerr"));

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
    dfuy = 0.2;
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

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, "#frac{#Lambda + #bar{#Lambda}}{2K^{0}_{S}}  0-100%");
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{2K^{0}_{S}}  0-100%%", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toKRatioData_Cent0.eps", sType.Data()));
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toKRatioData_Cent0.pdf", sType.Data()));
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toKRatioData_Cent0.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void PLRatio_pPb(TString sType = "Xi"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  auto l = (TList*)f->Get(Form("%s_toLRatio_0100", sType.Data()));
  f->Close();
  TH1D *h[4]; TGraph *gE[4];
  h[0] = (TH1D*)l->FindObject(Form("hInR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJCR"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JCRerr"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UERerr"));
  h[3] = (TH1D*)l->FindObject(Form("hJER"));   gE[3] = (TGraphErrors*)l->FindObject(Form("JERerr"));

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
    dfuy = 0.2;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / (#Lambda + #bar{#Lambda})";
  }
 
  SetStyle(kTRUE);
//=============================================================================

  auto can(MakeCanvas(Form("%s_LRatio", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Lambda + #bar{#Lambda}}  0-100%%", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toLRatioData_Cent0.eps", sType.Data()));
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toLRatioData_Cent0.pdf", sType.Data()));
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toLRatioData_Cent0.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void PXRatio_pPb(TString sType = "Omega"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  auto l = (TList*)f->Get(Form("%s_toXRatio_0100", sType.Data()));
  f->Close();
  TH1D *h[4]; TGraph *gE[4];
  h[0] = (TH1D*)l->FindObject(Form("hInR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJCR"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JCRerr"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UERerr"));
  h[3] = (TH1D*)l->FindObject(Form("hJER"));   gE[3] = (TGraphErrors*)l->FindObject(Form("JERerr"));

//=============================================================================
  auto dflx(0.), dfux(5.);
  auto dfly(0.), dfuy(1.0);
  
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
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Xi^{-} + #bar{#Xi}^{+}} in 0-100%%", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toXRatioData_Cent0.eps", sType.Data()));
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toXRatioData_Cent0.pdf", sType.Data()));
  can->SaveAs(Form("./figure/15_pPb5d02TeV_%s_toXRatioData_Cent0.png", sType.Data()));
  CanvasEnd(can);

  return;
}

