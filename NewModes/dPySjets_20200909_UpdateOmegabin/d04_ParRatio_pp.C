#include "inc/PyJetUtils.h"

void PKRatio_pp(TString sType = "Lambda_sum");
void PLRatio_pp(TString sType = "Xi");
void PXRatio_pp(TString sType = "Omega");

void d04_ParRatio_pp(){
  PKRatio_pp("Lambda_sum");
  PKRatio_pp("Xi");
  PKRatio_pp("Omega");
  PLRatio_pp("Xi");
  PLRatio_pp("Omega");
  PXRatio_pp("Omega");
}
void PKRatio_pp(TString sType = "Lambda_sum"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio", sType.Data()));
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
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, "#frac{#Lambda + #bar{#Lambda}}{2K^{0}_{S}}");
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{2K^{0}_{S}}", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toKRatioData.eps", sType.Data()));
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toKRatioData.pdf", sType.Data()));
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toKRatioData.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void PLRatio_pp(TString sType = "Xi"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toLRatio", sType.Data()));
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
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Lambda + #bar{#Lambda}}", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toLRatioData.eps", sType.Data()));
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toLRatioData.pdf", sType.Data()));
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toLRatioData.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void PXRatio_pp(TString sType = "Omega"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toXRatio", sType.Data()));
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
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Xi^{-} + #bar{#Xi}^{+}}", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toXRatioData.eps", sType.Data()));
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toXRatioData.pdf", sType.Data()));
  can->SaveAs(Form("./figure/04_pp13TeV_%s_toXRatioData.png", sType.Data()));
  CanvasEnd(can);

  return;
}

