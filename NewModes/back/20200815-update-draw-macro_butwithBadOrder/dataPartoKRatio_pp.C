#include "inc/PyJetUtils.h"

TH1D *h[4]; TGraph *gE[4];
auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(1.0);


auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("Ratio: (#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");

void  DrawCones(TString sType = "Lambda_sum");

void dataPartoKRatio_pp(TString sType = "Lambda_sum"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio", sType.Data()));
  f->Close();
  h[0] = (TH1D*)l->FindObject(Form("hInclR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InclRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJCR"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JCRerr"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UERerr"));
  h[3] = (TH1D*)l->FindObject(Form("hJER"));   gE[3] = (TGraphErrors*)l->FindObject(Form("JERerr"));

//=============================================================================
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
  DrawCones(sType);
//=============================================================================
  return;
}
//=============================================================================

void  DrawCones(TString sType = "Lambda_sum")
{
  auto can(MakeCanvas(Form("Ratio_%s", sType.Data())));
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
  
  can->SaveAs(Form("./figure/pp13TeV_%s_toKRatioData.eps", sType.Data()));
  can->SaveAs(Form("./figure/pp13TeV_%s_toKRatioData.pdf", sType.Data()));
  can->SaveAs(Form("./figure/pp13TeV_%s_toKRatioData.png", sType.Data()));
  CanvasEnd(can);

  return;
}

