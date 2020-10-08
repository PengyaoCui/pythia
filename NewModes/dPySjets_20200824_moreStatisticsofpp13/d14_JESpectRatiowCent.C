#include "inc/PyJetUtils.h"

void JESpectRatiowCent(const TString sType = "Kshort");
void d14_JESpectRatiowCent(){
  JESpectRatiowCent("Kshort");
  JESpectRatiowCent("Lambda");
  JESpectRatiowCent("Xi");
}
void JESpectRatiowCent(const TString sType = "Kshort"){

  TString sLatex(Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
  TList *l[4]; TH1D *h[4]; TH1D *e[4];
  auto f = TFile::Open("./data/Data/pPb.root", "read");
  l[0] = (TList*)f->Get(Form("%s_0100", sType.Data()));
  l[1] = (TList*)f->Get(Form("%s_010", sType.Data()));
  l[2] = (TList*)f->Get(Form("%s_1040", sType.Data()));
  l[3] = (TList*)f->Get(Form("%s_40100", sType.Data()));
  f->Close();
  for(Int_t i = 0; i < 4 ; i++) { h[i] = (TH1D*)l[i]->FindObject("JECen"); e[i] = (TH1D*)l[i]->FindObject("JEErr"); }

  TH1D *hE[4]; TH1D *hR[4]; TGraph *gE[4];
  for(Int_t i = 0; i<4; i++){
    hR[i] = (TH1D*)h[i]->Clone(Form("hR_%d", i));
    hR[i]->Divide(h[0]);
    TH1D* he[2] = {e[0], e[i]};
    hE[i] = (TH1D*)QuadraticSum(2, he);
    gE[i] = (TGraphErrors*)ConvHistogramToGraphErrors(hR[i], hE[i], hR[i]->GetNbinsX());
  }


  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(3.5);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
 
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("0-100% as Ref");

  if(sType == "Xi"){
    dfux = 8.;
  }
  if(sType == "Omega"){
    dfux = 5.;
  }

  SetStyle(kTRUE);
//=============================================================================

  auto can(MakeCanvas(Form("SpectRatio_%s", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(hR[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(hR[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(hR[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(hR[1], "0-10%",   "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[2], "10-40%",  "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[3], "40-100%", "LP")->SetTextSizePixels(24);
  //leg->AddEntry(gE[1], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, "Spectra centrality dependent");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.72, Form("K^{0}_{S} in jet"));
  if(sType == "Lambda") tex->DrawLatex(0.16, 0.72, Form("#Lambda + #bar{#Lambda} in jet"));
  if(sType == "Xi") tex->DrawLatex(0.16, 0.72, Form("#Xi^{-} + #bar{#Xi}^{+} in jet"));

  can->SaveAs(Form("./figure/14_pPb5d02TeV_%s_Spect_DoubleRatio_Data_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/14_pPb5d02TeV_%s_Spect_DoubleRatio_Data_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/14_pPb5d02TeV_%s_Spect_DoubleRatio_Data_Cone3.png", sType.Data()));
  CanvasEnd(can);
 
  return;
}
