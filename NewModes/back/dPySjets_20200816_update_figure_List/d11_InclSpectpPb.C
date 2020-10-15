#include "inc/PyJetUtils.h"

void InclSpectpPb(TString sType = "Kshort");
void d11_InclSpectpPb(){
  InclSpectpPb("Kshort");
  InclSpectpPb("Lambda_sum");
  InclSpectpPb("Xi");

}

void InclSpectpPb(TString sType = "Kshort"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[4];
  l[0]= (TList*)f->Get(Form("%s_0100", sType.Data()));
  l[1]= (TList*)f->Get(Form("%s_010", sType.Data()));
  l[2]= (TList*)f->Get(Form("%s_1040", sType.Data()));
  l[3]= (TList*)f->Get(Form("%s_40100", sType.Data()));
  f->Close();
  TH1D *h[4]; TGraph *gE[4];
  for(Int_t i = 0; i<4; i++){
    h[i] = (TH1D*)l[i]->FindObject(Form("InclCen")); gE[i] = (TGraphErrors*)l[i]->FindObject(Form("Inclerr"));
  }
//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(1e-6), dfuy(1e2);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("d#it{#rho}/d#it{p}_{T} (#it{c}/GeV)");
  if(sType == "Lambda_sum"){ 
    dfly = 1e-7;
    dfuy = 5e2;
  }
  if(sType == "Xi"){
    dfux = 8.;
    dfly = 1e-6;
    dfuy = 1e1;
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfly = 1e-6;
    dfuy = 1e-1;
  }
  
  SetStyle(kTRUE);
  
  auto can(MakeCanvas(Form("%s_Spect", sType.Data())));
  can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "0-100%",  "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "0-10%",   "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "10-40%",  "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "40-100%", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("Inclusive K^{0}_{S}"));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("Inclusive #Lambda + #bar{#Lambda}"));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("Inclusive #%s^{-} + #bar{#%s}^{+}", sType.Data(), sType.Data()));

  can->SaveAs(Form("./figure/11_pPb5d02TeV_%s_SpectData_Cone0.eps", sType.Data()));
  can->SaveAs(Form("./figure/11_pPb5d02TeV_%s_SpectData_Cone0.pdf", sType.Data()));
  can->SaveAs(Form("./figure/11_pPb5d02TeV_%s_SpectData_Cone0.png", sType.Data()));
  CanvasEnd(can);


  return;
}

