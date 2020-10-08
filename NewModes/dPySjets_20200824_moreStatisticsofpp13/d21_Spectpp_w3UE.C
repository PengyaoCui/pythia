#include "inc/PyJetUtils.h"

void Spectpp_w3UE(TString sType = "Kshort");
void d21_Spectpp_w3UE(){
  Spectpp_w3UE("Kshort");
  Spectpp_w3UE("Lambda_sum");
  Spectpp_w3UE("Xi");
  Spectpp_w3UE("Omega");

}
void Spectpp_w3UE(TString sType = "Kshort"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s", sType.Data()));
  f->Close();
  TH1D *h[5]; TGraph *gE[5];
  h[0] = (TH1D*)l->FindObject(Form("InclCen")); gE[0] = (TGraphErrors*)l->FindObject(Form("Inclerr"));
  h[1] = (TH1D*)l->FindObject(Form("JCCen"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JCerr"));
  h[2] = (TH1D*)l->FindObject(Form("UECen"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UEerr"));
  h[3] = (TH1D*)l->FindObject(Form("OCCen"));   gE[3] = (TGraphErrors*)l->FindObject(Form("OCerr"));
  h[4] = (TH1D*)l->FindObject(Form("NJCen"));   gE[4] = (TGraphErrors*)l->FindObject(Form("NJerr"));

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
  
  auto can(MakeCanvas(Form("Data_%s", sType.Data())));
  can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[4], wmk[4], "same"); DrawGraph(gE[3], wcl[4], "E2");
  DrawHisto(h[4], wcl[5], wmk[5], "same"); DrawGraph(gE[4], wcl[5], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "PC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "OC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[4], "NJ", "LP")->SetTextSizePixels(24);
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("K^{0}_{S}"));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("#Lambda + #bar{#Lambda}"));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#%s^{-} + #bar{#%s}^{-}", sType.Data(), sType.Data()));
  
  can->SaveAs(Form("./figure/21_pp13TeV_%s_SpectData_w3UE.eps", sType.Data()));
  can->SaveAs(Form("./figure/21_pp13TeV_%s_SpectData_w3UE.pdf", sType.Data()));
  can->SaveAs(Form("./figure/21_pp13TeV_%s_SpectData_w3UE.png", sType.Data()));
  CanvasEnd(can);

  return;
}

