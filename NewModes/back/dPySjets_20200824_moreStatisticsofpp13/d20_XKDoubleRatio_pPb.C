#include "inc/PyJetUtils.h"

void XKDoubleRatio_pPb();
void XLDoubleRatio_pPb();

void d20_XKDoubleRatio_pPb(){
  XKDoubleRatio_pPb();
  XLDoubleRatio_pPb();

}
void XKDoubleRatio_pPb(){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[4]; TH1D *h[4]; TH1D *hE[4]; TH1D *hR[4]; TGraph *gE[4];
  l[0]= (TList*)f->Get(Form("Xi_toKRatio_0100"));
  l[1]= (TList*)f->Get(Form("Xi_toKRatio_010"));
  l[2]= (TList*)f->Get(Form("Xi_toKRatio_1040"));
  l[3]= (TList*)f->Get(Form("Xi_toKRatio_40100"));
  f->Close();
  for(Int_t i = 0; i<4; i++){ h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));   hE[i] = (TH1D*)l[i]->FindObject(Form("heJER")); }
  TH1D *heSum[4];
  for(Int_t i = 0; i<4; i++){
    hR[i] = (TH1D*)h[i]->Clone(Form("hR_%d", i));
    hR[i]->Divide(h[0]);
    TH1D* he[2] = {hE[0], hE[i]};
    heSum[i] = (TH1D*)QuadraticSum(2, he);
    gE[i] = (TGraphErrors*)ConvHistogramToGraphErrors(hR[i], heSum[i], hR[i]->GetNbinsX()); 
  }

  //Draw hR[1 to 3] gE[1 to 3]
//=============================================================================
  auto dflx(0.), dfux(8.);
  auto dfly(0.), dfuy(3.5);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Double ratio : 0-100% as Ref");
   
  SetStyle(kTRUE);
//=============================================================================

  auto can(MakeCanvas(Form("DoubleRatio_XK")));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(hR[1], wcl[0], wmk[0], "same"); DrawGraph(gE[1], wcl[0], "E2");
  DrawHisto(hR[2], wcl[1], wmk[1], "same"); DrawGraph(gE[2], wcl[1], "E2");
  DrawHisto(hR[3], wcl[2], wmk[2], "same"); DrawGraph(gE[3], wcl[2], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(hR[1], "0-10%",   "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[2], "10-40%",  "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[3], "40-100%", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[1], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#Xi^{-} + #bar{#Xi}^{+}}{2 K^{0}_{S}} in jet"));
  
  can->SaveAs(Form("./figure/20_pPb5d02TeV_Xi_toKRatio_DoubleRatio_Data_Cone3.eps"));
  can->SaveAs(Form("./figure/20_pPb5d02TeV_Xi_toKRatio_DoubleRatio_Data_Cone3.pdf"));
  can->SaveAs(Form("./figure/20_pPb5d02TeV_Xi_toKRatio_DoubleRatio_Data_Cone3.png"));
  CanvasEnd(can);


  return;
}

void XLDoubleRatio_pPb(){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[4]; TH1D *h[4]; TH1D *hE[4]; TH1D *hR[4]; TGraph *gE[4];
  l[0]= (TList*)f->Get(Form("Xi_toLRatio_0100"));
  l[1]= (TList*)f->Get(Form("Xi_toLRatio_010"));
  l[2]= (TList*)f->Get(Form("Xi_toLRatio_1040"));
  l[3]= (TList*)f->Get(Form("Xi_toLRatio_40100"));
  f->Close();
  for(Int_t i = 0; i<4; i++){ h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));   hE[i] = (TH1D*)l[i]->FindObject(Form("heJER")); }
  TH1D *heSum[4];
  for(Int_t i = 0; i<4; i++){
    hR[i] = (TH1D*)h[i]->Clone(Form("hR_%d", i));
    hR[i]->Divide(h[0]);
    TH1D* he[2] = {hE[0], hE[i]};
    heSum[i] = (TH1D*)QuadraticSum(2, he);
    gE[i] = (TGraphErrors*)ConvHistogramToGraphErrors(hR[i], heSum[i], hR[i]->GetNbinsX()); 
  }

  //Draw hR[1 to 3] gE[1 to 3]
//=============================================================================
  auto dflx(0.), dfux(8.);
  auto dfly(0.), dfuy(3.5);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Double ratio : 0-100% as Ref");
   
  SetStyle(kTRUE);
//=============================================================================

  auto can(MakeCanvas(Form("DoubleRatio_XL")));
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
  tex->DrawLatex(0.16, 0.82, Form("#frac{#Xi^{-} + #bar{#Xi}^{+}}{#Lambda + #bar{#Lambda}} in jet"));
  
  can->SaveAs(Form("./figure/20_pPb5d02TeV_Xi_toLRatio_DoubleRatio_Data_Cone3.eps"));
  can->SaveAs(Form("./figure/20_pPb5d02TeV_Xi_toLRatio_DoubleRatio_Data_Cone3.pdf"));
  can->SaveAs(Form("./figure/20_pPb5d02TeV_Xi_toLRatio_DoubleRatio_Data_Cone3.png"));
  CanvasEnd(can);


  return;
}

