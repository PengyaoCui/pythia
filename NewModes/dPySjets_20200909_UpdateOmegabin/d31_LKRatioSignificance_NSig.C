#include "inc/PyJetUtils.h"


void d31_LKRatioSignificance_NSig(){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[3]; TH1D *h[3]; TH1D *hE[3]; TH1D *hR[3]; TGraph *gE[3];
  l[0]= (TList*)f->Get(Form("Lambda_sum_toKRatio_010"));
  l[1]= (TList*)f->Get(Form("Lambda_sum_toKRatio_1040"));
  l[2]= (TList*)f->Get(Form("Lambda_sum_toKRatio_40100"));
  f->Close();
  for(Int_t i = 0; i<3; i++){ h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));   hE[i] = (TH1D*)l[i]->FindObject(Form("heJER")); }
  Double_t dCent[3]; Double_t dStat[3]; Double_t dSyst[3]; Double_t dErr[3];
  auto hNSigma0 = (TH1D*)h[0]->Clone("hNSigma0"); hNSigma0->Reset(); 
  auto hNSigma2 = (TH1D*)h[0]->Clone("hNSigma2"); hNSigma2->Reset();
  for(Int_t i = 1; i<= h[0]->GetNbinsX(); i++){
    dCent[0] =  h[0]->GetBinContent(i); dCent[1] =  h[1]->GetBinContent(i); dCent[2] =  h[2]->GetBinContent(i);
    dStat[0] =  h[0]->GetBinError(i);   dStat[1] =  h[1]->GetBinError(i);   dStat[2] =  h[2]->GetBinError(i);
    dSyst[0] = hE[0]->GetBinContent(i); dSyst[1] = hE[1]->GetBinContent(i); dSyst[2] = hE[2]->GetBinContent(i);
    auto dNSigma0 = TMath::Abs(dCent[0] - dCent[1])/TMath::Sqrt(dStat[0]*dStat[0] + dSyst[0]*dSyst[0] + dStat[1]*dStat[1] + dSyst[1]*dSyst[1]);
    auto dNSigma2 = TMath::Abs(dCent[2] - dCent[1])/TMath::Sqrt(dStat[2]*dStat[2] + dSyst[2]*dSyst[2] + dStat[1]*dStat[1] + dSyst[1]*dSyst[1]);
    hNSigma0->SetBinContent(i, dNSigma0); hNSigma0->SetBinError(i, 0);
    hNSigma2->SetBinContent(i, dNSigma2); hNSigma2->SetBinError(i, 0);
  } 

  //Draw hR[1 to 3] gE[1 to 3]
//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(2.5);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("N#sigma : 10-40\% as Ref");
   
  SetStyle(kTRUE);
//=============================================================================

  auto can(MakeCanvas(Form("NSigma")));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(hNSigma0, wcl[1], wmk[1], "same");
  DrawHisto(hNSigma2, wcl[3], wmk[3], "same");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(hNSigma0, "0-10%",   "LP")->SetTextSizePixels(24);
  leg->AddEntry(hNSigma2, "40-100%", "LP")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#Lambda + #bar{#Lambda}}{2 K^{0}_{S}} JE"));
  
  can->SaveAs(Form("./figure/d31_pPb5d02TeV_LKRatioSignificance_NSig_Cone3.eps"));
  can->SaveAs(Form("./figure/d31_pPb5d02TeV_LKRatioSignificance_NSig_Cone3.pdf"));
  can->SaveAs(Form("./figure/d31_pPb5d02TeV_LKRatioSignificance_NSig_Cone3.png"));
  CanvasEnd(can);


  return;
}

