#include "inc/PyJetUtils.h"

TH1D *h[4][4]; TGraph *gE[4][4];
auto dflx(0.), dfux(12.);
auto dfly(1e-6), dfuy(1e2);
 

auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("d#it{#rho}/d#it{p}_{T} (#it{c}/GeV)");

void  DrawCones(TString sType = "Kshort", Int_t nCent = 0);
void  DrawCents(TString sType = "Kshort", Int_t nCone = 0);

void dataSpect_pPb(TString sType = "Kshort", Bool_t IsDiffCone = 1){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[4];
  l[0]= (TList*)f->Get(Form("%s_0100", sType.Data()));
  l[1]= (TList*)f->Get(Form("%s_010", sType.Data()));
  l[2]= (TList*)f->Get(Form("%s_1040", sType.Data()));
  l[3]= (TList*)f->Get(Form("%s_40100", sType.Data()));
  f->Close();
  for(Int_t i = 0; i<4; i++){
    h[i][0] = (TH1D*)l[i]->FindObject(Form("InclCen")); gE[i][0] = (TGraphErrors*)l[i]->FindObject(Form("Inclerr"));
    h[i][1] = (TH1D*)l[i]->FindObject(Form("JCCen"));   gE[i][1] = (TGraphErrors*)l[i]->FindObject(Form("JCerr"));
    h[i][2] = (TH1D*)l[i]->FindObject(Form("UECen"));   gE[i][2] = (TGraphErrors*)l[i]->FindObject(Form("UEerr"));
    h[i][3] = (TH1D*)l[i]->FindObject(Form("JECen"));   gE[i][3] = (TGraphErrors*)l[i]->FindObject(Form("JEerr"));
  }
//=============================================================================
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
  if (IsDiffCone){for(Int_t i = 0; i<4; i++) DrawCones(sType, i);} 
  if(!IsDiffCone){DrawCents(sType, 0); DrawCents(sType, 3);}// DrawCents(sType, 2); DrawCents(sType, 3);
//=============================================================================
  return;
}
//=============================================================================

void  DrawCones(TString sType = "Kshort", Int_t nCent = 0)
{
  TString sCent[] = {"0-100%", "0-10%", "10-40%", "40-100%"};
  auto can(MakeCanvas(Form("Data_%s_Cent%d", sType.Data(), nCent)));
  can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[nCent][0], wcl[0], wmk[0], "same"); DrawGraph(gE[nCent][0], wcl[0], "E2");
  DrawHisto(h[nCent][1], wcl[1], wmk[1], "same"); DrawGraph(gE[nCent][1], wcl[1], "E2");
  DrawHisto(h[nCent][2], wcl[2], wmk[2], "same"); DrawGraph(gE[nCent][2], wcl[2], "E2");
  DrawHisto(h[nCent][3], wcl[3], wmk[3], "same"); DrawGraph(gE[nCent][3], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[nCent][0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[nCent][1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[nCent][2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[nCent][3], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[nCent][0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("K^{0}_{S} %s", sCent[nCent].Data()));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("#Lambda + #bar{#Lambda} %s", sCent[nCent].Data()));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#%s^{-} + #bar{#%s}^{-} %s", sType.Data(), sType.Data(), sCent[nCent].Data()));
  
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_SpectData_Cent%d.eps", sType.Data(), nCent));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_SpectData_Cent%d.pdf", sType.Data(), nCent));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_SpectData_Cent%d.png", sType.Data(), nCent));
  CanvasEnd(can);

  return;
}

void  DrawCents(TString sType = "Kshort", Int_t nCone = 0)
{
  TString sCone[] = {"Inclusive", "JC", "UE", "JE"};
  auto can(MakeCanvas(Form("Data_%s_Cone%d", sType.Data(), nCone)));
  can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0][nCone], wcl[0], wmk[0], "same"); DrawGraph(gE[0][nCone], wcl[0], "E2");
  DrawHisto(h[1][nCone], wcl[1], wmk[1], "same"); DrawGraph(gE[1][nCone], wcl[1], "E2");
  DrawHisto(h[2][nCone], wcl[2], wmk[2], "same"); DrawGraph(gE[2][nCone], wcl[2], "E2");
  DrawHisto(h[3][nCone], wcl[3], wmk[3], "same"); DrawGraph(gE[3][nCone], wcl[3], "E2");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0][nCone], "0-100%",  "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1][nCone], "0-10%",   "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2][nCone], "10-40%",  "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3][nCone], "40-100%", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0][nCone], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("K^{0}_{S} %s", sCone[nCone].Data()));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("#Lambda + #bar{#Lambda} %s", sCone[nCone].Data()));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#%s^{-} + #bar{#%s}^{-} %s", sType.Data(), sType.Data(), sCone[nCone].Data()));

  can->SaveAs(Form("./figure/pPb5d02TeV_%s_SpectData_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_SpectData_Cone%d.pdf", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_SpectData_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);


  return;
}

