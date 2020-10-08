#include "inc/PyJetUtils.h"

TH1D *h[4][4]; TGraph *gE[4][4];

auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(1.0);


auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Xi^{-} + #bar{#Xi}^{+}");

 
void  DrawCones(TString sType = "Omega", Int_t nCent = 0);
void  DrawCents(TString sType = "Omega", Int_t nCone = 0);

void dataPartoXRatio_pPb(TString sType = "Omega", Bool_t IsDiffCone = kTRUE){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[4];
  l[0]= (TList*)f->Get(Form("%s_toXRatio_0100", sType.Data()));
  l[1]= (TList*)f->Get(Form("%s_toXRatio_010", sType.Data()));
  l[2]= (TList*)f->Get(Form("%s_toXRatio_1040", sType.Data()));
  l[3]= (TList*)f->Get(Form("%s_toXRatio_40100", sType.Data()));
  f->Close();
  for(Int_t i = 0; i<4; i++){
    h[i][0] = (TH1D*)l[i]->FindObject(Form("hInclR")); gE[i][0] = (TGraphErrors*)l[i]->FindObject(Form("InclRerr"));
    h[i][1] = (TH1D*)l[i]->FindObject(Form("hJCR"));   gE[i][1] = (TGraphErrors*)l[i]->FindObject(Form("JCRerr"));
    h[i][2] = (TH1D*)l[i]->FindObject(Form("hUER"));   gE[i][2] = (TGraphErrors*)l[i]->FindObject(Form("UERerr"));
    h[i][3] = (TH1D*)l[i]->FindObject(Form("hJER"));   gE[i][3] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }
//=============================================================================
  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 1.;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / (#Xi^{-} + #bar{#Xi}^{+})";
  }
 
  SetStyle(kTRUE);
  if(IsDiffCone)for(Int_t i = 0; i<4; i++) DrawCones(sType, i); 
  if(!IsDiffCone){DrawCents(sType, 0); DrawCents(sType, 3);}// DrawCents(sType, 2); DrawCents(sType, 3);
//=============================================================================
  return;
}
//=============================================================================

void  DrawCones(TString sType = "Omega", Int_t nCent = 0)
{
  TString sCent[] = {"0-100%", "0-10%", "10-40%", "40-100%"};
  auto can(MakeCanvas(Form("DataRatio_%s_Cent%d", sType.Data(), nCent)));
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
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Xi^{-} + #bar{#Xi}^{-}} %s", sType.Data(), sType.Data(), sCent[nCent].Data()));
  
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toXRatioData_Cent%d.eps", sType.Data(), nCent));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toXRatioData_Cent%d.pdf", sType.Data(), nCent));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toXRatioData_Cent%d.png", sType.Data(), nCent));
  CanvasEnd(can);

  return;
}

void  DrawCents(TString sType = "Omega", Int_t nCone = 0)
{
  TString sCone[] = {"Inclusive", "JC", "UE", "JE"};
  auto can(MakeCanvas(Form("DataRatio_%s_Cone%d", sType.Data(), nCone)));
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
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{-}}{#Xi^{-} + #bar{#Xi}^{+}} %s", sType.Data(), sType.Data(), sCone[nCone].Data()));
  
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toXRatioData_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toXRatioData_Cone%d.pdf", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toXRatioData_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);


  return;
}

