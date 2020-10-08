#include "inc/PyJetUtils.h"


TH1D *h[3]; TGraph *gE[3]; TGraph *g[3][nm];
auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(1.0);
 

auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Lambda + #bar{#Lambda}");

void  DrawMode(TString sType = "Omega", Int_t nCone = 0);
void  DrawCones(TString sType = "Omega",Int_t nMode = 0);

void PartoLRatio_pp(TString sType = "Omega", Bool_t IsDiffCone = kTRUE){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toLRatio", sType.Data()));
  f->Close();
  h[0] = (TH1D*)l->FindObject(Form("hInclR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InclRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJER"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JERerr"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UERerr"));


  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Xi"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioXL(i,sd));
      g[1][i] = new TGraph(RatioXL(i, sd, sj, "JC04", "PC04"));
      g[2][i] = new TGraph(RatioXL(i, sd, sj, "OC08"));
    }
  }

  if(sType == "Omega"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioOL(i,sd));
      g[1][i] = new TGraph(RatioOL(i, sd, sj, "JC04", "PC04"));
      g[2][i] = new TGraph(RatioOL(i, sd, sj, "OC08"));
    }
  }
  for (auto i=0; i<nm; ++i) for (auto j=0; j<3; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================

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
  if(!IsDiffCone)for(Int_t i = 0; i< 3; i++) DrawMode(sType, i);
  if(IsDiffCone)for(Int_t i = 0; i< 4; i++) DrawCones(sType, i);
//=============================================================================
  return;
}
//=============================================================================

void  DrawMode(TString sType = "Omega", Int_t nCone = 0)
{
  auto can(MakeCanvas(Form("PyDataCone_%d", nCone)));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[nCone], wcl[0], wmk[0], "same"); DrawGraph(gE[nCone], wcl[0], "E2");
  
  DrawGraph(g[nCone][1],  wcl[0], "C");
  DrawGraph(g[nCone][0],  wcl[1], "C");
  DrawGraph(g[nCone][2],  wcl[2], "C");
  DrawGraph(g[nCone][3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[nCone], "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(nCone == 0) tex->DrawLatex(0.16, 0.82, "Inclusive");
  if(nCone == 1) tex->DrawLatex(0.16, 0.82, "JE");
  if(nCone == 2) tex->DrawLatex(0.16, 0.82, "UE");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/pp13TeV_%s_toLRatio_DataPy_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_toLRatio_DataPy_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);
  
  return;  
}

void  DrawCones(TString sType = "Omega", Int_t nMode = 0)
{
  auto can(MakeCanvas(Form("PyDataMode_%d", nMode)));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");

  DrawGraph(g[0][nMode],  wcl[0], "C");
  DrawGraph(g[1][nMode],  wcl[1], "C");
  DrawGraph(g[2][nMode],  wcl[2], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->AddEntry(g[0][nMode], "PYTHIA", "L")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(nMode == 1)tex->DrawLatex(0.16, 0.82, "PYTHIA 8 Monash");
  if(nMode == 0)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 0");
  if(nMode == 2)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 2");
  if(nMode == 3)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 3");
  can->SaveAs(Form("./figure/pp13TeV_%s_toLRatio_DataPy_Mode%d.eps", sType.Data(), nMode));
  can->SaveAs(Form("./figure/pp13TeV_%s_toLRatio_DataPy_Mode%d.png", sType.Data(), nMode));
  CanvasEnd(can);

  return;
}

