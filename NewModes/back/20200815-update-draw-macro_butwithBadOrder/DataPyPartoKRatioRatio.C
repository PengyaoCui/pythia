#include "inc/PyJetUtils.h"

TH1D *h[4];
TH1D *hE[4];  
TH1D *hPy[4][nm];
TH1D *hR[4][nm];
TGraph *gER[4][nm]; 
auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(2.5);
 

auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("PYTHIA/DATA");
TString r("#frac{#Lambda + #bar{#Lambda}}{2 K^{0}_{S}}");

void  DrawMode(TString sType = "Lambda_sum", Int_t nCone = 0);

void DataPyPartoKRatioRatio(TString sType = "Lambda_sum"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio", sType.Data()));
  f->Close();
  h[0] = (TH1D*)l->FindObject(Form("hInclR")); hE[0] = (TH1D*)l->FindObject(Form("heInclR"));//Ratio relative sys error
  h[1] = (TH1D*)l->FindObject(Form("hJCR"));   hE[1] = (TH1D*)l->FindObject(Form("heJCR"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   hE[2] = (TH1D*)l->FindObject(Form("heUER"));
  h[3] = (TH1D*)l->FindObject(Form("hJER"));   hE[3] = (TH1D*)l->FindObject(Form("heJER"));


  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      hPy[0][i] = RatioLK(i, sd);
      hPy[1][i] = RatioLK(i, sd, sj, "JC04");
      hPy[2][i] = RatioLK(i, sd, sj, "OC08");
      hPy[3][i] = RatioLK(i, sd, sj, "JC04", "PC04");
    }
  }
  if(sType == "Xi"){
    for (auto i=0; i<nm; ++i){
      hPy[0][i] = RatioXK(i, sd);
      hPy[1][i] = RatioXK(i, sd, sj, "JC04");
      hPy[2][i] = RatioXK(i, sd, sj, "OC08");
      hPy[3][i] = RatioXK(i, sd, sj, "JC04", "PC04");
    }
  }
  if(sType == "Omega"){
    for (auto i=0; i<nm; ++i){
      hPy[0][i] = RatioOK(i, sd);
      hPy[1][i] = RatioOK(i, sd, sj, "JC04");
      hPy[2][i] = RatioOK(i, sd, sj, "OC08");
      hPy[3][i] = RatioOK(i, sd, sj, "JC04", "PC04");
    }
  }
  
  for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) {
    hR[j][i] = (TH1D*)hPy[j][i]->Clone(Form("Ratio_%d%d", j, i));
    
    DeNormBinningHistogram(hR[j][i]);
    hR[j][i] = MakeRebinTH1D(hR[j][i], h[j]);
    NormBinningHistogram(hR[j][i]);
    hR[j][i]->Divide(h[j]);
    gER[j][i] = (TGraphErrors*)ConvHistogramToGraphErrors(hR[j][i], hE[j], hE[j]->GetNbinsX()); gER[j][i]->SetName(Form("gER_%d%d", i, j));
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
   if(sType == "Xi"){
    dfux = 8.;
    r = "#frac{#Xi^{-} + #bar{#Xi}^{+}}{2 K^{0}_{S}}";
    //dfuy = 0.2;
  }
  if(sType == "Omega"){
    dfux = 5.;
    r = "#frac{#Omega^{-} + #bar{#Omega}^{+}}{2 K^{0}_{S}}";
    //dfuy = 0.05;
  }
 
  SetStyle(kTRUE);
  for(Int_t i = 0; i< 4; i++) DrawMode(sType, i);
  //DrawMode(sType, 0);
//=============================================================================
  return;
}
//=============================================================================

void  DrawMode(TString sType = "Lambda_sum", Int_t nCone = 0)
{
  auto can(MakeCanvas(Form("PyDataCone_%d", nCone)));
  //can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  //DrawHisto(h[nCone], wcl[0], wmk[0], "same"); DrawGraph(gE[nCone], wcl[0], "E2");
  DrawHisto(hR[nCone][1], wcl[0], wmk[0], "same"); DrawGraph(gER[nCone][1], wcl[0], "E2");
  DrawHisto(hR[nCone][0], wcl[1], wmk[1], "same"); DrawGraph(gER[nCone][0], wcl[1], "E2");
  DrawHisto(hR[nCone][2], wcl[2], wmk[2], "same"); DrawGraph(gER[nCone][2], wcl[2], "E2");
  DrawHisto(hR[nCone][3], wcl[3], wmk[3], "same"); DrawGraph(gER[nCone][3], wcl[3], "E2");
  
  //DrawGraph(g[nCone][1],  wcl[0], "C");
  //DrawGraph(g[nCone][0],  wcl[1], "C");
  //DrawGraph(g[nCone][2],  wcl[2], "C");
  //DrawGraph(g[nCone][3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(hR[nCone][1], "Monash", "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[nCone][0], "mode 0", "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[nCone][2], "mode 2", "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[nCone][3], "mode 3", "LP")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(nCone == 0) tex->DrawLatex(0.16, 0.82, Form("Inclusive %s", r.Data()));
  if(nCone == 1) tex->DrawLatex(0.16, 0.82, Form("JC %s", r.Data()));
  if(nCone == 2) tex->DrawLatex(0.16, 0.82, Form("UE %s", r.Data()));
  if(nCone == 3) tex->DrawLatex(0.16, 0.82, Form("JE %s", r.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPyRatioPartoKRatio_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPyRatioPartoKRatio_Cone%d.pdf", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPyRatioPartoKRatio_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);
  
  return;  
}

