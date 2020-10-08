#include "inc/PyJetUtils.h"

TH1D *h[4];
TH1D *hE[4];  
TH1D *hPy[4][nm];
TH1D *hR[4][nm];
TGraph *gER[4][nm]; 
auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(1.3);
 

auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("PYTHIA/DATA : K^{0}_{S}");

void  DrawMode(TString sType = "Kshort", Int_t nCone = 0);

void DataPySpectRatio(TString sType = "Kshort"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s", sType.Data()));
  f->Close();
  h[0] = (TH1D*)l->FindObject(Form("InclCen")); hE[0] = (TH1D*)l->FindObject(Form("InclErr"));//Spectra relative sys error
  h[1] = (TH1D*)l->FindObject(Form("JCCen"));   hE[1] = (TH1D*)l->FindObject(Form("JCErr"));
  h[2] = (TH1D*)l->FindObject(Form("UECen"));   hE[2] = (TH1D*)l->FindObject(Form("UEErr"));
  h[3] = (TH1D*)l->FindObject(Form("JECen"));   hE[3] = (TH1D*)l->FindObject(Form("JEErr"));


  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      hPy[0][i] = Spectrum(i, sd, "Lambda");
      hPy[1][i] = Spectrum(i, sd, "Lambda", sj, "JC04");
      hPy[2][i] = Spectrum(i, sd, "Lambda", sj, "OC08");
      hPy[3][i] = Spectrum(i, sd, "Lambda", sj, "JC04", "PC04");
    }
  }
  if(sType != "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      hPy[0][i] = Spectrum(i, sd, sType);
      hPy[1][i] = Spectrum(i, sd, sType, sj, "JC04");
      hPy[2][i] = Spectrum(i, sd, sType, sj, "OC08");
      hPy[3][i] = Spectrum(i, sd, sType, sj, "JC04", "PC04");
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
  if(sType == "Kshort"){ dfly = 0.2;}
  if(sType == "Lambda_sum"){ 
    //dfly = 1e-7;
    dfuy = 2.;
    stny = "PYTHIA/DATA : #Lambda + #bar{#Lambda}";
  }
  if(sType == "Xi"){
    dfux = 8.;
    //dfly = 1e-6;
    //dfuy = 1e1;
    stny = "PYTHIA/DATA : #Xi^{-} + #bar{#Xi}^{+}";
  }
  if(sType == "Omega"){
    dfux = 5.;
    //dfly = 1e-6;
    dfuy = 0.6;
    stny = "PYTHIA/DATA : #Omega^{-} + #bar{#Omega}^{+}";
  }
  
  SetStyle(kTRUE);
  //for(Int_t i = 0; i< 4; i++) DrawMode(sType, i);
  DrawMode(sType, 0);
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
  if(nCone == 0) tex->DrawLatex(0.16, 0.82, "Inclusive");
  if(nCone == 1) tex->DrawLatex(0.16, 0.82, "JC");
  if(nCone == 2) tex->DrawLatex(0.16, 0.82, "UE");
  if(nCone == 3) tex->DrawLatex(0.16, 0.82, "JE");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPyRatio_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPyRatio_Cone%d.pdf", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPyRatio_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);
  
  return;  
}

