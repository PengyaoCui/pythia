#include "inc/PyJetUtils.h"


void  InclSpectDataPyRatio(TString sType = "Kshort");

void d03_InclSpectDataPyRatio(){
  InclSpectDataPyRatio("Kshort");
  InclSpectDataPyRatio("Lambda_sum");
  InclSpectDataPyRatio("Xi");
  InclSpectDataPyRatio("Omega");
}

void InclSpectDataPyRatio(TString sType = "Kshort"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s", sType.Data()));
  f->Close();
  TH1D *hPy[nm]; TH1D *hR[nm]; TGraph *gER[nm]; 
  auto h = (TH1D*)l->FindObject(Form("InclCen")); auto hE = (TH1D*)l->FindObject(Form("InclErr"));//Spectra relative sys error

  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      hPy[i] = Spectrum(i, sd, "Lambda");
    }
  }
  if(sType != "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      hPy[i] = Spectrum(i, sd, sType);
    }
  }
  
  for (auto i=0; i<nm; ++i){
    hR[i] = (TH1D*)hPy[i]->Clone(Form("Ratio_%d", i));
    
    DeNormBinningHistogram(hR[i]);
    hR[i] = MakeRebinTH1D(hR[i], h);
    NormBinningHistogram(hR[i]);
    hR[i]->Divide(h);
    gER[i] = (TGraphErrors*)ConvHistogramToGraphErrors(hR[i], hE, hE->GetNbinsX()); gER[i]->SetName(Form("gER_%d", i));
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(1.3);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("PYTHIA/DATA");
  if(sType == "Kshort"){ dfly = 0.2;}
  if(sType == "Lambda_sum"){ 
    //dfly = 1e-7;
    dfuy = 2.;
  }
  if(sType == "Xi"){
    dfux = 8.;
    //dfly = 1e-6;
    //dfuy = 1e1;
  }
  if(sType == "Omega"){
    dfux = 5.;
    //dfly = 1e-6;
    dfuy = 0.6;
  }
  
  SetStyle(kTRUE);
  
  auto can(MakeCanvas(Form("PyData_%s", sType.Data())));
  //can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(hR[1], wcl[0], wmk[0], "same"); DrawGraph(gER[1], wcl[0], "E2");
  DrawHisto(hR[0], wcl[1], wmk[1], "same"); DrawGraph(gER[0], wcl[1], "E2");
  DrawHisto(hR[2], wcl[2], wmk[2], "same"); DrawGraph(gER[2], wcl[2], "E2");
  DrawHisto(hR[3], wcl[3], wmk[3], "same"); DrawGraph(gER[3], wcl[3], "E2");
  

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(hR[1], "Monash", "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[0], "mode 0", "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[2], "mode 2", "LP")->SetTextSizePixels(24);
  leg->AddEntry(hR[3], "mode 3", "LP")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("Inclusive K^{0}_{S}"));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("Inclusive #Lambda + #bar{#Lambda}"));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("Inclusive #%s^{-} + #bar{#%s}^{-}", sType.Data(), sType.Data()));

  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/03_pp13TeV_%s_DataPyRatio_Cone0.eps", sType.Data()));
  can->SaveAs(Form("./figure/03_pp13TeV_%s_DataPyRatio_Cone0.pdf", sType.Data()));
  can->SaveAs(Form("./figure/03_pp13TeV_%s_DataPyRatio_Cone0.png", sType.Data()));
  CanvasEnd(can);
  
  return;  
}

