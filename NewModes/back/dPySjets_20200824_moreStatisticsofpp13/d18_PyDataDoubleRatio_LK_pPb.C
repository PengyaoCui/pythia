#include "inc/PyJetUtils.h"


void d18_PyDataDoubleRatio_LK_pPb(){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  auto l = (TList*)f->Get(Form("Lambda_sum_toKRatio_0100"));
  f->Close();
  auto h = (TH1D*)l->FindObject(Form("hJER")); auto hE = (TH1D*)l->FindObject(Form("heJER"));//Ratio relative sys error
  
  const TString sj("Jet10");
  const TString sd("pp05d02TeVrs");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  TH1D *hPy[nm]; TH1D *hR[nm]; TGraph *gER[nm];
  for (auto i=0; i<nm; ++i) hPy[i] = RatioLK(i, sd, sj, "JC04", "PC04");
  
  for (auto i=0; i<nm; ++i) {
    hR[i] = (TH1D*)hPy[i]->Clone(Form("Ratio_%d", i));
    
    DeNormBinningHistogram(hR[i]);
    
    hR[i] = MakeRebinTH1D(hR[i], h);
    NormBinningHistogram(hR[i]);

    hR[i]->Divide(h);
    gER[i] = (TGraphErrors*)ConvHistogramToGraphErrors(hR[i], hE, hE->GetNbinsX()); gER[i]->SetName(Form("gER_%d", i));
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);


auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(3.5);


auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("PYTHIA/DATA");
TString r("#frac{#Lambda + #bar{#Lambda}}{2 K^{0}_{S}}");


//=============================================================================
  SetStyle(kTRUE);

  auto can(MakeCanvas(Form("DataPyDoubleRatio")));
  //can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  //DrawHisto(h[nCone], wcl[0], wmk[0], "same"); DrawGraph(gE[nCone], wcl[0], "E2");
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
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, Form("JE %s", r.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/18_pPb5d02TeV_DataPyDoubleRatio_Lambda_sum_toKRatio_Cone3.eps"));
  can->SaveAs(Form("./figure/18_pPb5d02TeV_DataPyDoubleRatio_Lambda_sum_toKRatio_Cone3.pdf"));
  can->SaveAs(Form("./figure/18_pPb5d02TeV_DataPyDoubleRatio_Lambda_sum_toKRatio_Cone3.png"));
  CanvasEnd(can);
  
  return;  
}

