#include "inc/PyJetUtils.h"

void d24_InclJELKRatio_DataPy_pp(){
//=============================================================================
  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("Lambda_sum_toKRatio"));
  f->Close();
  TH1D* h[2]; TGraph *gE[2];
  h[0] = (TH1D*)l->FindObject(Form("hInR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJER"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JERerr"));

  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  TGraph *gI = new TGraph(RatioLK(3, sd));
  TGraph *gJ = new TGraph(RatioLK(3, sd, sj, "JC04", "PC04"));
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(1.0);
   
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");
  
  SetStyle(kTRUE);

  auto can(MakeCanvas(Form("Ratio_PyData")));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[3], wmk[3], "same"); DrawGraph(gE[1], wcl[3], "E2");
  
  DrawGraph(gI,  wcl[0], "C");
  DrawGraph(gJ,  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gI, "PYTHIA", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, "#frac{#Lambda + #bar{#Lambda}}{ 2 K^{0}_{S}} in jet");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8 BLC 3");
  can->SaveAs(Form("./figure/24_pp13TeV_Lambda_sum_toKRatio_DataPy_Cone0_3.eps"));
  can->SaveAs(Form("./figure/24_pp13TeV_Lambda_sum_toKRatio_DataPy_Cone0_3.pdf"));
  can->SaveAs(Form("./figure/24_pp13TeV_Lambda_sum_toKRatio_DataPy_Cone0_3.png"));
  CanvasEnd(can);

  return;
}

