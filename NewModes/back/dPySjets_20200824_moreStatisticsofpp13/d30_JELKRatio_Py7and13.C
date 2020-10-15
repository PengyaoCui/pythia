#include "inc/PyJetUtils.h"

void d30_JELKRatio_Py7and13(){
//=============================================================================
  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("Lambda_sum_toKRatio"));
  f->Close();
  
  auto h = (TH1D*)l->FindObject(Form("hJER")); auto gE = (TGraphErrors*)l->FindObject(Form("JERerr"));

  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  TGraph *g[2];
   g[0] = new TGraph(RatioLK(3, sd, sj, "JC04", "PC04"));
   g[1] = new TGraph(RatioLK(3, "pp07d00TeV", sj, "JC04", "PC04"));
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

  DrawHisto(h, wcl[0], wmk[0], "same"); DrawGraph(gE, wcl[0], "E2");
  
  DrawGraph(g[0],  wcl[0], "C");
  DrawGraph(g[1],  wcl[1], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "13 TeV", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "7 TeV", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, "#frac{#Lambda + #bar{#Lambda}}{ 2 K^{0}_{S}} in jet");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8 BLC mode 3");
  can->SaveAs(Form("./figure/30_pp13TeV_pp7TeV_Lambda_sum_toKRatio_DataPy_Cone3.eps"));
  can->SaveAs(Form("./figure/30_pp13TeV_pp7TeV_Lambda_sum_toKRatio_DataPy_Cone3.pdf"));
  can->SaveAs(Form("./figure/30_pp13TeV_pp7TeV_Lambda_sum_toKRatio_DataPy_Cone3.png"));
  CanvasEnd(can);

  return;
}

