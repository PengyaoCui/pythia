#include "inc/PyJetUtils.h"

void PartoKRatio(TString sType = "Xi"){
//=============================================================================

  TH1D *h[3]; TGraph *gE[3];
  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio", sType.Data()));
  f->Close();
  h[0] = (TH1D*)l->FindObject(Form("hInclR")); gE[0] = (TGraphErrors*)l->FindObject(Form("InclRerr"));
  h[1] = (TH1D*)l->FindObject(Form("hJER"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JERerr"));
  h[2] = (TH1D*)l->FindObject(Form("hUER"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UERerr"));


  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
                                    //               used to fit p-Pb acceptance
  TGraph *g[3][nm];
  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioLK(i,sd));
      g[1][i] = new TGraph(RatioLK(i, sd, sj, "JC04"));
      g[2][i] = new TGraph(RatioLK(i, sd, sj, "JC04", "PC04"));
    }
  }
  if(sType == "Xi"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioXK(i,sd));
      g[1][i] = new TGraph(RatioXK(i, sd, sj, "JC04"));
      g[2][i] = new TGraph(RatioXK(i, sd, sj, "JC04", "PC04"));
    }
  }
  if(sType == "Omega"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioOK(i,sd));
      g[1][i] = new TGraph(RatioOK(i, sd, sj, "JC04"));
      g[2][i] = new TGraph(RatioOK(i, sd, sj, "JC04", "PC04"));
    }
  }
  for (auto i=0; i<nm; ++i) for (auto j=0; j<3; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================

  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(1.0);
   

  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);

  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");

  if(sType == "Xi"){
    dfux = 8.;
    dfuy = 0.2;
    stny = "Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / 2K_{S}^{0}";
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 0.05;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / 2K_{S}^{0}";
  }
//=============================================================================

  SetStyle(kTRUE);
  auto can(MakeCanvas("PyMode"));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);
  
  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2"); 

  DrawGraph(g[0][1],  wcl[0], "C");
  DrawGraph(g[0][0],  wcl[1], "C");
  DrawGraph(g[0][2],  wcl[2], "C");
  DrawGraph(g[0][3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Data", "LP")->SetTextSizePixels(24);
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->AddEntry(g[0][1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0][0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0][2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0][3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "PYTHIA 8");
  CanvasEnd(can);

  return;
}
