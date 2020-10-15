#include "inc/PyJetUtils.h"

void JE_XKRatio_pPb();
void JE_XLRatio_pPb();

void d19_JE_XPRatio_DataPy_pPb(){
  JE_XKRatio_pPb();
  JE_XLRatio_pPb();

}

void JE_XKRatio_pPb(){

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TH1D *h[4]; TGraph *gE[4]; 
  TList* l[4];
  l[0] = (TList*)f->Get(Form("Xi_toKRatio_0100"));
  l[1] = (TList*)f->Get(Form("Xi_toKRatio_010"));
  l[2] = (TList*)f->Get(Form("Xi_toKRatio_1040"));
  l[3] = (TList*)f->Get(Form("Xi_toKRatio_40100"));
  f->Close();
  for(Int_t i = 0; i<4; i++){
    h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));   gE[i] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }
  
  const TString sj("Jet10");
  const TString sd("pp05d02TeVrs");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance
  TGraph *g[nm];
  for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioXK(i, sd, sj, "JC04", "PC04"));
  
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(8.);
  auto dfly(0.), dfuy(0.2);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / 2K_{S}^{0}");
  
  SetStyle(kTRUE);
  
  auto can(MakeCanvas(Form("Ratio_XK")));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");
  
  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  DrawGraph(g[2],  wcl[2], "C");
  DrawGraph(g[3],  wcl[3], "C");

  auto leg(new TLegend(0.52, 0.50, 0.98, 0.92)); SetupLegend(leg);
  leg -> SetNColumns(2); 
 
  leg->AddEntry(h[0], "DATA", "");                             leg->AddEntry(g[1], "PYTHIA", "");
  leg->AddEntry(h[0], "0-100%", "LP")->SetTextSizePixels(24);  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "0-10%", "LP")->SetTextSizePixels(24);   leg->AddEntry(g[0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "10-40%", "LP")->SetTextSizePixels(24);  leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "40-100%", "LP")->SetTextSizePixels(24); leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, "#frac{#Xi^{-} + #bar{#Xi}^{+}}{2K_{S}^{0}} in jet");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/19_pPb5d02TeV_Xi_toKRatio_DataPy_Cone3.eps"));
  can->SaveAs(Form("./figure/19_pPb5d02TeV_Xi_toKRatio_DataPy_Cone3.pdf"));
  can->SaveAs(Form("./figure/19_pPb5d02TeV_Xi_toKRatio_DataPy_Cone3.png"));
  CanvasEnd(can);
  
  return;  
}

void JE_XLRatio_pPb(){

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TH1D *h[4]; TGraph *gE[4];
  TList* l[4];
  l[0] = (TList*)f->Get(Form("Xi_toLRatio_0100"));
  l[1] = (TList*)f->Get(Form("Xi_toLRatio_010"));
  l[2] = (TList*)f->Get(Form("Xi_toLRatio_1040"));
  l[3] = (TList*)f->Get(Form("Xi_toLRatio_40100"));
  f->Close();
  for(Int_t i = 0; i<4; i++){
    h[i] = (TH1D*)l[i]->FindObject(Form("hJER"));   gE[i] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }

  const TString sj("Jet10");
  const TString sd("pp05d02TeVrs");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
                                    //               used to fit p-Pb acceptance
  TGraph *g[nm];
  for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioXL(i, sd, sj, "JC04", "PC04"));

  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(8.);
  auto dfly(0.), dfuy(0.5);

  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);

  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / #Lambda + #bar{#Lambda}");

  SetStyle(kTRUE);

  auto can(MakeCanvas(Form("Ratio_XL")));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  DrawGraph(g[2],  wcl[2], "C");
  DrawGraph(g[3],  wcl[3], "C");

  auto leg(new TLegend(0.52, 0.50, 0.98, 0.92)); SetupLegend(leg);
  leg -> SetNColumns(2);
  leg->AddEntry(h[0], "DATA", "");                             leg->AddEntry(g[1], "PYTHIA", "");
  leg->AddEntry(h[0], "0-100%", "LP")->SetTextSizePixels(24);  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "0-10%", "LP")->SetTextSizePixels(24);   leg->AddEntry(g[0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "10-40%", "LP")->SetTextSizePixels(24);  leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "40-100%", "LP")->SetTextSizePixels(24); leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV");
  tex->DrawLatex(0.16, 0.82, "#frac{#Xi^{-} + #bar{#Xi}^{+}}{#Lambda + #bar{#Lambda}} in jet");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/19_pPb5d02TeV_Xi_toLRatio_DataPy_Cone3.eps"));
  can->SaveAs(Form("./figure/19_pPb5d02TeV_Xi_toLRatio_DataPy_Cone3.pdf"));
  can->SaveAs(Form("./figure/19_pPb5d02TeV_Xi_toLRatio_DataPy_Cone3.png"));
  CanvasEnd(can);

  return;
}

