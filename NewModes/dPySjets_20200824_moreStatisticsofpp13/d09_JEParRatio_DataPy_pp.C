#include "inc/PyJetUtils.h"

void JEPKRatio_DataPy_pp(TString sType = "Xi");
void JEPLRatio_DataPy_pp(TString sType = "Xi");
void JEPXRatio_DataPy_pp(TString sType = "Omega");
void d09_JEParRatio_DataPy_pp(){
  JEPKRatio_DataPy_pp("Xi");
  JEPKRatio_DataPy_pp("Omega");
  JEPLRatio_DataPy_pp("Xi");
  JEPLRatio_DataPy_pp("Omega");
  JEPXRatio_DataPy_pp("Omega");
}

void JEPKRatio_DataPy_pp(TString sType = "Xi"){
//=============================================================================
  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio", sType.Data()));
  f->Close();
  
  auto h = (TH1D*)l->FindObject(Form("hJER")); auto gE = (TGraphErrors*)l->FindObject(Form("JERerr"));

  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  TGraph *g[nm];
  if(sType == "Xi"){for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioXK(i, sd, sj, "JC04", "PC04"));}
  if(sType == "Omega"){for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioOK(i, sd, sj, "JC04", "PC04"));}
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(8.);
  auto dfly(0.), dfuy(0.2);
   
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / 2K_{S}^{0}");
 
  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 0.03;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / 2K_{S}^{0}";
  }
  
  SetStyle(kTRUE);

  auto can(MakeCanvas(Form("%s_toKRatio_PyData", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h, wcl[0], wmk[0], "same"); DrawGraph(gE, wcl[0], "E2");
  
  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  DrawGraph(g[2],  wcl[2], "C");
  DrawGraph(g[3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{+}}{ 2 K^{0}_{S}} in jet", sType.Data(), sType.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toKRatio_DataPy_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toKRatio_DataPy_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toKRatio_DataPy_Cone3.png", sType.Data()));
  CanvasEnd(can);

  return;
}

void JEPLRatio_DataPy_pp(TString sType = "Xi"){


  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toLRatio", sType.Data()));
  f->Close();
  auto h = (TH1D*)l->FindObject(Form("hJER"));   auto gE = (TGraphErrors*)l->FindObject(Form("JERerr"));


  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  TGraph *g[nm];
  if(sType == "Xi"){
    for (auto i=0; i<nm; ++i){
      g[i] = new TGraph(RatioXL(i, sd, sj, "JC04", "PC04"));
    }
  }

  if(sType == "Omega"){
    for (auto i=0; i<nm; ++i){
      g[i] = new TGraph(RatioOL(i, sd, sj, "JC04", "PC04"));
    }
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(8.);
  auto dfly(0.), dfuy(0.5);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);

  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Xi^{+} + #bar{#Xi}^{-}) / #Lambda + #bar{#Lambda}");

  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 0.1;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / (#Lambda + #bar{#Lambda})";
  }
  
  SetStyle(kTRUE);
  auto can(MakeCanvas(Form("%s_toLRatio_PyData", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h, wcl[0], wmk[0], "same"); DrawGraph(gE, wcl[0], "E2");
  
  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  DrawGraph(g[2],  wcl[2], "C");
  DrawGraph(g[3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{+}}{#Lambda + #bar{#Lambda}} in jet", sType.Data(), sType.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toLRatio_DataPy_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toLRatio_DataPy_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toLRatio_DataPy_Cone3.png", sType.Data()));
  CanvasEnd(can);

  return;
}


void JEPXRatio_DataPy_pp(TString sType = "Omega"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s_toXRatio", sType.Data()));
  f->Close();
  auto h = (TH1D*)l->FindObject(Form("hJER"));   auto gE = (TGraphErrors*)l->FindObject(Form("JERerr"));

  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  TGraph *g[nm];
  for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioOX(i, sd, sj, "JC04", "PC04"));
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(5.);
  auto dfly(0.), dfuy(0.6);
   
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Xi^{-} + #bar{#Xi}^{+}");

  SetStyle(kTRUE);
  auto can(MakeCanvas(Form("%s_toXRatio_PyData", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h, wcl[0], wmk[0], "same"); DrawGraph(gE, wcl[0], "E2");
  
  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  DrawGraph(g[2],  wcl[2], "C");
  DrawGraph(g[3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#Omega^{+} + #bar{#Omega}^{-}}{#Xi^{-} + #bar{#Xi}^{+}} in jet"));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toXRatio_DataPy_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toXRatio_DataPy_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toXRatio_DataPy_Cone3.png", sType.Data()));
  CanvasEnd(can);
  
  return;  
}


