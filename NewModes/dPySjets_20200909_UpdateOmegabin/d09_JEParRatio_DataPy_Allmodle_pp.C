#include "inc/PyJetUtils.h"

void JEPKRatio_DataPy_pp(TString sType = "Xi");
void JEPLRatio_DataPy_pp(TString sType = "Xi");
void JEPXRatio_DataPy_pp(TString sType = "Omega");
void d09_JEParRatio_DataPy_Allmodle_pp(){
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

  TGraph *g[nm+2];
  if(sType == "Xi"){for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioXK(i, sd, sj, "JC04", "PC04"));}
  if(sType == "Omega"){for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioOK(i, sd, sj, "JC04", "PC04"));}

  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  const Double_t dPtBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5.0, 8.0, 12., 15. };
  const auto nPtBin(sizeof(dPtBin) / sizeof(Double_t) - 1);

  auto rf = TFile::Open("./data/Rope/Hard_J.root", "read");
  auto rl = (TList*)rf->Get(Form("list_results"));
  rf->Close();
  TH1D* rh[2]; TH1D* hh[2];
  rh[0] = (TH1D*)rl->FindObject(Form("hKshort_Jet10_C04")); rh[0]->Scale(2.);
  hh[0]=(TH1D*)rh[0]->Rebin(nPtBin, "hh", dPtBin);
  if(sType == "Lambda_sum"){
    rh[1] = (TH1D*)rl->FindObject(Form("hLambda_Jet10_C04"));
    auto ha = (TH1D*)rl->FindObject(Form("hAntiLa_Jet10_C04"));
    rh[1]->Add(ha);
  }
  if(sType == "Xi"){
    rh[1] = (TH1D*)rl->FindObject(Form("hXiNeg_Jet10_C04"));
    auto ha = (TH1D*)rl->FindObject(Form("hXiPos_Jet10_C04"));
    rh[1]->Add(ha);
  }
  if(sType == "Omega"){
    rh[1] = (TH1D*)rl->FindObject(Form("hOmegaNeg_Jet10_C04"));
    auto ha = (TH1D*)rl->FindObject(Form("hOmegaPos_Jet10_C04"));
    rh[1]->Add(ha);
  }
  hh[1]=(TH1D*)rh[1]->Rebin(nPtBin, "h1", dPtBin);

  hh[1]->Divide(hh[0]);
  g[nm]= new TGraph(hh[1]);

  auto rF = TFile::Open("./data/Rope/Soft.root", "read");
  auto rL = (TList*)rF->Get(Form("list_results"));
  rF->Close();
  TH1D* rH[2]; TH1D* H[2];
  rH[0] = (TH1D*)rL->FindObject(Form("hKshort_Jet10_C04")); rH[0]->Scale(2.);
  H[0]=(TH1D*)rH[0]->Rebin(nPtBin, "H", dPtBin);
  if(sType == "Lambda_sum"){
    rH[1] = (TH1D*)rL->FindObject(Form("hLambda_Jet10_C04"));
    auto haa = (TH1D*)rL->FindObject(Form("hAntiLa_Jet10_C04"));
    rH[1]->Add(haa);
  }
  if(sType == "Xi"){
    rH[1] = (TH1D*)rL->FindObject(Form("hXiNeg_Jet10_C04"));
    auto haa = (TH1D*)rL->FindObject(Form("hXiPos_Jet10_C04"));
    rH[1]->Add(haa);
  }
  if(sType == "Omega"){
    rH[1] = (TH1D*)rL->FindObject(Form("hOmegaNeg_Jet10_C04"));
    auto haa = (TH1D*)rL->FindObject(Form("hOmegaPos_Jet10_C04"));
    rH[1]->Add(haa);
  }
  H[1]=(TH1D*)rH[1]->Rebin(nPtBin, "H1", dPtBin);

  H[1]->Divide(H[0]);
  g[nm+1]= new TGraph(hh[1]);
  g[nm]= new TGraph(H[1]);

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
  gStyle->SetErrorX(0);

  auto can(MakeCanvas(Form("%s_toKRatio_PyData", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h, wcl[1], wmk[0], "same"); DrawGraph(gE, wcl[1], "E2");
  
  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  //DrawGraph(g[2],  wcl[2], "C");
  //DrawGraph(g[3],  wcl[3], "C");
  DrawGraph(g[4],  wcl[4], "C");
  DrawGraph(g[5],  wcl[5], "C");

  auto leg(new TLegend(0.52, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "BLC mode 0", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[4], "Rope w/ soft QCD", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[5], "Rope w/ hard QCD", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{+}}{ 2 K^{0}_{S}} in jet", sType.Data(), sType.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toKRatio_DataPy_AllModle_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toKRatio_DataPy_AllModle_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toKRatio_DataPy_AllModle_Cone3.png", sType.Data()));
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

  TGraph *g[nm+2];
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
  const Double_t dPtBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5.0, 8.0, 12., 15. };
  const auto nPtBin(sizeof(dPtBin) / sizeof(Double_t) - 1);

  auto rf = TFile::Open("./data/Rope/Hard_J.root", "read");
  auto rl = (TList*)rf->Get(Form("list_results"));
  rf->Close();
  TH1D* rh[2]; TH1D* hh[2];

  rh[0] = (TH1D*)rl->FindObject(Form("hLambda_Jet10_C04"));
  auto hA = (TH1D*)rl->FindObject(Form("hAntiLa_Jet10_C04"));
  rh[0]->Add(hA);
  hh[0]=(TH1D*)rh[0]->Rebin(nPtBin, "hh", dPtBin);

  if(sType == "Xi"){
    rh[1] = (TH1D*)rl->FindObject(Form("hXiNeg_Jet10_C04"));
    auto ha = (TH1D*)rl->FindObject(Form("hXiPos_Jet10_C04"));
    rh[1]->Add(ha);
  }
  if(sType == "Omega"){
    rh[1] = (TH1D*)rl->FindObject(Form("hOmegaNeg_Jet10_C04"));
    auto ha = (TH1D*)rl->FindObject(Form("hOmegaPos_Jet10_C04"));
    rh[1]->Add(ha);
  }
  hh[1]=(TH1D*)rh[1]->Rebin(nPtBin, "h1", dPtBin);

  hh[1]->Divide(hh[0]);
  g[nm]= new TGraph(hh[1]);

  auto rF = TFile::Open("./data/Rope/Soft.root", "read");
  auto rL = (TList*)rF->Get(Form("list_results"));
  rF->Close();
  TH1D* rH[2]; TH1D* H[2];
  rH[0] = (TH1D*)rL->FindObject(Form("hLambda_Jet10_C04"));
  auto hAA = (TH1D*)rL->FindObject(Form("hAntiLa_Jet10_C04"));
  rH[0]->Add(hAA);
  H[0]=(TH1D*)rH[0]->Rebin(nPtBin, "H0", dPtBin);

  if(sType == "Xi"){
    rH[1] = (TH1D*)rL->FindObject(Form("hXiNeg_Jet10_C04"));
    auto haa = (TH1D*)rL->FindObject(Form("hXiPos_Jet10_C04"));
    rH[1]->Add(haa);
  }
  if(sType == "Omega"){
    rH[1] = (TH1D*)rL->FindObject(Form("hOmegaNeg_Jet10_C04"));
    auto haa = (TH1D*)rL->FindObject(Form("hOmegaPos_Jet10_C04"));
    rH[1]->Add(haa);
  }
  H[1]=(TH1D*)rH[1]->Rebin(nPtBin, "H1", dPtBin);

  H[1]->Divide(H[0]);
  g[nm+1]= new TGraph(hh[1]);
  g[nm]= new TGraph(H[1]);

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
  gStyle->SetErrorX(0);
  auto can(MakeCanvas(Form("%s_toLRatio_PyData", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h, wcl[1], wmk[0], "same"); DrawGraph(gE, wcl[1], "E2");
  
  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  //DrawGraph(g[2],  wcl[2], "C");
  //DrawGraph(g[3],  wcl[3], "C");
  DrawGraph(g[4],  wcl[4], "C");
  DrawGraph(g[5],  wcl[5], "C");

  auto leg(new TLegend(0.52, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "BLC mode 0", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[4], "Rope w/ soft QCD", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[5], "Rope w/ hard QCD", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#%s^{-} + #bar{#%s}^{+}}{#Lambda + #bar{#Lambda}} in jet", sType.Data(), sType.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toLRatio_DataPy_AllModle_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toLRatio_DataPy_AllModle_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toLRatio_DataPy_AllModle_Cone3.png", sType.Data()));
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

  TGraph *g[nm+2];
  for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioOX(i, sd, sj, "JC04", "PC04"));
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

  const Double_t dPtBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5.0, 8.0, 12., 15. };
  const auto nPtBin(sizeof(dPtBin) / sizeof(Double_t) - 1);
  
  auto rf = TFile::Open("./data/Rope/Hard_J.root", "read");
  auto rl = (TList*)rf->Get(Form("list_results"));
  rf->Close();
  TH1D* rh[2]; TH1D* hh[2];
  rh[0] = (TH1D*)rl->FindObject(Form("hXiNeg_Jet10_C04"));
  auto hA = (TH1D*)rl->FindObject(Form("hXiPos_Jet10_C04"));
  rh[0]->Add(hA);
  hh[0]=(TH1D*)rh[0]->Rebin(nPtBin, "hh", dPtBin);
 
  rh[1] = (TH1D*)rl->FindObject(Form("hOmegaNeg_Jet10_C04"));
  auto ha = (TH1D*)rl->FindObject(Form("hOmegaPos_Jet10_C04"));
  rh[1]->Add(ha); 
  hh[1]=(TH1D*)rh[1]->Rebin(nPtBin, "h1", dPtBin);

  hh[1]->Divide(hh[0]);
  g[nm]= new TGraph(hh[1]);

  auto rF = TFile::Open("./data/Rope/Soft.root", "read");
  auto rL = (TList*)rF->Get(Form("list_results"));
  rF->Close();
  TH1D* rH[2]; TH1D* H[2];
  rH[0] = (TH1D*)rL->FindObject(Form("hXiNeg_Jet10_C04"));
  auto hAA = (TH1D*)rL->FindObject(Form("hXiPos_Jet10_C04"));
  rH[0]->Add(hAA);
  H[0]=(TH1D*)rH[0]->Rebin(nPtBin, "H0", dPtBin);

  rH[1] = (TH1D*)rL->FindObject(Form("hOmegaNeg_Jet10_C04"));
  auto haa = (TH1D*)rL->FindObject(Form("hOmegaPos_Jet10_C04"));
  rH[1]->Add(haa);
  H[1]=(TH1D*)rH[1]->Rebin(nPtBin, "H1", dPtBin);

  H[1]->Divide(H[0]);
  g[nm+1]= new TGraph(hh[1]);
  g[nm]= new TGraph(H[1]);


//=============================================================================
  auto dflx(0.), dfux(5.);
  auto dfly(0.), dfuy(0.6);
   
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Xi^{-} + #bar{#Xi}^{+}");

  SetStyle(kTRUE);
  gStyle->SetErrorX(0);
  auto can(MakeCanvas(Form("%s_toXRatio_PyData", sType.Data())));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h, wcl[1], wmk[0], "same"); DrawGraph(gE, wcl[1], "E2");
  
  DrawGraph(g[1],  wcl[0], "L");
  DrawGraph(g[0],  wcl[1], "L");
  //DrawGraph(g[2],  wcl[2], "L");
  //DrawGraph(g[3],  wcl[3], "L");
  DrawGraph(g[4],  wcl[4], "L");
  DrawGraph(g[5],  wcl[5], "L");

  auto leg(new TLegend(0.52, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "BLC mode 0", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[4], "Rope w/ soft QCD", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[5], "Rope w/ hard QCD", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.16, 0.82, Form("#frac{#Omega^{+} + #bar{#Omega}^{-}}{#Xi^{-} + #bar{#Xi}^{+}} in jet"));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toXRatio_DataPy_Allmodle_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toXRatio_DataPy_Allmodle_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/09_pp13TeV_%s_toXRatio_DataPy_Allmodle_Cone3.png", sType.Data()));
  CanvasEnd(can);
  
  return;  
}


