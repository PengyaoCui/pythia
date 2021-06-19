#include "inc/PyJetUtils.h"

void d07_JELKRatio_DataPy_Allmodle_pp(){
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

  TGraph *g[nm+2];
  for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioLK(i, sd, sj, "JC04", "PC04"));
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);
  auto rf = TFile::Open("./data/Rope/Hard_J.root", "read");
  auto rl = (TList*)rf->Get(Form("list_results"));
  rf->Close();

  const Double_t dPtBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5.0, 8.0, 12., 15. };
  const auto nPtBin(sizeof(dPtBin) / sizeof(Double_t) - 1);

  TH1D* rh[2];
  rh[0] = (TH1D*)rl->FindObject(Form("hKshort_Jet10_C04")); rh[0]->Scale(2.); auto hK(rh[0]->Rebin(nPtBin, "rh0", dPtBin));
  rh[1] = (TH1D*)rl->FindObject(Form("hLambda_Jet10_C04"));
  auto hA = (TH1D*)rl->FindObject(Form("hAntiLa_Jet10_C04"));
  rh[1]->Add(hA);
  
  //rh[1]->Rebin(5);
  auto hL(rh[1]->Rebin(nPtBin, "rh1", dPtBin));

  hL->Divide(hK);
  //rh[1]->Divide(rh[0]);
  auto rF = TFile::Open("./data/Rope/Soft.root", "read");
  auto rL = (TList*)rF->Get(Form("list_results"));
  rF->Close();
  
  TH1D* rH[2];
  rH[0] = (TH1D*)rL->FindObject(Form("hKshort_Jet10_C04")); rH[0]->Scale(2.); auto hk(rH[0]->Rebin(nPtBin, "rH0", dPtBin));
  rH[1] = (TH1D*)rL->FindObject(Form("hLambda_Jet10_C04"));
  auto ha = (TH1D*)rL->FindObject(Form("hAntiLa_Jet10_C04"));
  rH[1]->Add(ha);
  
  auto hl(rH[1]->Rebin(nPtBin, "rH1", dPtBin));
  hl->Divide(hk);


  g[nm]= new TGraph(hl);
  g[nm+1]= new TGraph(hL);

//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(1.0);
   
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");
  
  SetStyle(kTRUE);
  gStyle->SetErrorX(0); 

  auto can(MakeCanvas(Form("Ratio_PyData")));
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
  tex->DrawLatex(0.16, 0.82, "#frac{#Lambda + #bar{#Lambda}}{ 2 K^{0}_{S}} in jet");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/07_pp13TeV_Lambda_sum_toKRatio_DataPy_Allmodle_Cone3.eps"));
  can->SaveAs(Form("./figure/07_pp13TeV_Lambda_sum_toKRatio_DataPy_Allmodle_Cone3.pdf"));
  can->SaveAs(Form("./figure/07_pp13TeV_Lambda_sum_toKRatio_DataPy_Allmodle_Cone3.png"));
  CanvasEnd(can);

  return;
}

