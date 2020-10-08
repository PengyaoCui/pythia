#include "inc/PyJetUtils.h"

TH1D *h[4]; TGraph *gE[4]; TGraph *g[4][nm];
auto dflx(0.), dfux(12.);
auto dfly(1e-6), dfuy(1e2);
 

auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("d#it{#rho}/d#it{p}_{T}: K^{0}_{S}");

void  DrawMode(TString sType = "Lambda_sum", Int_t nCone = 0);
void  DrawCones(TString sType = "Lambda_sum",Int_t nMode = 0);

void PtSpect_DataPy(TString sType = "Kshort", Bool_t IsDiffCone = kFALSE){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s", sType.Data()));
  f->Close();
  h[0] = (TH1D*)l->FindObject(Form("InclCen")); gE[0] = (TGraphErrors*)l->FindObject(Form("Inclerr"));
  h[1] = (TH1D*)l->FindObject(Form("JCCen"));   gE[1] = (TGraphErrors*)l->FindObject(Form("JCerr"));
  h[2] = (TH1D*)l->FindObject(Form("UECen"));   gE[2] = (TGraphErrors*)l->FindObject(Form("UEerr"));
  h[3] = (TH1D*)l->FindObject(Form("JECen"));   gE[3] = (TGraphErrors*)l->FindObject(Form("JEerr"));


  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(Spectrum(i, sd, "Lambda"));
      g[1][i] = new TGraph(Spectrum(i, sd, "Lambda", sj, "JC04"));
      g[2][i] = new TGraph(Spectrum(i, sd, "Lambda", sj, "OC08"));
      g[3][i] = new TGraph(Spectrum(i, sd, "Lambda", sj, "JC04", "PC04"));
    }
  }
  if(sType != "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(Spectrum(i, sd, sType));
      g[1][i] = new TGraph(Spectrum(i, sd, sType, sj, "JC04"));
      g[2][i] = new TGraph(Spectrum(i, sd, sType, sj, "OC08"));
      g[3][i] = new TGraph(Spectrum(i, sd, sType, sj, "JC04", "PC04"));
    }
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  if(sType == "Lambda_sum"){ 
    dfly = 1e-7;
    dfuy = 5e2;
    stny = "d#it{#rho}/d#it{p}_{T}: #Lambda + #bar{#Lambda}";
  }
  if(sType == "Xi"){
    dfux = 8.;
    dfly = 1e-6;
    dfuy = 1e1;
    stny = "d#it{#rho}/d#it{p}_{T}: #Xi^{-} + #bar{#Xi}^{+}";
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfly = 1e-6;
    dfuy = 1e-1;
    stny = "d#it{#rho}/d#it{p}_{T}: #Omega^{-} + #bar{#Omega}^{+}";
  }
  
  SetStyle(kTRUE);
  if(!IsDiffCone)for(Int_t i = 0; i< 4; i++) DrawMode(sType, i);
  if(IsDiffCone)for(Int_t i = 0; i< 4; i++) DrawCones(sType, i);
//=============================================================================
  return;
}
//=============================================================================

void  DrawMode(TString sType = "Lambda_sum", Int_t nCone = 0)
{
  auto can(MakeCanvas(Form("PyDataCone_%d", nCone)));
  can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[nCone], wcl[0], wmk[0], "same"); DrawGraph(gE[nCone], wcl[0], "E2");
  
  DrawGraph(g[nCone][1],  wcl[0], "C");
  DrawGraph(g[nCone][0],  wcl[1], "C");
  DrawGraph(g[nCone][2],  wcl[2], "C");
  DrawGraph(g[nCone][3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[nCone], "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[nCone][3], "mode 3", "L")->SetTextSizePixels(24);
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
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPy_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPy_Cone%d.pdf", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPy_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);
  
  return;  
}

void  DrawCones(TString sType = "Lambda_sum", Int_t nMode = 0)
{
  auto can(MakeCanvas(Form("PyDataMode_%d", nMode)));
  can->SetLogy();
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0], wcl[0], wmk[0], "same"); DrawGraph(gE[0], wcl[0], "E2");
  DrawHisto(h[1], wcl[1], wmk[1], "same"); DrawGraph(gE[1], wcl[1], "E2");
  DrawHisto(h[2], wcl[2], wmk[2], "same"); DrawGraph(gE[2], wcl[2], "E2");
  DrawHisto(h[3], wcl[3], wmk[3], "same"); DrawGraph(gE[3], wcl[3], "E2");

  DrawGraph(g[0][nMode],  wcl[0], "C");
  DrawGraph(g[1][nMode],  wcl[1], "C");
  DrawGraph(g[2][nMode],  wcl[2], "C");
  DrawGraph(g[3][nMode],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[3], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->AddEntry(g[0][nMode], "PYTHIA", "L")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(nMode == 1)tex->DrawLatex(0.16, 0.82, "PYTHIA 8 Monash");
  if(nMode == 0)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 0");
  if(nMode == 2)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 2");
  if(nMode == 3)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 3");
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPy_Mode%d.eps", sType.Data(), nMode));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPy_Mode%d.eps", sType.Data(), nMode));
  can->SaveAs(Form("./figure/pp13TeV_%s_DataPy_Mode%d.png", sType.Data(), nMode));
  CanvasEnd(can);

  return;
}

