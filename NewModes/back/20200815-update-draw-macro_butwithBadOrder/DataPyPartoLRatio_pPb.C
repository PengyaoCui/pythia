#include "inc/PyJetUtils.h"


TH1D *h[4][4]; TGraph *gE[4][4]; TGraph *g[4][nm];
auto dflx(0.), dfux(12.);
auto dfly(0.), dfuy(1.0);
 

auto dlsx(0.05), dlsy(0.05);
auto dtsx(0.05), dtsy(0.05);
auto dtox(1.30), dtoy(1.10);

TString stnx("#it{p}_{T} (GeV/#it{c})");
TString stny("Ratio: (#Omega^{+} + #bar{#Omega}^{-}) / #Lambda + #bar{#Lambda}");

void  DrawMode(TString sType = "Omega", Int_t nCone = 0);
void  DrawCones(TString sType = "Omega",Int_t nMode = 0);

void DataPyPartoLRatio_pPb(TString sType = "Omega", Bool_t IsDiffCone = kTRUE){
//=============================================================================

  auto f = TFile::Open("./data/Data/pPb.root", "read");
  TList* l[4];
  l[0] = (TList*)f->Get(Form("%s_toLRatio_0100", sType.Data()));
  l[1] = (TList*)f->Get(Form("%s_toLRatio_010", sType.Data()));
  l[2] = (TList*)f->Get(Form("%s_toLRatio_1040", sType.Data()));
  l[3] = (TList*)f->Get(Form("%s_toLRatio_40100", sType.Data()));
  f->Close();
  for(Int_t i = 0; i<4; i++){
    h[i][0] = (TH1D*)l[i]->FindObject(Form("hInclR")); gE[i][0] = (TGraphErrors*)l[i]->FindObject(Form("InclRerr"));
    h[i][1] = (TH1D*)l[i]->FindObject(Form("hJCR"));   gE[i][1] = (TGraphErrors*)l[i]->FindObject(Form("JCRerr"));
    h[i][2] = (TH1D*)l[i]->FindObject(Form("hUER"));   gE[i][2] = (TGraphErrors*)l[i]->FindObject(Form("UERerr"));
    h[i][3] = (TH1D*)l[i]->FindObject(Form("hJER"));   gE[i][3] = (TGraphErrors*)l[i]->FindObject(Form("JERerr"));
  }


  const TString sj("Jet10");
  const TString sd("pp05d02TeVrs");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Xi"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioXL(i,sd));
      g[1][i] = new TGraph(RatioXL(i, sd, sj, "JC04"));
      g[2][i] = new TGraph(RatioXL(i, sd, sj, "OC08"));
      g[3][i] = new TGraph(RatioXL(i, sd, sj, "JC04", "PC04"));
    }
  }

  if(sType == "Omega"){
    for (auto i=0; i<nm; ++i){
      g[0][i] = new TGraph(RatioOL(i,sd));
      g[1][i] = new TGraph(RatioOL(i, sd, sj, "JC04"));
      g[2][i] = new TGraph(RatioOL(i, sd, sj, "OC08"));
      g[3][i] = new TGraph(RatioOL(i, sd, sj, "JC04", "PC04"));
    }
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================

  if(sType == "Xi"){
    dfux = 8.;
    dfuy = 0.5;
    stny = "Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / (#Lambda + #bar{#Lambda})";
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfuy = 0.2;
    stny = "Ratio: (#Omega^{-} + #bar{#Omega}^{+}) / (#Lambda + #bar{#Lambda})";
  }
  
  SetStyle(kTRUE);
  if(!IsDiffCone)for(Int_t i = 0; i< 4; i++) DrawMode(sType, i);
  if(IsDiffCone)for(Int_t i = 0; i< 4; i++) DrawCones(sType, i);
//=============================================================================
  return;
}
//=============================================================================

void  DrawMode(TString sType = "Omega", Int_t nCone = 0)
{
  auto can(MakeCanvas(Form("PyDataCone_%d", nCone)));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);


  DrawHisto(h[0][nCone], wcl[0], wmk[0], "same"); DrawGraph(gE[0][nCone], wcl[0], "E2");
  DrawHisto(h[1][nCone], wcl[1], wmk[1], "same"); DrawGraph(gE[1][nCone], wcl[1], "E2");
  DrawHisto(h[2][nCone], wcl[2], wmk[2], "same"); DrawGraph(gE[2][nCone], wcl[2], "E2");
  DrawHisto(h[3][nCone], wcl[3], wmk[3], "same"); DrawGraph(gE[3][nCone], wcl[3], "E2");

  DrawGraph(g[nCone][1],  wcl[0], "C");
  DrawGraph(g[nCone][0],  wcl[1], "C");
  DrawGraph(g[nCone][2],  wcl[2], "C");
  DrawGraph(g[nCone][3],  wcl[3], "C");

  auto leg(new TLegend(0.52, 0.50, 0.98, 0.92)); SetupLegend(leg);
  leg -> SetNColumns(2);

  leg->AddEntry(h[0][nCone], "DATA", "");                            leg->AddEntry(g[nCone][1], "PYTHIA", "");
  leg->AddEntry(h[0][nCone], "0-100%", "LP")->SetTextSizePixels(24);  leg->AddEntry(g[nCone][1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[1][nCone], "0-10%", "LP")->SetTextSizePixels(24);   leg->AddEntry(g[nCone][0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[2][nCone], "10-40%", "LP")->SetTextSizePixels(24);  leg->AddEntry(g[nCone][2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(h[3][nCone], "40-100%", "LP")->SetTextSizePixels(24); leg->AddEntry(g[nCone][3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV,  V0A : 0-100%");
  if(nCone == 0) tex->DrawLatex(0.16, 0.82, "Inclusive");
  if(nCone == 1) tex->DrawLatex(0.16, 0.82, "JC");
  if(nCone == 2) tex->DrawLatex(0.16, 0.82, "UE");
  if(nCone == 3) tex->DrawLatex(0.16, 0.82, "JE");
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toLRatio_DataPy_Cone%d.eps", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toLRatio_DataPy_Cone%d.pdf", sType.Data(), nCone));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toLRatio_DataPy_Cone%d.png", sType.Data(), nCone));
  CanvasEnd(can);
  
  return;  
}

void  DrawCones(TString sType = "Omega", Int_t nMode = 0)
{
  auto can(MakeCanvas(Form("PyDataMode_%d", nMode)));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawHisto(h[0][0], wcl[0], wmk[0], "same"); DrawGraph(gE[0][0], wcl[0], "E2");
  DrawHisto(h[0][1], wcl[1], wmk[1], "same"); DrawGraph(gE[0][1], wcl[1], "E2");
  DrawHisto(h[0][2], wcl[2], wmk[2], "same"); DrawGraph(gE[0][2], wcl[2], "E2");
  DrawHisto(h[0][3], wcl[3], wmk[3], "same"); DrawGraph(gE[0][3], wcl[3], "E2");

  DrawGraph(g[0][nMode],  wcl[0], "C");
  DrawGraph(g[1][nMode],  wcl[1], "C");
  DrawGraph(g[2][nMode],  wcl[2], "C");
  DrawGraph(g[3][nMode],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h[0][0], "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[0][1], "JC", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[0][2], "UE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(h[0][3], "JE", "LP")->SetTextSizePixels(24);
  leg->AddEntry(gE[0][0], "Sys. Error", "f")->SetTextSizePixels(24);
  leg->AddEntry(g[0][nMode], "PYTHIA", "L")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV,  V0A : 0-100%");
  if(nMode == 1)tex->DrawLatex(0.16, 0.82, "PYTHIA 8 Monash");
  if(nMode == 0)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 0");
  if(nMode == 2)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 2");
  if(nMode == 3)tex->DrawLatex(0.16, 0.82, "PYTHIA 8, BLC Mode 3");
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toLRatio_DataPy_Mode%d.eps", sType.Data(), nMode));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toLRatio_DataPy_Mode%d.pdf", sType.Data(), nMode));
  can->SaveAs(Form("./figure/pPb5d02TeV_%s_toLRatio_DataPy_Mode%d.png", sType.Data(), nMode));
  CanvasEnd(can);

  return;
}

