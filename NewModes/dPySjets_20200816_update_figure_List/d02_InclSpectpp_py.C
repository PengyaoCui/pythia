#include "inc/PyJetUtils.h"



void InclSpectpp_py(TString sType = "Kshort");
void d02_InclSpectpp_py(){
  InclSpectpp_py("Kshort");
  InclSpectpp_py("Lambda_sum");
  InclSpectpp_py("Xi");
  InclSpectpp_py("Omega");

}
void InclSpectpp_py(TString sType = "Kshort"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s", sType.Data()));
  f->Close();
  auto h = (TH1D*)l->FindObject(Form("InclCen")); auto gE = (TGraphErrors*)l->FindObject(Form("Inclerr"));

  TGraph *g[nm];
  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[i] = new TGraph(Spectrum(i, sd, "Lambda"));
    }
  }
  if(sType != "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[i] = new TGraph(Spectrum(i, sd, sType));
    }
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(1e-6), dfuy(1e2);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("d#it{#rho}/d#it{p}_{T}: K^{0}_{S}");

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
  
  auto can(MakeCanvas(Form("PyData_%s", sType.Data())));
  can->SetLogy();
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
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("Inclusive K^{0}_{S}"));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("Inclusive #Lambda + #bar{#Lambda}"));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("Inclusive #%s^{-} + #bar{#%s}^{-}", sType.Data(), sType.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/02_pp13TeV_%s_DataPy_Cone0.eps", sType.Data()));
  can->SaveAs(Form("./figure/02_pp13TeV_%s_DataPy_Cone0.pdf", sType.Data()));
  can->SaveAs(Form("./figure/02_pp13TeV_%s_DataPy_Cone0.png", sType.Data()));
  CanvasEnd(can);
  
  return;  
}
