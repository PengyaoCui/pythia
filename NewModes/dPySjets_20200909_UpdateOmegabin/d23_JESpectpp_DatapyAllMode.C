#include "inc/PyJetUtils.h"



void JESpectpp_py(TString sType = "Kshort");
void d23_JESpectpp_DatapyAllMode(){
  JESpectpp_py("Kshort");
  JESpectpp_py("Lambda_sum");
  JESpectpp_py("Xi");
  JESpectpp_py("Omega");

}
void JESpectpp_py(TString sType = "Kshort"){
//=============================================================================

  auto f = TFile::Open("./data/Data/pp.root", "read");
  auto l = (TList*)f->Get(Form("%s", sType.Data()));
  f->Close();
  auto h = (TH1D*)l->FindObject(Form("JECen")); auto gE = (TGraphErrors*)l->FindObject(Form("JEerr"));

  TGraph *g[nm+2];
  const TString sj("Jet10");
  const TString sd("pp13d00TeV");  // pp13d00TeV:   pp at 13   TeV
                                    // pp05d02TeVrs: pp at 5.02 TeV
                                    //               with rapidity shift
				    //               used to fit p-Pb acceptance

  if(sType == "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[i] = new TGraph(Spectrum(i, sd, "Lambda", sj, "JC04", "PC04"));
    }
  }
  if(sType != "Lambda_sum"){
    for (auto i=0; i<nm; ++i){
      g[i] = new TGraph(Spectrum(i, sd, sType, sj, "JC04", "PC04"));
    }
  }
  //for (auto i=0; i<nm; ++i) for (auto j=0; j<4; ++j) g[j][i]->SetLineStyle(9);

  auto rf = TFile::Open("./data/Rope/Soft.root", "read");
  auto rl = (TList*)rf->Get(Form("list_results"));
  rf->Close();
  TH1D* rh[2];
  rh[0] = (TH1D*)rl->FindObject(Form("h%s_Jet10_C04", sType.Data()));
  if(sType == "Lambda_sum"){rh[0] = (TH1D*)rl->FindObject(Form("hLambda_Jet10_C04"));
                            auto hA = (TH1D*)rl->FindObject(Form("hAntiLa_Jet10_C04"));
                            rh[0]->Add(hA);}
  if(sType == "Xi"){rh[0] = (TH1D*)rl->FindObject(Form("hXiNeg_Jet10_C04"));
                            auto hA = (TH1D*)rl->FindObject(Form("hXiPos_Jet10_C04"));
                            rh[0]->Add(hA);}
  if(sType == "Omega"){rh[0] = (TH1D*)rl->FindObject(Form("hOmegaNeg_Jet10_C04"));
                            auto hA = (TH1D*)rl->FindObject(Form("hOmegaPos_Jet10_C04"));
                            rh[0]->Add(hA);
			    rh[0]->Rebin(2);}

  rh[0]->Scale(1./(137565.*0.06*0.75*TMath::TwoPi()));
  rh[0]->Rebin(5);
  NormBinningHistogram(rh[0]);

  auto hf = TFile::Open("./data/Rope/Hard.root", "read");
  auto hl = (TList*)hf->Get(Form("list_results"));
  hf->Close();
  rh[1] = (TH1D*)hl->FindObject(Form("h%s_Jet10_C04", sType.Data()));
  if(sType == "Lambda_sum"){rh[1] = (TH1D*)hl->FindObject(Form("hLambda_Jet10_C04"));
                            auto hA = (TH1D*)hl->FindObject(Form("hAntiLa_Jet10_C04"));
                            rh[1]->Add(hA);}
  if(sType == "Xi"){rh[1] = (TH1D*)hl->FindObject(Form("hXiNeg_Jet10_C04"));
                            auto hA = (TH1D*)hl->FindObject(Form("hXiPos_Jet10_C04"));
                            rh[1]->Add(hA);}
  if(sType == "Omega"){rh[1] = (TH1D*)hl->FindObject(Form("hOmegaNeg_Jet10_C04"));
                            auto hA = (TH1D*)hl->FindObject(Form("hOmegaPos_Jet10_C04"));
                            rh[1]->Add(hA);
                            rh[1]->Rebin(2);}
  //rh[1]->Scale(1./(2.*0.75*TMath::TwoPi()*70.));
  rh[1]->Rebin(5);
  
  rh[1]->Scale(1./(137565.*0.06*0.75*TMath::TwoPi()));
  NormBinningHistogram(rh[1]);

  g[nm] = new TGraph(rh[0]);
  g[nm+1] = new TGraph(rh[1]);



//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(1e-4), dfuy(1e1);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("d#it{#rho}/d#it{p}_{T}: K^{0}_{S}");

  if(sType == "Lambda_sum"){ 
    dfly = 1e-4;
    dfuy = 5e0;
    stny = "d#it{#rho}/d#it{p}_{T}: #Lambda + #bar{#Lambda}";
  }
  if(sType == "Xi"){
    dfux = 8.;
    dfly = 1e-4;
    dfuy = 1e-1;
    stny = "d#it{#rho}/d#it{p}_{T}: #Xi^{-} + #bar{#Xi}^{+}";
  }
  if(sType == "Omega"){
    dfux = 5.;
    dfly = 5e-5;
    dfuy = 1e-2;
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
  //DrawGraph(g[2],  wcl[2], "C");
  //DrawGraph(g[3],  wcl[3], "C");
  DrawGraph(g[4],  wcl[4], "C");
  DrawGraph(g[5],  wcl[5], "C");

  auto leg(new TLegend(0.52, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(h, "Data", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "BLC mode 0", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[2], "BLC mode 2", "L")->SetTextSizePixels(24);
  //leg->AddEntry(g[3], "BLC mode 3", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[4], "Rope w/ soft QCD", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[5], "Rope w/ hard QCD", "L")->SetTextSizePixels(24);
  leg->Draw();
  //leg->AddEntry(gE[0], "Sys. Error", "f")->SetTextSizePixels(24);

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  if(sType == "Kshort") tex->DrawLatex(0.16, 0.82, Form("K^{0}_{S} in jet"));
  if(sType == "Lambda_sum") tex->DrawLatex(0.16, 0.82, Form("#Lambda + #bar{#Lambda} in jet"));
  if(sType == "Xi" || sType == "Omega") tex->DrawLatex(0.16, 0.82, Form("#%s^{-} + #bar{#%s}^{-} in jet", sType.Data(), sType.Data()));
  tex->DrawLatex(0.16, 0.72, "PYTHIA 8");
  can->SaveAs(Form("./figure/23_pp13TeV_%s_DataPy_AllModel_Cone3.eps", sType.Data()));
  can->SaveAs(Form("./figure/23_pp13TeV_%s_DataPy_AllModel_Cone3.pdf", sType.Data()));
  can->SaveAs(Form("./figure/23_pp13TeV_%s_DataPy_AllModel_Cone3.png", sType.Data()));
  CanvasEnd(can);
  
  return;  
}
