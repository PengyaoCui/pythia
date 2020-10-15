#include "./SourceFun.h"

void draw_DataPy(const TString sType = "Xi"){


  auto sLatex("pp #sqrt{#it{s}} = 13 TeV");
  //if(sType != "Lambda_sum") sLatex = Form("#%s to K^{0}_{S} ratio", sType.Data());

  auto f = TFile::Open("./result/FinalSpect_ThisAna.root", "read");
  auto l = (TList*)f->Get(Form("%s_toKRatio", sType.Data()));
  f->Close();
  auto h0 = (TH1D*)l->FindObject(Form("hInclR")); h0->SetName("h0");
  auto h1 = (TH1D*)l->FindObject(Form("hJER")); h1->SetName("h1");
  auto h2 = (TH1D*)l->FindObject(Form("hUER")); h2->SetName("h2");
  auto g0 = (TGraphErrors*)l->FindObject(Form("InclRerr")); g0->SetName("g0");
  auto g1 = (TGraphErrors*)l->FindObject(Form("JERerr"));   g1->SetName("g1");
  auto g2 = (TGraphErrors*)l->FindObject(Form("UERerr"));   g2->SetName("g2");
  
  //auto f0 = TFile::Open("./result/ParToKRatio_PYTHIA_pp_SoftQCD.root", "read");
  auto f0 = TFile::Open("./result/ParToKRatio_PYTHIA_pp_HardQCD.root", "read");
  auto h00 = (TH1D*)f0->Get("LKRatio_Incl"); h00->SetName("h00");
  auto h01 = (TH1D*)f0->Get("LKRatio_JE");   h01->SetName("h01");
  auto h02 = (TH1D*)f0->Get("LKRatio_OC");   h02->SetName("h02");
  if(sType == "Xi"){
    h00 = (TH1D*)f0->Get(Form("XKRatio_Incl")); h00->SetName("h00");
    h01 = (TH1D*)f0->Get(Form("XKRatio_JE"));   h01->SetName("h01");
    h02 = (TH1D*)f0->Get(Form("XKRatio_OC"));   h02->SetName("h02");
  }

  if(sType == "Omega"){
    h00 = (TH1D*)f0->Get(Form("OKRatio_Incl")); h00->SetName("h00");
    h01 = (TH1D*)f0->Get(Form("OKRatio_JE"));   h01->SetName("h01");
    h02 = (TH1D*)f0->Get(Form("OKRatio_OC"));   h02->SetName("h02");
  }

 
  //=================================================================
  auto g00 = new TGraph(h00); g00->SetLineStyle(9); g00->SetName("g00");
  auto g01 = new TGraph(h01); g01->SetLineStyle(9); g01->SetName("g01");
  auto g02 = new TGraph(h02); g02->SetLineStyle(9); g02->SetName("g02");

  TCanvas *can = nullptr;
  can = MakeCanvas("can");
  can->SetTickx(); can->SetTicky();
  TLegend *leg = nullptr;
  leg = new TLegend(0.7,0.9,1.0,0.6); SetLegend(leg);

  Double_t YMax = 1.2; if(sType == "Xi") YMax = 0.2; if(sType == "Omega") YMax = 0.05;

  h0->GetYaxis()->SetRangeUser(0., YMax);
  if(sType == "Lambda_sum")h0->GetYaxis()->SetTitle("#Lambda/K^{0}_{S}");
  if(sType != "Lambda_sum")h0->GetYaxis()->SetTitle(Form("#%s/K^{0}_{S}", sType.Data()));
  
  DrawHisto(h0, cLine[0], sMark[0], "same"); DrawGraph(g0, cLine[0], "E2"); DrawGraph(g00, cLine[0], "C"); leg->AddEntry(h0, "Inclusive");
  DrawHisto(h1, cLine[1], sMark[1], "same"); DrawGraph(g1, cLine[1], "E2"); leg->AddEntry(h1, "JE"); DrawGraph(g01, cLine[1], "C"); 
  DrawHisto(h2, cLine[2], sMark[2], "same"); DrawGraph(g2, cLine[2], "E2"); leg->AddEntry(h2, "UE"); //DrawGraph(g02, cLine[2], "C"); 
  leg->AddEntry(g00, "PYTHIA8", "LP")->SetTextSizePixels(24);
  //leg->AddEntry(g0, "Sys.Error", "f"); 
  leg->Draw();

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.15, 0.9, sLatex);
  tex->DrawLatex(0.15, 0.83, "HardQCD + CR1 + Rope");
  if(sType == "Lambda_sum")tex->DrawLatex(0.15, 0.75, Form("#frac{#Lambda + #bar{#Lambda}}{2K^{0}_{S}}"));
  if(sType == "Xi" || sType == "Omega")tex->DrawLatex(0.15, 0.75, Form("#frac{#%s + #bar{#%s}}{2K^{0}_{S}}", sType.Data(), sType.Data()));

  gStyle->SetOptStat("");
 
  can->SaveAs(Form("./figure/%s_toKRatio_DataPy_Rope_hQCD.eps", sType.Data()));
  can->SaveAs(Form("./figure/%s_toKRatio_DataPy_Rope_hQCD.png", sType.Data()));
  CanvasEnd(can);
  
  f0->Close();
  return;
}
