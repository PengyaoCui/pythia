#include "./SourceFun.h"

void SpectRatiowCent(const TString sType = "Kshort"){

  TString sCone = "Incl";

  TString sLatex(Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
 
  auto f = TFile::Open("./result/FinalSpect_ThisAna.root", "read");
  auto l0 = (TList*)f->Get(Form("%s_0100", sType.Data()));
  auto l1 = (TList*)f->Get(Form("%s_010", sType.Data()));
  auto l2 = (TList*)f->Get(Form("%s_1040", sType.Data()));
  auto l3 = (TList*)f->Get(Form("%s_40100", sType.Data()));
  f->Close();

  auto h0 = (TH1D*)l0->FindObject(Form("%sCen", sCone.Data())); h0->SetName("h0");
  auto h1 = (TH1D*)l1->FindObject(Form("%sCen", sCone.Data())); h1->SetName("h1");
  auto h2 = (TH1D*)l2->FindObject(Form("%sCen", sCone.Data())); h2->SetName("h2");
  auto h3 = (TH1D*)l3->FindObject(Form("%sCen", sCone.Data())); h3->SetName("h3");

  auto e0 = (TH1D*)l0->FindObject(Form("%sErr", sCone.Data())); e0->SetName("e0");
  auto e1 = (TH1D*)l1->FindObject(Form("%sErr", sCone.Data())); e1->SetName("e1");
  auto e2 = (TH1D*)l2->FindObject(Form("%sErr", sCone.Data())); e2->SetName("e2");
  auto e3 = (TH1D*)l3->FindObject(Form("%sErr", sCone.Data())); e3->SetName("e3");

  sCone = "UE";
  auto hUE0= (TH1D*)l0->FindObject(Form("%sCen", sCone.Data())); hUE0->SetName("hUE0");
  auto hUE1= (TH1D*)l1->FindObject(Form("%sCen", sCone.Data())); hUE1->SetName("hUE1");
  auto hUE2= (TH1D*)l2->FindObject(Form("%sCen", sCone.Data())); hUE2->SetName("hUE2");
  auto hUE3= (TH1D*)l3->FindObject(Form("%sCen", sCone.Data())); hUE3->SetName("hUE3");
  
  auto eUE0 = (TH1D*)l0->FindObject(Form("%sErr", sCone.Data())); eUE0->SetName("eUE0");
  auto eUE1 = (TH1D*)l1->FindObject(Form("%sErr", sCone.Data())); eUE1->SetName("eUE1");
  auto eUE2 = (TH1D*)l2->FindObject(Form("%sErr", sCone.Data())); eUE2->SetName("eUE2");
  auto eUE3 = (TH1D*)l3->FindObject(Form("%sErr", sCone.Data())); eUE3->SetName("eUE3");
  
  sCone = "JE";
  auto hJE0= (TH1D*)l0->FindObject(Form("%sCen", sCone.Data())); hJE0->SetName("hJE0");
  auto hJE1= (TH1D*)l1->FindObject(Form("%sCen", sCone.Data())); hJE1->SetName("hJE1");
  auto hJE2= (TH1D*)l2->FindObject(Form("%sCen", sCone.Data())); hJE2->SetName("hJE2");
  auto hJE3= (TH1D*)l3->FindObject(Form("%sCen", sCone.Data())); hJE3->SetName("hJE3");
  
  auto eJE0 = (TH1D*)l0->FindObject(Form("%sErr", sCone.Data())); eJE0->SetName("eJE0");
  auto eJE1 = (TH1D*)l1->FindObject(Form("%sErr", sCone.Data())); eJE1->SetName("eJE1");
  auto eJE2 = (TH1D*)l2->FindObject(Form("%sErr", sCone.Data())); eJE2->SetName("eJE2");
  auto eJE3 = (TH1D*)l3->FindObject(Form("%sErr", sCone.Data())); eJE3->SetName("eJE3");

  //-----------------------------------
  TCanvas *can = nullptr;
  TLatex* tex = nullptr;
  TLegend *leg = nullptr;
  can = MakeCanvas("can");
  //can->SetLogy();
  //can->SetGridx(); can->SetGridy();
  can->SetTickx(); can->SetTicky();
  leg = new TLegend(0.7,0.9,0.9,0.6); SetLegend(leg);
  //leg = new TLegend(0.5,0.9,0.9,0.6); SetLegend(leg);
  //leg -> SetNColumns(3);
  auto XMin = 0.; auto XMax = 0.; auto YMin = 0.001; auto YMax = 3.;
  if(sType == "Kshort"){ XMax = 12.;}
  if(sType == "Lambda_sum"){ XMax = 12.;}
  if(sType == "Xi"){ XMax = 8.;}
  if(sType == "Omega"){ XMax = 5.;}

  TH1D* h = new TH1D("h0", "", 100, 0, 20);

  h->GetXaxis()->SetRangeUser(XMin, XMax);
  h->GetYaxis()->SetRangeUser(YMin, YMax);
  SetFrame(h, "#it{p}_{T}(GeV/c)", "0-100% as ref");
  
  DrawHisto(h, cLine[0], sMark[0], "same"); 
  //-----------------------------------

  h1->Divide(h0); hUE1->Divide(hUE0); hJE1->Divide(hJE0);
  h2->Divide(h0); hUE2->Divide(hUE0); hJE2->Divide(hJE0);
  h3->Divide(h0); hUE3->Divide(hUE0); hJE3->Divide(hJE0);
  
  TH1D* hE1[2] = {e0, e1};       auto hEincl1 = (TH1D*)QuadraticSum(2,hE1); 
  TH1D* hEUE1[2] = {eUE0, eUE1}; auto hEue1 = (TH1D*)QuadraticSum(2,hEUE1);
  TH1D* hEJE1[2] = {eJE0, eJE1}; auto hEje1 = (TH1D*)QuadraticSum(2,hEJE1);
  
  TH1D* hE2[2] = {e0, e2};       auto hEincl2 = (TH1D*)QuadraticSum(2,hE2); 
  TH1D* hEUE2[2] = {eUE0, eUE2}; auto hEue2 = (TH1D*)QuadraticSum(2,hEUE2);
  TH1D* hEJE2[2] = {eJE0, eJE2}; auto hEje2 = (TH1D*)QuadraticSum(2,hEJE2);

  TH1D* hE3[2] = {e0, e3};       auto hEincl3 = (TH1D*)QuadraticSum(2,hE3); 
  TH1D* hEUE3[2] = {eUE0, eUE3}; auto hEue3 = (TH1D*)QuadraticSum(2,hEUE3);
  TH1D* hEJE3[2] = {eJE0, eJE3}; auto hEje3 = (TH1D*)QuadraticSum(2,hEJE3);
  
  //DrawHisto(h1, cLine[0], sMark[0], "same"); DrawHisto(hUE1, cLine[3], sMark[1], "same"); DrawHisto(hJE1, cLine[6], sMark[2], "same"); 
  //DrawHisto(h2, cLine[1], sMark[0], "same"); DrawHisto(hUE2, cLine[4], sMark[1], "same"); DrawHisto(hJE2, cLine[7], sMark[2], "same"); 
  //DrawHisto(h3, cLine[2], sMark[0], "same"); DrawHisto(hUE3, cLine[5], sMark[1], "same"); DrawHisto(hJE3, cLine[8], sMark[2], "same"); 
  
  DrawHisto(hJE1, cLine[0], sMark[0], "same"); 
  DrawHisto(hJE2, cLine[1], sMark[1], "same"); 
  DrawHisto(hJE3, cLine[2], sMark[2], "same"); 

  //auto gIncl1 = (TGraphErrors*)ConvHistogramToGraphErrors(h1, hEincl1, h1->GetNbinsX()); DrawGraph(gIncl1, cLine[0], "E2");
  //auto gIncl2 = (TGraphErrors*)ConvHistogramToGraphErrors(h2, hEincl2, h2->GetNbinsX()); DrawGraph(gIncl2, cLine[1], "E2");
  //auto gIncl3 = (TGraphErrors*)ConvHistogramToGraphErrors(h3, hEincl3, h3->GetNbinsX()); DrawGraph(gIncl3, cLine[2], "E2");
  //auto gUE1 = (TGraphErrors*)ConvHistogramToGraphErrors(hUE1, hEue1, hUE1->GetNbinsX()); DrawGraph(gUE1, cLine[3], "E2");
  //auto gUE2 = (TGraphErrors*)ConvHistogramToGraphErrors(hUE2, hEue2, hUE2->GetNbinsX()); DrawGraph(gUE2, cLine[4], "E2");
  //auto gUE3 = (TGraphErrors*)ConvHistogramToGraphErrors(hUE3, hEue3, hUE3->GetNbinsX()); DrawGraph(gUE3, cLine[5], "E2");
  //auto gJE1 = (TGraphErrors*)ConvHistogramToGraphErrors(hJE1, hEje1, hJE1->GetNbinsX()); DrawGraph(gJE1, cLine[6], "E2");
  //auto gJE2 = (TGraphErrors*)ConvHistogramToGraphErrors(hJE2, hEje2, hJE2->GetNbinsX()); DrawGraph(gJE2, cLine[7], "E2");
  //auto gJE3 = (TGraphErrors*)ConvHistogramToGraphErrors(hJE3, hEje3, hJE3->GetNbinsX()); DrawGraph(gJE3, cLine[8], "E2");
  auto gJE1 = (TGraphErrors*)ConvHistogramToGraphErrors(hJE1, hEje1, hJE1->GetNbinsX()); DrawGraph(gJE1, cLine[0], "E2");
  auto gJE2 = (TGraphErrors*)ConvHistogramToGraphErrors(hJE2, hEje2, hJE2->GetNbinsX()); DrawGraph(gJE2, cLine[1], "E2");
  auto gJE3 = (TGraphErrors*)ConvHistogramToGraphErrors(hJE3, hEje3, hJE3->GetNbinsX()); DrawGraph(gJE3, cLine[2], "E2");
  
  //leg->AddEntry(h1, "Incl", "p"); //leg->AddEntry(hUE1, "UE", "p"); leg->AddEntry(hJE1, "JE", "p");
  //leg->AddEntry(h1, "0-10%");     //leg->AddEntry(hUE1, "0-10%");   leg->AddEntry(hJE1, "0-10%");
  //leg->AddEntry(h2, "10-40%");    //leg->AddEntry(hUE2, "10-40%");  leg->AddEntry(hJE2, "10-40%");
  //leg->AddEntry(h3, "40-100%");   //leg->AddEntry(hUE3, "40-100%"); leg->AddEntry(hJE3, "40-100%");
  leg->AddEntry(hJE1, "0-10%");
  leg->AddEntry(hJE2, "10-40%");
  leg->AddEntry(hJE3, "40-100%");
  leg->AddEntry(gJE1, "Sys. Error", "f");
  leg->Draw();

  TLine* l = new TLine(XMin, 1., XMax, 1.);
  l->SetLineColor(kRed);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  l->Draw("same");
  
  //-----------------------------------
  tex =  new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.15, 0.9, "Spectra centrality dependent");
  tex->DrawLatex(0.15, 0.8, sLatex);
 
  if(sType == "Kshort"){
      tex->DrawLatex(0.15, 0.7, Form("K^{0}_{S} in jets"));
  }
  if(sType == "Lambda_sum"){
      tex->DrawLatex(0.15, 0.7, Form("#Lambda + #bar{#Lambda} in jets"));
  }
  if(sType == "Xi" || sType == "Omega"){
    tex->DrawLatex(0.15, 0.7, Form("#%s in jets", sType.Data()));
  }
  gStyle->SetOptStat("");
  can->SaveAs(Form("./figure/ThisAna/JE_%s_SpectRatio_0100.eps", sType.Data()));
  //DrawAliLogo(0.65, 0.90, 24, kTRUE);
  CanvasEnd(can);
  return;
}
