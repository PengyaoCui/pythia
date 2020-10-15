#include "./SourceFun.h"

void Draw(){
  auto Path = "";
  auto File = "";
  TString sPath[] = {"./Data", "./Data"};
  TString sFile[] = {"ParToKRatio_DATA_pPb.root", "ParToKRatio_PYTHIA_pPb.root"};
  TString sList[] = {"", ""};
  TString sHist[] = {"OKRatio", "OKRatio"};

  TString sLeg[] = {"pp","p-Pb"}; 

  auto sLatex("#Omega to K^{0}_{S} ratio");
  //-----------------------------------
  TCanvas *can = nullptr;

  can = MakeCanvas("can");
  //can->SetLogy();
  //can->SetGridx(); can->SetGridy();
  
  Double_t dBin[] = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.3, 3.6, 3.9, 4.2, 4.6, 5, 5.4, 5.9, 6.5, 7, 7.5, 8, 8.5, 9.2, 10, 11, 12, 13.5, 15 };
  
  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  //-----------------------------------

  auto leg = new TLegend(0.8,0.9,1.0,0.65);    
  //-----------------------------------
  auto nHist = sizeof(sHist)/sizeof(TString);
  
  for (Int_t i = 0; i< nHist; i++){
    //TString sMyPath = Form("%s_%s", Path.Data(), sPath[i].Data()); 
    TString sMyPath = sPath[i]; 
    TString sMyFile = sFile[i]; 
    TString sMyList = sList[i]; 
    TString sMyHist = sHist[i]; 
    
    TH1D* h1 = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHist.Data());
    
    //DeNormHistBinWidth(h1); 
    //h1 = (TH1D*) h1-> Rebin(nBin, "h1", dBin);
    //NormHistBinWidth(h1); 
    h1->SetName(Form("hist_%d", i));
    
    //for(int j =1; j<=h1->GetNbinsX(); j++){cout<<h1->GetBinLowEdge(j)<<endl;}
    auto dMini = h1->GetMinimum(); 
    auto dMaxi = h1->GetMaximum(); 
    //h1->GetXaxis()->SetRangeUser(0.6, 12.);
    //h1->GetYaxis()->SetRangeUser(0.1*dMini, 3.*dMaxi);
    h1->GetYaxis()->SetRangeUser(0.001*dMini, 0.02);
    h1->SetTitle(""); 
    SetFrame(h1, "#it{p}_{T}(GeV/#it{c})", "#Omega/K^{0}_{S}");
    DrawHisto(h1, cLine[i], sMark[i], "same");
    leg->AddEntry(h1,sLeg[i],"lp");
  }

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  SetLegend(leg);
  tex->DrawLatex(0.16, 0.91, sLatex);
  tex->DrawLatex(0.16, 0.8, "|#eta| < 0.75");
  tex->DrawLatex(0.16, 0.7, "40-100%");
  //tex->DrawLatex(0.16, 0.7, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  leg->Draw();
  gStyle->SetOptStat("");

  //can->SaveAs(Form("./figure/%s_pp_pPb.eps", sHist[0].Data()));
  CanvasEnd(can);
  return;
}
