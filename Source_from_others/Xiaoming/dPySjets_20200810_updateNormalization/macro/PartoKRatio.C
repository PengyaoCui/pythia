#include "./SourceFun.h"

void PartoKRatio(TString sType = "Lambda_sum"){

  TString sCone = "JC";
  TString sPath[] = {"../", "../", "../", "../"};
  TString sFile[] = {"AnalysisResults_Monash.root", "AnalysisResults_BLCmode0.root", "AnalysisResults_BLCmode2.root", "AnalysisResults_BLCmode3.root"};
  TString sList[] = {"list_results", "list_results", "list_results", "list_results"};
  
  TString sK[] = {"hKshort_Jet10_JC04"};
  TString sP[] = {"hLambda_Jet10_JC04"};
  TString sA[] = {"hAntiLa_Jet10_JC04"};
  
  TString sName = "LKRatio";
  Double_t XMin = 0.; Double_t XMax = 12; Double_t YMin = 0.; Double_t YMax = 1.;
  Double_t dBin[] = {0.6, 1.0, 1.6, 2.4, 3.2, 5.0, 8, 12., 16};
  auto sLatex(Form("%s #Lambda to K^{0}_{S} ratio in pp at #sqrt{s} = 13 TeV", sCone.Data()));

  if(sType == "Xi"){
    sP[0] = "hXiNeg_Jet10_JC04";
    sA[0] = "hXiPos_Jet10_JC04";
    sName = "XKRatio";
    XMax = 8; YMax = 0.2;
    sLatex = (Form("%s #Xi to K^{0}_{S} ratio in pp at #sqrt{s} = 13 TeV", sCone.Data()));
  } 
  
  if(sType == "Omega"){
    sP[0] = "hOmegaNeg_Jet10_JC04";
    sA[0] = "hOmegaPos_Jet10_JC04";
    sName = "OKRatio";
    XMax = 5; YMax = 0.04;
    sLatex = (Form("%s #Omega to K^{0}_{S} ratio in pp at #sqrt{s} = 13 TeV", sCone.Data()));
  } 


  TString sLeg[] = {"Data", "Monash", "Mode 0", "Mode 2", "Mode 3"};

  //-----------------------------------
 
  TCanvas *can = nullptr;
  TLatex* tex = nullptr;
  TLegend *leg = nullptr;

  //TFile *f = TFile::Open("./result/ParToKRatio_PYTHIA_pp.root", "UPDATE");
  can = MakeCanvas("can");

  leg = new TLegend(0.7,0.9,1.0,0.6); SetLegend(leg);
  //-----------------------------------
  auto nHist = sizeof(sFile)/sizeof(TString);
  
  TH1D* hP = nullptr;
  TH1D* hA = nullptr;
  TH1D* hK = nullptr;
  TH1D* hRatio = nullptr;

  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  auto *h = new TH1D("h", "", nBin, dBin);
  h->GetXaxis()->SetRangeUser(XMin, XMax);
  h->GetYaxis()->SetRangeUser(YMin, YMax);
  SetFrame(h, "#it{p}_{T} (GeV/#it{c})", "#Lambda/K^{0}_{S}");
  if(sType == "Xi")SetFrame(h, "#it{p}_{T} (GeV/#it{c})", "#Xi/K^{0}_{S}");
  if(sType == "Omega")SetFrame(h, "#it{p}_{T} (GeV/#it{c})", "#Omega/K^{0}_{S}");
  DrawHisto(h, cLine[0], sMark[0], "same");
  
  auto f0 = TFile::Open("./result/FinalSpect_ThisAna.root", "read");
  auto l = (TList*)f0->Get(Form("%s_toKRatio", sType.Data()));
  f0->Close();
  auto h0 = (TH1D*)l->FindObject(Form("hInclR")); h0->SetName("h0");
  auto h1 = (TH1D*)l->FindObject(Form("hJER")); h1->SetName("h1");
  auto h2 = (TH1D*)l->FindObject(Form("hUER")); h2->SetName("h2");
  auto g0 = (TGraphErrors*)l->FindObject(Form("InclRerr")); g0->SetName("g0");
  auto g1 = (TGraphErrors*)l->FindObject(Form("JERerr"));   g1->SetName("g1");
  auto g2 = (TGraphErrors*)l->FindObject(Form("UERerr"));   g2->SetName("g2");
  //DrawHisto(h0, cLine[0], sMark[0], "same"); DrawGraph(g0, cLine[0], "E2"); 
  DrawHisto(h1, cLine[0], sMark[0], "same"); DrawGraph(g1, cLine[0], "E2");//if(sCone == "JC" || sCone == "JE")
  //DrawHisto(h2, cLine[0], sMark[0], "same"); DrawGraph(g2, cLine[0], "E2");//if(sCone == "UE" || sCone == "PC")
  leg->AddEntry(h0, sLeg[0], "LP");

  for (Int_t i = 0; i< nHist; ++i){
    TString sMyPath = sPath[i]; 
    TString sMyFile = sFile[i]; 
    TString sMyList = sList[i]; 
    TString sMyHistP = sP[0]; 
    TString sMyHistA = sA[0]; 
    TString sMyHistK = sK[0]; 

    hK = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistK.Data());
    hK->SetName(Form("Kshort_%d", i)); 
    //DeNormHistBinWidth(hK); 
    hK=RebinTH1D(hK, h);
    //NormHistBinWidth(hK);
    hK->Scale(2.); 

    hP = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistP.Data());
    hP->SetName(Form("Lambda_%d", i)); 
    hP=RebinTH1D(hP, hK);
    
    hA = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistA.Data());
    hA->SetName(Form("AntiLa_%d", i)); 
    hA=RebinTH1D(hA, hK);
    hP->Add(hA);
    
    hRatio = (TH1D*) hK->Clone(Form("hRatio_%d", i)); hRatio->Reset();
    hRatio->Divide(hP, hK);
    hRatio->SetName(Form("hratio_%d", i));
    auto g00 = new TGraph(hRatio); g00->SetLineStyle(9); g00->SetName(Form("g%d", i));
    DrawGraph(g00, cLine[i+1], "C");
    leg->AddEntry(g00, sLeg[i+1], "L");
  }

  //f->Close();
  tex =  new TLatex();
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.15, 0.92, sLatex);
  leg->Draw();
  gStyle->SetOptStat("");
  CanvasEnd(can);
  can->SaveAs(Form("./figure/%s_toKRatio_wDiffCRMode_%s.eps", sType.Data(), sCone.Data()));
  can->SaveAs(Form("./figure/%s_toKRatio_wDiffCRMode_%s.png", sType.Data(), sCone.Data()));
  return;
}
