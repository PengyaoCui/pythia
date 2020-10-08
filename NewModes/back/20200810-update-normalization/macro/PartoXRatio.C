#include "./SourceFun.h"

void PartoXRatio(TString sType = "Omega"){
  TString sCone = "PC";
  TString sPath[] = {"../", "../", "../", "../"};
  TString sFile[] = {"AnalysisResults_Monash.root", "AnalysisResults_BLCmode0.root", "AnalysisResults_BLCmode2.root", "AnalysisResults_BLCmode3.root"};
  TString sList[] = {"list_results", "list_results", "list_results", "list_results"};
  
  TString sL[] = {Form("hXiNeg_Jet10_%s04", sCone.Data())};
  TString sL0[] = {Form("hXiPos_Jet10_%s04", sCone.Data())};
  TString sP[] = {Form("h%sNeg_Jet10_%s04", sType.Data(), sCone.Data())};
  TString sA[] = {Form("h%sPos_Jet10_%s04", sType.Data(), sCone.Data())};
  //TString sL[] = {Form("hXiNeg")};
  //TString sL0[] = {Form("hXiPos")};
  //TString sP[] = {Form("h%sNeg", sType.Data())};
  //TString sA[] = {Form("h%sPos", sType.Data())};
  
  TString sName = "OXRatio";
  Double_t XMin = 0.6; Double_t XMax = 8; Double_t YMin = 0.; Double_t YMax = 0.5;
  Double_t dBin[] = {0.6, 1.0, 1.6, 2.4, 3.2, 5.0, 8, 12., 16};
  auto sLatex = (Form("%s #%s to #Xi ratio in pp at #sqrt{s} = 13 TeV", sCone.Data(), sType.Data()));
  
  if(sType == "Omega"){
    sName = "OXRatio";
    XMax = 5; YMax = 0.6;
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
  TH1D* hL = nullptr;
  TH1D* hL0 = nullptr;
  TH1D* hRatio = nullptr;

  Int_t nBin = sizeof(dBin)/sizeof(Double_t) - 1;
  auto *h = new TH1D("h", "", nBin, dBin);
  h->GetXaxis()->SetRangeUser(XMin, XMax);
  h->GetYaxis()->SetRangeUser(YMin, YMax);
  SetFrame(h, "#it{p}_{T} (GeV/#it{c})", Form("#%s/#Xi", sType.Data()));
  DrawHisto(h, cLine[0], sMark[0], "same");
  
  auto f0 = TFile::Open("./result/FinalSpect_ThisAna.root", "read");
  auto l = (TList*)f0->Get(Form("%s_toXRatio", sType.Data()));
  f0->Close();
  auto h0 = (TH1D*)l->FindObject(Form("hInclR")); h0->SetName("h0");
  auto h1 = (TH1D*)l->FindObject(Form("hJER")); h1->SetName("h1");
  auto h2 = (TH1D*)l->FindObject(Form("hUER")); h2->SetName("h2");
  auto g0 = (TGraphErrors*)l->FindObject(Form("InclRerr")); g0->SetName("g0");
  auto g1 = (TGraphErrors*)l->FindObject(Form("JERerr"));   g1->SetName("g1");
  auto g2 = (TGraphErrors*)l->FindObject(Form("UERerr"));   g2->SetName("g2");
  //DrawHisto(h0, cLine[0], sMark[0], "same"); DrawGraph(g0, cLine[0], "E2"); 
  //DrawHisto(h1, cLine[0], sMark[0], "same"); DrawGraph(g1, cLine[0], "E2");//if(sCone == "JC" || sCone == "JE")
  DrawHisto(h2, cLine[0], sMark[0], "same"); DrawGraph(g2, cLine[0], "E2");//if(sCone == "UE" || sCone == "PC")
  leg->AddEntry(h0, sLeg[0], "LP");

  for (Int_t i = 0; i< nHist; ++i){
    TString sMyPath = sPath[i]; 
    TString sMyFile = sFile[i]; 
    TString sMyList = sList[i]; 
    TString sMyHistP = sP[0]; 
    TString sMyHistA = sA[0]; 
    TString sMyHistL = sL[0]; 
    TString sMyHistL0 = sL0[0]; 

    hL = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistL.Data());
    hL0 = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistL0.Data());
    hL->SetName(Form("XiNeg_%d", i)); 
    hL0->SetName(Form("XiPos_%d", i)); 
    //DeNormHistBinWidth(hL); 
    hL->Add(hL0);
    hL=RebinTH1D(hL, h);
    //NormHistBinWidth(hL);

    hP = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistP.Data());
    hP->SetName(Form("Particle_%d", i)); 
    hP=RebinTH1D(hP, hL);
    
    hA = (TH1D*)GetTH1D(sMyPath.Data(), sMyFile.Data(), sMyList.Data(), sMyHistA.Data());
    hA->SetName(Form("AntiPar_%d", i)); 
    hA=RebinTH1D(hA, hL);
    hP->Add(hA);
    
    hRatio = (TH1D*) hL->Clone(Form("hRatio_%d", i)); hRatio->Reset();
    hRatio->Divide(hP, hL);
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
  can->SaveAs(Form("./figure/%s_toXRatio_wDiffCRMode_%s.eps", sType.Data(), sCone.Data()));
  can->SaveAs(Form("./figure/%s_toXRatio_wDiffCRMode_%s.png", sType.Data(), sCone.Data()));
  return;
}
