#include "TPlotStd.h"
#include "TAliFigs.h"
//=============================================================================

void Plotd()
{
  const TString  sPtHat[] = { "PtHat_005_011", "PtHat_011_021", "PtHat_021_036", "PtHat_036_057", "PtHat_057_INF" };
  const Double_t dPtBin[] = { 0.6, 1.6, 2.2, 2.8, 3.7, 5.0, 8.0, 12., 15. };
  const auto nPtBin(sizeof(dPtBin) / sizeof(Double_t) - 1);
//=============================================================================

  TH1D *hKshortF(nullptr);
  TH1D *hLambdaF(nullptr);
  TH1D *hAntiLaF(nullptr);
  for (const auto s : sPtHat) {
    auto file(TFile::Open(Form("data/AnalysisResults_%s.root",s.Data()),"READ"));
    auto lNor(static_cast<TList*>(file->Get("list_pyxsect")));
    auto lRes(static_cast<TList*>(file->Get("list_results")));
    file->Close();

    auto hWgtEv(static_cast<TH1D*>(lNor->FindObject("hTrials")));
    hWgtEv->SetName(Form("hWgtEv_%s",s.Data()));

    auto hXsect(static_cast<TProfile*>(lNor->FindObject("hXsect")));
    hXsect->SetName(Form("hXsect_%s",s.Data()));

    const auto dWgt(hXsect->GetBinContent(1) / hWgtEv->GetBinContent(1));
//=============================================================================

    auto hK(static_cast<TH1D*>(lRes->FindObject("hKshort_Jet10_C04")));
    hK->SetName(Form("hK_%s",s.Data()));
    hK->Scale(dWgt);
    if (hKshortF) {
      hKshortF->Add(hK);
    } else {
      hKshortF = static_cast<TH1D*>(hK->Clone("hKshortF"));
    }

    auto hL(static_cast<TH1D*>(lRes->FindObject("hLambda_Jet10_C04")));
    hL->SetName(Form("hL_%s",s.Data()));
    hL->Scale(dWgt);
    if (hLambdaF) {
      hLambdaF->Add(hL);
    } else {
      hLambdaF = static_cast<TH1D*>(hL->Clone("hLambdaF"));
    }

    auto hA(static_cast<TH1D*>(lRes->FindObject("hAntiLa_Jet10_C04")));
    hA->SetName(Form("hA_%s",s.Data()));
    hA->Scale(dWgt);
    if (hAntiLaF) {
      hAntiLaF->Add(hA);
    } else {
      hAntiLaF = static_cast<TH1D*>(hA->Clone("hAntiLaF"));
    }
  }
//=============================================================================

  auto hKshortR(hKshortF->Rebin(nPtBin,"hKshortR",dPtBin));
  auto hLambdaR(hLambdaF->Rebin(nPtBin,"hLambdaR",dPtBin));
  auto hAntiLaR(hAntiLaF->Rebin(nPtBin,"hAntiLaR",dPtBin));
//=============================================================================

  auto h(static_cast<TH1D*>(hLambdaR->Clone("h")));
  h->Add(hAntiLaR);
  h->Divide(hKshortR);
  h->Scale(0.25);
//=============================================================================


  SetStyle();
  auto can(MakeCanvas("Incl"));
  auto frm(can->DrawFrame(0.,0.,12.,1.));
  SetupFrame(frm, "#it{p}_{T}", "Ratio");

  auto g(new TGraph(h));
  DrawGraph(g, wcl[0], "L");

  CanvasEnd(can);
//=============================================================================

  return;
}
