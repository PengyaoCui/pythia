#include "inc/PyJetUtils.h"
//=============================================================================

const auto m(0);
const TString sj("Jet10");
const TString sJC("JC04"), sUE("OC08");

const TString sd("pp05d02TeVrs");  // pp13d00TeV:   pp at 13   TeV
                                   // pp05d02TeVrs: pp at 5.02 TeV
                                   //               with rapidity shift
                                   //               used to fit p-Pb acceptance
//=============================================================================

void PlotSrcs()
{
  auto gIN(new TGraph(RatioLK(m,sd)));
  auto gJC(new TGraph(RatioLK(m,sd,sj,sJC)));
  auto gUE(new TGraph(RatioLK(m,sd,sj,sUE)));
  auto gHD(new TGraph(RatioLK(m,sd,sj,sJC,sUE)));
//=============================================================================

  const auto dflx(0.), dfux(12.);
  const auto dfly(0.), dfuy(1.0);

  const auto dlsx(0.05), dlsy(0.05);
  const auto dtsx(0.05), dtsy(0.05);
  const auto dtox(1.30), dtoy(1.10);

  const TString stnx("#it{p}_{T} (GeV/#it{c})");
  const TString stny("Ratio: (#Lambda + #bar{#Lambda}) / 2K_{S}^{0}");
//=============================================================================

  SetStyle(kTRUE);
  auto can(MakeCanvas("PySrcs"));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawGraph(gIN,  wcl[0], "C");
  DrawGraph(gJC,  wcl[1], "C");
  DrawGraph(gUE,  wcl[2], "C");
  DrawGraph(gHD,  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(gIN, "Inclusive", "L")->SetTextSizePixels(24);
  leg->AddEntry(gJC, "In jets", "L")->SetTextSizePixels(24);
  leg->AddEntry(gUE, "UE", "L")->SetTextSizePixels(24);
  leg->AddEntry(gHD, "UE sub", "L")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "PYTHIA 8");
  CanvasEnd(can);
//=============================================================================

  return;
}
