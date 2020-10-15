#include "inc/PyJetUtils.h"
//=============================================================================

const TString sj("Jet10");
const TString sd("pp05d02TeVrs");  // pp13d00TeV:   pp at 13   TeV
                                   // pp05d02TeVrs: pp at 5.02 TeV
                                   //               with rapidity shift
                                   //               used to fit p-Pb acceptance
//=============================================================================

void PlotMode()
{
  TGraph *g[nm];
  for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioLK(i,sd));
//for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioLK(i, sd, sj, "JC04"));
//for (auto i=0; i<nm; ++i) g[i] = new TGraph(RatioLK(i, sd, sj, "JC04", "PC04"));
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
  auto can(MakeCanvas("PyMode"));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawGraph(g[1],  wcl[0], "C");
  DrawGraph(g[0],  wcl[1], "C");
  DrawGraph(g[2],  wcl[2], "C");
  DrawGraph(g[3],  wcl[3], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(g[1], "Monash", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[0], "mode 0", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[2], "mode 2", "L")->SetTextSizePixels(24);
  leg->AddEntry(g[3], "mode 3", "L")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "PYTHIA 8");
  CanvasEnd(can);
//=============================================================================

  return;
}
