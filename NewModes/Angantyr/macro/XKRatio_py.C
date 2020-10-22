#include "inc/PyJetUtils.h"

void XKRatio_py(){
//=============================================================================

  //auto f = TFile::Open("./result/ParToKRatio_PYTHIA_pp_HardQCD.root", "read");
  auto f = TFile::Open("./result/ParToKRatio_PYTHIA_pp_SoftQCD.root", "read");
  auto h0 = (TH1D*)f->Get("XKRatio_Incl"); h0->SetName("h0");
  auto h1 = (TH1D*)f->Get("XKRatio_JE");   h1->SetName("h1");
  auto h2 = (TH1D*)f->Get("XKRatio_OC");   h2->SetName("h2");
  auto g0 = new TGraph(h0); g0->SetName("g0");
  auto g1 = new TGraph(h1); g1->SetName("g1");


//=============================================================================
  auto dflx(0.), dfux(12.);
  auto dfly(0.), dfuy(0.2);
  
  auto dlsx(0.05), dlsy(0.05);
  auto dtsx(0.05), dtsy(0.05);
  auto dtox(1.30), dtoy(1.10);
  
  TString stnx("#it{p}_{T} (GeV/#it{c})");
  TString stny("Ratio: (#Xi^{-} + #bar{#Xi}^{+}) / 2K_{S}^{0}");
 
  SetStyle(kTRUE);

  auto can(MakeCanvas(Form("XKRatio_py")));
  auto hfm(can->DrawFrame(dflx, dfly, dfux, dfuy));
  SetupFrame(hfm, stnx, stny, dlsx, dlsy, dtsx, dtsy, dtox, dtoy);
  hfm->GetXaxis()->SetNdivisions(510);
  hfm->GetYaxis()->SetNdivisions(510);

  DrawGraph(g0, wcl[0], "C");
  DrawGraph(g1, wcl[1], "C");

  auto leg(new TLegend(0.72, 0.60, 0.98, 0.92)); SetupLegend(leg);
  leg->AddEntry(g0, "Inclusive", "LP")->SetTextSizePixels(24);
  leg->AddEntry(g1, "in jet", "LP")->SetTextSizePixels(24);
  leg->Draw();

  auto tex(new TLatex());
  tex->SetNDC();
  tex->SetTextSizePixels(24);
  tex->DrawLatex(0.16, 0.92, "pp at #sqrt{#it{s}} = 13 TeV");
  tex->DrawLatex(0.15, 0.83, "SoftQCD + CR1 + Rope");
  tex->DrawLatex(0.16, 0.72, Form("#frac{#Xi^{-} + #bar{#Xi}^{+}}{2K^{0}_{S}}"));
  
  can->SaveAs(Form("./figure/XitoKRatioPy_SoftQCD.eps"));
  can->SaveAs(Form("./figure/XitoKRatioPy_SoftQCD.pdf"));
  can->SaveAs(Form("./figure/XitoKRatioPy_SoftQCD.png"));
  CanvasEnd(can);

  return;
}
