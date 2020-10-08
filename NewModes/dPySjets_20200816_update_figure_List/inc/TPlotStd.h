#if !defined(__CINT__) || defined(__MAKECINT__) || !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TString.h>
#endif

TCanvas *MakeCanvas(const TString s)
{
  auto c(new TCanvas(Form("c%s",s.Data()), s.Data(), 700, 500));
  c->Range(0., 0., 1., 1.);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.15);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);

  return c;
}

//_____________________________________________________________________________
void CanvasEnd(TCanvas* const c)
{
  if (!c) return;

  c->Modified();
  c->cd();
  c->SetSelected(c);

  return;
}
