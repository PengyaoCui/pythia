#if !defined(__CINT__) || defined(__MAKECINT__) || !defined(__CLING__) || defined(__ROOTCLING__)
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#endif

//_____________________________________________________________________________
TH1D *MakeRebinTH1D(TH1D const *hRaw, TH1D const *hRef)
{
  if ((!hRaw) || (!hRef)) return nullptr;

  const auto nRef(hRef->GetNbinsX());

  const auto dLower(hRef->GetXaxis()->GetBinLowEdge(1));
  const auto dUpper(hRef->GetXaxis()->GetBinUpEdge(nRef));

  auto hNew(static_cast<TH1D*>(hRef->Clone(Form("h_%s_%s", hRaw->GetName(),
                                                           hRef->GetName()))));
  hNew->Reset();
//=============================================================================

  for (auto k=1; k<=hRaw->GetNbinsX(); ++k) {
    const auto dXvar(hRaw->GetBinCenter(k));
    if ((dXvar<dLower) || (dXvar>=dUpper)) continue;

    const auto iBinX(hNew->FindBin(dXvar));
    const auto dYvar(hRaw->GetBinContent(k));
    const auto dYerr(hRaw->GetBinError(k));
    const auto dYsw2(hNew->GetBinError(iBinX));

    hNew->SetBinContent(iBinX, hNew->GetBinContent(iBinX) + dYvar);
    hNew->SetBinError(iBinX, TMath::Sqrt(dYsw2*dYsw2 + dYerr*dYerr));
  }
//=============================================================================

  return hNew;
}

//_____________________________________________________________________________
void NormBinningHistogram(TH1D* const h, const Double_t dw=1.)
{
  if (!h) return;
//=============================================================================

  for (auto k=1; k<=h->GetNbinsX(); ++k) {
    const auto dBW(h->GetBinWidth(k) * dw);
    h->SetBinContent(k, h->GetBinContent(k)/dBW);
    h->SetBinError(k, h->GetBinError(k)/dBW);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void DeNormBinningHistogram(TH1D* const h)
{
  if (!h) return;
//=============================================================================

  for (auto k=1; k<=h->GetNbinsX(); ++k) {
    const auto dBW(h->GetBinWidth(k));
    h->SetBinContent(k, h->GetBinContent(k) * dBW);
    h->SetBinError(k, h->GetBinError(k) * dBW);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void CalcRelDeviation(TH1D* const hi, TH1D* const hr)
{
  if ((!hi) || (!hr)) return;
//=============================================================================

  for (auto k=1; k<=hr->GetNbinsX(); ++k) {
    auto dVarI = hi->GetBinContent(k); if (TMath::Abs(dVarI)<1e-12) dVarI = 1e-12;
    auto dErrI = hi->GetBinError(k) / dVarI;

    auto dVarR = hr->GetBinContent(k); if (TMath::Abs(dVarR)<1e-12) dVarR = 1e-12;
    auto dErrR = hr->GetBinError(k) / dVarR;

    const auto dVarD = TMath::Abs(dVarI/dVarR - 1.);
    const auto dErrD = dVarD * TMath::Sqrt(dErrI*dErrI + dErrR*dErrR);
    hi->SetBinContent(k, dVarD); hi->SetBinError(k, dErrD);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
TH1D *MaxHistograms(const Int_t n, TH1D **h)
{
  if ((n<=0) || (!h[0])) return nullptr;
//=============================================================================

  auto d(new Double_t[n]);
  auto hMax(static_cast<TH1D*>(h[0]->Clone("hMax")));

  for (auto k=1; k<=hMax->GetNbinsX(); ++k) {
    for (auto j=0; j<n; ++j) d[j] = h[j]->GetBinContent(k);
    hMax->SetBinContent(k, TMath::MaxElement(n,d));
    hMax->SetBinError(k, TMath::RMS(n,d));
  }
//=============================================================================

  return hMax;
}

//_____________________________________________________________________________
TH1D *QuadraticSum(const Int_t n, TH1D **h)
{

  if ((n<=0) || (!h)) return nullptr;
//=============================================================================

  auto hSum(static_cast<TH1D*>(h[0]->Clone("hSum")));

  for (auto k=1; k<=hSum->GetNbinsX(); ++k) {
    auto dSum(0.);
    for (auto j=0; j<n; ++j) {
      const auto d(h[j]->GetBinContent(k));
      dSum += (d*d);
    }

    hSum->SetBinContent(k, TMath::Sqrt(dSum));
  }
//=============================================================================

  return hSum;
}

//_____________________________________________________________________________
TH1D *ConvRelError(TH1D *hi, TH1D *hr)
{
  auto hc(static_cast<TH1D*>(hi->Clone(Form("h_%s_%s",hi->GetName(),hr->GetName()))));
  for (auto k=1; k<=hc->GetNbinsX(); ++k) hc->SetBinError(k, hi->GetBinContent(k) * hr->GetBinContent(k));
//=============================================================================

  return hc;
}

//_____________________________________________________________________________
TGraphErrors *ConvHistogramToGraphErrors(TH1D* const hVar, TH1D* const hErr, const Int_t n)
{
//auto dtm(new Double_t[n]);

  Double_t dvx[n], dvy[n];
  Double_t dex[n], dey[n];
  for (auto i=0, k=1; i<n; ++i, ++k) {
    dvx[i] = hVar->GetBinCenter(k);
    dex[i] = hVar->GetBinWidth(k) / 2.;

    if (dex[i]>0.2) dex[i] = 0.1;

    dvy[i] = hVar->GetBinContent(k);
    if (hErr) {
      dey[i] = hErr->GetBinContent(k) * dvy[i];
    } else {
      dey[i] = hVar->GetBinError(k);
    }
  }
//=============================================================================

  return (new TGraphErrors(n,dvx,dvy,dex,dey));
}

//_____________________________________________________________________________
TGraphErrors *ConvHistogramToGraphErrorL(TH1D* const hVar, TH1D* const hErr, const Int_t n)
{
//auto dtm(new Double_t[n]);

  Double_t dvx[n], dvy[n];
  Double_t dex[n], dey[n];
  for (auto i=0, k=1; i<n; ++i, ++k) {
    dvx[i] = hVar->GetBinCenter(k);
    dex[i] = hVar->GetBinWidth(k) / 2.;

    if (dex[i]>0.3) dvx[i] -= 0.2;
    if (dex[i]>0.2) dex[i]  = 0.1;

    dvy[i] = hVar->GetBinContent(k);

    if (hErr) {
      dey[i] = hErr->GetBinContent(k) * dvy[i];
    } else {
      dey[i] = hVar->GetBinError(k);
    }
  }
//=============================================================================

  return (new TGraphErrors(n,dvx,dvy,dex,dey));
}

//_____________________________________________________________________________
TGraphErrors *ConvHistogramToGraphErrorR(TH1D* const hVar, TH1D* const hErr, const Int_t n)
{
//auto dtm(new Double_t[n]);

  Double_t dvx[n], dvy[n];
  Double_t dex[n], dey[n];
  for (auto i=0, k=1; i<n; ++i, ++k) {
    dvx[i] = hVar->GetBinCenter(k);
    dex[i] = hVar->GetBinWidth(k) / 2.;

    if (dex[i]>0.3) dvx[i] += 0.2;
    if (dex[i]>0.2) dex[i]  = 0.1;

    dvy[i] = hVar->GetBinContent(k);

    if (hErr) {
      dey[i] = hErr->GetBinContent(k) * dvy[i];
    } else {
      dey[i] = hVar->GetBinError(k);
    }
  }
//=============================================================================

  return (new TGraphErrors(n,dvx,dvy,dex,dey));
}

//_____________________________________________________________________________
TGraphAsymmErrors *ConvHistogramToGraphErrors(TH1D* const hVar, TH1D* const hErl, TH1D* const hErh)
{
//const auto m(n);
//auto dtm(new Double_t[m]);

  const auto n(hVar->GetNbinsX());
  Double_t dvx[n], dlx[n], dhx[n];
  Double_t dvy[n], dly[n], dhy[n];
  for (auto i=0, k=1; i<n; ++i, ++k) {
    dvx[i] = hVar->GetBinCenter(k);
    dlx[i] = 0.05; // hVar->GetBinWidth(k) / 2.;
    dhx[i] = 0.05; // hVar->GetBinWidth(k) / 2.;

    dvy[i] = hVar->GetBinContent(k);
    dly[i] = hErl->GetBinContent(k) * dvy[i];
    dhy[i] = hErh->GetBinContent(k) * dvy[i];
  }
//=============================================================================

  return (new TGraphAsymmErrors(n,dvx,dvy,dlx,dhx,dly,dhy));
}

//_____________________________________________________________________________
TGraphErrors *ConvHistogramToGraphErrors(TH1D* const hVar, const Double_t dd=0.05)
{
  const auto n(hVar->GetNbinsX());

  Double_t dvx[n], dex[n];
  Double_t dvy[n], dey[n];
  for (auto i=0, k=1; i<n; ++i, ++k) {
    dvx[i] = hVar->GetBinCenter(k);
    dex[i] = dd;

    dvy[i] = hVar->GetBinContent(k);
    dey[i] = hVar->GetBinError(k);
  }
//=============================================================================

  return (new TGraphErrors(n,dvx,dvy,dex,dey));
}
