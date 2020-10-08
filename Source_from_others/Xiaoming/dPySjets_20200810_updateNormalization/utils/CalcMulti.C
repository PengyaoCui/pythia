const TString sPref("PyMulti");
const TString sPath(Form("../data/%s",sPref.Data()));

const TString sTune[] = { "Monash", "BLCmode0", "BLCmode2", "BLCmode3"};
const auto nTune(sizeof(sTune) / sizeof(TString));
//=============================================================================

void MeanMulti(const TString &, const TString &, Double_t &, Double_t &);
//=============================================================================

void CalcMulti(const TString se="pp13d00TeV")
{
  const TString sp(Form("%s/%s_%s",sPath.Data(),sPref.Data(),se.Data()));
  if (gSystem->AccessPathName(sp)) {
    ::Error("CalcMulti::MeanMulti", "No path: %s", sp.Data());
    exit(-1);
  }
//=============================================================================

  const TString so(Form("%s/AnalysisOutputs.root",sp.Data()));
  if (!gSystem->AccessPathName(so)) {
    ::Warning("CalcMulti::MeanMulti", "Has output file: %s", so.Data());
    return;
  }
//=============================================================================

  auto list(new TList());
  auto hMulti(new TH1D("hMulti", "", nTune, -0.5, -0.5+nTune));
  list->Add(hMulti);

  auto hMult1(new TH1D("hMult1", "", nTune, -0.5, -0.5+nTune));
  list->Add(hMult1);
//=============================================================================

  for (auto i=0, k=1; i<nTune; ++i, ++k) {
    const auto &st(sTune[i]);
    auto dMulti(0.), dMult1(0.);
    MeanMulti(st, sp, dMulti, dMult1);

    hMulti->GetXaxis()->SetBinLabel(k,st);
    hMulti->SetBinContent(k, dMulti);

    hMult1->GetXaxis()->SetBinLabel(k,st);
    hMult1->SetBinContent(k, dMult1);
  }

  list->Print();
//=============================================================================

  auto file(TFile::Open(so, "NEW"));
  list->Write("list_results", TObject::kSingleKey);
  file->Close();
//=============================================================================

  return;
}

//_____________________________________________________________________________
void MeanMulti(const TString &st, const TString &sp, Double_t &dMulti, Double_t &dMult1)
{
  if (st.IsNull() || sp.IsNull()) {
    ::Warning("CalcMulti::MeanMulti", "No input");
    return;
  }
//=============================================================================

  const TString sz(Form("%s/AnalysisResults_root.zip",sp.Data()));
  if (gSystem->AccessPathName(sz)) {
    ::Error("CalcMulti::MeanMulti", "No file: %s", sz.Data());
    exit(-1);
  }
//=============================================================================

  auto file(TFile::Open(Form("%s#AnalysisResults_%s.root",sz.Data(),st.Data()), "READ"));
  auto list(static_cast<TList*>(file->Get("list_results")));
  file->Close();
//=============================================================================

  auto hist(static_cast<TH1D*>(list->FindObject("hMulti")));
//=============================================================================

  auto dM0(0.), dW0(0.);
  auto dM1(0.), dW1(0.);
  for (auto k=1; k<=hist->GetNbinsX(); ++k) {
    const auto dm(hist->GetBinCenter(k));
    const auto dw(hist->GetBinContent(k));
    const auto ds(dm * dw);

    dW0 += dw;
    dM0 += ds;

    if (dm>0.5) {
      dW1 += dw;
      dM1 += ds;
    }
  }

  dM0 /= dW0;
  dM1 /= dW1;
//=============================================================================

  dMulti = dM0;
  dMult1 = dM1;
//=============================================================================

  return;
}
