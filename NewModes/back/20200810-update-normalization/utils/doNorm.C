void doNorm(
  const TString sm = "",
  const TString ss = "",
  const TString sj = "",
  const TString sd = "pp13d00TeV")
{
  if (sm.IsNull() || ss.IsNull() || sj.IsNull()) {
    ::Error("doNorm::doNorm", "No inputs!!!");
    return;
  }
//=============================================================================

  const TString sType(Form("%s_%s",ss.Data(),sj.Data()));
  const TString sOutF(Form("../data/%s/AnalysisOutputs_%s.root",sd.Data(), sm.Data()));
  const auto b1st(gSystem->AccessPathName(sOutF.Data()));

  if (!b1st) {
    auto file(TFile::Open(sOutF.Data(),"READ"));
    auto list(static_cast<TList*>(file->Get(Form("list_%s",sType.Data()))));
    if (list) {
      ::Warning("doNorm::doNorm", "List was created!!");
      return;
    }
  }
//=============================================================================

  const TString sp(Form("../data/%s/AnalysisResults_root.zip",sd.Data()));

  if (gSystem->AccessPathName(sp)) {
    ::Error("doNorm::doNorm", "No file: %s", sp.Data());
    exit(-1);
  }
//=============================================================================

  auto file(TFile::Open(Form("%s#AnalysisResults_%s.root",sp.Data(),sm.Data()),"READ"));
  auto list_evnInfo(static_cast<TList*>(file->Get("list_evInfo")));
  auto list_weights(static_cast<TList*>(file->Get("list_weights")));
  auto list_results(static_cast<TList*>(file->Get(Form("list_%s",ss.Data()))));
  file->Close();
//=============================================================================

  const auto dA(1.5*TMath::TwoPi());
  auto hEv(static_cast<TH1D*>(list_evnInfo->FindObject("hEvent")));
  auto hWgt(static_cast<TH1D*>(list_weights->FindObject(Form("h%s",sj.Data()))));
  auto aWgt(hWgt->GetXaxis());

  TH1D *histo(nullptr);
  auto list_outputs(new TList());
  const auto kEv(hEv->GetXaxis()->FindBin("ALL"));
  histo = static_cast<TH1D*>(list_results->FindObject(Form("h%s",ss.Data())));
  histo->Scale(1. /  hEv->GetBinContent(kEv) / dA);
  list_outputs->Add(histo);

  const auto kej(hEv->GetXaxis()->FindBin(sj));
  const auto dej(hEv->GetBinContent(kej));
  histo = static_cast<TH1D*>(list_results->FindObject(Form("h%s_%s",ss.Data(),sj.Data())));
  histo->SetName(Form("h%s_JE",ss.Data()));
  histo->Scale(1. /  dej / dA);
  list_outputs->Add(histo);

  const auto dd(hWgt->GetBinContent(1));
  for (auto k=2; k<=aWgt->GetNbins(); ++k) {
    const TString s(aWgt->GetBinLabel(k));
    const auto    d(hWgt->GetBinContent(k) / dd);
    histo = static_cast<TH1D*>(list_results->FindObject(Form("h%s_%s_%s",ss.Data(),sj.Data(),s.Data())));
    histo->SetName(Form("h%s_%s",ss.Data(),s.Data()));
    histo->Scale(1. /  dej / dA / d);
    list_outputs->Add(histo);
  }

  const auto kNJ(hEv->GetXaxis()->FindBin("NJ"));
  histo = static_cast<TH1D*>(list_results->FindObject(Form("h%s_NJ",ss.Data())));
  histo->Scale(1. /  hEv->GetBinContent(kNJ) / dA);
  list_outputs->Add(histo);
//=============================================================================

  file = TFile::Open(sOutF.Data(), b1st ? "NEW" : "UPDATE");
  list_outputs->Write(Form("list_%s",sType.Data()), TObject::kSingleKey);
  file->Close();
//=============================================================================

  return;
}
