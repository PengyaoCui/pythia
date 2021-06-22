const TString gksp("hard");
const TString gksf("AnalysisResults");

const Int_t gklb[] = { 5, 11, 21, 36, 57, -1 };
const auto  gknb(sizeof(gklb) / sizeof(Int_t) - 1);
//=============================================================================

TString  GetFileName(const Int_t);
Double_t GetNormFactor(TList *, const Int_t);
Double_t GetJetEventWeight();
void     Normalization(TList *, const Int_t);
//=============================================================================

void MergePyHard(const Bool_t b=kFALSE)
{
  auto list(new TList());
  for (auto i=0; i<gknb; ++i) Normalization(list, i);
  list->Print();
//=============================================================================

  if (b) {
    auto file(TFile::Open(Form("%s_hard.root",gksf.Data()),"NEW"));
    list->Write("list_results", TObject::kSingleKey);
    file->Close();
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void Normalization(TList *lr, const Int_t l)
{
  const auto sf(GetFileName(l));
  ::Info("MergePyHard.C::Normalization", "File name = %s", sf.Data());
//=============================================================================

  auto f(TFile::Open(sf.Data(),"READ"));
  auto ln(static_cast<TList*>(f->Get("list_pyxsect")));
  auto lh(static_cast<TList*>(f->Get("list_results")));
  f->Close();
//=============================================================================

  const auto dn(GetNormFactor(ln,l));
  const auto dr(GetJetEventWeight());
//=============================================================================

  TIter next(lh);
  TH1D *hb(nullptr);
  while ((hb = static_cast<TH1D*>(next()))) {
    const TString sh(hb->GetName());
    hb->SetName(Form("%s_%d",sh.Data(),l));
    hb->Scale(dn);

    auto ha(static_cast<TH1D*>(lr->FindObject(sh.Data())));
    if (ha) {
      ha->Add(hb);
      //ha->Scale(1./dr);
    } else {
      lr->Add(static_cast<TH1D*>(hb->Clone(sh.Data())));
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
Double_t GetNormFactor(TList *ln, const Int_t l)
{
  auto hn(static_cast<TH1D*>(ln->FindObject("hTrials")));
  hn->SetName(Form("%s_%d",hn->GetName(),l));
  
  //auto hn(static_cast<TH1D*>(lh->FindObject("hJEvent")));
  //hn->SetName(Form("%s_%d",hn->GetName(),l));

  auto hx(static_cast<TProfile*>(ln->FindObject("hXsect")));
  hx->SetName(Form("%s_%d",hx->GetName(),l));

  return (hx->GetBinContent(1) / hn->GetBinContent(1));
}

//_____________________________________________________________________________
Double_t GetJetEventWeight()
{
  auto r(0.);
  for (auto i=0; i<gknb; ++i){
    auto sf(GetFileName(i));
//=============================================================================

    auto f(TFile::Open(sf.Data(),"READ"));
    auto ln(static_cast<TList*>(f->Get("list_pyxsect")));
    f->Close();
    auto hn(static_cast<TH1D*>(ln->FindObject("hTrials")));
    hn->SetName(Form("%s_%d",hn->GetName(),i));
    
    auto hj(static_cast<TH1D*>(ln->FindObject("hJEvent")));
    hj->SetName(Form("%s_%d",hj->GetName(),i));

    auto hx(static_cast<TProfile*>(ln->FindObject("hXsect")));
    hx->SetName(Form("%s_%d",hx->GetName(),i));
    r = r + ((hx->GetBinContent(1) * hj->GetEntries()) / hn->GetBinContent(1));
  }

  return (r);
}
//_____________________________________________________________________________
TString GetFileName(const Int_t l)
{
  const auto k(l+1);
  const TString s((gklb[k]==-1) ? Form("%03d_INF",  gklb[l]) :
                                  Form("%03d_%03d", gklb[l], gklb[k]));

  return Form("%s/%s_PtHat_%s.root",gksp.Data(),gksf.Data(),s.Data());
}
