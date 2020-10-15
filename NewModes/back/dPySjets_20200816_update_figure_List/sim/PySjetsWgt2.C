#ifdef CMODE
#include "utils.h"

int main(int argc, char *argv[])
{
  TString sf(__FILE__);
  sf.ReplaceAll(".C", "");
  const TString ksp(Form("%s::%s",sf.Data(),__FUNCTION__));
  for (auto i=0; i<argc; ++i) ::Info(ksp.Data(), "argv[%d] = %s", i, argv[i]);
//=============================================================================

#if CMODE == 0
  ::Info(ksp.Data(), "Local dev / test mode");
#elif CMODE == 1
  ::Info(ksp.Data(), "Batch mode on CERN LXPLUS");
  assert(argc==3);
#elif CMODE == 2
  ::Info(ksp.Data(), "Batch mode on CCNU PC farm");
  assert(argc==6);
#else
  ::Error(ksp.Data(), "Compiled with unknown mode");
  exit(0);
#endif
//=============================================================================

#if CMODE == 0
  const auto modeBLC(-1);
  const auto djR(0.4), dEs(0.);
#else
  const auto djR(TString(argv[1]).Atof());
  const auto dEs(TString(argv[2]).Atof());
  const auto modeBLC(TString(argv[3]).Atoi());
  const auto kClusID(TString(argv[4]).Atoi());
  const auto kProcID(TString(argv[5]).Atoi());
#endif

  const auto bLC((modeBLC==0) || (modeBLC==2) || (modeBLC==3));
  const TString kTune(bLC ? Form("BLCmode%d",modeBLC) : "Monash");

  if (bLC) {
    ::Info(ksp.Data(), "Tune:pp = BLC mode %d", modeBLC);
  } else {
    ::Info(ksp.Data(), "Tune:pp = Monash 2013");
  }
//=============================================================================

  ::Info(ksp.Data(), "Jet radius = %f", djR);
  ::Info(ksp.Data(), "Rapidity shift = %f", dEs);

#if CMODE != 0
  ::Info(ksp.Data(), "Cluster ID = %d", kClusID);
  ::Info(ksp.Data(), "Process ID = %d", kProcID);
#endif
//=============================================================================

  const auto kSeed(RandomSeed());
  ::Info(ksp.Data(), "Random seed = %d", kSeed);

#if CMODE == 0
  const auto kSSeed(kSeed);
#else
  const auto kSSeed(kSeed - kClusID + 9*kProcID);
#endif
//=============================================================================

  Pythia8::Pythia pythia;
  auto &pyReco(pythia.event);

#if CMODE != 0
  pythia.readString("Next:numberCount = 0");
#endif

  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");

  pythia.readString("333:mayDecay = off");
  pythia.readString("310:mayDecay = off");
  pythia.readString("3122:mayDecay = off");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d", kSSeed));

  pythia.readFile(Form("cfgs/%s.cmnd",sf.Data()));
  if (bLC) pythia.readFile(Form("cfgs/Py%s.cmnd",kTune.Data()));

  pythia.init();
//=============================================================================

  const auto djPtMin (1.00);
  const auto djEtaMin(-0.35 + dEs);
  const auto djEtaMax( 0.35 + dEs);
  const auto djAcut  (0.6 * TMath::Pi() * djR *djR);

  const auto dcPtMin (0.15);
  const auto dcEtaMin(-0.90 + dEs);
  const auto dcEtaMax( 0.90 + dEs);

  const auto dsEtaMin(-0.75 + dEs);
  const auto dsEtaMax( 0.75 + dEs);

  const fastjet::GhostedAreaSpec aGhostSpec(1.5);
  const fastjet::AreaDefinition  aAreaDef(fastjet::active_area, aGhostSpec);
  const fastjet::JetDefinition   aJetDef(fastjet::antikt_algorithm, djR, fastjet::BIpt_scheme, fastjet::Best);
  const auto aSelEta(fastjet::SelectorEtaRange(djEtaMin,djEtaMax));

  std::vector<fastjet::PseudoJet> vConstis, vStrgs;
//=============================================================================

  auto file(TFile::Open(Form("wgt/AccptWgt_%s.root",kTune.Data()),"NEW"));
  auto list_weights(static_cast<TList*>(file->Get("list_weights")));
  file->Close();
//=============================================================================

#if CMODE == 0
  file = TFile::Open(Form("AnalysisResults_%d.root",kSeed),"NEW");
#else
  const TString kFileID(Form("%d_%d",kClusID,kProcID));

#if CMODE == 1
  file = TFile::Open(Form("AnalysisResults_%s.root",kFileID.Data()),"NEW");
#elif CMODE == 2
  const TString kPathID(Form("out/%s_%s",sf.Data(),kTune.Data()));
  file = TFile::Open(Form("%s/AnalysisResults_%s_%d.root",kPathID.Data(),kFileID.Data(),kSeed),"NEW");
#endif
#endif
//=============================================================================

  auto list_weight2(new TList());
  for (const auto &ss : gksStrgs) for (const auto &sj : gksJets) {
    const auto n(1+gknJCs+gknJCs+gknOCs);
    auto h(new TH1D(Form("hw%s_%s",ss.Data(),sj.Data()), "", n, -0.5, -0.5+n));
    list_weight2->Add(h);

    auto k(1);
    auto aX(h->GetXaxis());
    aX->SetBinLabel(k, "ALL");
    for (const auto &sc : gksJCs) aX->SetBinLabel(++k, Form("J%s",sc.Data()));
    for (const auto &sc : gksJCs) aX->SetBinLabel(++k, Form("P%s",sc.Data()));
    for (const auto &sc : gksOCs) aX->SetBinLabel(++k, Form("O%s",sc.Data()));
  }

  CallSumw2(list_weight2);
//=============================================================================

  if (gRandom) {
    delete gRandom;
    gRandom = nullptr;
  }

  gRandom = new TRandom3(kSSeed);
//=============================================================================

  TStopwatch timer;
  timer.Start();

  ReducedDict acEv;
  for (auto iEvent=0; iEvent<pythia.mode("Main:numberOfEvents"); ++iEvent) if (pythia.next()) {
    vConstis.resize(0);
    vStrgs.resize(0);

    for (auto i=0; i<pyReco.size(); ++i) {
      const auto &ap(pyReco[i]);
      const auto de(ap.eta());

      if (ap.isFinal()     &&
          ap.isVisible()   &&
          ap.isCharged()   &&
         (ap.pT()>dcPtMin) && (de>dcEtaMin) && (de<dcEtaMax)) {
        vConstis.emplace_back(ap.px(),ap.py(),ap.pz(),ap.p().pAbs());
        vConstis.back().set_user_index(i);
      }
//=============================================================================

      if (!(ap.isHadron())) continue;
      if ((de<dsEtaMin) || (de>=dsEtaMax)) continue;

      const auto ks(StrgType(ap,pyReco));
      if (ks==EStrg::Undef) continue;

      vStrgs.emplace_back(ap);
      vStrgs.back().set_user_index(i);
      vStrgs.back().set_user_info(new StrgInfo(ks));
    }
//=============================================================================

    fastjet::ClusterSequenceArea acs(vConstis, aJetDef, aAreaDef);
    const auto &vJets(aSelEta(acs.inclusive_jets(djPtMin)));

    auto djMax(-1.);
    if (vJets.size()>0) for (const auto &aj : vJets) if (aj.area()>djAcut) {
      const auto dj(aj.pt());
      if (dj>djMax) djMax = dj;
    }

    acEv.Reset();
    if (djMax>0.) for (const auto &sj : gksJets) {
      const auto dj(JetVal(sj));
      if (djMax>dj) acEv.SetPair(sj,dj);
    }

    if (acEv.Empty()) continue;
//=============================================================================

    TVector3 vg;
    for (const auto &as : vStrgs) {
      auto ps(as.user_info_shared_ptr());
      const auto es((static_cast<StrgInfo*>(ps.get()))->GetType());
      const auto ks(static_cast<unsigned int>(es));
      const auto &ss(gksStrgs[ks]);
      const auto dd(as.pt());
//=============================================================================

      for (const auto &sj : gksJets) if (acEv.HasKey(sj)) {
        auto hWgt(static_cast<TH1D*>(list_weight2->FindObject(Form("hw%s_%s",ss.Data(),sj.Data()))));

        Double_t de(0.), df(0.);
        (static_cast<TH2D*>(list_weights->FindObject(Form("h2%s_%s",ss.Data(),sj.Data()))))->GetRandom2(de,df);
        vg.SetPtEtaPhi(dd, de, df);

        TVector3 vj;
        auto doc(-1.);
        ReducedDict agJC, agPC;
        const auto djCut(acEv.GetVal(sj));
        for (const auto &aj : vJets) if (aj.area()>djAcut) {
          const auto dj(aj.pt());
          if (dj<djCut) continue;

          vj.SetPtEtaPhi(dj, aj.eta(), aj.phi());
          const auto ddj(vj.DeltaR(vg));

          for (const auto &sc : gksJCs) {
            if (agJC.NoKey(sc)) {
              const auto dc(ConeVal(sc));
              if (ddj < dc) agJC.SetPair(sc,dc);
            }

            if (agPC.NoKey(sc)) {
              const auto dc(ConeVal(sc));
              if (PerpDeltaR(vj,vg) < dc) agPC.SetPair(sc,dc);
            }
          }

          if ((1./ddj) > (1./doc)) doc = ddj;
        }
//=============================================================================

        auto Bname([] (const TString st, const TString &sc) {
          return (Form("%s%s",st.Data(),sc.Data()));
        });

        for (const auto &sc : gksJCs) {
          if (agJC.HasKey(sc)) hWgt->AddBinContent(hWgt->GetXaxis()->FindBin(Bname("J",sc)), 1.);
          if (agPC.HasKey(sc)) hWgt->AddBinContent(hWgt->GetXaxis()->FindBin(Bname("P",sc)), 1.);
        }

        for (const auto &sc : gksOCs) if (doc>ConeVal(sc)) {
          hWgt->AddBinContent(hWgt->GetXaxis()->FindBin(Bname("O",sc)), 1.);
        }
      }
    }
  }
//=============================================================================

  timer.Stop();
  timer.Print();
//=============================================================================

  file->cd();
  list_weight2->Write("list_weight2", TObject::kSingleKey);
  file->Close();
//=============================================================================

#if CMODE != 0
  DumpPyCfgs(modeBLC,sf,pythia);
#endif

  ::Info(ksp.Data(), "DONE");
//=============================================================================

  return 0;
}

#endif
