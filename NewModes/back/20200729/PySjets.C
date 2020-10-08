#include "utils.h"

int main(int argc, char *argv[])
{
  assert(argc>4);
//=============================================================================

  TString sf(__FILE__);
  sf.ReplaceAll(".C", "");
  const TString ksp(Form("%s::%s",sf.Data(),__FUNCTION__));
  for (auto i=0; i<argc; ++i) ::Info(ksp.Data(), "argv[%d] = %s", i, argv[i]);
//=============================================================================

  const auto dJetR(TString(argv[1]).Atof());
  ::Info(ksp.Data(), "Jet radius = %f", dJetR);

  const auto dEshift(TString(argv[2]).Atof());
  ::Info(ksp.Data(), "Rapidity shift = %f", dEshift);

  const auto kClusID(TString(argv[3]).Atoi());
  const auto kProcID(TString(argv[4]).Atoi());
  ::Info(ksp.Data(), "Cluster ID = %d", kClusID);
  ::Info(ksp.Data(), "Process ID = %d", kProcID);

  const auto kSeed(RandomSeed());
  ::Info(ksp.Data(), "Random seed = %d", kSeed);
//=============================================================================

  Pythia8::Pythia pythia;
  auto &pyReco(pythia.event);
  auto &pyInfo(pythia.info);

  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");

  pythia.readString("333:mayDecay = off");
  pythia.readString("310:mayDecay = off");
  pythia.readString("3122:mayDecay = off");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d", kSeed - kClusID + 9*kProcID));

  pythia.readFile(Form("%s.cmnd",sf.Data()));
  const auto mode(pythia.mode("Main:spareMode1"));

  if ((mode==0) || (mode==2) || (mode==3)) {
    pythia.readFile(Form("PyBLCmode%d.cmnd",mode));
  }

  pythia.init();
//=============================================================================

  const auto dJetPtMin(1.00);
  const auto dJetEtaMin(-0.35 + dEshift);
  const auto dJetEtaMax( 0.35 + dEshift);
  const auto dJetAreaMin(0.6 * TMath::Pi() * dJetR *dJetR);

  const auto dConstiPtMin(0.15);
  const auto dConstiEtaMin(-0.90 + dEshift);
  const auto dConstiEtaMax( 0.90 + dEshift);

  const auto dStrgEtaMin(-0.75 + dEshift);
  const auto dStrgEtaMax( 0.75 + dEshift);

  const fastjet::GhostedAreaSpec aGhostSpec(1.5);
  const fastjet::AreaDefinition  aAreaDef(fastjet::active_area, aGhostSpec);
  const fastjet::JetDefinition   aJetDef(fastjet::antikt_algorithm, dJetR, fastjet::BIpt_scheme, fastjet::Best);
  const auto aSelEta(fastjet::SelectorEtaRange(dJetEtaMin,dJetEtaMax));

  std::vector<fastjet::PseudoJet> vConstis, vStrgs;
//=============================================================================

  //auto file(TFile::Open(Form("AnalysisResults_%d_%d_%d.root",kSeed,kClusID,kProcID),"NEW"));
  auto file(TFile::Open(Form("AnalysisResults_%d_%d.root", kClusID, kProcID),"NEW"));

  auto list_pyxsect(new TList());
  auto hTrials(new TH1D("hTrials",    "", 1, 0., 1.));
  list_pyxsect->Add(hTrials);

  auto hXsect(new TProfile("hXsect",  "", 1, 0., 1.));
  list_pyxsect->Add(hXsect);
//=============================================================================

  auto list_results(new TList());
  auto hPtHat(new TH1D("hPtHat", "", 1000, 0., 1000.));
  list_results->Add(hPtHat);

  auto hEvents(new TH1D("hEvents", "", gknJets+2, -0.5, gknJets+1.5));
  list_results->Add(hEvents);

  auto aX(hEvents->GetXaxis());
  aX->SetBinLabel(1, "ALL");
  aX->SetBinLabel(2, "NJ");
  for (unsigned long i=0, k=3; i<gknJets; ++i, ++k) aX->SetBinLabel(k, gksJets[i]);
//=============================================================================

  auto hJet(new TH1D("hJet", "", 500, 0., 500.));
  list_results->Add(hJet);

  auto hConsti(new TH1D("hConsti", "", 1000, 0., 100.));
  list_results->Add(hConsti);

  for (const auto &ss : gksStrgs) {
    list_results->Add(new TH1D(Form("h%s",   ss.Data()), "", 1000, 0., 100.));
    list_results->Add(new TH1D(Form("h%s_NJ",ss.Data()), "", 1000, 0., 100.));

    for (const auto &sj : gksJets) {
      list_results->Add(new TH1D(Form("h%s_%s",ss.Data(),sj.Data()), "", 1000, 0., 100.));

      for (const auto &sc : gksJCs) {
        list_results->Add(new TH1D(Form("h%s_%s_J%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
        list_results->Add(new TH1D(Form("h%s_%s_P%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
      }

      for (const auto &sc : gksOCs) {
        list_results->Add(new TH1D(Form("h%s_%s_O%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
      }
    }
  }
//=============================================================================

  TObject *p(nullptr);
  TIter next(list_results);
  while ((p = next())) {
    auto h(dynamic_cast<TH1*>(p));
    if (h) h->Sumw2();
  }
//=============================================================================

  TStopwatch timer;
  timer.Start();

  ReducedDict aCountEv, aCountJC, aCountPC;
  for (auto iEvent=0; iEvent<pythia.mode("Main:numberOfEvents"); ++iEvent) if (pythia.next()) {
    vConstis.resize(0);
    vStrgs.resize(0);

    for (auto i=0; i<pyReco.size(); ++i) {
      const auto &ap(pyReco[i]);
      const auto dpPt(ap.pT()), dpEta(ap.eta());

      if (ap.isFinal()       &&
          ap.isVisible()     &&
          ap.isCharged()     &&
         (dpPt>dConstiPtMin) && (dpEta>dConstiEtaMin) && (dpEta<dConstiEtaMax)) {
        vConstis.emplace_back(ap.px(),ap.py(),ap.pz(),ap.p().pAbs());
        vConstis.back().set_user_index(i);
        hConsti->Fill(dpPt);
      }
//=============================================================================

      if (!(ap.isHadron())) continue;
      if ((dpEta<dStrgEtaMin) || (dpEta>dStrgEtaMax)) continue;

      const auto ks(StrgType(ap,pyReco));
      if (ks==EStrg::Undef) continue;

      vStrgs.emplace_back(ap);
      vStrgs.back().set_user_index(i);
      vStrgs.back().set_user_info(new StrgInfo(ks));
    }
//=============================================================================

    fastjet::ClusterSequenceArea acs(vConstis, aJetDef, aAreaDef);
    const auto &vJets(aSelEta(acs.inclusive_jets(dJetPtMin)));

    auto djMax(-1.);
    for (const auto &aj : vJets) if (aj.area()>dJetAreaMin) {
      const auto dj(aj.pt());

      hJet->Fill(dj);
      if (dj>djMax) djMax = dj;
    }

    aCountEv.Reset();
    if (djMax<5.) { // Hard code: jet pT < 5 GeV/c for NJ
      aCountEv.SetPair("NJ",5.);
    } else {
      for (const auto &sj : gksJets) {
        const auto dj(JetVal(sj));
        if (djMax>dj) aCountEv.SetPair(sj,dj);
      }
    }
//=============================================================================

    TVector3 vs, vj;
    for (const auto &as : vStrgs) {
      auto ps(as.user_info_shared_ptr());
      const auto ks((static_cast<StrgInfo*>(ps.get()))->GetType());
      const auto ss(StrgName(ks));

      const auto ds(as.pt());
      vs.SetPtEtaPhi(ds, as.eta(), as.phi());
      (static_cast<TH1D*>(list_results->FindObject(Form("h%s",ss.Data()))))->Fill(ds);

      if (aCountEv.HasKey("NJ")) {
        (static_cast<TH1D*>(list_results->FindObject(Form("h%s_NJ",ss.Data()))))->Fill(ds);
        continue;
      }
//=============================================================================

      for (const auto &sj : gksJets) if (aCountEv.HasKey(sj)) {
        const auto djCut(aCountEv.GetVal(sj));
        (static_cast<TH1D*>(list_results->FindObject(Form("h%s_%s",ss.Data(),sj.Data()))))->Fill(ds);

        auto doc(-1.);
        aCountJC.Reset();
        aCountPC.Reset();
        for (const auto &aj : vJets) if (aj.area()>dJetAreaMin) {
          const auto dj(aj.pt());
          if (dj<djCut) continue;

          vj.SetPtEtaPhi(dj, aj.eta(), aj.phi());
          const auto dsj(vj.DeltaR(vs));

          for (const auto &sc : gksJCs) {
            if (aCountJC.NoKey(sc)) {
              const auto dc(ConeVal(sc));
              if (dsj < dc) aCountJC.SetPair(sc,dc);
            }

            if (aCountPC.NoKey(sc)) {
              const auto dc(ConeVal(sc));
              if (PerpDeltaR(vj,vs) < dc) aCountPC.SetPair(sc,dc);
            }
          }

          if ((1./dsj) > (1./doc)) doc = dsj;
        }
//=============================================================================

        auto Hname([&ss,&sj] (const TString st, const TString &sc) {
          return (Form("h%s_%s_%s%s",ss.Data(),sj.Data(),st.Data(),sc.Data()));
        });

        for (const auto &sc : gksJCs) {
          if (aCountJC.HasKey(sc)) (static_cast<TH1D*>(list_results->FindObject(Hname("J",sc))))->Fill(ds);
          if (aCountPC.HasKey(sc)) (static_cast<TH1D*>(list_results->FindObject(Hname("P",sc))))->Fill(ds);
        }

        for (const auto &sc : gksOCs) if (doc>ConeVal(sc)) {
          (static_cast<TH1D*>(list_results->FindObject(Hname("O",sc))))->Fill(ds);
        }
      }
    }
//=============================================================================

    hEvents->AddBinContent(1,1.);
    for (auto k=2; k<=hEvents->GetNbinsX(); ++k)
      if (aCountEv.HasKey(aX->GetBinLabel(k))) hEvents->AddBinContent(k,1.);

    hPtHat->Fill(pyInfo.pTHat());
  }
//=============================================================================

  timer.Stop();
  timer.Print();
//=============================================================================

  hXsect ->Fill(0.5, pyInfo.sigmaGen());
  hTrials->Fill(0.5, pyInfo.weightSum());

  list_pyxsect->Print();
  list_results->Print();

  file->cd();
  list_pyxsect->Write("list_pyxsect", TObject::kSingleKey);
  list_results->Write("list_results", TObject::kSingleKey);
  file->Close();
//=============================================================================

  DumpPyCfgs(pythia, sf);
  ::Info(ksp.Data(), "DONE");
//=============================================================================

  return 0;
}
