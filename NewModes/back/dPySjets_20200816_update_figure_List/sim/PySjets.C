#ifdef CMODE
#include "utils.h"

int main(int argc, char *argv[])
{
  TString seID(__FILE__);
  seID.ReplaceAll(".C", "");
  const TString ksp(Form("%s::%s",seID.Data(),__FUNCTION__));
  for (auto i=0; i<argc; ++i) ::Info(ksp.Data(), "argv[%d] = %s", i, argv[i]);
//=============================================================================

#if CMODE == 0
  ::Info(ksp.Data(), "Local dev / test mode");
#elif CMODE == 1
  ::Info(ksp.Data(), "Batch mode on CERN LXPLUS");
  assert(argc==3);
#elif CMODE == 2
  ::Info(ksp.Data(), "Batch mode on CCNU PC farm");
  assert(argc==5);
#else
  ::Error(ksp.Data(), "Compiled with unknown mode");
  exit(0);
#endif
//=============================================================================

#if CMODE == 0
  const auto ktID(-1);
#else
  const auto ktID(TString(argv[1]).Atoi());
#endif

  const auto bLC((ktID==0) || (ktID==2) || (ktID==3));
  const TString stID(bLC ? Form("BLCmode%d",ktID) : "Monash");

  if (bLC) {
    ::Info(ksp.Data(), "Tune:pp = BLC mode %d", ktID);
  } else {
    ::Info(ksp.Data(), "Tune:pp = Monash 2013");
  }
//=============================================================================

#if CMODE != 0
  const TString sjID(argv[4]);
  ::Info(ksp.Data(), "Job name = %s", sjID.Data());

  const auto kcID(TString(argv[2]).Atoi());
  const auto kpID(TString(argv[3]).Atoi());
  ::Info(ksp.Data(), "Cluster ID = %d", kcID);
  ::Info(ksp.Data(), "Process ID = %d", kpID);
#endif
//=============================================================================

  const auto kSeed(RandomSeed());
  ::Info(ksp.Data(), "Random seed = %d", kSeed);

#if CMODE == 0
  const auto kSSeed(kSeed);
#else
  const auto kSSeed(kSeed - kcID + 9*kpID);
#endif
//=============================================================================

  Pythia8::Pythia pythia;
  auto &pyReco(pythia.event);
  auto &pyInfo(pythia.info);

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

  pythia.readString("Tune:pp = 14");
  pythia.readString("SoftQCD:inelastic = on");

#if CMODE == 0
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Main:spareParm1 = 0.4");
  pythia.readString("Main:spareParm2 = 0.0");
  pythia.readString("Main:spareFlag1 = off");
#else
  auto LoadCfgFile([&ksp,&pythia](const TString sCfg) {
    if (gSystem->AccessPathName(sCfg.Data())) {
      ::Error(ksp.Data(), "No file %s",sCfg.Data());
      exit(-1);
    } else {
      pythia.readFile(sCfg.Data());
      ::Info(ksp.Data(), "Loaded file %s",sCfg.Data());
    }
  });

  LoadCfgFile(Form("etc/cfgs/%s.cmnd",seID.Data()));
  if (sjID != seID) LoadCfgFile(Form("etc/cfgs/cJobID/%s.cmnd",sjID.Data()));
  if (bLC) LoadCfgFile(Form("etc/cfgs/cTunes/Py%s.cmnd",stID.Data()));
#endif

  pythia.init();
//=============================================================================

  ::Info(ksp.Data(), "Beams:eCM = %f", pythia.parm("Beams:eCM"));

  const auto nEv(pythia.mode("Main:numberOfEvents"));
  ::Info(ksp.Data(), "N events  = %d", nEv);

  const auto djR(pythia.parm("Main:spareParm1"));
  ::Info(ksp.Data(), "Jet radius = %f", djR);

  const auto drs(pythia.parm("Main:spareParm2"));
  ::Info(ksp.Data(), "Rapidity shift = %f", drs);

  const auto bw2(pythia.flag("Main:spareFlag1"));
  ::Info(ksp.Data(), "Weigth 2D outputs = %s", bw2 ? "true" : "flase");
//=============================================================================

  const auto djPtMin (1.00);
  const auto djEtaMin(-0.35 + drs);
  const auto djEtaMax( 0.35 + drs);
  const auto djAcut  (0.6 * TMath::Pi() * djR *djR);

  const auto dcPtMin (0.15);
  const auto dcEtaMin(-0.90 + drs);
  const auto dcEtaMax( 0.90 + drs);

  const auto dsEtaMin(-0.75 + drs);
  const auto dsEtaMax( 0.75 + drs);

  const fastjet::GhostedAreaSpec aGhostSpec(1.5);
  const fastjet::AreaDefinition  aAreaDef(fastjet::active_area, aGhostSpec);
  const fastjet::JetDefinition   aJetDef(fastjet::antikt_algorithm, djR, fastjet::BIpt_scheme, fastjet::Best);
  const auto aSelEta(fastjet::SelectorEtaRange(djEtaMin,djEtaMax));

  std::vector<fastjet::PseudoJet> vConstis, vStrgs;
//=============================================================================

  auto list_pyxsect(new TList());
  auto hTrials(new TH1D("hTrials",     "", 1, 0., 1.));
  list_pyxsect->Add(hTrials);

  auto hXsect(new TProfile("hXsect",  "", 1, 0., 1.));
  list_pyxsect->Add(hXsect);
//=============================================================================

  auto list_evInfo(new TList());
  auto hPtHat(new TH1D("hPtHat", "", 1000, 0., 1000.));
  list_evInfo->Add(hPtHat);

  auto hEvent(new TH1D("hEvent", "", gknJets+2, -0.5, gknJets+1.5));
  list_evInfo->Add(hEvent);

  auto aEventX(hEvent->GetXaxis());
  aEventX->SetBinLabel(1, "ALL");
  aEventX->SetBinLabel(2, "NJ");
  for (unsigned long i=0, k=3; i<gknJets; ++i, ++k) aEventX->SetBinLabel(k, gksJets[i]);

  CallSumw2(list_evInfo);
//=============================================================================

  auto list_Jets(new TList());
  auto hJet(new TH1D("hJet", "", 500, 0., 500.));
  list_Jets->Add(hJet);

  auto hConsti(new TH1D("hConsti", "", 1000, 0., 100.));
  list_Jets->Add(hConsti);

  CallSumw2(list_Jets);
//=============================================================================

  auto list_weights(new TList());
  for (const auto &sj : gksJets) {
    const auto n(1+gknJCs+gknJCs+gknOCs);
    auto h(new TH1D(Form("h%s",sj.Data()), "", n, -0.5, -0.5+n));
    list_weights->Add(h);

    auto k(1);
    auto aX(h->GetXaxis());
    aX->SetBinLabel(k, "ALL");
    for (const auto &sc : gksJCs) aX->SetBinLabel(++k, Form("J%s",sc.Data()));
    for (const auto &sc : gksJCs) aX->SetBinLabel(++k, Form("P%s",sc.Data()));
    for (const auto &sc : gksOCs) aX->SetBinLabel(++k, Form("O%s",sc.Data()));
  }

  CallSumw2(list_weights);
//=============================================================================

  const auto ns(static_cast<unsigned int>(EStrg::Undef));

  TList *list_Strgs[ns];
  for (unsigned int i=0; i<ns; ++i) {
    const auto &ss(gksStrgs[i]);

    auto list(list_Strgs[i] = new TList());
    list->Add(new TH1D(Form("h%s", ss.Data()), "", 1000, 0., 100.));

    for (const auto &sj : gksJets) {
      list->Add(new TH1D(Form("h%s_%s",ss.Data(),sj.Data()), "", 1000, 0., 100.));

      for (const auto &sc : gksJCs) {
        list->Add(new TH1D(Form("h%s_%s_J%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
      }

       for (const auto &sc : gksJCs) {
        list->Add(new TH1D(Form("h%s_%s_P%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
      }

      for (const auto &sc : gksOCs) {
        list->Add(new TH1D(Form("h%s_%s_O%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
      }
    }

    list->Add(new TH1D(Form("h%s_NJ",ss.Data()), "", 1000, 0., 100.));
    CallSumw2(list);
  }
//=============================================================================

  TList *list_Strg2(nullptr);

  if (bw2) {
    list_Strg2 = new TList();
    for (const auto &ss : gksStrgs) for (const auto &sj : gksJets) {
      list_Strg2->Add(new TH2D(Form("h2%s_%s",ss.Data(),sj.Data()), "", 1500, dsEtaMin, dsEtaMax, 360, 0., d2PI));
    }

    CallSumw2(list_Strg2);
  }
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
      const auto dd(ap.pT()), de(ap.eta());

      if (ap.isFinal()   &&
          ap.isVisible() &&
          ap.isCharged() &&
         (dd>dcPtMin) && (de>dcEtaMin) && (de<dcEtaMax)) {
        vConstis.emplace_back(ap.px(),ap.py(),ap.pz(),ap.p().pAbs());
        vConstis.back().set_user_index(i);
        hConsti->Fill(dd);
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
      hJet->Fill(dj);
    }

    acEv.Reset();
    if (djMax>0.) {
      if (djMax<5.) { // Hard code: jet pT < 5 GeV/c for NJ
        acEv.SetPair("NJ",5.);
      } else {
        for (const auto &sj : gksJets) {
          const auto dj(JetVal(sj));
          if (djMax>dj) acEv.SetPair(sj,dj);
        }
      }
    }
//=============================================================================

    TVector3 vs, vg;
    for (const auto &as : vStrgs) {
      auto ps(as.user_info_shared_ptr());
      const auto es((static_cast<StrgInfo*>(ps.get()))->GetType());
      const auto ks(static_cast<unsigned int>(es));
      const auto &ss(gksStrgs[ks]);
      auto list(list_Strgs[ks]);

      const auto dd(as.pt());
      (static_cast<TH1D*>(list->FindObject(Form("h%s",ss.Data()))))->Fill(dd);

      if (acEv.Empty()) continue; // No jet event
//=============================================================================

      if (acEv.HasKey("NJ")) {
        (static_cast<TH1D*>(list->FindObject(Form("h%s_NJ",ss.Data()))))->Fill(dd);
        continue;
      }
//=============================================================================

      const auto de(as.eta()), df(as.phi());

      vs.SetPtEtaPhi(dd, de, df);
      vg.SetPtEtaPhi(dd, gRandom->Uniform(dsEtaMin,dsEtaMax),
                         gRandom->Uniform(0.,d2PI));
//=============================================================================

      for (const auto &sj : gksJets) if (acEv.HasKey(sj)) {
        (static_cast<TH1D*>(list->FindObject(Form("h%s_%s",ss.Data(),sj.Data()))))->Fill(dd);
        if (bw2) (static_cast<TH2D*>(list_Strg2->FindObject(Form("h2%s_%s",ss.Data(),sj.Data()))))->Fill(de,df);
        auto hWgt(static_cast<TH1D*>(list_weights->FindObject(Form("h%s",sj.Data()))));
        hWgt->AddBinContent(1, 1.);

        TVector3 vj;
        ReducedDict asJC, asPC;
        ReducedDict agJC, agPC;
        auto dso(-1.), dgo(-1.);
        const auto djCut(acEv.GetVal(sj));
        for (const auto &aj : vJets) if (aj.area()>djAcut) {
          const auto dj(aj.pt());
          if (dj<djCut) continue;

          vj.SetPtEtaPhi(dj, aj.eta(), aj.phi());
          auto MatchJets([&vj](const TVector3 &v, ReducedDict &aJC, ReducedDict &aPC, double &doc) {
            const auto ddj(vj.DeltaR(v));

            for (const auto &sc : gksJCs) {
              if (aJC.NoKey(sc)) {
                const auto dc(ConeVal(sc));
                if (ddj<dc) aJC.SetPair(sc,dc);
              }

              if (aPC.NoKey(sc)) {
                const auto dc(ConeVal(sc));
                if (PerpDeltaR(vj,v)<dc) aPC.SetPair(sc,dc);
              }
            }

            if ((1./ddj) > (1./doc)) doc = ddj;
          });

          MatchJets(vs, asJC, asPC, dso);
          MatchJets(vg, agJC, agPC, dgo);
        }
//=============================================================================

        auto Bname([] (const TString st, const TString &sc) {
          return (Form("%s%s",st.Data(),sc.Data()));
        });

        auto Hname([&ss,&sj,&Bname] (const TString st, const TString &sc) {
          return (Form("h%s_%s_%s",ss.Data(),sj.Data(),Bname(st,sc)));
        });

        for (const auto &sc : gksJCs) {
          if (asJC.HasKey(sc)) (static_cast<TH1D*>(list->FindObject(Hname("J",sc))))->Fill(dd);
          if (asPC.HasKey(sc)) (static_cast<TH1D*>(list->FindObject(Hname("P",sc))))->Fill(dd);

          if (agJC.HasKey(sc)) hWgt->AddBinContent(hWgt->GetXaxis()->FindBin(Bname("J",sc)), 1.);
          if (agPC.HasKey(sc)) hWgt->AddBinContent(hWgt->GetXaxis()->FindBin(Bname("P",sc)), 1.);
        }

        for (const auto &sc : gksOCs) if (dso>ConeVal(sc)) {
          (static_cast<TH1D*>(list->FindObject(Hname("O",sc))))->Fill(dd);
        }

        for (const auto &sc : gksOCs) if (dgo>ConeVal(sc)) {
          hWgt->AddBinContent(hWgt->GetXaxis()->FindBin(Bname("O",sc)), 1.);
        }
      }
    }
//=============================================================================

    hPtHat->Fill(pyInfo.pTHat());
    hEvent->AddBinContent(1,1.);
    for (auto k=2; k<=hEvent->GetNbinsX(); ++k) if (acEv.HasKey(aEventX->GetBinLabel(k))) hEvent->AddBinContent(k,1.);
  }
//=============================================================================

  hXsect ->Fill(0.5, pyInfo.sigmaGen());
  hTrials->Fill(0.5, pyInfo.weightSum());

  timer.Stop();
  timer.Print();
//=============================================================================

#if CMODE == 0
  list_pyxsect -> Print();
  list_evInfo  -> Print();
  list_weights -> Print();
  for (auto list : list_Strgs) list->Print();
  if (bw2) list_Strg2->Print();
#endif
//=============================================================================

#if CMODE == 0
  auto file(TFile::Open(Form("AnalysisResults_%d.root",kSeed),"NEW"));
#else
  const TString sPref((sjID == seID) ? seID : Form("%s_%s",seID.Data(),sjID.Data()));
  const TString sFile(Form("%s_%s_%d_%d",sPref.Data(),stID.Data(),kcID,kpID));

#if CMODE == 1
  auto file(TFile::Open(Form("%s.root",sFile.Data()),"NEW"));
#elif CMODE == 2
  auto file(TFile::Open(Form("out/%s_%d.root",sFile.Data(),kSeed),"NEW"));
#endif
#endif

  list_pyxsect -> Write("list_pyxsect", TObject::kSingleKey);
  list_evInfo  -> Write("list_evInfo", TObject::kSingleKey);

  for (unsigned int i=0; i<ns; ++i) {
    list_Strgs[i]->Write(Form("list_%s",gksStrgs[i].Data()), TObject::kSingleKey);
  }

  list_weights -> Write("list_weights", TObject::kSingleKey);
  file->Close();
//=============================================================================

  if (bw2) {
#if CMODE == 0
    file = TFile::Open(Form("AccptWgt_%d.root",kSeed),"NEW");
#else
#if CMODE == 1
    file = TFile::Open(Form("AccptWgt_%s.root",sFile.Data()),"NEW");
#elif CMODE == 2
    file = TFile::Open(Form("out/AccptWgt_%s_%d.root",sFile.Data(),kSeed),"NEW");
#endif
#endif

    list_Strg2->Write("list_weights", TObject::kSingleKey);
    file->Close();
  }
//=============================================================================

#if CMODE == 0
  DumpPyCfgs(ktID,seID,pythia);
#endif

  ::Info(ksp.Data(), "DONE");
//============================================================================

  return 0;
}

#endif
