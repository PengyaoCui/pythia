#include "utils.h"
using namespace Pythia8;

//=============================================================================

int main(int argc, char *argv[])
{
  assert(argc>1);
//=============================================================================

  TString sf(__FILE__);
  sf.ReplaceAll(".C", "");
  const TString gksp(Form("%s::%s",sf.Data(),__FUNCTION__));
  for (auto i=0; i<argc; ++i) ::Info(gksp.Data(), "argv[%d] = %s", i, argv[i]);
//=============================================================================

  const auto dJetR(TString(argv[1]).Atof());
  ::Info(gksp.Data(), "Jet radius = %f", dJetR);

  const auto dEtaShift(TString(argv[2]).Atof());
  ::Info(gksp.Data(), "Rapidity shift = %f", dEtaShift);

  const auto kClusID(TString(argv[3]).Atoi());
  const auto kProcID(TString(argv[4]).Atoi());
  ::Info(gksp.Data(), "Cluster ID = %d", kClusID);
  ::Info(gksp.Data(), "Process ID = %d", kProcID);

  const auto kSeed(RandomSeed());
  ::Info(gksp.Data(), "Random seed = %d", kSeed);
//=============================================================================

//TApplication app(gksp.Data(), &argc, argv);
//=============================================================================

  Pythia8::Pythia pythia;
  auto &pyReco(pythia.event);
  auto &pyInfo(pythia.info);

  pythia.readFile(Form("%s.cmnd",sf.Data()));

  //pythia.readString("Tune:pp = 14");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  //pythia.readString("ParticleDecays:limitTau0 = 0");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("333:mayDecay = off");

  pythia.readString("310:mayDecay = off");// \KShort will not decay
  pythia.readString("3122:mayDecay = off");// \Lambda will not decay

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d", kSeed - kClusID + 9*kProcID));


//============================================================================
  // There will be eight centrality bins based on the sum transverse
  // emergy in a rapidity interval between -4.9 and -3.2. The borders
  // between the classes have been read off the plot in the paper:
  double explim[] = {90.0, 66.0, 53.0, 41.0, 32.0, 24.0, 13.0, 6.0};
  // Alternatively we can obtain the borders from the generated
  // transverse energy spectrum. The default settings should give
  // approximately the following:
  double genlim[] = {77.5, 54.5, 44.4, 34.1, 27.6, 22.5, 14.2, 3.5};
  // If you change any parameters these should also be changed.

  // The upper edge of the correponding percentiles:
  double pclim[] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9};

  // Book the pseudorapidity histograms and get counter for sum of
  // event weights:
  typedef map<double,int,std::greater<double> > MapIdx;
  MapIdx expetaidx, genetaidx;
  vector<Hist*> expetadist(8), genetadist(8);
  string expetaname("EtadistCexp"), genetaname("EtadistCgen");
  vector<double> expsumw(8, 0.0), gensumw(8, 0.0);
  for ( int i = 0; i < 8; ++i ) {
    expetaidx[explim[i]] = i;
    expetadist[i] = new Hist(expetaname + char('0' + i), 54, -2.7, 2.7);
    genetaidx[genlim[i]] = i;
    genetadist[i] = new Hist(genetaname + char('0' + i), 54, -2.7, 2.7);
  }

  // Book histogram for the centrality measure.
  Hist sumet("SumETfwd", 100, 0.0, 200.0);
  // Also make a map of all weight to check the generated centrality
  // classes.
  multimap<double,double> gencent;

  // Book a histogram for the distribution of number of wounded
  // nucleons.
  Hist wounded("Nwounded", 60, -0.5, 59.5);

  // Sum up the weights of all generated events.
  double sumw = 0.0;

  // Initialise Pythia.
//============================================================================
  pythia.init();
//=============================================================================

  const auto dJetPtMin(1.00);
  const auto dJetEtaMin(-0.35 + dEtaShift);
  const auto dJetEtaMax( 0.35 + dEtaShift);
  const auto dJetAreaMin(0.6 * TMath::Pi() * dJetR *dJetR);

  const auto dConstiPtMin(0.15);
  const auto dConstiEtaMin(-0.90 + dEtaShift);
  const auto dConstiEtaMax( 0.90 + dEtaShift);

  const auto dStrgEtaMin(-0.75 + dEtaShift);
  const auto dStrgEtaMax( 0.75 + dEtaShift);

  const fastjet::GhostedAreaSpec aGhostSpec(1.5);
  const fastjet::AreaDefinition  aAreaDef(fastjet::active_area, aGhostSpec);
  const fastjet::JetDefinition   aJetDef(fastjet::antikt_algorithm, dJetR, fastjet::BIpt_scheme, fastjet::Best);

  const auto aSelEta(fastjet::SelectorEtaRange(dJetEtaMin,dJetEtaMax));

  std::vector<fastjet::PseudoJet> vConstis, vStrgs;
//=============================================================================

  auto file(TFile::Open(Form("AnalysisResults_%d_%d.root",kClusID,kProcID),"NEW"));
  //auto file(TFile::Open(Form("AnalysisResults_%d_%d_%d.root",kSeed,kClusID,kProcID),"NEW"));

  auto list_pyxsect(new TList());
  auto hTrials(new TH1D("hTrials",     "", 1, 0., 1.));
  list_pyxsect->Add(hTrials);

  auto hXsect (new TProfile("hXsect",  "", 1, 0., 1.));
  list_pyxsect->Add(hXsect);

  auto list_results(new TList());
  auto hPtHat(new TH1D("hPtHat", "", 1000, 0., 1000.));
  list_results->Add(hPtHat);

  auto hJet(new TH1D("hJet", "", 500, 0., 500.));
  list_results->Add(hJet);

  auto hConsti(new TH1D("hConsti", "", 1000, 0., 100.));
  list_results->Add(hConsti);

  for (const auto &ss : gksStrgs) {
    list_results->Add(new TH1D(Form("h%s",ss.Data()), "", 1000, 0., 100.));

    for (const auto &sj : gksJets) for (const auto &sc : gksStrgJCs) {
      list_results->Add(new TH1D(Form("h%s_%s_%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));
      list_results->Add(new TH1D(Form("h%s_PC%s_%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));//NEW
      list_results->Add(new TH1D(Form("h%s_OC%s_%s",ss.Data(),sj.Data(),sc.Data()), "", 1000, 0., 100.));//NEW
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

  TStopwatch timer; timer.Start();
  for (auto iEvent=0; iEvent<pythia.mode("Main:numberOfEvents"); ++iEvent) if (pythia.next()) {

    vConstis.resize(0);
    vStrgs.resize(0);
    double etfwd = 0.0;
    bool trigfwd = false;
    bool trigbwd = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() ) {
        double eta = p.eta();
        if ( p.isCharged() && p.pT() > 0.1 && eta < -2.09 && eta > -3.84 )
          trigfwd = true;
        if ( p.isCharged() && p.pT() > 0.1 && eta > 2.09 && eta < 3.84 )
          trigbwd = true;
        if ( p.pT() > 0.1 && eta < -3.2 && eta > -4.9 )
          etfwd += p.eT();
      }
    }
    // Skip if not triggered
    if ( !(trigfwd && trigbwd) ) continue;

    // Keep track of the sum of waights
    double weight = pythia.info.weight();
    sumw += weight;

    // Histogram and save the summed Et.
    sumet.fill(etfwd, weight);
    gencent.insert(make_pair(etfwd, weight));
    // Also fill the number of (absorptively and diffractively)
    // wounded nucleaons.
    int nw = pythia.info.hiinfo->nAbsTarg() +
             pythia.info.hiinfo->nDiffTarg();
    wounded.fill(nw, weight);

    // Find the correct centrality histograms.
    MapIdx::iterator expit =  expetaidx.upper_bound(etfwd);
    int expidx = expit== expetaidx.end()? -1: expit->second;
    MapIdx::iterator genit = genetaidx.upper_bound(etfwd);
    int genidx = genit== genetaidx.end()? -1: genit->second;

    // Sum the weights in the centrality classes, skip if not in a class.
    if ( expidx < 0 && genidx < 0 ) continue;
    if ( expidx >= 0 ) expsumw[expidx] += weight;
    if ( genidx >= 0 ) gensumw[genidx] += weight;

    // Go through the event again and fill the eta distributions.
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle & p = pythia.event[i];
      if ( p.isFinal() && p.isCharged() &&
           abs(p.eta()) < 2.7 && p.pT() > 0.1 ) {
        if ( expidx >= 0 ) expetadist[expidx]->fill(p.eta(), weight);
        if ( genidx >= 0 ) genetadist[genidx]->fill(p.eta(), weight);
      }
    }

//=============================================================================

    for (auto i=0; i<pyReco.size(); ++i) {
      const auto &ap(pyReco[i]);
      const auto dpPt(ap.pT()), dpEta(ap.eta());

      if (ap.isFinal()       &&
          ap.isVisible()     &&
          ap.isCharged()     &&
         (dpPt>dConstiPtMin) && (dpEta>dConstiEtaMin)  && (dpEta<dConstiEtaMax)) {
        vConstis.emplace_back(ap.px(),ap.py(),ap.pz(),ap.p().pAbs());
        vConstis.back().set_user_index(i);
        hConsti->Fill(dpPt);
      }
//=============================================================================

      if (!(ap.isHadron())) continue;
      if ((dpEta<dStrgEtaMin) || (dpEta>dStrgEtaMax)) continue;
//=============================================================================

      auto ks(EStrg::Undef);
      const auto id(ap.id());
      if (id==310) ks = EStrg::Kshort;

      if (TMath::Abs(id)==3122) {
        const auto m(ap.mother1());

        if (id>0) {
          if (m>0 ? (pyReco[m].id()==3312) || (pyReco[m].id()==3322) : false)
            ks = EStrg::LambdaFd;
          else
            ks = EStrg::Lambda;
        } else {
          if (m>0 ? (pyReco[m].id()==-3312) || (pyReco[m].id()==-3322) : false)
            ks = EStrg::AntiLaFd;
          else
            ks = EStrg::AntiLa;
        }
      }

      if (id==3312) ks = EStrg::XiNeg;

      if (id==-3312) ks = EStrg::XiPos;

      if (id==3334) ks = EStrg::OmegaNeg;

      if (id==-3334) ks = EStrg::OmegaPos;

      if (ks==EStrg::Undef) continue;
//=============================================================================

      vStrgs.emplace_back(ap);
      vStrgs.back().set_user_index(i);
      vStrgs.back().set_user_info(new StrgInfo(ks));

      const auto ss(StrgName(ks));
      if (!ss.IsNull()) (static_cast<TH1D*>(list_results->FindObject(Form("h%s",ss.Data()))))->Fill(dpPt);
    }
//=============================================================================

    fastjet::ClusterSequenceArea acs(vConstis, aJetDef, aAreaDef);
    const auto &vJets(aSelEta(acs.inclusive_jets(dJetPtMin)));
    for (const auto &aj : vJets) if (aj.area()>dJetAreaMin) hJet->Fill(aj.pt());
//=============================================================================

    for (const auto &av : vStrgs) {
      auto ps(av.user_info_shared_ptr());
      const auto ks((static_cast<StrgInfo*>(ps.get()))->GetType());
      const auto ss(StrgName(ks));

      TVector3 strg, vj, vl1, vl2, vu1, vu2;
      strg.SetPtEtaPhi(av.pt(), av.eta(), av.phi());

      bool bJC[gknStrgJCs];
      bool bPC[gknStrgJCs];
      for (const auto &sj : gksJets) {
        for (auto &b : bJC) b = false;
        
	for (auto &b : bPC) b = false;

        const auto dJetPtCut(JetVal(sj));
        for (const auto &aj : vJets) if (aj.area()>dJetAreaMin) {
          const auto dj(aj.pt());
          if (dj<dJetPtCut) continue;
          vj.SetPtEtaPhi(dj, aj.eta(), aj.phi());
          vl1.SetPtEtaPhi(dj, aj.eta(), aj.phi()); vl1.RotateZ(TMath::PiOver2());
          vl2.SetPtEtaPhi(dj, aj.eta(), aj.phi()); vl2.RotateZ(-1.*TMath::PiOver2());
          vu1.SetPtEtaPhi(dj, -1.*aj.eta(), aj.phi()); vu1.RotateZ(TMath::PiOver2());
          vu2.SetPtEtaPhi(dj, -1.*aj.eta(), aj.phi()); vu2.RotateZ(-1.*TMath::PiOver2());
          
	  double d(vj.DeltaR(strg));
	  double dl1(vl1.DeltaR(strg));
	  double dl2(vl2.DeltaR(strg));
	  double du1(vu1.DeltaR(strg));
	  double du2(vu2.DeltaR(strg));
          for (unsigned long i=0; i<gknStrgJCs; ++i) {
            if (d<StrgJC(gksStrgJCs[i])) bJC[i] = true;
	    if (dl1<StrgJC(gksStrgJCs[i]) || dl2<StrgJC(gksStrgJCs[i]) || du1<StrgJC(gksStrgJCs[i]) || du2<StrgJC(gksStrgJCs[i])) bPC[i] = true; 
	  }
	}

        for (unsigned long i=0; i<gknStrgJCs; ++i) {
	  if (bJC[i]) {
            const TString s(Form("h%s_%s_%s",ss.Data(),sj.Data(),gksStrgJCs[i].Data()));
            (static_cast<TH1D*>(list_results->FindObject(s.Data())))->Fill(av.pt());
	  }else{
            const TString s(Form("h%s_OC%s_%s",ss.Data(),sj.Data(),gksStrgJCs[i].Data()));
            (static_cast<TH1D*>(list_results->FindObject(s.Data())))->Fill(av.pt());
	  }
	  if(bPC[i]){
            const TString s(Form("h%s_PC%s_%s",ss.Data(),sj.Data(),gksStrgJCs[i].Data()));
            (static_cast<TH1D*>(list_results->FindObject(s.Data())))->Fill(av.pt());
	  }
	} 
      }
    }
//=============================================================================

    hPtHat->Fill(pyInfo.pTHat());
  }

  timer.Stop();
  timer.Print();
//=============================================================================

  hXsect->Fill(0.5, pyInfo.sigmaGen());
  hTrials->Fill(0.5, pyInfo.weightSum());
//=============================================================================

  list_pyxsect->Print();
  list_results->Print();
//=============================================================================

  file->cd();
  list_pyxsect->Write("list_pyxsect", TObject::kSingleKey);
  list_results->Write("list_results", TObject::kSingleKey);
  file->Close();
//=============================================================================

  ::Info(gksp.Data(), "DONE");
//=============================================================================

  return 0;
}
