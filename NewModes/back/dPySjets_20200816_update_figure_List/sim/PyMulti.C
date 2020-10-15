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
//=============================================================================

  const auto dcMin(-1.), dcMax(1.);

  auto list_pyxsect(new TList());
  auto hTrials(new TH1D("hTrials",     "", 1, 0., 1.));
  list_pyxsect->Add(hTrials);

  auto hXsect(new TProfile("hXsect",  "", 1, 0., 1.));
  list_pyxsect->Add(hXsect);
//=============================================================================

  auto list_results(new TList());
  auto hPtHat(new TH1D("hPtHat", "", 1000, 0., 1000.));
  list_results->Add(hPtHat);

  auto hMulti(new TH1D("hMulti", "", 1000, -0.5, 999.5));
  list_results->Add(hMulti);

  CallSumw2(list_results);
//=============================================================================

  TStopwatch timer;
  timer.Start();

  for (auto iEvent=0; iEvent<pythia.mode("Main:numberOfEvents"); ++iEvent) if (pythia.next()) {

    auto dNch(0.);
    for (auto i=0; i<pyReco.size(); ++i) {
      const auto &ap(pyReco[i]);

      if (ap.isFinal() && ap.isVisible() && ap.isCharged()) {
        const auto de(ap.eta());
        if ((de>dcMin) && (de<=dcMax)) dNch += 1.;
      }
    }

    hMulti->Fill(dNch);
    hPtHat->Fill(pyInfo.pTHat());
  }
//=============================================================================

  hXsect ->Fill(0.5, pyInfo.sigmaGen());
  hTrials->Fill(0.5, pyInfo.weightSum());

  timer.Stop();
  timer.Print();
//=============================================================================

#if CMODE == 0
  list_pyxsect->Print();
  list_results->Print();
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

  list_pyxsect->Write("list_pyxsect", TObject::kSingleKey);
  list_results->Write("list_results", TObject::kSingleKey);
  file->Close();
//=============================================================================

#if CMODE == 0
  DumpPyCfgs(ktID,seID,pythia);
#endif

  ::Info(ksp.Data(), "DONE");
//============================================================================

  return 0;
}

#endif
