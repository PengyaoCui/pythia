#ifndef ROOT_PY_FJ_UTILS_H
#define ROOT_PY_FJ_UTILS_H

#include <cassert>

#include "TSystem.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "TMath.h"
#include "TVector3.h"
#include "TString.h"

#include "TFile.h"
#include "TList.h"

#include "TH2D.h"
#include "TProfile.h"

#include "Pythia8/Pythia.h"
//nclude "Pythia8/ColourReconnection.h"

#include "Pythia8Plugins/FastJet3.h"
//nclude "Pythia8Plugins/ColourReconnectionHooks.h"

#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

//nclude "fastjet/tools/JetMedianBackgroundEstimator.hh"
//nclude "fastjet/tools/Subtractor.hh"
//nclude "fastjet/tools/Filter.hh"
//=============================================================================
//=============================================================================


const TString gksJets[] = {
  "Jet08", "Jet09", "Jet10", "Jet11", "Jet12",
  "Jet15", "Jet18", "Jet20", "Jet22" };
const auto gknJets(sizeof(gksJets) / sizeof(TString));

//_____________________________________________________________________________
float JetVal(TString s)
{
  s.ReplaceAll("Jet", "");
  return s.Atof();
}
//=============================================================================

const TString gksJCs[] = { "C02", "C03", "C04" };
const auto gknJCs(sizeof(gksJCs) / sizeof(TString));

const TString gksOCs[] = { "C04", "C06", "C08" };
const auto gknOCs(sizeof(gksOCs) / sizeof(TString));

//_____________________________________________________________________________
float ConeVal(TString s)
{
  s.ReplaceAll("C0", "");
  return (s.Atof() / 10.);
}
//=============================================================================
//=============================================================================

const TString gksStrgs[] = { "Kshort", "Lambda", "LambdaFd", "Xi", "Omega" };
enum struct EStrg : unsigned int { Kshort=0, Lambda, LambdaFd, Xi, Omega, Undef };
//=============================================================================

//_____________________________________________________________________________
EStrg StrgType(const Pythia8::Particle &particle, Pythia8::Event &event)
{
  const auto id(particle.idAbs());
  if (id==310) return EStrg::Kshort;
//=============================================================================

  if (id==3122) {
    const auto m(particle.mother1());
    if (m>0 ? (event[m].idAbs()==3312) || (event[m].idAbs()==3322) : false) {
      return EStrg::LambdaFd;
    } else {
      return EStrg::Lambda;
    }
  }
//=============================================================================

  if (id==3312) return EStrg::Xi;
  if (id==3334) return EStrg::Omega;
//=============================================================================

  return EStrg::Undef;
}

//_____________________________________________________________________________
class StrgInfo : public fastjet::PseudoJet::UserInfoBase {

 public :

  StrgInfo() = default;
  StrgInfo(EStrg e) : UserInfoBase(), fType(e) { }

  StrgInfo(const StrgInfo &src) : UserInfoBase(src), fType(src.fType) { }

  StrgInfo &operator=(const StrgInfo &src) {
    if(&src==this) return *this;

    fType = src.fType;
    return *this;
  }
//=============================================================================

  EStrg GetType() const { return fType; }
//=============================================================================

 private :

  EStrg fType = EStrg::Undef;
};
//=============================================================================
//=============================================================================

const auto d2PI(TMath::TwoPi());
const auto dPIm(-1.*TMath::Pi());
const auto dPI2(TMath::PiOver2());
Double_t PerpDeltaR(const TVector3 &vj, const TVector3 &vp)
{
  auto dp(-1.);
  TVector3 vt(vj);

  auto InverseDeltaR([&] {
    const auto dt(1. / vt.DeltaR(vp));
    if (dp<dt) dp = dt; });

  vt.RotateZ(dPI2);
  InverseDeltaR();

  vt.RotateZ(dPIm);
  InverseDeltaR();

  vt.SetPtEtaPhi(vj.Pt(), -1.*vj.Eta(), vj.Phi());
  vt.RotateZ(dPI2);
  InverseDeltaR();

  vt.RotateZ(dPIm);
  InverseDeltaR();

  return (1./dp);
}

//_____________________________________________________________________________
class ReducedDict {

 public :

  ReducedDict() = default;
//=============================================================================

  bool NoKey (const TString &s) { return (fMap.find(s)==fMap.end()); }
  bool HasKey(const TString &s) { return !(NoKey(s)); }

  void SetPair(const TString &s, const float d) {
    if (NoKey(s)) fMap.insert(std::pair<TString,float>(s,d));
  }

  float GetVal(const TString &s) { return fMap[s]; }

  void Reset() { fMap.clear(); }
  bool Empty() { return fMap.empty(); }
//=============================================================================

 private :

  std::map<TString, float> fMap;
};
//=============================================================================
//=============================================================================

void CallSumw2(TList *list)
{
  TObject *p(nullptr);
  TIter next(list);
  while ((p = next())) {
    auto h(dynamic_cast<TH1*>(p));
    if (h) h->Sumw2();
  }

  return;
}

//_____________________________________________________________________________
int RandomSeed()
{
  TDatime adt;
  if (!gSystem) new TSystem();
  const auto kTime((UInt_t)adt.Get());
  const auto kProc((UInt_t)gSystem->GetPid());

  const auto kInit(kTime - kProc);
  return ((int)kInit - (int)(((int)(kInit/1e8))*1e8));
}

//_____________________________________________________________________________
void DumpPyCfgs(const int mode, const TString &sf, Pythia8::Pythia &pythia)
{
  const TString ksp(Form("%s::%s",sf.Data(),__FUNCTION__));
//=============================================================================

  auto PrintNull([&ksp] { ::Info(ksp, "=============================="); });
//=============================================================================

  PrintNull();
  const auto bLC((mode==0) || (mode==2) || (mode==3));

  if (bLC) {
    ::Info(ksp.Data(), "Tune:pp = BLC mode %d", mode);
  } else {
    ::Info(ksp.Data(), "Tune:pp = Monash 2013");
  }
//=============================================================================

  auto PrintNone([&ksp] (const TString s) {
    ::Info(ksp, "%s = None", s.Data());
  });

  auto PrintFlag([&ksp,&pythia] (const TString s) {
    const auto b(pythia.flag(s.Data()));
    ::Info(ksp, "%s = %s", s.Data(), b ? "on" : "off");
  });

  auto PrintMode([&ksp,&pythia] (const TString s) {
    ::Info(ksp, "%s = %d", s.Data(), pythia.mode(s.Data()));
  });

  auto PrintParm([&ksp,&pythia] (const TString s) {
    ::Info(ksp, "%s = %f", s.Data(), pythia.parm(s.Data()));
  });

  auto PrintPvec([&ksp,&pythia] (const TString s) {
    const auto v(pythia.settings.pvec(s.Data()));
    std::cout << "Info in <" << ksp.Data() << ">: " << s.Data() << " =";
    for (const auto &e : v) std::cout << " " << e;
    std::cout << std::endl;
  });
//=============================================================================

  PrintNull();
  PrintParm("StringPT:sigma");
  PrintParm("StringZ:aLund");
  PrintParm("StringZ:bLund");
  PrintParm("StringFlav:probQQtoQ");
  PrintParm("StringFlav:ProbStoUD");
  PrintPvec("StringFlav:probQQ1toQQ0join");

  PrintNull();
  PrintParm("MultiPartonInteractions:pT0Ref");

  PrintNull();
  PrintMode("BeamRemnants:remnantMode");
  if (bLC) {
    PrintParm("BeamRemnants:saturation");
  } else {
    PrintNone("BeamRemnants:saturation");
  }

  PrintNull();
  PrintMode("ColourReconnection:mode");
  PrintFlag("ColourReconnection:allowDoubleJunRem");

  if (bLC) {
    PrintParm("ColourReconnection:m0");
    PrintFlag("ColourReconnection:allowJunctions");
    PrintParm("ColourReconnection:junctionCorrection");
    PrintMode("ColourReconnection:timeDilationMode");
    if (mode>0) {
      PrintParm("ColourReconnection:timeDilationPar");
    } else {
      PrintNone("ColourReconnection:timeDilationPar");
    }
  } else {
    PrintNone("ColourReconnection:m0");
    PrintNone("ColourReconnection:allowJunctions");
    PrintNone("ColourReconnection:junctionCorrection");
    PrintNone("ColourReconnection:timeDilationMode");
    PrintNone("ColourReconnection:timeDilationPar");
  }

  PrintNull();
  PrintNull();
//=============================================================================

  return;
}
//=============================================================================
//=============================================================================

#endif
