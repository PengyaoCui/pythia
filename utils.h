#include "TSystem.h"
#include "TDatime.h"
#include "TMath.h"
#include "TVector3.h"
//=============================================================================

enum { kR02=0, kR03, kR04 };
enum { kC02=0, kC03, kC04 };
enum { kKshortIn=0, kLambdaDr, kLambdaFd, kAntiLaDr, kAntiLaFd };
const double dCutArea[] = { 0.6*TMath::Pi()*0.2*0.2, 0.6*TMath::Pi()*0.3*0.3, 0.6*TMath::Pi()*0.4*0.4 };
//=============================================================================

class RecoInfo : public fastjet::PseudoJet::UserInfoBase {

 public :

  RecoInfo() : UserInfoBase(), fType(-1) {
    for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) fIsMatch[j][i] = false;
  }

  RecoInfo(int l) : UserInfoBase(), fType(l) {
    for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) fIsMatch[j][i] = false;
  }

  RecoInfo(const RecoInfo &src) : UserInfoBase(src), fType(src.fType) {
    for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) fIsMatch[j][i] = src.fIsMatch[j][i];
  }

  RecoInfo &operator=(const RecoInfo &src) {
    if(&src==this) return *this;

    fType = src.fType;
    for (int j=0; j<3; ++j) for (int i=0; i<3; ++i) fIsMatch[j][i] = src.fIsMatch[j][i];

    return *this;
  }
//=============================================================================

  int GetType() const { return fType; }

  bool GetIsMatch(int j, int i) const {
    if ((j<0) || (j>=3) || (i<0) || (i>=3)) return false;
    return fIsMatch[j][i];
  }
//=============================================================================

  void SetIsMatch(fastjet::PseudoJet &aRec, const std::vector<fastjet::PseudoJet> &aJets, int lj) {
    if ((lj<0) || (lj>=3)) return;
    TVector3 vRec; vRec.SetPtEtaPhi(aRec.pt(), aRec.eta(), aRec.phi());

    TVector3 vJet;
    for (unsigned long j=0; j<aJets.size(); ++j) if ((aJets[j].area()>dCutArea[lj]) &&
                                                     (aJets[j].pt()>10.) &&
                                                     (!(fIsMatch[lj][0]  &&
                                                        fIsMatch[lj][1]  &&
                                                        fIsMatch[lj][2]))) {
      vJet.SetPtEtaPhi(aJets[j].pt(), aJets[j].eta(), aJets[j].phi());

      double d = vJet.DeltaR(vRec);
      if (d<0.2) fIsMatch[lj][0] = true;
      if (d<0.3) fIsMatch[lj][1] = true;
      if (d<0.4) fIsMatch[lj][2] = true;
    }

    return;
  }
//=============================================================================

 private :

  int  fType;
  bool fIsMatch[3][3];
};
//=============================================================================

int GetRandomSeed()
{
  TDatime adt;
  if(!gSystem) new TSystem();

  UInt_t lTime = (UInt_t)adt.Get();
  UInt_t lProc = (UInt_t)gSystem->GetPid();
  UInt_t lInit = lTime - lProc;
  int lSeed = (int)lInit-(int)(((int)(lInit/1e8))*1e8);
  return lSeed;
}
