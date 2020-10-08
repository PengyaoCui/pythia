#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "fastjet/ClusterSequence.hh"

#include "Pythia8/Pythia.h"
#include "Pythia8/Analysis.h"
#include "Pythia8/ColourReconnection.h"
#include "Pythia8Plugins/ColourReconnectionHooks.h"

#include "TDatime.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

//=============================================================================

using namespace std;
using namespace Pythia8;
using namespace fastjet;
//=============================================================================

void PrintInfo();
bool CheckQCDFlag(TString);
bool CheckCRFlag (TString);
int  CheckCRMode (TString);

int GetRandomSeed();

float GetSphericity(vector<TVector3> tracks);

double GetGraphVal(TProfile *hStrange, TProfile *hPi, int i=0);
double GetGraphError(TProfile *hStrange, TProfile *hPi, int i=0);

TH1D *getLogXBin(TString histname, TString histaxis, int bin_num, double bin_min, double bin_max);

const TString srcName = "AnaPythiaRope";
//=============================================================================

int main(int argc, char* argv[])
{
  const TString sMethod = Form("%s::main", srcName.Data());
//=============================================================================

  TApplication theApp(srcName.Data(), &argc, argv);
  for (int i=0; i<argc; i++) ::Info(sMethod.Data(), "argv[%d] = %s", i, argv[i]);

//=============================================================================

  int kSeed = GetRandomSeed();
//=============================================================================

  Pythia8::Pythia pythia;
  Pythia8::Event& event = pythia.event;
  Pythia8::Info&  info  = pythia.info;

  pythia.readString("Main:numberOfEvents = 1000000");
  pythia.readString("Beams:eCM = 7000.");
  pythia.readString("SoftQCD:nonDiffractive = on");
  // Enabling setting of vertex information.
  pythia.readString("PartonVertex:setVertex = on");
  pythia.readString("PartonVertex:modeVertex = 2");
  pythia.readString("PartonVertex:ProtonRadius = 0.7");
  pythia.readString("PartonVertex:EmissionWidth = 0.1");

  // Enabling flavour ropes, setting model parameters.
  // The model is still untuned. These parameter values
  // are choosen for illustrative purposes.
  pythia.readString("Ropewalk:RopeHadronization = on");
  //pythia.readString("Ropewalk:doShoving = on");
  pythia.readString("Ropewalk:doFlavour = on");
  pythia.readString("Ropewalk:r0 = 0.5");
  pythia.readString("Ropewalk:m0 = 0.2");
  pythia.readString("Ropewalk:beta = 0.2");
  //pythia.readString("Ropewalk:setFixedKappa = on");
  //pythia.readString("Ropewalk:presetKappa = 2");


  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10.");
  pythia.readString("333:mayDecay = off");
  pythia.readString("3334:mayDecay = off");

  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",kSeed));

  pythia.init();

//=============================================================================

  //TFile *file = TFile::Open(Form("%s_%d_%s_%s.root",srcName.Data(),kSeed,sID.Data(),sQCD.Data()),"NEW"); TList *list = new TList();
  TList *list = new TList();

  TH1D     *hTrials = new     TH1D("hTrials", "", 1, 0., 1.); list->Add(hTrials);
  TProfile  *hXsect = new TProfile("hXsect",  "", 1, 0., 1.); list->Add(hXsect);

  TH1D *hKappa = new TH1D("hKappa", ";#kappa", 50, 0.5, 1.5); list->Add(hKappa);

  //event-wise variables
  TH1D *hMultInel = new TH1D("hMultInel", "charged multiplicity |eta|<1 for INEL>0;N_{ch};P(N_{ch})", 100, 0, 100); hMultInel->Sumw2(); list->Add(hMultInel);

  TH2D *hEvent = new TH2D("hEvent", "", 1000, 0., 1000., 2000, -0.5, 1999.5); hEvent->Sumw2(); list->Add(hEvent);

  TH2D *hFwdVsMid = new TH2D("hFwdVsMid", ";N_{trk}^{Mid} |#eta|<0.5; N_{trk}^{Fwd} V0 region", 500, -0.5, 499.5, 500, -0.5, 499.5); hFwdVsMid->Sumw2(); list->Add(hFwdVsMid);

  TH1D *hNtrkV0 = new TH1D("hNtrkV0", ";N_{trk}; N_{evt}", 500, -0.5, 499.5); hNtrkV0->Sumw2(); list->Add(hNtrkV0);
  
  TH2D *hSphericity = new TH2D("hSphericity", ";N_{ch}^{fwd}; S", 100, 0, 200, 50, 0, 1);
  hSphericity->Sumw2(); list->Add(hSphericity);
  
  TH2D *hSphericity_Soft = new TH2D("hSphericity_Soft", ";N_{ch}^{fwd}; S", 100, 0, 200, 50, 0, 1);
  hSphericity_Soft->Sumw2(); list->Add(hSphericity_Soft);
  
  TH2D *hSphericity_Hard = new TH2D("hSphericity_Hard", ";N_{ch}^{fwd}; S", 100, 0, 200, 50, 0, 1);
  hSphericity_Hard->Sumw2(); list->Add(hSphericity_Hard);
  
  TH3D *hMPINumSphNch = new TH3D("hMPINumSphNch", ";N_{ch}^{fwd};S_{T};nMPI", 100, 0, 200, 50, 0, 1, 50, -0.5, 49.5);
  hMPINumSphNch->Sumw2(); list->Add(hMPINumSphNch);
  
  //inclusive distribution
  TH1D *hChEta = new TH1D("hChEta", ";#eta; N_{trig}", 50, -5, 5); list->Add(hChEta);
  TH1D *hChpT = new TH1D("hChpT", ";p_{T} [GeV]; N_{ch}", 100, 0, 20); list->Add(hChpT);
  TH1D *hChEta_Inel = new TH1D("hChEta_Inel", ";#eta; N_{trig}", 50, -5, 5); list->Add(hChEta_Inel);
  TH1D *hChpT_Inel = new TH1D("hChpT_Inel", ";p_{T} [GeV]; N_{ch}", 100, 0, 20); list->Add(hChpT_Inel);
  TH1D *hKS_pT = new TH1D("hKS_pT", ";p_{T} [GeV]; N_{ch}", 100, 0, 20); list->Add(hKS_pT);
  TH1D *hLambda_pT = new TH1D("hLambda_pT", ";p_{T} [GeV]; N_{ch}", 100, 0, 20); list->Add(hLambda_pT);
  TH1D *hXi_pT = new TH1D("hXi_pT", ";p_{T} [GeV]; N_{ch}", 100, 0, 20); list->Add(hXi_pT);
  
  TH3D *hKCh = new TH3D("hKCh", ";p_{T} [GeV];N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hKCh->Sumw2(); list->Add(hKCh);
  
  TH3D *hPiCh = new TH3D("hPiCh", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hPiCh->Sumw2(); list->Add(hPiCh);
  
  TH3D *hKshort = new TH3D("hKshort", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hKshort->Sumw2(); list->Add(hKshort);
  
  TH3D *hLambda = new TH3D("hLambda", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hLambda->Sumw2(); list->Add(hLambda);
  
  TH3D *hProton = new TH3D("hProton", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hProton->Sumw2(); list->Add(hProton);
  
  TH3D *hPhi = new TH3D("hPhi", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hPhi->Sumw2(); list->Add(hPhi);
  
  TH3D *hOmega = new TH3D("hOmega", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hOmega->Sumw2(); list->Add(hOmega);
  
  TH3D *hCascade = new TH3D("hCascade", ";p_{T} [GeV] ;N_{ch}^{fwd}; S", 100, 0, 20, 500, -0.5, 499.5, 50, 0, 1);
  hCascade->Sumw2(); list->Add(hCascade);


  //double dNtrkEdge[11]={0 , 20, 31, 40, 53, 68, 79, 92, 112, 147 ,244 };
  //event class for inclusive
  //double dNtrkEdge[11]={0 , 20, 31, 40, 53, 68, 79, 92, 112, 147 ,400 };
  //double dNtrkEdge[11]={ 0 , 23, 36, 44, 55, 67, 76, 85, 100, 128 ,400 }; //mode 2 INEL>0
  double dNtrkEdge[11]={0 , 22, 35, 45, 57, 72, 84, 97, 116, 150 ,244 }; //mode 3 INEL>0
													
  TProfile*hdNchdEta = new TProfile("hdNchdEta", ";Event Class ; <dN_{ch}/d#eta>;", 11, 0.5, 11.5 );
  hdNchdEta->Sumw2(); list->Add(hdNchdEta);
  TProfile *hPiChClass = new TProfile("hPiChClass", ";Event Class ; <N^{#pi}>", 11, 0.5, 11.5);
  hPiChClass->Sumw2(); list->Add(hPiChClass);
  TProfile *hKshortClass = new TProfile("hKshortClass", ";Event Class ; <N^{K_{S}}>", 11, 0.5, 11.5);
  hKshortClass->Sumw2(); list->Add(hKshortClass);
  TProfile *hLambdaClass = new TProfile("hLambdaClass", ";Event Class ; <N^{#Lambda}>", 11, 0.5, 11.5);
  hLambdaClass->Sumw2(); list->Add(hLambdaClass);
  TProfile *hPhiClass = new TProfile("hPhiClass", ";Event Class ; <N^{#phi}>", 11, 0.5, 11.5);
  hPhiClass->Sumw2(); list->Add(hPhiClass);
  TProfile *hCascadeClass = new TProfile("hCascadeClass", ";Event Class ; <N^{#Xi}>", 11, 0.5, 11.5);
  hCascadeClass->Sumw2(); list->Add(hCascadeClass);
  TProfile *hOmegaClass = new TProfile("hOmegaClass", ";Event Class ; <N^{#Omega}>", 11, 0.5, 11.5);
  hOmegaClass->Sumw2(); list->Add(hOmegaClass);
  
  TH1D *hNLambdaOverNfwd = getLogXBin("hNLambdaOverNfwd", ";N_{fwd}; N(#Lambda)", 20, 1, 200);
  hNLambdaOverNfwd->Sumw2(); list->Add(hNLambdaOverNfwd);
  TH1D *hNXiOverNfwd = getLogXBin("hNXiOverNfwd", ";N_{fwd}; N(#Xi)", 20, 1, 200);
  hNXiOverNfwd->Sumw2(); list->Add(hNXiOverNfwd);
  TH1D *hNOmegaOverNfwd = getLogXBin("hNOmegaOverNfwd", ";N_{fwd}; N(#Omega)", 20, 1, 200);
  hNOmegaOverNfwd->Sumw2(); list->Add(hNOmegaOverNfwd);
  TH1D *hNKShortOverNfwd = getLogXBin("hNKShortOverNfwd", ";N_{fwd}; N(K_{S})", 20, 1, 200);
  hNKShortOverNfwd->Sumw2(); list->Add(hNKShortOverNfwd);
  TH1D *hNPiOverNfwd = getLogXBin("hNPiOverNfwd", ";N_{fwd}; N(#pi^{#pm})", 20, 1, 200);
  hNPiOverNfwd->Sumw2(); list->Add(hNPiOverNfwd);
  
  //jet related
  TH1D *hJetPt = new TH1D("hJetPt", ";p_{T}^{jet} [GeV]; counts", 100, 0, 50);
  hJetPt->Sumw2(); list->Add(hJetPt);
  
  TH1D *hJetPtLeading = new TH1D("hJetPtLeading", ";p_{T}^{jet} [GeV]; counts", 100, 0, 50);
  hJetPtLeading->Sumw2(); list->Add(hJetPtLeading);
  
  //charged pion with jet
  TH3D *hInJet_PiCh = new TH3D("hInJet_PiCh", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hInJet_PiCh->Sumw2(); list->Add(hInJet_PiCh);
  
  TH3D *hPerpJet_PiCh = new TH3D("hPerpJet_PiCh", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hPerpJet_PiCh->Sumw2(); list->Add(hPerpJet_PiCh);
  
  TH3D *hNoJet_PiCh = new TH3D("hNoJet_PiCh", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hNoJet_PiCh->Sumw2(); list->Add(hNoJet_PiCh);
  
  
  //Kshort with jet
  TH3D *hInJet_Kshort = new TH3D("hInJet_Kshort", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hInJet_Kshort->Sumw2(); list->Add(hInJet_Kshort);
  
  TH3D *hPerpJet_Kshort = new TH3D("hPerpJet_Kshort", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hPerpJet_Kshort->Sumw2(); list->Add(hPerpJet_Kshort);
  
  TH2D *hInJet_Kshort_test = new TH2D("hInJet_Kshort_test", ";p_{T}^{jet} [GeV]; N_{ch}^{fwd}", 
  		100, 0, 30, 100, 0, 200);
  hInJet_Kshort_test->Sumw2(); list->Add(hInJet_Kshort_test);
  
  TH3D *hNoJet_Kshort = new TH3D("hNoJet_Kshort", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hNoJet_Kshort->Sumw2(); list->Add(hNoJet_Kshort);
  
  //Lambda with jet
  TH3D *hInJet_Lambda = new TH3D("hInJet_Lambda", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hInJet_Lambda->Sumw2(); list->Add(hInJet_Lambda);
  
  TH3D *hPerpJet_Lambda = new TH3D("hPerpJet_Lambda", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hPerpJet_Lambda->Sumw2(); list->Add(hPerpJet_Lambda);
  
  TH2D *hInJet_Lambda_test = new TH2D("hInJet_Lambda_test", ";p_{T}^{jet} [GeV]; N_{ch}^{fwd}", 
  		100, 0, 20, 100, 0, 200);
  hInJet_Lambda_test->Sumw2(); list->Add(hInJet_Lambda_test);
  
  TH3D *hNoJet_Lambda = new TH3D("hNoJet_Lambda", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hNoJet_Lambda->Sumw2(); list->Add(hNoJet_Lambda);
  
  //Cascade with jet
  TH3D *hInJet_Cascade = new TH3D("hInJet_Cascade", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hInJet_Cascade->Sumw2(); list->Add(hInJet_Cascade);
  
  TH3D *hPerpJet_Cascade = new TH3D("hPerpJet_Cascade", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hPerpJet_Cascade->Sumw2(); list->Add(hPerpJet_Cascade);
  
  TH3D *hNoJet_Cascade = new TH3D("hNoJet_Cascade", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hNoJet_Cascade->Sumw2(); list->Add(hNoJet_Cascade);
  
  //Omega with jet
  TH3D *hInJet_Omega = new TH3D("hInJet_Omega", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hInJet_Omega->Sumw2(); list->Add(hInJet_Omega);
  
  TH3D *hPerpJet_Omega = new TH3D("hPerpJet_Omega", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hPerpJet_Omega->Sumw2(); list->Add(hPerpJet_Omega);
  
  TH3D *hNoJet_Omega = new TH3D("hNoJet_Omega", ";p_{T} [GeV]; N_{ch}^{fwd}; S", 
  		100, 0, 20, 100, 0, 200, 50, 0, 1);
  hNoJet_Omega->Sumw2(); list->Add(hNoJet_Omega);


//=============================================================================

  //cuts
  //alice cut:
  //INEL>0 |eta|<1
  //mid-rap ch density |eta|<0.5
  //inclusive strange |y|<0.5
  //cms cut:
  //inclusive strange |y|<2
  const double dcEtaCut  = 1;
  const double dfEtaMin  = 2.;
  const double dfEtaMax  = 5.;
  const double dPtCut    = 0.15;

  const double dPtTrMax = 3;
  const double dPtTrMin = 1;
  const double dPtAsMax = 3;
  const double dPtAsMin = 1;
  
  const double dHighMult = 80;
  const double dLowMult = 20;

  //temperary containers
  TVector3 vec_tmp;
  TVector3 vec_LeadJet;
  vector<TVector3> sph_container;
  
  //jet related
  vector<PseudoJet> jetParticles;
  double R=0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  cout << "Clustering with " << jet_def.description() << endl;

//=============================================================================
  //event loop starts
  for (int iEvent=0; iEvent<pythia.mode("Main:numberOfEvents"); iEvent++) if (pythia.next()) { 
    //event wise initializations
    double dFwdCh = 0.;
    double dMidCh = 0.;
    double dNch05 = 0.;
    double dNtrkV0 = 0.;
    sph_container.clear();
    jetParticles.clear();

//=============================================================================
    //determine charged track number in mid/fwd rapidity in this event 
    for (int i=0; i<event.size(); i++) if (event[i].isFinal()   &&
                                           event[i].isVisible() &&
//                                           event[i].isCharged() &&
                                          (event[i].pT()>dPtCut)) {
      if( ((event[i].eta()>2.8&&event[i].eta()<5.1) || 
			      (event[i].eta()>-3.7&&event[i].eta()<-1.7)) 
		      && event[i].isCharged()) dNtrkV0++;

      double dEtaAbs = TMath::Abs(event[i].eta());
      double dRapAbs = TMath::Abs(event[i].y());
      vec_tmp.SetXYZ(event[i].px(), event[i].py(), event[i].pz());
  
      if( event[i].isCharged() ){
        hChEta->Fill(event[i].eta());
      
	if ((dEtaAbs>dfEtaMin) && (dEtaAbs<=dfEtaMax)) {
	  dFwdCh += 1.; 
	  //sph_container.push_back(vec_tmp); //too small sphericity
	} 

	if( dEtaAbs<dcEtaCut ) {
	  hChpT->Fill(event[i].pT(), 1./event[i].pT());
	  dMidCh++;
	  jetParticles.push_back( PseudoJet(event[i].px(), event[i].py(), event[i].pz(), event[i].e()) );
	  sph_container.push_back(vec_tmp);
	  if( dEtaAbs<0.5 ) dNch05++;
       	} 
      
      }//fi 
   
    } 
   
    if(dMidCh<1) continue; 
    hKappa->Fill(pythia.getKappa()); 
    hNtrkV0->Fill(dNtrkV0); 
    hMultInel->Fill(dMidCh); 
    //hFwdVsMid->Fill(dMidCh, dFwdCh);
    hFwdVsMid->Fill(dNch05, dNtrkV0);
    //sphericity analysis
    double dSph = GetSphericity(sph_container);
    if(dSph>1.||dSph<0.) cout<<"dsph="<<dSph<<" dMidCh="<<dMidCh<<endl;
    if(dMidCh>0) 
      hSphericity->Fill(dFwdCh, dSph); 
      hMPINumSphNch->Fill(dFwdCh, dSph, info.nMPI()); 
    
    
    //jet analysis
    //run the clustering, extract the jets
    ClusterSequence cs(jetParticles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());	
    vector<TVector3> jets_selected;
    bool bTagSoft = true;
    if(jets.size()>0){
      hJetPtLeading->Fill(jets[0].pt());
      vec_LeadJet.SetXYZ(jets[0].px(), jets[0].py(), jets[0].pz());
    } 
    for(int i=0; i<jets.size(); i++){
      if(fabs(jets[i].eta())<(dcEtaCut-R)){
          hJetPt->Fill(jets[i].pt());
	  if(jets[i].pt()>5){
	    bTagSoft=false;
	    if(jets[i].pt()>10){ 
      	      vec_tmp.SetXYZ(jets[i].px(), jets[i].py(), jets[i].pz());
	      jets_selected.push_back(vec_tmp);
	    }
	  }
      }
    }//for jets 
    if(bTagSoft) 
      hSphericity_Soft->Fill(dFwdCh, dSph);
    else if(jets_selected.size()>0) 
      hSphericity_Hard->Fill(dFwdCh, dSph);
	

//=============================================================================
    double dNPiCh=0;
    double dNKshort=0;
    double dNLambda=0;
    double dNCascade=0;
    double dNOmega=0;
    double dNPhi=0;

    //connect the inclusive particle properties with the mid/fwd charged
    //track number
    for (int i=0; i<event.size(); i++) if (event[i].isFinal() &&
                                           event[i].isVisible() &&
					   event[i].pT()>dPtCut) {

      double dEtaAbs = TMath::Abs(event[i].eta()); 
      double dRapAbs = TMath::Abs(event[i].y());
      if(dRapAbs<2){
        if(event[i].idAbs()==310) hKS_pT->Fill(event[i].pT()); 
	if(event[i].idAbs()==3122) hLambda_pT->Fill(event[i].pT()); 
	if(event[i].idAbs()==3312) hXi_pT->Fill(event[i].pT());
      } 
      if(event[i].isCharged()) hChEta_Inel->Fill(event[i].eta()); 
      if (dEtaAbs>dcEtaCut) continue; 
      if(event[i].isCharged()){
        hChpT_Inel->Fill(event[i].pT(), 1./event[i].pT());
      }
      
      int id  = event[i].idAbs();
      bool bKshort = (id==310);
      bool bKCh = (id==321);
      bool bLambda = (id==3122);
      bool bProton = (id==2212);
      bool bPiCh = (id==211);
      bool bPhi = (id==333);
      bool bOmega = (id==3334);
      bool bCascade = (id==3312);
      
      if ( (!bKshort) && (!bKCh) && (!bLambda) && (!bProton) && (!bPiCh) && (!bPhi) && (!bOmega) && (!bCascade)	) continue;

      double dPt = event[i].pT();
      vec_tmp.SetXYZ(event[i].px(), event[i].py(), event[i].pz());
      if (bKshort) { 
	dNKshort++;
	//hKshort->Fill(dPt, dFwdCh,dSph);
	if(dRapAbs<0.5) hKshort->Fill(dPt, dNtrkV0, dSph);
       	if(jets_selected.size()>0 && fabs(fabs(TVector2::Phi_mpi_pi(vec_tmp.DeltaPhi(vec_LeadJet)))-pi/2.)<pi/4.){
       	  hPerpJet_Kshort->Fill(dPt, dFwdCh, dSph);
       	} 
	for(int j_jets=0; j_jets<jets_selected.size(); j_jets++){
	  if(jets_selected.at(j_jets).DeltaR(vec_tmp)<R) {
	    hInJet_Kshort->Fill(dPt, dFwdCh, dSph);
	    hInJet_Kshort_test->Fill(jets_selected.at(j_jets).Pt(), dFwdCh);
	  }
       	}//for jets
       	if(bTagSoft){
	  hNoJet_Kshort->Fill(dPt, dFwdCh, dSph);
       	}
      }

      if (bKCh) { 
	//hKCh->Fill(dPt, dFwdCh,dSph);
	if(dRapAbs<0.5) hKCh->Fill(dPt, dNtrkV0, dSph);
      }

      if (bLambda) { 
	dNLambda++;
       	//hLambda->Fill(dPt, dFwdCh,dSph);
	if(dRapAbs<0.5) hLambda->Fill(dPt, dNtrkV0, dSph);
       	if(jets_selected.size()>0&&fabs(fabs(TVector2::Phi_mpi_pi(vec_tmp.DeltaPhi(vec_LeadJet)))-pi/2.)<pi/4.){
	  hPerpJet_Lambda->Fill(dPt, dFwdCh, dSph);
       	} 
	for(int j_jets=0; j_jets<jets_selected.size(); j_jets++){
	  if(jets_selected.at(j_jets).DeltaR(vec_tmp)<R) {
	    hInJet_Lambda->Fill(dPt, dFwdCh, dSph);
	    hInJet_Lambda_test->Fill(jets_selected.at(j_jets).Pt(), dFwdCh);
	  }
       	}//for jets
       	if(bTagSoft){
	  hNoJet_Lambda->Fill(dPt, dFwdCh, dSph);
       	}
      }

      if (bProton) { 
        //hProton->Fill(dPt, dFwdCh, dSph);
	if(dRapAbs<0.5) hProton->Fill(dPt, dNtrkV0, dSph);
      }

      if (bPiCh) { 
        dNPiCh++;
       	//hPiCh->Fill(dPt, dFwdCh,dSph);
	if(dRapAbs<0.5) hPiCh->Fill(dPt, dNtrkV0, dSph);
       	if(jets_selected.size()>0&&fabs(fabs(TVector2::Phi_mpi_pi(vec_tmp.DeltaPhi(vec_LeadJet)))-pi/2.)<pi/4.){
       	  hPerpJet_PiCh->Fill(dPt, dFwdCh, dSph);
       	} 
	for(int j_jets=0; j_jets<jets_selected.size(); j_jets++){
	  if(jets_selected.at(j_jets).DeltaR(vec_tmp)<R) {
	    hInJet_PiCh->Fill(dPt, dFwdCh, dSph);
	  }
       	}//for jets
       	if(bTagSoft){
	  hNoJet_PiCh->Fill(dPt, dFwdCh, dSph);
       	}
      }

      if (bPhi) { 
	dNPhi++;
       	//hPhi->Fill(dPt, dFwdCh, dSph);
	if(dRapAbs<0.5) hPhi->Fill(dPt, dNtrkV0, dSph);
      }
     
      if (bOmega) {
        dNOmega++;
       	//hOmega->Fill(dPt, dFwdCh, dSph);
	if(dRapAbs<0.5) hOmega->Fill(dPt, dNtrkV0, dSph);
       	if(jets_selected.size()>0&&fabs(fabs(TVector2::Phi_mpi_pi(vec_tmp.DeltaPhi(vec_LeadJet)))-pi/2.)<pi/4.){
	  hPerpJet_Omega->Fill(dPt, dFwdCh, dSph);
       	} 
	for(int j_jets=0; j_jets<jets_selected.size(); j_jets++){
	  if(jets_selected.at(j_jets).DeltaR(vec_tmp)<R) {
	    hInJet_Omega->Fill(dPt, dFwdCh, dSph);
	  }
       	}//for jets
       	if(bTagSoft){
	  hNoJet_Omega->Fill(dPt, dFwdCh, dSph);
       	}
      } 
      if (bCascade) {
        dNCascade++;
       	//hCascade->Fill(dPt, dFwdCh, dSph);
	if(dRapAbs<0.5) hCascade->Fill(dPt, dNtrkV0, dSph);
       	if(jets_selected.size()>0&&fabs(fabs(TVector2::Phi_mpi_pi(vec_tmp.DeltaPhi(vec_LeadJet)))-pi/2.)<pi/4.){
	  hPerpJet_Cascade->Fill(dPt, dFwdCh, dSph);
       	} 
	
	for(int j_jets=0; j_jets<jets_selected.size(); j_jets++){
	  if(jets_selected.at(j_jets).DeltaR(vec_tmp)<R) {
	    hInJet_Cascade->Fill(dPt, dFwdCh, dSph);
	  }
       	}//for jets
       	if(bTagSoft){
	  hNoJet_Cascade->Fill(dPt, dFwdCh, dSph);
       	}
      } 
    }//track loop end 
   
    for(int i_event(0); i_event<11; i_event++){
      double LowerEdge=-1E5;
      double UpperEdge=1E5;
      if(i_event<10){
        LowerEdge=dNtrkEdge[i_event];
       	UpperEdge=dNtrkEdge[i_event+1];
      }
      if(dNtrkV0>LowerEdge&&dNtrkV0<UpperEdge){
	hdNchdEta->Fill(i_event+1, dMidCh);
	hPiChClass->Fill(i_event+1, dNPiCh);
	hKshortClass->Fill(i_event+1, dNKshort);
	hLambdaClass->Fill(i_event+1, dNLambda);
	hPhiClass->Fill(i_event+1, dNPhi);
	hCascadeClass->Fill(i_event+1, dNCascade);
	hOmegaClass->Fill(i_event+1, dNOmega);
      }
    } 
    hNKShortOverNfwd->Fill(dFwdCh, dNKshort);
    hNLambdaOverNfwd->Fill(dFwdCh, dNLambda);
    hNXiOverNfwd->Fill(dFwdCh, dNCascade);
    hNOmegaOverNfwd->Fill(dFwdCh, dNOmega);
    hNPiOverNfwd->Fill(dFwdCh, dNPiCh);

    hEvent->Fill(info.pTHat(), dFwdCh);
  }
//=============================================================================

	
  hXsect->Fill(0.5, info.sigmaGen());
  hTrials->Fill(0.5, info.weightSum());

  //get normalized TGraph for strangeness
  TGraphErrors *grKshortClass = new TGraphErrors();
  TGraphErrors *grLambdaClass = new TGraphErrors();
  TGraphErrors *grCascadeClass = new TGraphErrors();
  TGraphErrors *grOmegaClass = new TGraphErrors();
  TGraphErrors *grPhiClass = new TGraphErrors();
  grKshortClass->SetName("grKshortClass");
  list->Add(grKshortClass);
  grLambdaClass->SetName("grLambdaClass");
  list->Add(grLambdaClass);
  grPhiClass->SetName("grPhiClass");
  list->Add(grPhiClass);
  grCascadeClass->SetName("grCascadeClass");
  list->Add(grCascadeClass);
  grOmegaClass->SetName("grOmegaClass");
  list->Add(grOmegaClass);
  
  cout<<"now ready to get graph"<<endl;
  for(int i_point(0); i_point<10; i_point++){
    grKshortClass->SetPoint(i_point, hdNchdEta->GetBinContent(i_point+1),
    					GetGraphVal(hKshortClass, hPiChClass, i_point+1) );
    grKshortClass->SetPointError(i_point, 0., GetGraphError(hKshortClass, hPiChClass, i_point+1) );
    
    grOmegaClass->SetPoint(i_point, hdNchdEta->GetBinContent(i_point+1),
    					GetGraphVal(hOmegaClass, hPiChClass, i_point+1) );
    grOmegaClass->SetPointError(i_point, 0., GetGraphError(hOmegaClass, hPiChClass, i_point+1) );
    
    
    grCascadeClass->SetPoint(i_point, hdNchdEta->GetBinContent(i_point+1),
    					GetGraphVal(hCascadeClass, hPiChClass, i_point+1) );
    grCascadeClass->SetPointError(i_point, 0., GetGraphError(hCascadeClass, hPiChClass, i_point+1) );
    
    grPhiClass->SetPoint(i_point, hdNchdEta->GetBinContent(i_point+1),
    					GetGraphVal(hPhiClass, hPiChClass, i_point+1) );
    grPhiClass->SetPointError(i_point, 0., GetGraphError(hPhiClass, hPiChClass, i_point+1) );
    
    grLambdaClass->SetPoint(i_point, hdNchdEta->GetBinContent(i_point+1),
    					GetGraphVal(hLambdaClass, hPiChClass, i_point+1) );
    grLambdaClass->SetPointError(i_point, 0., GetGraphError(hLambdaClass, hPiChClass, i_point+1) );
  }



  //TFile *file = new TFile(Form("%s.root",srcName.Data()),"recreate"); 
  TFile *file = new TFile(Form("%s_default.root",srcName.Data()),"recreate"); 
  file->cd(); 
  list->Write(); 

  file->Close();
//=============================================================================

  ::Info(sMethod.Data(), "DONE");
  return 0;
}

//_____________________________________________________________________________
int GetRandomSeed()
{
  const TString sMethod = Form("%s::GetRandomSeed", srcName.Data());
//=============================================================================

  TDatime adt;
  new TSystem();
  UInt_t kTime = (UInt_t)adt.Get();
  UInt_t kProc = (UInt_t)gSystem->GetPid();//Get process id.
  UInt_t kInit = kTime - kProc*10;
  int kSeed = (int)kInit - (int)(((int)(kInit/1e8))*1e8);

  ::Info(sMethod.Data(), "Proc ID     = %d", kProc);
  ::Info(sMethod.Data(), "System time = %d", kTime);
  ::Info(sMethod.Data(), "Random number seed = %d", kSeed);

  return kSeed;
}

//_____________________________________________________________________________
bool CheckQCDFlag(TString s)
{
  const TString sMethod = Form("%s::CheckQCDFlag", srcName.Data());
//=============================================================================

  if (s.IsNull()) {
    ::Fatal(sMethod.Data(), "QCD processes has not yet been set!!");
  } else {
    if ((s!="Soft") && (s!="Hard")) {
      ::Fatal(sMethod.Data(), "QCD processe (= %s) is neither \'Soft\' nor \'Hard\'!!", s.Data());
      return true;
    } else {
      if (s=="Soft") ::Info(sMethod.Data(), "SoftQCD:nonDiffractive = on");
      if (s=="Hard") ::Info(sMethod.Data(), "HardQCD:all = on");
    }
  }

  return false;
}

//_____________________________________________________________________________
bool CheckCRFlag(TString s)
{
  const TString sMethod = Form("%s::CheckCRFlag", srcName.Data());
//=============================================================================

	if (s!="ON"&&s!="OFF") {
		::Warning(sMethod.Data(), "CR flag: %s not defined, end in CheckCRFlag!", s.Data());
		exit(-1);
	}

  ::Info(sMethod.Data(), "CR status: %s", s.Data());

	if (s=="ON") {
		::Info(sMethod.Data(), "ColourReconnection:reconnect = on");
		return true;
	} else {
		::Info(sMethod.Data(), "ColourReconnection:reconnect = off");
		return false;
	}

}

//_____________________________________________________________________________
int CheckCRMode(TString s)
{
  const TString sMethod = Form("%s::CheckCRMode", srcName.Data());
//=============================================================================

  int kMode = 1;
  if (s.IsNull()) {
    ::Warning(sMethod.Data(), "CR mode has not yet been set, using defaul (= 1)");
  } else {
    kMode = s.Atoi();
    if ((kMode<0) || (kMode>4)) {
      ::Fatal(sMethod.Data(), "CR mode (= %d) is invalidated!!", kMode);
    } else {
      ::Info(sMethod.Data(), "ColourReconnection:mode = %d", kMode);
    }
  }

  return kMode;
}

//_____________________________________________________________________________
void PrintInfo()
{
  const TString sMethod = Form("%s::PrintInfo", srcName.Data());
//=============================================================================

  ::Info(sMethod.Data(), "./%s.exe arg1 [arg2] [arg3]", srcName.Data());

  ::Info(sMethod.Data(), "arg1: QCD processes");
  ::Info(sMethod.Data(), "  = Soft: SoftQCD:nonDiffractive  = on");
  ::Info(sMethod.Data(), "  = Hard: HardQCD:all = on -- pT hat setting is required");
  ::Info(sMethod.Data(), "arg2: optional, CR flag (default = off)");
  ::Info(sMethod.Data(), "  = wCR: ColourReconnection:reconnect = on");
  ::Info(sMethod.Data(), "arg3: optional, CR mode (default = 1)");
  ::Info(sMethod.Data(), "  =  0: The MPI-based original Pythia 8 scheme");
  ::Info(sMethod.Data(), "  =  1: The new more QCD based scheme");
  ::Info(sMethod.Data(), "  =  2: The new gluon-move model");
  ::Info(sMethod.Data(), "  =  3: The SK  I e^{+} e^{-} CR model");
  ::Info(sMethod.Data(), "  =  4: The SK II e^{+} e^{-} CR model");

  return;
}

float GetSphericity(vector<TVector3> tracks)
{
  // a critical check
  if(tracks.size()<=1) return 0; 
  double sumPt = 0;
  // first fill the momentum tensor
  // TMatrixDSym MomentumTensor(3);
  TMatrixDSym MomentumTensor(2);
  for(vector<TVector3>::iterator itTrack = tracks.begin(); itTrack!=tracks.end(); ++itTrack) {
  //std::vector<double> momentum(3);
  std::vector<double> momentum(2);
  momentum[0] = itTrack->X();
  momentum[1] = itTrack->Y();
  double trackPt = itTrack->Pt();
  sumPt+=trackPt;
  //momentum[2] = itTrack->Z();
  //for(unsigned int i=0;i<3;i++)
  for(unsigned int i=0;i<2;i++)
    for(unsigned int j=0;j<=i;j++) {
      MomentumTensor[i][j] += momentum[i]*momentum[j]/trackPt;
    }
  }
  //MomentumTensor*=1/(MomentumTensor[0][0]+MomentumTensor[1][1]+MomentumTensor[2][2]);
  MomentumTensor*=1/(sumPt);
  // find the eigen values
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();
  //vector<float> eigenvaluess(3);
  vector<float> eigenvaluess(2);
  eigenvaluess[0] = eigenvals[0];
  eigenvaluess[1] = eigenvals[1];
  //eigenvaluess[2] = eigenvals[2];
  //sort the eigen value from low to high
  sort(eigenvaluess.begin(),eigenvaluess.end());
  //compute spericity
  //float sph = ( 1.5*(1-eigenvaluess[2]));
  float sph = 2*eigenvaluess[0]/(eigenvaluess[0]+eigenvaluess[1]);
  return sph;
}

double GetGraphVal(TProfile *hStrange, TProfile *hPi, int i){
	double content = -9;
	if(hStrange->GetBinContent(11)<1E-6||hPi->GetBinContent(i)<1E-6) return content;
	content = hStrange->GetBinContent(i)/hPi->GetBinContent(i)/(hStrange->GetBinContent(11)/hPi->GetBinContent(11));
	return content; 
}

double GetGraphError(TProfile *hStrange, TProfile *hPi, int i){
	double content = GetGraphVal(hStrange, hPi, i);
	if(content<0) return -9;
	double err_str = hStrange->GetBinError(i)/hStrange->GetBinContent(i);
	double err_pi = hPi->GetBinError(i)/hStrange->GetBinContent(i);	
	double err_str_tot = hStrange->GetBinError(11)/hStrange->GetBinContent(11);
	double err_pi_tot = hPi->GetBinError(11)/hPi->GetBinContent(11);
	return content*sqrt(err_str*err_str+ err_pi*err_pi + err_str_tot*err_str_tot + err_pi_tot*err_pi_tot);
}

TH1D *getLogXBin(TString histname, TString histaxis, int bin_num, double bin_min, double bin_max){
	double x_edge[bin_num+1];
       	for(int i(0); i<=bin_num; i++){
	      double step = (log10(bin_max)-log10(bin_min))/bin_num;
	      x_edge[i]=pow(10, log10(bin_min)+i*step);
       	}
       	TH1D *htemp = new TH1D(histname, histaxis, bin_num, x_edge);
       	return htemp;
}
