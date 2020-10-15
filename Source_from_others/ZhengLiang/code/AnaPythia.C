#include <iostream>
#include <vector>

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
#include "TProfile.h"
#include "TVector3.h"
//=============================================================================

using namespace std;
using namespace Pythia8;
//=============================================================================

void PrintInfo();
bool CheckQCDFlag(TString);
bool CheckCRFlag (TString);
int  CheckCRMode (TString);

int GetRandomSeed();//得到种子
bool findIn(int, int, int, std::vector<int>);

const TString srcName = "AnaPythiaCR";
//=============================================================================

int main(int argc, char* argv[])
{
  const TString sMethod = Form("%s::main", srcName.Data());
//=============================================================================

  TApplication theApp(srcName.Data(), &argc, argv);
  for (int i=0; i<argc; i++) ::Info(sMethod.Data(), "argv[%d] = %s", i, argv[i]);

  if (argc<3) {
	  ::Error(sMethod.Data(), "Number of arguments (= %d) < required value", argc);
	  PrintInfo();
	  return -1;
  }
//=============================================================================

  const TString sID  = argv[1];
  const TString sQCD = argv[2]; CheckQCDFlag(sQCD);
//=============================================================================

  bool bCRf = CheckCRFlag(argv[3]);
  int  kCRm = (bCRf ? CheckCRMode(argv[4]) : 1);
  int kSeed = GetRandomSeed();
//=============================================================================

  Pythia8::Pythia pythia;//定义一个pythia
  Pythia8::Event& event = pythia.event;//固定的写法
  Pythia8::Info&  info  = pythia.info;//固定的写法

  pythia.readString("Main:numberOfEvents = 10000");//对一组参数做了设定

  //mode  Tune:pp (default = 5; minimum = -1; maximum = 14)
  //pythia.readString("Tune:pp = 5"); //tune 4C
  //pythia.readString("Tune:pp = 14"); //Monash 2013 tune
  pythia.readString("Tune:pp = 18"); //CMS Monashstar

  pythia.readString("Beams:eCM = 13000.");//13TeV

  if (sQCD=="Soft") {
	  pythia.readString("SoftQCD:nonDiffractive = on");//对应实验上的 minimum-bias trigger
  }

  if (sQCD=="Hard") {
	  pythia.readString("HardQCD:all = on");//Common switch for the group of all hard QCD processes
	  pythia.readFile(Form("%s.cmnd",srcName.Data()));//????
  }

  pythia.readString("Next:numberShowInfo = 0");//The number of events to list the Info information for, where relevant. default = 1
  pythia.readString("Next:numberShowProcess = 0");//The number of events to list the process record for, where relevant. default = 1
  pythia.readString("Next:numberShowEvent = 0");//The number of events to list the event record for, where relevant. default = 1 

  pythia.readString("ParticleDecays:limitTau0 = on");//only particles with tau0 < tau0Max are decayed.
  pythia.readString("ParticleDecays:tau0Max = 10.");//only particles with tau0 < tau0Max are decayed.
  pythia.readString("333:mayDecay = off");// \phi will not decay
  pythia.readString("3334:mayDecay = off");// \Omega^- will not decay

  pythia.readString("Random:setSeed = on");//用自己定义的种子
  pythia.readString(Form("Random:seed = %d",kSeed));

  if (bCRf) {
	  pythia.readString(Form("ColourReconnection:mode  = %d",kCRm));//0-4
	  pythia.readString(Form("BeamRemnants:remnantMode = %d",kCRm));// if ColourReconnection:mode == 1, then must BeamRemnants:remnantMode ==1
  } else {
	  pythia.readString("ColourReconnection:reconnect = off");
  }

  pythia.init();

  if (bCRf) {
	  ::Info(sMethod.Data(), "ColourReconnection:mode  = %d", pythia.mode("ColourReconnection:mode"));
	  ::Info(sMethod.Data(), "BeamRemnants:remnantMode = %d", pythia.mode("BeamRemnants:remnantMode"));
  }
//=============================================================================
  //Define out put file
  TFile *file = TFile::Open(Form("%s_%d_%s_%s.root",srcName.Data(),kSeed,sID.Data(),sQCD.Data()),"NEW"); TList *list = new TList();
  
  //Define Histograms....
  TH1D *hTrials    = new TH1D("hTrials", "", 1, 0., 1.); list->Add(hTrials);
  TProfile *hXsect = new TProfile("hXsect",  "", 1, 0., 1.); list->Add(hXsect);

  TH2D *hEvent = new TH2D("hEvent", "", 1000, 0., 1000., 2000, -0.5, 1999.5); hEvent->Sumw2(); list->Add(hEvent);

  TH2D *hKshort    = new TH2D("hKshort", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hKshort->Sumw2(); list->Add(hKshort);
  TH2D *hKCh       = new TH2D("hKCh", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hKCh->Sumw2(); list->Add(hKCh);
  TH2D *hLambda    = new TH2D("hLambda", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hLambda->Sumw2(); list->Add(hLambda);
  TH2D *hProton    = new TH2D("hProton", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hProton->Sumw2(); list->Add(hProton);
  TH2D *hPiCh      = new TH2D("hPiCh", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hPiCh->Sumw2(); list->Add(hPiCh);
  TH2D *hResPhi    = new TH2D("hResPhi", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hResPhi->Sumw2(); list->Add(hResPhi);
  TH2D *hLambdaMid = new TH2D("hLambdaMid", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hLambdaMid->Sumw2(); list->Add(hLambdaMid);
  TH2D *hKshortMid = new TH2D("hKshortMid", "N_{trk} vs p_{T}; p_{T} [GeV]; N_{trk}", 200, 0., 20., 2000, -0.5, 1999.5); hKshortMid->Sumw2(); list->Add(hKshortMid);

  TH2D *hFwdVsMid    = new TH2D("hFwdVsMid", ";N_{trk}^{Mid}; N_{trk}^{Fwd}", 2000, -0.5, 1999.5, 2000, -0.5, 1999.5); hFwdVsMid->Sumw2(); list->Add(hFwdVsMid);
  TH2D *hFwdVsMidTag = new TH2D("hFwdVsMidTag", ";N_{trk}^{Mid}; N_{trk}^{Fwd}", 2000, -0.5, 1999.5, 2000, -0.5, 1999.5); hFwdVsMidTag->Sumw2(); list->Add(hFwdVsMidTag);

  TProfile *hProfKshort = new TProfile("hProfKshort", ";N_{trk}; <p_{T}> [GeV]", 2000, -0.5, 1999.5); list->Add(hProfKshort);
  TProfile *hProfKCh    = new TProfile("hProfKCh", ";N_{trk}; <p_{T}> [GeV]", 2000, -0.5, 1999.5); list->Add(hProfKCh);
  TProfile *hProfLambda = new TProfile("hProfLambda", ";N_{trk}; <p_{T}> [GeV]", 2000, -0.5, 1999.5); list->Add(hProfLambda);
  TProfile *hProfProton = new TProfile("hProfProton", ";N_{trk}; <p_{T}> [GeV]", 2000, -0.5, 1999.5); list->Add(hProfProton);
  TProfile *hProfPiCh   = new TProfile("hProfPiCh", ";N_{trk}; <p_{T}> [GeV]", 2000, -0.5, 1999.5); list->Add(hProfPiCh);
  TProfile *hProfResPhi = new TProfile("hProfResPhi", ";N_{trk}; <p_{T}> [GeV]", 2000, -0.5, 1999.5); list->Add(hProfResPhi);

  TH1D *hTrEta = new TH1D("hTrEta", ";#eta; N_{trig}", 50, -2, 2); list->Add(hTrEta);
  TH1D *hTrPt  = new TH1D("hTrPt", ";p_{T}; N_{trig}", 50, 0, 10); list->Add(hTrPt);
  TH1D *hAsPt  = new TH1D("hAsPt", ";p_{T}; N_{assc}", 50, 0, 10); list->Add(hAsPt);

  const double pi = 3.14159;// why not TMath::pi()
  TH2D *hDEtaDPhiSameEvent = new TH2D("hDEtaDPhiSameEvent", "; #Delta#phi; #Delta#eta", 40, -0.5*pi, 1.5*pi, 40, -4, 4); hDEtaDPhiSameEvent->Sumw2(); list->Add(hDEtaDPhiSameEvent);
  TH2D *hDEtaDPhiMixEvent  = new TH2D("hDEtaDPhiMixEvent", "; #Delta#phi; #Delta#eta", 40, -0.5*pi, 1.5*pi, 40, -4, 4); hDEtaDPhiMixEvent->Sumw2(); list->Add(hDEtaDPhiMixEvent);
  
  TH2D *hDEtaDPhiSameEventLowMid = new TH2D("hDEtaDPhiSameEventLowMid", "; #Delta#phi; #Delta#eta", 40, -0.5*pi, 1.5*pi, 40, -4, 4); hDEtaDPhiSameEventLowMid->Sumw2(); list->Add(hDEtaDPhiSameEventLowMid);
  TH2D *hDEtaDPhiMixEventLowMid  = new TH2D("hDEtaDPhiMixEventLowMid", "; #Delta#phi; #Delta#eta", 40, -0.5*pi, 1.5*pi, 40, -4, 4); hDEtaDPhiMixEventLowMid->Sumw2(); list->Add(hDEtaDPhiMixEventLowMid);
  
  TH2D *hDEtaDPhiSameEventHighMid = new TH2D("hDEtaDPhiSameEventHighMid", "; #Delta#phi; #Delta#eta", 40, -0.5*pi, 1.5*pi, 40, -4, 4); hDEtaDPhiSameEventHighMid->Sumw2(); list->Add(hDEtaDPhiSameEventHighMid);
  TH2D *hDEtaDPhiMixEventHighMid  = new TH2D("hDEtaDPhiMixEventHighMid", "; #Delta#phi; #Delta#eta", 40, -0.5*pi, 1.5*pi, 40, -4, 4); hDEtaDPhiMixEventHighMid->Sumw2(); list->Add(hDEtaDPhiMixEventHighMid);
  
  TH3D *hNchDEtaDPhiSameEvent = new TH3D("hNchDEtaDPhiSameEvent", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiSameEvent->Sumw2(); list->Add(hNchDEtaDPhiSameEvent);
  TH3D *hNchDEtaDPhiMixEvent  = new TH3D("hNchDEtaDPhiMixEvent", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiMixEvent->Sumw2(); list->Add(hNchDEtaDPhiMixEvent);
  
  TH3D *hNchDEtaDPhiSameEventKshort = new TH3D("hNchDEtaDPhiSameEventKshort", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiSameEventKshort->Sumw2(); list->Add(hNchDEtaDPhiSameEventKshort);
  TH3D *hNchDEtaDPhiMixEventKshort  = new TH3D("hNchDEtaDPhiMixEventKshort", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiMixEventKshort->Sumw2(); list->Add(hNchDEtaDPhiMixEventKshort);
  
  TH3D *hNchDEtaDPhiSameEventLambda = new TH3D("hNchDEtaDPhiSameEventLambda", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiSameEventLambda->Sumw2(); list->Add(hNchDEtaDPhiSameEventLambda);
  TH3D *hNchDEtaDPhiMixEventLambda  = new TH3D("hNchDEtaDPhiMixEventLambda", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiMixEventLambda->Sumw2(); list->Add(hNchDEtaDPhiMixEventLambda);
  
  TH3D *hNchDEtaDPhiSameEventOmega = new TH3D("hNchDEtaDPhiSameEventOmega", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiSameEventOmega->Sumw2(); list->Add(hNchDEtaDPhiSameEventOmega);
  TH3D *hNchDEtaDPhiMixEventOmega  = new TH3D("hNchDEtaDPhiMixEventOmega", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiMixEventOmega->Sumw2(); list->Add(hNchDEtaDPhiMixEventOmega);
  
  TH3D *hNchDEtaDPhiSameEventCascade = new TH3D("hNchDEtaDPhiSameEventCascade", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiSameEventCascade->Sumw2(); list->Add(hNchDEtaDPhiSameEventCascade);
  TH3D *hNchDEtaDPhiMixEventCascade  = new TH3D("hNchDEtaDPhiMixEventCascade", "; #Delta#phi; #Delta#eta; N_{trk}^{mid}", 40, -0.5*pi, 1.5*pi, 40, -4, 4, 300, -0.5, 299.5); hNchDEtaDPhiMixEventCascade->Sumw2(); list->Add(hNchDEtaDPhiMixEventCascade);
  
  TH2D *hSphericity = new TH2D("hSphericity", ";S;N_{ch}^{mid}", 50, 0, 1, 20, 0, 200); hSphericity->Sumw2(); list->Add(hSphericity);


//=============================================================================

  //cuts
  const double dcEtaCut  = 2;//Mid rapidity eta cut
  const double dfEtaMin  = 2.;//forward rapidity eta min
  const double dfEtaMax  = 5.;//forward rapidity eta max
  const double dPtCut    = 0.15;// all track pT cut

  const double dPtTrMax = 50;//trigger track pT max
  const double dPtTrMin = 3;//trigger track pT min
  const double dPtAsMax = 2.5;//associate track pT max
  const double dPtAsMin = 1;//associate track pT min
  
  //temperary containers
  TVector3 vec_tmp;
  vector<int> iTrIndex;
  vector<int> iAsIndex;
  vector<TVector3> vTrArray;
  const int iMixSize=20;
  vector< vector<TVector3> > vTrArrayPool;

//=============================================================================

//event loop starts
  for (int iEvent=0; iEvent<pythia.mode("Main:numberOfEvents"); iEvent++) if (pythia.next()) { //pythia.next() Generate an(other) event
	  //event wise initializations 
	  double dFwdCh = 0.;
	  double dMidCh = 0.;
	  iTrIndex.clear();
	  iAsIndex.clear();
	  vTrArray.clear();

//=============================================================================
         //determine charged track number in mid/fwd rapidity in this event
	 
	 //track && particles loop
	 for (int i=0; i<event.size(); i++) if (event[i].isFinal()   &&
                                           event[i].isVisible() &&
                                           //event[i].isCharged() &&
                                           (event[i].pT()>dPtCut)) { 
		 double dEtaAbs = TMath::Abs(event[i].eta());
		 if( event[i].isCharged() ){
			 if ((dEtaAbs>dfEtaMin) && (dEtaAbs<=dfEtaMax)) dFwdCh += 1.;//2<|\eta|<5
			 if( dEtaAbs<dcEtaCut ) dMidCh++;//|\eta|<2

		       	 //find trig
			 if( dEtaAbs<dcEtaCut && (event[i].pT()>dPtTrMin && event[i].pT()<dPtTrMax) ) {//trigger track 3<pT<50
				 iTrIndex.push_back(i);
				 vec_tmp.SetXYZ(event[i].px(), event[i].py(), event[i].pz());
				 vTrArray.push_back(vec_tmp);
				 hTrEta->Fill(vec_tmp.Eta());
				 hTrPt->Fill(vec_tmp.Pt());
			 }
		 } 
		 //find associate
		 if( dEtaAbs<dcEtaCut && (event[i].pT()>dPtAsMin && event[i].pT()<dPtAsMax) ) iAsIndex.push_back(i); 
	 } 
	 hFwdVsMid->Fill(dMidCh, dFwdCh);
//=============================================================================
         //connect the inclusive particle properties with the forward charged
         //track number
         for (int i=0; i<event.size(); i++) if (event[i].isFinal() &&
                                           event[i].isVisible()) { 
		 double dEtaAbs = TMath::Abs(event[i].eta()); if (dEtaAbs>=dcEtaCut) continue;

                 int  id      = event[i].idAbs();
                 bool bKshort = (id==310);
                 bool bKCh    = (id==321);
                 bool bLambda = (id==3122);
                 bool bProton = (id==2212);
                 bool bPiCh   = (id==211);
                 bool bResPhi = (id==333);

		 bool bOmega   = (id==3334);
		 bool bCascade = (id==3312); 

		 if ( (!bKshort) && (!bLambda) && (!bProton) && (!bResPhi) && (!bPiCh) && (!bKCh)) continue; 

		 double dPt = event[i].pT();
		 if (bKshort) { hKshort->Fill(dPt,dFwdCh); hProfKshort->Fill(dFwdCh,dPt); hKshortMid->Fill(dPt, dMidCh); }
		 if (bKCh) { hKshort->Fill(dPt,dFwdCh); hProfKCh->Fill(dFwdCh,dPt); }
		 if (bLambda) { hLambda->Fill(dPt,dFwdCh); hProfLambda->Fill(dFwdCh,dPt); hLambdaMid->Fill(dPt, dMidCh); }
		 if (bProton) { hProton->Fill(dPt,dFwdCh); hProfProton->Fill(dFwdCh,dPt); }
		 if (bPiCh) { hPiCh->Fill(dPt,dFwdCh); hProfPiCh->Fill(dFwdCh,dPt); }
		 if (bResPhi) { hResPhi->Fill(dPt,dFwdCh); hProfResPhi->Fill(dFwdCh,dPt); }
	 }

//=============================================================================

	//do correlation in the same event
	//findIn removes the possibility to double count the same pair
	for ( int iTr=0; iTr<iTrIndex.size(); iTr++ ) if(iTrIndex.size()>=1) {
	       	for ( int iAs=iTr; iAs<iAsIndex.size(); iAs++ ) {
		       	if( iAsIndex[iAs]!=iTrIndex[iTr] && !findIn(0, iTr, iAsIndex.at(iAs), iTrIndex) ) {
			       	int iIdxAs = iAsIndex[iAs];
				vec_tmp.SetXYZ(event[iIdxAs].px(), event[iIdxAs].py(), event[iIdxAs].pz());
				double dDPhi = vec_tmp.DeltaPhi(vTrArray.at(iTr));
				if(dDPhi<-0.5*pi) dDPhi=dDPhi + 2*pi;
				if(dDPhi>1.5*pi) dDPhi=dDPhi - 2*pi;
				double dDEta = vec_tmp.Eta() - vTrArray.at(iTr).Eta();
				if(event[iIdxAs].isCharged()) {
				       	hDEtaDPhiSameEvent->Fill(dDPhi, dDEta);
					hNchDEtaDPhiSameEvent->Fill(dDPhi, dDEta, dMidCh);
					hAsPt->Fill(vec_tmp.Pt());

					if(dMidCh<20) 
						hDEtaDPhiSameEventLowMid->Fill(dDPhi, dDEta);
						else if(dMidCh>80)
						       	hDEtaDPhiSameEventHighMid->Fill(dDPhi, dDEta);
					}
			       	if(event[iIdxAs].idAbs()==310)
					hNchDEtaDPhiSameEventKshort->Fill(dDPhi, dDEta, dMidCh);
				else if(event[iIdxAs].idAbs()==3122)
					hNchDEtaDPhiSameEventLambda->Fill(dDPhi, dDEta, dMidCh);
				else if(event[iIdxAs].idAbs()==3334)
					hNchDEtaDPhiSameEventOmega->Fill(dDPhi, dDEta, dMidCh);
				else if(event[iIdxAs].idAbs()==3312)
					hNchDEtaDPhiSameEventCascade->Fill(dDPhi, dDEta, dMidCh); 
			
			}//fi associate indx different from trig
	       	}//for associate
       	}//for trig

	//if mix pool exceeds limit erases the 0th element
	if(vTrArrayPool.size()>iMixSize) 
		vTrArrayPool.erase(vTrArrayPool.begin());

		//do mix event correlation
		//only mix when accumulated pool size reaches the
		//requirement
		for ( int iTrEvent=0; iTrEvent<vTrArrayPool.size(); iTrEvent++ ) if(vTrArrayPool.size()==iMixSize) {
			for ( int iTr=0; iTr<vTrArrayPool.at(iTrEvent).size(); iTr++ ) {
				for ( int iAs=0; iAs<iAsIndex.size(); iAs++ ) {
					int iIdxAs = iAsIndex[iAs];
					vec_tmp.SetXYZ(event[iIdxAs].px(), event[iIdxAs].py(), event[iIdxAs].pz());
					double dDPhi = vec_tmp.DeltaPhi(vTrArrayPool.at(iTrEvent).at(iTr));
					if(dDPhi<-0.5*pi) dDPhi=dDPhi + 2*pi;
					if(dDPhi>1.5*pi) dDPhi=dDPhi - 2*pi;
					double dDEta = vec_tmp.Eta() - vTrArrayPool.at(iTrEvent).at(iTr).Eta();
					if(event[iIdxAs].isCharged()) {
						hDEtaDPhiMixEvent->Fill(dDPhi, dDEta);
						hNchDEtaDPhiMixEvent->Fill(dDPhi, dDEta, dMidCh);

						if(dMidCh<20) 
							hDEtaDPhiMixEventLowMid->Fill(dDPhi, dDEta);
						else if(dMidCh>80)
							hDEtaDPhiMixEventHighMid->Fill(dDPhi, dDEta);
					}

					if(event[iIdxAs].idAbs()==310)
						hNchDEtaDPhiMixEventKshort->Fill(dDPhi, dDEta, dMidCh);
					else if(event[iIdxAs].idAbs()==3122)
						hNchDEtaDPhiMixEventLambda->Fill(dDPhi, dDEta, dMidCh);
					else if(event[iIdxAs].idAbs()==3334)
						hNchDEtaDPhiMixEventOmega->Fill(dDPhi, dDEta, dMidCh);
					else if(event[iIdxAs].idAbs()==3312)
						hNchDEtaDPhiMixEventCascade->Fill(dDPhi, dDEta, dMidCh);

				}//for as
			}//for trig in one loop
		}//for trip pool loop

		//found 1 accept trig list, save to the
		//event mixing pool
		if(iTrIndex.size()>=1) {
			vTrArrayPool.push_back(vTrArray);
			hFwdVsMidTag->Fill(dMidCh, dFwdCh);
		}

    hEvent->Fill(info.pTHat(), dFwdCh);
  }
//=============================================================================

  hXsect->Fill(0.5, info.sigmaGen());
  hTrials->Fill(0.5, info.weightSum());
  file->cd(); list->Write(); file->Close();
//=============================================================================

  //pythia.stat();//obtain statistics on the number of events generated of the different kinds, and the estimated cross sections.
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
  UInt_t kTime = (UInt_t)adt.Get();//系统时间
  UInt_t kProc = (UInt_t)gSystem->GetPid();//Get process id.
  UInt_t kInit = kTime - kProc;
  int kSeed = (int)kInit - (int)(((int)(kInit/1e8))*1e8);//种子最大是8位数,对数位大的部分进行裁剪

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

  ::Info(sMethod.Data(), "CR status: %s", s.Data());

  if (s.IsNull()) {
    ::Warning(sMethod.Data(), "CR flag has not yet been set, using defaul (= off)");
  } else {
    if (s=="ON") {
      ::Info(sMethod.Data(), "ColourReconnection:reconnect = on");
      return true;
    } else {
      ::Info(sMethod.Data(), "ColourReconnection:reconnect = off");
      return false;
    }
  }

  ::Info(sMethod.Data(), "ColourReconnection:reconnect = off");
  return false;
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

bool findIn(int start, int end, int object, std::vector<int> index){
	for(int iter=start; iter<end; iter++){
		if(index.at(iter)==object) 
			return true;
	} 
	return false;
}
  

