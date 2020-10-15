//
//TH1D *GetTH1D(char *file,char *list , char * name)
//TH1D *GetTH2D(char *file,char *list , char * name)
//void DrawHisto(TH1D *h, Color_t wc, Style_t ws, Option_t *opt)
//TH1D* RebinTH1D(TH1D const *hRaw, TH1D const *hRef)
//TH2D * RebinTH2D(TH2D const *hF, TH2D const *hR)
//void NormHistBinWidth(TH1D *h, Double_t dw=1.)
//void DeNormHistBinWidth(TH1D* h)
//TH1D* MaxHistograms(const Int_t n, TH1D **h) //找到所有直方图中每个Bin最大的值组成hMax
//TCanvas* MakeCanvas(TString sName)
//void CanvasEnd(TCanvas *c)
//void SetLegend(TLegend *l)
//void SetFrame(TH1 *f, const TString stnx="",   const TString stny="",...........
//void DrawAliLogo(Double_t dXmin, Double_t dYmin, Bool_t bPre=kTRUE, Int_t ip=24)
//TH1D* NormIncl(TH1D* hIncl, TH1D* hEvent)
//TH1D* NormJC(TH1D* hJC, TH1D* hEvent, Double_t R = 0.4)
//TH1D* NormPC(TH1D* hPCL, TH1D* hPCU, TH1D* hEvent, Double_t R = 0.4)
//
//
//

//_____________________________________________________________________________
const Color_t cLine[] = { kBlack,  kRed+1,  kBlue+1, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2,  kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
//const Color_t wcf[] = { kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7 };
const Style_t sMark[] = { kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kFullSquare,  kOpenSquare, kFullCross, kOpenCross, kFullStar, kOpenStar};

//_____________________________________________________________________________
TH1D *GetTH1D(const TString path, const TString file, const TString list , const TString name)
{
  TString sFile = Form("%s/%s", path.Data(), file.Data());
  TFile *f1 = TFile::Open(sFile.Data(),"read");

  if(!(list.IsNull())){
    TList *flist = (TList *)f1->Get(list);
    cout<<(TH1D*)flist->FindObject(name)<<endl;
    return (TH1D*)flist->FindObject(name);
  }else{
    return (TH1D*)f1->Get(name);
  }

  cout<<(TH1D*)f1->Get(name)<<endl;
}

//_____________________________________________________________________________
TH2D *GetTH2D(const TString path, const TString file, const TString list , const TString name)
{
  TString sFile = Form("%s/%s", path.Data(), file.Data());
  TFile *f1 = TFile::Open(sFile.Data(),"read");

  if(!(list.IsNull())){
    TList *flist = (TList *)f1->Get(list);
    cout<<(TH2D*)flist->FindObject(name)<<endl;
    return (TH2D*)flist->FindObject(name);
  }else{
    return (TH2D*)f1->Get(name);
  }

  cout<<(TH2D*)f1->Get(name)<<endl;
}

//_____________________________________________________________________________
void DrawHisto(TH1D *h, Color_t wc, Style_t ws, Option_t *opt)
{
  h->SetLineWidth(2);
  h->SetLineColor(wc);
  h->SetMarkerStyle(ws);
  h->SetMarkerColor(wc);
  h->SetFillStyle(0);
  //h->DrawCopy(opt);
  h->Draw(opt);

  return;
}

//_____________________________________________________________________________
TH1D* RebinTH1D(TH1D const *hRaw, TH1D const *hRef)
{
  if ((!hRaw) || (!hRef)) return 0x0;

  const Int_t nRef = hRef->GetNbinsX();

  const Double_t dLower = hRef->GetXaxis()->GetBinLowEdge(1);
  const Double_t dUpper = hRef->GetXaxis()->GetBinUpEdge(nRef);

  auto hNew =(TH1D*)hRef->Clone(Form("h_%s_%s", hRaw->GetName(),
                                           hRef->GetName()));
  hNew->Reset();
  //=============================================================================

  for (Int_t k=1; k<=hRaw->GetNbinsX(); k++) {
    Double_t dXvar = hRaw->GetBinCenter(k);
    if ((dXvar<dLower) || (dXvar>=dUpper)) continue;

    Int_t iBinX = hNew->FindBin(dXvar);
    Double_t dYvar = hRaw->GetBinContent(k);
    Double_t dYerr = hRaw->GetBinError(k);
    Double_t dYsw2 = hNew->GetBinError(iBinX);

    hNew->SetBinContent(iBinX, hNew->GetBinContent(iBinX) + dYvar);
    hNew->SetBinError(iBinX, TMath::Sqrt(dYsw2*dYsw2 + dYerr*dYerr));
  }

  return hNew;
}


//_____________________________________________________________________________
void NormHistBinWidth(TH1D *h, Double_t dw=1.)
{
  if (!h) return;
  //=============================================================================

  for (Int_t k=1; k<=h->GetNbinsX(); k++) {
    Double_t dBW = h->GetBinWidth(k) * dw;
    h->SetBinContent(k, h->GetBinContent(k)/dBW);
    h->SetBinError(k, h->GetBinError(k)/dBW);
  }

  return;
}

//_____________________________________________________________________________
void DeNormHistBinWidth(TH1D* h)
{
  if (!h) return;
  //=============================================================================

  for (Int_t k=1; k<=h->GetNbinsX(); k++) {
    Double_t dBW = h->GetBinWidth(k);
    h->SetBinContent(k, h->GetBinContent(k) * dBW);
    h->SetBinError(k, h->GetBinError(k) * dBW);
  }

  return;
}

//____________________________________________________________________
TH2D * RebinTH2D(TH2D const *hF, TH2D const *hR)
{
  if ((!hF) || (!hR)) return 0x0;
  const Int_t nx1 = hR->GetNbinsX();
  const Int_t ny1 = hR->GetNbinsY();
  const Int_t nx = nx1;
  const Int_t ny = ny1;

  Double_t dx[nx+1];
  Double_t dy[ny+1];

  dx[nx] = hR->GetXaxis()->GetBinUpEdge(nx);
  dy[ny] = hR->GetYaxis()->GetBinUpEdge(ny);
  for (Int_t i=0, k=1; i<nx; i++, k++) dx[i] = hR->GetXaxis()->GetBinLowEdge(k);
  for (Int_t j=0, l=1; j<ny; j++, l++) dy[j] = hR->GetYaxis()->GetBinLowEdge(l);

  TH2D *hN = new TH2D(Form("hN_%s_%s",hF->GetName(),hR->GetName()), "", nx, dx, ny, dy); hN->Sumw2();
  hN->Reset();

  for (Int_t k=1; k<=hF->GetNbinsX(); k++) {
    Double_t dXvar = hF->GetXaxis()->GetBinCenter(k);
    if ((dXvar<dx[0]) || (dXvar>=dx[nx])) continue;

    for (Int_t l=1; l<=hF->GetNbinsY(); l++) {
      Double_t dYvar = hF->GetYaxis()->GetBinCenter(l);
      if ((dYvar<dy[0]) || (dYvar>=dy[ny])) continue;

      Double_t dvar = hF->GetBinContent(k, l);
      Double_t derr = hF->GetBinError(k, l);

      Int_t nbx = hN->GetXaxis()->FindBin(dXvar);
      Int_t nby = hN->GetYaxis()->FindBin(dYvar);
      Double_t dsw2 = hN->GetBinError(nbx, nby);

      hN->SetBinContent(nbx, nby, dvar + hN->GetBinContent(nbx, nby));
      hN->SetBinError(nbx, nby, TMath::Sqrt(derr*derr + dsw2*dsw2));
    }
  }

  return hN;
}

//_____________________________________________________________________________
TH1D* MaxHistograms(const Int_t n, TH1D **h)
{
  if ((n<=0) || (!h[0])) return 0x0;
  //=============================================================================

  Double_t *d = new Double_t[n];
  TH1D *hMax = (TH1D*)h[0]->Clone("hMax");
  for (Int_t k=1; k<=hMax->GetNbinsX(); k++) {
    for (Int_t j=0; j<n; j++) d[j] = h[j]->GetBinContent(k);
    hMax->SetBinContent(k, TMath::MaxElement(n,d));
    hMax->SetBinError(k, TMath::RMS(n,d));
  }

  return hMax;
}

//_____________________________________________________________________________
TH1D* QuadraticSum(const Int_t n, TH1D **h)
{

  if ((n<=0) || (!h)) return 0x0;
//=============================================================================

  TH1D *hSum = (TH1D*)h[0]->Clone("hSum");
  for (Int_t k=1; k<=hSum->GetNbinsX(); k++) {
    Double_t dSum = 0.;
    for (Int_t j=0; j<n; j++) {
      Double_t d = h[j]->GetBinContent(k);
      dSum += (d*d);
    }

    hSum->SetBinContent(k, TMath::Sqrt(dSum));
  }

  return hSum;
}

//_____________________________________________________________________________
TH1D* CancleMaterBudget(TH1D *h)
{

  if (!h) return 0x0;
//=============================================================================

  TH1D *h0 = (TH1D*)h->Clone("h0");
  for (Int_t k=1; k<=h0->GetNbinsX(); k++) {
    Double_t d = h->GetBinContent(k);
      h0->SetBinContent(k, TMath::Sqrt(TMath::Abs(d*d - 0.04*0.04))); 
    }
  return h0;
}

//_____________________________________________________________________________
void SetFrame(TH1 *f, const TString stnx="",   const TString stny="",
                      const Float_t dlsx=0.05, const Float_t dlsy=0.05,
                      const Float_t dtsx=0.06, const Float_t dtsy=0.06,
                      const Float_t dtox=1.10, const Float_t dtoy=1.00)
{
  f->GetXaxis()->SetLabelSize(dlsx);
  f->GetYaxis()->SetLabelSize(dlsy);

  f->GetXaxis()->SetTitleSize(dtsx);
  f->GetYaxis()->SetTitleSize(dtsy);

  f->GetXaxis()->SetTitleOffset(dtox);
  f->GetYaxis()->SetTitleOffset(dtoy);

  f->SetXTitle(stnx.Data());
  f->SetYTitle(stny.Data());

  f->GetYaxis()->SetNdivisions(505);

  return;
}

//_____________________________________________________________________________
TCanvas* MakeCanvas(TString sName)
{
  TCanvas *c = new TCanvas(Form("c%s",sName.Data()), sName.Data(), 700, 500);
  c->Range(0., 0., 1., 1.);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.13);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.15);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);

  return c;
}

//_____________________________________________________________________________
void CanvasEnd(TCanvas *c)
{
  if (!c) return;

  c->Modified();
  c->cd();
  c->SetSelected(c);

  return;
}

//_____________________________________________________________________________
TPad* MakePadT(TString sName)
{
  TPad *p = new TPad(Form("c%s_U",sName.Data()), Form("%s_U",sName.Data()), 0., 0.35, 1., 1.);
  p->Range(0., 0., 1., 1.);
  p->SetFillColor(0);
  p->SetBorderMode(0);
  p->SetBorderSize(0);
  p->SetRightMargin(0.03);
  p->SetLeftMargin(0.18);
  p->SetTopMargin(0.02);
  p->SetBottomMargin(0.);
  p->SetFrameFillStyle(0);
  p->SetFrameBorderMode(0);
  p->Draw();
  p->cd();

  return p;
}

//_____________________________________________________________________________
TPad* MakePadB(TString sName)
{
  TPad *p = new TPad(Form("c%s_D",sName.Data()), Form("%s_D",sName.Data()), 0., 0. , 1., 0.35);
  p->Range(0., 0., 1., 1.);
  p->SetFillColor(0);
  p->SetBorderMode(0);
  p->SetBorderSize(0);
  p->SetRightMargin(0.03);
  p->SetLeftMargin(0.18);
  p->SetTopMargin(0.);
  p->SetBottomMargin(0.25);
  p->SetFrameFillStyle(0);
  p->SetFrameBorderMode(0);
  p->Draw();
  p->cd();

  return p;
}

//_____________________________________________________________________________
TPad* MakePadL(TString sName)
{
  TPad *p = new TPad(Form("c%s_L",sName.Data()), Form("%s_L",sName.Data()), 0., 0., 0.52, 1.);
  p->Range(0., 0., 1., 1.);
  p->SetFillColor(0);
  p->SetBorderMode(0);
  p->SetBorderSize(0);
  p->SetRightMargin(0.);
  p->SetLeftMargin(0.18);
  p->SetTopMargin(0.02);
  p->SetBottomMargin(0.15);
  p->SetFrameFillStyle(0);
  p->SetFrameBorderMode(0);
  p->Draw();
  p->cd();

  return p;
}

//_____________________________________________________________________________
TPad* MakePadR(TString sName)
{
  TPad *p = new TPad(Form("c%s_R",sName.Data()), Form("%s_R",sName.Data()), 0.52, 0., 1., 1.);
  p->Range(0., 0., 1., 1.);
  p->SetFillColor(0);
  p->SetBorderMode(0);
  p->SetBorderSize(0);
  p->SetRightMargin(0.03);
  p->SetLeftMargin(0.);
  p->SetTopMargin(0.02);
  p->SetBottomMargin(0.15);
  p->SetFrameFillStyle(0);
  p->SetFrameBorderMode(0);
  p->Draw();
  p->cd();

  return p;
}


//_____________________________________________________________________________
void SetLegend(TLegend *l)
{
  l->SetFillColor(0);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.8*gStyle->GetTextSize());

  return;
}


//_____________________________________________________________________________
TH1D* NormIncl(TH1D* hIncl, TH1D* hEvent)
{

  //Double_t dEvent = hEvent->GetEntries();
  Double_t dEvent = 0.;
  for(Int_t i=1; i<= hEvent->GetNbinsX(); i++) {dEvent += hEvent->GetBinContent(i);}
  cout<<"..............."<<"Incl Events = "<<dEvent<<"................"<<endl;
  //hIncl->Scale(1./(dEvent*2*0.75*TMath::TwoPi()));
  hIncl->Scale(1./(dEvent*0.5*TMath::TwoPi()));
  //hIncl->Scale(1./(dEvent*0.5));
  return hIncl;
}

//_____________________________________________________________________________
TH1D* NormJC(TH1D* hJC, TH1D* hEvent, Double_t R = 0.4)
{
  //Double_t dJetEvent = hEvent->GetEntries();
  Double_t dJetEvent = 0.;
  for(Int_t i=1; i<= hEvent->GetNbinsX(); i++) {dJetEvent += hEvent->GetBinContent(i);}
  cout<<"..............."<<"Jet Events = "<<dJetEvent<<"................"<<endl;

  hJC->Scale(1./(dJetEvent*2*0.75*TMath::TwoPi()*0.06));
  return hJC;
}

//_____________________________________________________________________________
TH1D* NormPC(TH1D* hPCL, TH1D* hPCU, TH1D* hEvent, Double_t R = 0.4)
{
  //Double_t dJetEvent = hEvent->GetEntries();
  Double_t dJetEvent = 0.;
  for(Int_t i=1; i<= hEvent->GetNbinsX(); i++) {dJetEvent += hEvent->GetBinContent(i);}
  cout<<"..............."<<"Jet Events = "<<dJetEvent<<"................"<<endl;
  
  auto hPC = (TH1D*)hPCL->Clone("hPC"); hPC->Reset();
  hPC->Add(hPCL, hPCU, 1, 1);
  hPC->Scale(1./(dJetEvent*2*0.75*TMath::TwoPi()*0.24));
  return hPC;
}

//_____________________________________________________________________________
void CalcRelDeviation(TH1D *hi, TH1D *hr)
{
  if ((!hi) || (!hr)) return;
//=============================================================================

  for (Int_t k=1; k<=hr->GetNbinsX(); k++) {
    Double_t dVarI = hi->GetBinContent(k); if (TMath::Abs(dVarI)<1e-12) dVarI = 1e-12;
    Double_t dErrI = hi->GetBinError(k) / dVarI;

    Double_t dVarR = hr->GetBinContent(k); if (TMath::Abs(dVarR)<1e-12) dVarR = 1e-12;
    Double_t dErrR = hr->GetBinError(k) / dVarR;

    Double_t dVarD = TMath::Abs(dVarI/dVarR - 1.);
    Double_t dErrD = dVarD * TMath::Sqrt(dErrI*dErrI + dErrR*dErrR);
    hi->SetBinContent(k, dVarD); hi->SetBinError(k, dErrD);
  }

  return;
}

//_____________________________________________________________________________
TGraphErrors* ConvHistogramToGraphErrors(TH1D *hVar, TH1D *hErr, const Int_t n, const Bool_t bSmall=kFALSE)
{
  Double_t *dtm = new Double_t[n];

  Double_t dvx[n], dvy[n];
  Double_t dex[n], dey[n];
  for (Int_t i=0, k=1; i<n; i++, k++) {
    dvx[i] = hVar->GetBinCenter(k);
    dex[i] = hVar->GetBinWidth(k) / 2.;

    if (bSmall) {
      dex[i]  = 0.01;
    } else {
      if (dex[i]>0.2) dex[i] = 0.1;
    }

    dvy[i] = hVar->GetBinContent(k);
    if (hErr) {
      dey[i] = hErr->GetBinContent(k) * dvy[i];
    } else {
      dey[i] = hVar->GetBinError(k);
    }
  }

  TGraphErrors *g = new TGraphErrors(n, dvx, dvy, dex, dey);

  return g;
}

//_____________________________________________________________________________
void DrawAliLogo (Double_t dXmin, Double_t dYmin, Int_t ip=24, Bool_t bPre=kFALSE)
{
  TLatex *t = new TLatex(dXmin, dYmin, (bPre ? "ALICE Preliminary" : "ALICE"));
  t->SetNDC();
  t->SetTextFont(42);
  t->SetTextSizePixels(ip);
  t->Draw();

  return;
}

//_____________________________________________________________________________
void DrawGraph(TGraph *g, Color_t wc, Option_t *opt)
{
  g->SetLineWidth(2);
  g->SetLineColor(wc);
  g->SetFillStyle(0);
  g->SetFillColor(g->GetLineColor());
  g->Draw(opt);

  return;
}


//_____________________________________________________________________________
void DrawGraph(TGraph *g, Color_t wc, Style_t ws, Option_t *opt)
{
  g->SetLineWidth(2);
  g->SetLineColor(wc);
  g->SetMarkerStyle(ws);
  g->SetMarkerColor(wc);
  g->SetFillStyle(0);
  g->Draw(opt);

  return;
}

//_____________________________________________________________________________
void DrawGraph(TGraphErrors *g, Color_t wc, Option_t *opt)
{
  g->SetLineWidth(2);
  g->SetLineColor(wc);
  g->SetFillStyle(0);
  g->SetFillColor(g->GetLineColor());
  g->Draw(opt);

  return;
}

//_____________________________________________________________________________
void DrawGraph(TGraphErrors *g, Color_t wc, Style_t ws, Option_t *opt)
{
  g->SetLineWidth(2);
  g->SetLineColor(wc);
  g->SetMarkerStyle(ws);
  g->SetMarkerColor(wc);
  g->SetFillStyle(0);
  g->Draw(opt);

  return;
}
//_____________________________________________________________________________
void DrawGraph(TGraphAsymmErrors *g, Color_t wc, Option_t *opt)
{
  g->SetLineWidth(2);
  g->SetLineColor(wc);
  g->SetFillStyle(0);
  g->SetFillColor(g->GetLineColor());
  g->Draw(opt);

  return;
}

//_____________________________________________________________________________
void DrawGraph(TGraphAsymmErrors*g, Color_t wc, Style_t ws, Option_t *opt)
{
  g->SetLineWidth(2);
  g->SetLineColor(wc);
  g->SetMarkerStyle(ws);
  g->SetMarkerColor(wc);
  g->SetFillStyle(0);
  g->Draw(opt);

  return;
}

