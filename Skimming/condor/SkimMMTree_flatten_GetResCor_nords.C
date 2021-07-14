#include <ctime>

#include <TLorentzVector.h>
#include "../../HeaderFiles/commonUtility.h"
#include "../../HeaderFiles/HiEvtPlaneList.h"
#include "../../HeaderFiles/cutsAndBinUpsilonV2.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"

const Double_t pi = 3.141592653589;

Double_t getDPHI_Jared( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;
 
  if ( dphi > pi ) dphi = dphi - 2*pi;
  if ( dphi <= -pi ) dphi = dphi + 2*pi;
  if ( fabs(dphi) > pi ) {
    //cout << " getDPHI error!!! dphi=" << phi1 << "-" << phi2 << " is bigger than 3.141592653589 " << endl;
    dphi = -10;
  }

  return dphi;
}

Double_t restrictToPiOver2(Double_t phi) {

  Double_t newphi = phi;
  if (phi>pi/2) newphi = phi-pi;
  else if (phi<-pi/2) newphi = phi+pi;
  return newphi;
}

bool isAbout(float num1=0.0, float num2=0.0) {

  float numavg = (num1+num2)/2;
  float percdiff = abs(num1-num2)/numavg;
  if (percdiff>0.01) return kFALSE;
  else return kTRUE;

}

static const long MAXTREESIZE = 1000000000000;

void SkimMMTree_flatten_GetResCor_nords(int nevt=-1,
      int dateStr=20210712,
      bool flattenBinByBin=kTRUE) 
{

  const int flatOrder = 21;

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  float ptbins[5] = {0,3,6,10,50};
  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  float ybins[4] = {0.0,0.8,1.6,2.4};
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  Double_t cbins[4] = {10,30,50,90};
  const int numcbins = sizeof(cbins)/sizeof(double)-1;

  TH1D* hpt = new TH1D("hpt","hist vs pt",numptbins,ptbins);
  TH1D* hy = new TH1D("hy","hist vs y",numybins,ybins);
  TH1D* hcent = new TH1D("hcent","hist vs cent",numcbins,cbins);

  //gStyle->SetOptStat(0);
  //TH1D* hCosAC = new TH1D("hCosAC","cos(2*(psiA-psiC))",50,-1.2,1.2);
  TH1D* hRpt = new TH1D("hRpt","EP Resolution factor vs pt",numptbins,ptbins);
  TH1D* hRy = new TH1D("hRy","EP Resolution factor vs y",numybins,ybins);
  TH1D* hRc = new TH1D("hRc","EP Resolution factor vs cent",numcbins,cbins);

  //Event planes
  TH1D* hEpHF2old = new TH1D("hEpHF2old","hEpHF2old",50,-2,2);
  TH1D* hEpHF2new = new TH1D("hEpHF2new","hEpHF2new",50,-2,2);
  TH1D* hEpHFm2old = new TH1D("hEpHFm2old","hEpHFm2old",50,-2,2);
  TH1D* hEpHFm2new = new TH1D("hEpHFm2new","hEpHFm2new",50,-2,2);
  TH1D* hEpHFp2old = new TH1D("hEpHFp2old","hEpHFp2old",50,-2,2);
  TH1D* hEpHFp2new = new TH1D("hEpHFp2new","hEpHFp2new",50,-2,2);
  TH1D* hEptrackmid2old = new TH1D("hEptrackmid2old","hEptrackmid2old",50,-2,2);
  TH1D* hEptrackmid2new = new TH1D("hEptrackmid2new","hEptrackmid2new",50,-2,2);

  TH1D* hEpHF2 = new TH1D("hEpHF2","hEpHF2",50,-2,2);
  TH1D* hEpHFm2 = new TH1D("hEpHFm2","hEpHFm2",50,-2,2);
  TH1D* hEpHFp2 = new TH1D("hEpHFp2","hEpHFp2",50,-2,2);
  TH1D* hEptrackmid2 = new TH1D("hEptrackmid2","hEptrackmid2",50,-2,2);

  TH1D* hphi = new TH1D("hphi","hphi",50,-2,2);
  TH1D* hdphi = new TH1D("hdphi","hdphi",50,-2,2);
  //TH1D* hdphiw = new TH1D("hdphiw","hdphiw",50,-2,2);

  TH1F* massHistoUpsilonCandidates = new TH1F( "massHistoUpsilonCandidates",  "massHistoUpsilonCandidates;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", 60, 8, 14);

  TString inFileName = Form("skims/newOniaTree_Skim_UpsTrig_MM_RECENTERED_nords_n%i_20201118.root",nevt);
  if (nevt==-1) inFileName = "skims/newOniaTree_Skim_UpsTrig_MM_RECENTERED_nords_n-1_20201118.root";
  TFile* inFile = TFile::Open(inFileName,"READ");
  TTree* mytree = (TTree*)inFile->Get("mmep");
  TBranch* mm = (TBranch*)mytree->GetBranch("mm");

  TString binByBinStr = "";
  if (flattenBinByBin) binByBinStr = "BinByBin";

  TFile* newfile;
  newfile = new TFile(Form("skims/newOniaTree_Skim_UpsTrig_MM_flattened%s_order%i_n%i_%i.root",binByBinStr.Data(),flatOrder,nevt,dateStr),"recreate");

  DiMuon dm;
  TTree* mmtree = new TTree("mmep","dimuonAndEventPlanes");
  mmtree->SetMaxTreeSize(MAXTREESIZE);
  mmtree->Branch("mm",&dm.run,branchString.Data());

  int ALLPASS = 0;

  if(nevt == -1) nevt = 56114317;

  int nevtReal = mytree->GetEntries();

  //averages for resolution correction
  Double_t avgCosAB = 0;
  Double_t avgCosAC = 0;
  Double_t avgCosBC = 0;
  Double_t avgCosABpt[numptbins] = {0};
  Double_t avgCosACpt[numptbins] = {0};
  Double_t avgCosBCpt[numptbins] = {0};
  Double_t avgCosABy[numybins] = {0};
  Double_t avgCosACy[numybins] = {0};
  Double_t avgCosBCy[numybins] = {0};
  Double_t avgCosABcent[numcbins] = {0};
  Double_t avgCosACcent[numcbins] = {0};
  Double_t avgCosBCcent[numcbins] = {0};
  //errors on averages
  Double_t sumsqrsCosAB = 0;
  Double_t sumsqrsCosAC = 0;
  Double_t sumsqrsCosBC = 0;
  Double_t sumsqrsCosABpt[numptbins] = {0};
  Double_t sumsqrsCosACpt[numptbins] = {0};
  Double_t sumsqrsCosBCpt[numptbins] = {0};
  Double_t sumsqrsCosABy[numybins] = {0};
  Double_t sumsqrsCosACy[numybins] = {0};
  Double_t sumsqrsCosBCy[numybins] = {0};
  Double_t sumsqrsCosABcent[numcbins] = {0};
  Double_t sumsqrsCosACcent[numcbins] = {0};
  Double_t sumsqrsCosBCcent[numcbins] = {0};

  //get values for flattening
  TFile* avgFile = TFile::Open(Form("averages/avgEpFile%i.root",nevt),"READ");
  TH1D* havgCosEp = (TH1D*)avgFile->Get("havgCosEp;1");
  TH1D* havgSinEp = (TH1D*)avgFile->Get("havgSinEp;1");
  TH1D* havgCosEpHFm2 = (TH1D*)avgFile->Get("havgCosEpHFm2;1");
  TH1D* havgSinEpHFm2 = (TH1D*)avgFile->Get("havgSinEpHFm2;1");
  TH1D* havgCosEpHFp2 = (TH1D*)avgFile->Get("havgCosEpHFp2;1");
  TH1D* havgSinEpHFp2 = (TH1D*)avgFile->Get("havgSinEpHFp2;1");
  TH1D* havgCosEptrackmid2 = (TH1D*)avgFile->Get("havgCosEptrackmid2;1");
  TH1D* havgSinEptrackmid2 = (TH1D*)avgFile->Get("havgSinEptrackmid2;1");
  Double_t avgCosEp[flatOrder] = {0};
  Double_t avgSinEp[flatOrder] = {0};
  Double_t avgCosEpHFm2[flatOrder] = {0};
  Double_t avgSinEpHFm2[flatOrder] = {0};
  Double_t avgCosEpHFp2[flatOrder] = {0};
  Double_t avgSinEpHFp2[flatOrder] = {0};
  Double_t avgCosEptrackmid2[flatOrder] = {0};
  Double_t avgSinEptrackmid2[flatOrder] = {0};
  for (int n=1; n<=flatOrder; n++) {
    avgCosEp[n-1] = havgCosEp->GetBinContent(n);
    avgSinEp[n-1] = havgSinEp->GetBinContent(n);
    avgCosEpHFm2[n-1] = havgCosEpHFm2->GetBinContent(n);
    avgSinEpHFm2[n-1] = havgSinEpHFm2->GetBinContent(n);
    avgCosEpHFp2[n-1] = havgCosEpHFp2->GetBinContent(n);
    avgSinEpHFp2[n-1] = havgSinEpHFp2->GetBinContent(n);
    avgCosEptrackmid2[n-1] = havgCosEptrackmid2->GetBinContent(n);
    avgSinEptrackmid2[n-1] = havgSinEptrackmid2->GetBinContent(n);
    cout << "avg(Cos(2*" << n << "*Psi)) = " << avgCosEp[n-1] << endl;
    cout << "avg(Sin(2*" << n << "*Psi)) = " << avgSinEp[n-1] << endl;
  }

  TH1D* havgCosEppt[numptbins];
  TH1D* havgSinEppt[numptbins];
  TH1D* havgCosEpHFm2pt[numptbins];
  TH1D* havgSinEpHFm2pt[numptbins];
  TH1D* havgCosEpHFp2pt[numptbins];
  TH1D* havgSinEpHFp2pt[numptbins];
  TH1D* havgCosEptrackmid2pt[numptbins];
  TH1D* havgSinEptrackmid2pt[numptbins];

  TH1D* havgCosEpcent[numcbins];
  TH1D* havgSinEpcent[numcbins];
  TH1D* havgCosEpHFm2cent[numcbins];
  TH1D* havgSinEpHFm2cent[numcbins];
  TH1D* havgCosEpHFp2cent[numcbins];
  TH1D* havgSinEpHFp2cent[numcbins];
  TH1D* havgCosEptrackmid2cent[numcbins];
  TH1D* havgSinEptrackmid2cent[numcbins];

  Double_t avgCosEppt[numptbins][flatOrder] = {0};
  Double_t avgSinEppt[numptbins][flatOrder] = {0};
  Double_t avgCosEpHFm2pt[numptbins][flatOrder] = {0};
  Double_t avgSinEpHFm2pt[numptbins][flatOrder] = {0};
  Double_t avgCosEpHFp2pt[numptbins][flatOrder] = {0};
  Double_t avgSinEpHFp2pt[numptbins][flatOrder] = {0};
  Double_t avgCosEptrackmid2pt[numptbins][flatOrder] = {0};
  Double_t avgSinEptrackmid2pt[numptbins][flatOrder] = {0};

  Double_t avgCosEpy[numybins][flatOrder] = {0};
  Double_t avgSinEpy[numybins][flatOrder] = {0};
  Double_t avgCosEpHFm2y[numybins][flatOrder] = {0};
  Double_t avgSinEpHFm2y[numybins][flatOrder] = {0};
  Double_t avgCosEpHFp2y[numybins][flatOrder] = {0};
  Double_t avgSinEpHFp2y[numybins][flatOrder] = {0};
  Double_t avgCosEptrackmid2y[numybins][flatOrder] = {0};
  Double_t avgSinEptrackmid2y[numybins][flatOrder] = {0};

  Double_t avgCosEpcent[numcbins][flatOrder] = {0};
  Double_t avgSinEpcent[numcbins][flatOrder] = {0};
  Double_t avgCosEpHFm2cent[numcbins][flatOrder] = {0};
  Double_t avgSinEpHFm2cent[numcbins][flatOrder] = {0};
  Double_t avgCosEpHFp2cent[numcbins][flatOrder] = {0};
  Double_t avgSinEpHFp2cent[numcbins][flatOrder] = {0};
  Double_t avgCosEptrackmid2cent[numcbins][flatOrder] = {0};
  Double_t avgSinEptrackmid2cent[numcbins][flatOrder] = {0};

  for (int i=0; i<numptbins; i++) {
    havgCosEppt[i] = (TH1D*)avgFile->Get(Form("havgCosEppt[%i];1",i));
    havgSinEppt[i] = (TH1D*)avgFile->Get(Form("havgSinEppt[%i];1",i));
    havgCosEpHFm2pt[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFm2pt[%i];1",i));
    havgSinEpHFm2pt[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFm2pt[%i];1",i));
    havgCosEpHFp2pt[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFp2pt[%i];1",i));
    havgSinEpHFp2pt[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFp2pt[%i];1",i));
    havgCosEptrackmid2pt[i] = (TH1D*)avgFile->Get(Form("havgCosEptrackmid2pt[%i];1",i));
    havgSinEptrackmid2pt[i] = (TH1D*)avgFile->Get(Form("havgSinEptrackmid2pt[%i];1",i));
    for (int n=1; n<=flatOrder; n++) {
      avgCosEppt[i][n-1] = havgCosEppt[i]->GetBinContent(n);
      avgSinEppt[i][n-1] = havgSinEppt[i]->GetBinContent(n);
      avgCosEpHFm2pt[i][n-1] = havgCosEpHFm2pt[i]->GetBinContent(n);
      avgSinEpHFm2pt[i][n-1] = havgSinEpHFm2pt[i]->GetBinContent(n);
      avgCosEpHFp2pt[i][n-1] = havgCosEpHFp2pt[i]->GetBinContent(n);
      avgSinEpHFp2pt[i][n-1] = havgSinEpHFp2pt[i]->GetBinContent(n);
      avgCosEptrackmid2pt[i][n-1] = havgCosEptrackmid2pt[i]->GetBinContent(n);
      avgSinEptrackmid2pt[i][n-1] = havgSinEptrackmid2pt[i]->GetBinContent(n);
    }
  }

  for (int i=0; i<numcbins; i++) {
    havgCosEpcent[i] = (TH1D*)avgFile->Get(Form("havgCosEpcent[%i];1",i));
    havgSinEpcent[i] = (TH1D*)avgFile->Get(Form("havgSinEpcent[%i];1",i));
    havgCosEpHFm2cent[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFm2cent[%i];1",i));
    havgSinEpHFm2cent[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFm2cent[%i];1",i));
    havgCosEpHFp2cent[i] = (TH1D*)avgFile->Get(Form("havgCosEpHFp2cent[%i];1",i));
    havgSinEpHFp2cent[i] = (TH1D*)avgFile->Get(Form("havgSinEpHFp2cent[%i];1",i));
    havgCosEptrackmid2cent[i] = (TH1D*)avgFile->Get(Form("havgCosEptrackmid2cent[%i];1",i));
    havgSinEptrackmid2cent[i] = (TH1D*)avgFile->Get(Form("havgSinEptrackmid2cent[%i];1",i));
    for (int n=1; n<=flatOrder; n++) {
      avgCosEpcent[i][n-1] = havgCosEpcent[i]->GetBinContent(n);
      avgSinEpcent[i][n-1] = havgSinEpcent[i]->GetBinContent(n);
      avgCosEpHFm2cent[i][n-1] = havgCosEpHFm2cent[i]->GetBinContent(n);
      avgSinEpHFm2cent[i][n-1] = havgSinEpHFm2cent[i]->GetBinContent(n);
      avgCosEpHFp2cent[i][n-1] = havgCosEpHFp2cent[i]->GetBinContent(n);
      avgSinEpHFp2cent[i][n-1] = havgSinEpHFp2cent[i]->GetBinContent(n);
      avgCosEptrackmid2cent[i][n-1] = havgCosEptrackmid2cent[i]->GetBinContent(n);
      avgSinEptrackmid2cent[i][n-1] = havgSinEptrackmid2cent[i]->GetBinContent(n);
    }
  }

  avgFile->Close();
  delete avgFile;

  cout << "Done." << endl;

  //check the recentered averages, they should be zero.
  Double_t avgqxrec = 0;
  Double_t avgqyrec = 0;
  Double_t avgqxHFm2rec = 0;
  Double_t avgqyHFm2rec = 0;
  Double_t avgqxHFp2rec = 0;
  Double_t avgqyHFp2rec = 0;
  Double_t avgqxtrackmid2rec = 0;
  Double_t avgqytrackmid2rec = 0;

  //Calculate average q-vector products for the scalar-product method.
  Double_t avgqAqB = 0;//averaged over all events
  Double_t avgqAqC = 0;
  Double_t avgqBqC = 0;
  Double_t avgqqAUps = 0;//averaged over upsilon 1s candidates

  newfile->cd();

  int ptPASS[numptbins] = {0};
  int yPASS[numybins] = {0};
  int centPASS[numcbins] = {0};

  cout << "Total events = " << nevtReal << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevtReal ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/nevtReal) << "%)" << endl;

    mytree->GetEntry(iev);
    TLeaf *runLeaf = mm->GetLeaf("run");
    TLeaf *lumiLeaf = mm->GetLeaf("lumi");
    TLeaf *eventLeaf = mm->GetLeaf("event");
    TLeaf *vzLeaf = mm->GetLeaf("vz");
    TLeaf *massLeaf = mm->GetLeaf("mass");
    TLeaf *ptLeaf = mm->GetLeaf("pt");
    TLeaf *yLeaf = mm->GetLeaf("y");
    TLeaf *phiLeaf = mm->GetLeaf("phi");
    TLeaf *etaLeaf = mm->GetLeaf("eta");
    TLeaf *pt1Leaf = mm->GetLeaf("pt1");
    TLeaf *phi1Leaf = mm->GetLeaf("phi1");
    TLeaf *eta1Leaf = mm->GetLeaf("eta1");
    TLeaf *pt2Leaf = mm->GetLeaf("pt2");
    TLeaf *phi2Leaf = mm->GetLeaf("phi2");
    TLeaf *eta2Leaf = mm->GetLeaf("eta2");
    TLeaf *ep2Leaf = mm->GetLeaf("ep2");
    TLeaf *dphiEp2Leaf = mm->GetLeaf("dphiEp2");
    TLeaf *cBinLeaf = mm->GetLeaf("cBin");
    TLeaf *qxaLeaf = mm->GetLeaf("qxa");
    TLeaf *qyaLeaf = mm->GetLeaf("qya");
    TLeaf *qxbLeaf = mm->GetLeaf("qxb");
    TLeaf *qybLeaf = mm->GetLeaf("qyb");
    TLeaf *qxcLeaf = mm->GetLeaf("qxc");
    TLeaf *qycLeaf = mm->GetLeaf("qyc");
    TLeaf *qxdLeaf = mm->GetLeaf("qxdimu");
    TLeaf *qydLeaf = mm->GetLeaf("qydimu");

    Double_t qx[21] = {0};
    Double_t qy[21] = {0};

    qx[6] = qxaLeaf->GetValue();
    qy[6] = qyaLeaf->GetValue();
    qx[7] = qxbLeaf->GetValue();
    qy[7] = qybLeaf->GetValue();
    qx[8] = qxcLeaf->GetValue();
    qy[8] = qycLeaf->GetValue();
    qx[9] = qxdLeaf->GetValue();
    qy[9] = qydLeaf->GetValue();

    //Before any cuts, find averages of q-vector products.
    //Nominal: HF2 (8)
    //A: HFm2 (6), B: HFp2 (7), C: trackmid2 (9)
    Double_t qAx = qx[6]; Double_t qAy = qy[6];
    Double_t qBx = qx[7]; Double_t qBy = qy[7];
    Double_t qCx = qx[9]; Double_t qCy = qy[9];
    Double_t qAqB = qAx*qBx-qAy*qBy;
    Double_t qAqC = qAx*qCx-qAy*qCy;
    Double_t qBqC = qBx*qCx-qBy*qCy;
    avgqAqB += qAqB;
    avgqAqC += qAqC;
    avgqBqC += qBqC;
   
      Double_t qxrec = qx[8];
      Double_t qyrec = qy[8];
      Double_t qxHFm2rec = qx[6];
      Double_t qyHFm2rec = qy[6];
      Double_t qxHFp2rec = qx[7];
      Double_t qyHFp2rec = qy[7];
      Double_t qxtrackmid2rec = qx[9];
      Double_t qytrackmid2rec = qy[9];

      Double_t epHFm2rec = atan2(qyHFm2rec,qxHFm2rec)/2;
      Double_t epHFp2rec = atan2(qyHFp2rec,qxHFp2rec)/2;
      Double_t epHF2rec = atan2(qyrec,qxrec)/2;
      Double_t eptrackmid2rec = atan2(qytrackmid2rec,qxtrackmid2rec)/2;

      Double_t epHFm2 = epHFm2rec;
      Double_t epHFp2 = epHFp2rec;
      Double_t epHF2 = epHF2rec;
      Double_t eptrackmid2 = eptrackmid2rec;

      if (abs(epHF2)>5)
        continue;
      if (abs(epHFm2)>5)
        continue;
      if (abs(epHFp2)>5)
        continue;
      if (abs(eptrackmid2)>5)
        continue;

      ALLPASS++;

      //cout << "************epHF2 = " << epHF2 << endl;
      //cout << "************epHFm2 = " << epHFm2 << endl;
      //cout << "************epHFp2 = " << epHFp2 << endl;
      //cout << "************eptrackmid2 = " << eptrackmid2 << endl;

      dm.clear();      // clear the output tree: 
      dm.run = runLeaf->GetValue();
      dm.lumi = lumiLeaf->GetValue();
      dm.event = eventLeaf->GetValue();
      dm.vz = vzLeaf->GetValue();
      dm.cBin = cBinLeaf->GetValue();

      dm.phi = phiLeaf->GetValue();
      dm.ep2 = epHF2;
      dm.dphiEp2 = getDPHI_Jared( dm.phi, dm.ep2);

      dm.qxa = qxHFm2rec;
      dm.qya = qyHFm2rec;
      dm.qxb = qxHFp2rec;
      dm.qyb = qyHFp2rec;
      dm.qxc = qxrec;
      dm.qyc = qyrec;
      dm.qxdimu = qxtrackmid2rec;
      dm.qydimu = qytrackmid2rec;

      dm.mass = massLeaf->GetValue();
      dm.pt = ptLeaf->GetValue();

      dm.y = yLeaf->GetValue();
      dm.eta = etaLeaf->GetValue();
      dm.pt1 = pt1Leaf->GetValue();
      dm.eta1 = eta1Leaf->GetValue();
      dm.phi1 = phi1Leaf->GetValue();
      dm.pt2 = pt2Leaf->GetValue();
      dm.eta2 = eta2Leaf->GetValue();
      dm.phi2 = phi2Leaf->GetValue();
      dm.weight = 1.;

      massHistoUpsilonCandidates->Fill(dm.mass);

      //Find averages of q-vector products.
      Double_t q0x = qx[8]; Double_t q0y = qy[8];
      Double_t qqA = q0x*qAx-q0y*qAy;
      avgqqAUps += qqA;

      //Double_t epHFm2 = epang[6];
      //Double_t epHFp2 = epang[7];
      //Double_t eptrackmid2 = epang[9];

      //flatten event plane angles:
      //epang[0] = HFm1
      //epang[1] = HFp1
      //epang[2] = HF1
      //epang[3] = trackm1
      //epang[4] = trackp1
      //epang[5] = Castor1
      //epHFm2 = HFm2*  -5<eta<-3
      //epHFp2 = HFp2*   3<eta<5
      //epang[8] = HF2
      //eptrackmid2 = trackmid2*   -0.8<eta<0.8
      //epang[10] = trackm2
      //epang[11] = trackp2
      //epang[12] = Castor2
      Double_t DeltaEp2 = 0.0;
      Double_t DeltaEpHFm2 = 0.0;
      Double_t DeltaEpHFp2 = 0.0;
      Double_t DeltaEptrackmid2 = 0.0;
      for (int n=1; n<=flatOrder; n++) {
        DeltaEp2 += -cos(2*n*epHF2)*avgSinEp[n-1]/n;
        DeltaEp2 += sin(2*n*epHF2)*avgCosEp[n-1]/n;
        DeltaEpHFm2 += -cos(2*n*epHFm2)*avgSinEpHFm2[n-1]/n;
        DeltaEpHFm2 += sin(2*n*epHFm2)*avgCosEpHFm2[n-1]/n;
        DeltaEpHFp2 += -cos(2*n*epHFp2)*avgSinEpHFp2[n-1]/n;
        DeltaEpHFp2 += sin(2*n*epHFp2)*avgCosEpHFp2[n-1]/n;
        DeltaEptrackmid2 += -cos(2*n*eptrackmid2)*avgSinEptrackmid2[n-1]/n;
        DeltaEptrackmid2 += sin(2*n*eptrackmid2)*avgCosEptrackmid2[n-1]/n;
      }
      if (!flattenBinByBin) {
        epHF2 = epHF2 + DeltaEp2;
        epHFm2 = epHFm2 + DeltaEpHFm2;
        epHFp2 = epHFp2 + DeltaEpHFp2;
        eptrackmid2 = eptrackmid2 + DeltaEptrackmid2;
      }
      if (epHF2<-pi/2) epHF2 = epHF2+pi;
      if (epHF2>pi/2) epHF2 = epHF2-pi;
      if (epHFm2<-pi/2) epHFm2 = epHFm2+pi;
      if (epHFm2>pi/2) epHFm2 = epHFm2-pi;
      if (epHFp2<-pi/2) epHFp2 = epHFp2+pi;
      if (epHFp2>pi/2) epHFp2 = epHFp2-pi;
      if (eptrackmid2<-pi/2) eptrackmid2 = eptrackmid2+pi;
      if (eptrackmid2>pi/2) eptrackmid2 = eptrackmid2-pi;

      //flatten event plane angles bin by bin in centrality:
      int whichptbin = hpt->FindBin(dm.pt)-1;
      int whichybin = hy->FindBin(dm.y)-1;
      int whichcbin = hcent->FindBin(dm.cBin/2)-1;

      Double_t epHF2cent = epHF2;
      Double_t epHFm2cent = epHFm2;
      Double_t epHFp2cent = epHFp2;
      Double_t eptrackmid2cent = eptrackmid2;
      Double_t DeltaEp2cent = 0.0;
      Double_t DeltaEpHFm2cent = 0.0;
      Double_t DeltaEpHFp2cent = 0.0;
      Double_t DeltaEptrackmid2cent = 0.0;
      for (int n=1; n<=flatOrder; n++) {
        DeltaEp2cent += -cos(2*n*epHF2cent)*avgSinEpcent[whichcbin][n-1]/n;
        DeltaEp2cent += sin(2*n*epHF2cent)*avgCosEpcent[whichcbin][n-1]/n;
        DeltaEpHFm2cent += -cos(2*n*epHFm2cent)*avgSinEpHFm2cent[whichcbin][n-1]/n;
        DeltaEpHFm2cent += sin(2*n*epHFm2cent)*avgCosEpHFm2cent[whichcbin][n-1]/n;
        DeltaEpHFp2cent += -cos(2*n*epHFp2cent)*avgSinEpHFp2cent[whichcbin][n-1]/n;
        DeltaEpHFp2cent += sin(2*n*epHFp2cent)*avgCosEpHFp2cent[whichcbin][n-1]/n;
        DeltaEptrackmid2cent += -cos(2*n*eptrackmid2cent)*avgSinEptrackmid2cent[whichcbin][n-1]/n;
        DeltaEptrackmid2cent += sin(2*n*eptrackmid2cent)*avgCosEptrackmid2cent[whichcbin][n-1]/n;
      }
      epHF2cent = epHF2cent + DeltaEp2cent;
      epHFm2cent = epHFm2cent + DeltaEpHFm2cent;
      epHFp2cent = epHFp2cent + DeltaEpHFp2cent;
      eptrackmid2cent = eptrackmid2cent + DeltaEptrackmid2cent;
      if (epHF2cent<-pi/2) epHF2cent = epHF2cent+pi;
      if (epHF2cent>pi/2) epHF2cent = epHF2cent-pi;
      if (epHFm2cent<-pi/2) epHFm2cent = epHFm2cent+pi;
      if (epHFm2cent>pi/2) epHFm2cent = epHFm2cent-pi;
      if (epHFp2cent<-pi/2) epHFp2cent = epHFp2cent+pi;
      if (epHFp2cent>pi/2) epHFp2cent = epHFp2cent-pi;
      if (eptrackmid2cent<-pi/2) eptrackmid2cent = eptrackmid2cent+pi;
      if (eptrackmid2cent>pi/2) eptrackmid2cent = eptrackmid2cent-pi;

      //Choose flattening method.
      if (flattenBinByBin) {
        epHF2 = epHF2cent;
        epHFm2 = epHFm2cent;
        epHFp2 = epHFp2cent;
        eptrackmid2 = eptrackmid2cent;
      }

      //Fill event plane histograms.
      hEpHF2->Fill(epHF2);
      hEpHFm2->Fill(epHFm2);
      hEpHFp2->Fill(epHFp2);
      hEptrackmid2->Fill(eptrackmid2);

      dm.ep2 = epHF2;
      dm.dphiEp2 = restrictToPiOver2(dm.phi-dm.ep2);

      hphi->Fill(dm.phi);
      hdphi->Fill(dm.dphiEp2);

      Double_t epHF2old = atan2(qy[8],qx[8])/2;
      Double_t epHFm2old = atan2(qy[6],qx[6])/2;
      Double_t epHFp2old = atan2(qy[7],qx[7])/2;
      Double_t eptrackmid2old = atan2(qy[9],qx[9])/2;
      hEpHF2old->Fill(epHF2old);
      hEpHF2new->Fill(epHF2rec);
      hEpHFm2old->Fill(epHFm2old);
      hEpHFm2new->Fill(epHFm2rec);
      hEpHFp2old->Fill(epHFp2old);
      hEpHFp2new->Fill(epHFp2rec);
      hEptrackmid2old->Fill(eptrackmid2old);
      hEptrackmid2new->Fill(eptrackmid2rec);

      //cout << "ALLPASS=" << ALLPASS << " : qxrec = " << qxrec << " : avgqxrec = " << avgqxrec << endl;
      avgqxrec += qxrec;
      //cout << "ALLPASS=" << ALLPASS << " : qxrec = " << qxrec << " : avgqxrec = " << avgqxrec << endl;
      avgqyrec += qyrec;
      avgqxHFm2rec += qxHFm2rec;
      avgqyHFm2rec += qyHFm2rec;
      avgqxHFp2rec += qxHFp2rec;
      avgqyHFp2rec += qyHFp2rec;
      avgqxtrackmid2rec += qxtrackmid2rec;
      avgqytrackmid2rec += qytrackmid2rec;

      //get values to use for event plane resolution.
      Double_t psiA = 0;
      Double_t psiB = 0;
      Double_t psiC = 0;
      if (dm.eta<0.0) {
        psiA = epHFp2;
        psiB = epHFm2;
        psiC = eptrackmid2;
      }
      else {
        psiA = epHFm2;
        psiB = epHFp2;
        psiC = eptrackmid2;
      }

      Double_t deltaAB = psiA-psiB;
      Double_t deltaAC = psiA-psiC;
      Double_t deltaBC = psiB-psiC;
//cout << "deltaAB" << deltaAB << endl;
      avgCosAB += cos(2*(deltaAB));
      avgCosAC += cos(2*(deltaAC));
      avgCosBC += cos(2*(deltaBC));

      //cout << "CosAB = " << cos(2*(deltaAB)) << endl;
      //hCosAC->Fill(cos(2*(deltaAC)));

      sumsqrsCosAB += pow(cos(2*(deltaAB)),2);
      sumsqrsCosAC += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBC += pow(cos(2*(deltaBC)),2);

      //binned averages:
      avgCosABpt[whichptbin] += cos(2*(deltaAB));
      avgCosACpt[whichptbin] += cos(2*(deltaAC));
      avgCosBCpt[whichptbin] += cos(2*(deltaBC));
      sumsqrsCosABpt[whichptbin] += pow(cos(2*(deltaAB)),2);
      sumsqrsCosACpt[whichptbin] += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBCpt[whichptbin] += pow(cos(2*(deltaBC)),2);
      ptPASS[whichptbin] += 1;

      avgCosABy[whichybin] += cos(2*(deltaAB));
      avgCosACy[whichybin] += cos(2*(deltaAC));
      avgCosBCy[whichybin] += cos(2*(deltaBC));
      sumsqrsCosABy[whichybin] += pow(cos(2*(deltaAB)),2);
      sumsqrsCosACy[whichybin] += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBCy[whichybin] += pow(cos(2*(deltaBC)),2);
      yPASS[whichybin] += 1;

      avgCosABcent[whichcbin] += cos(2*(deltaAB));
      avgCosACcent[whichcbin] += cos(2*(deltaAC));
      avgCosBCcent[whichcbin] += cos(2*(deltaBC));
      sumsqrsCosABcent[whichcbin] += pow(cos(2*(deltaAB)),2);
      sumsqrsCosACcent[whichcbin] += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBCcent[whichcbin] += pow(cos(2*(deltaBC)),2);
      centPASS[whichcbin] += 1;

      mmtree->Fill();

  } //end of event loop

  mmtree->Write();  // Don't need to call Write() for trees
  newfile->Close();

  //Averaged for scalar product method:
  TFile* avgSPFile = new TFile(Form("averages/avgSPFile_n%i_%i.root",nevt,dateStr),"RECREATE");
  TH1D* havgSP = new TH1D("havgSP","avgSP",4,0,4);
  cout << "sumqqAUps = " << avgqqAUps << endl;
  cout << "sumqAqB = " << avgqAqB << endl;
  cout << "sumqAqC = " << avgqAqC << endl;
  cout << "sumqBqC = " << avgqBqC << endl;
  avgqAqB = avgqAqB/nevt;
  avgqAqC = avgqAqC/nevt;
  avgqBqC = avgqBqC/nevt;
  avgqqAUps = avgqqAUps/ALLPASS;
  cout << "avgqqAUps = " << avgqqAUps << endl;
  cout << "avgqAqB = " << avgqAqB << endl;
  cout << "avgqAqC = " << avgqAqC << endl;
  cout << "avgqBqC = " << avgqBqC << endl;
  havgSP->SetBinContent(1,avgqqAUps);
  havgSP->SetBinContent(2,avgqAqB);
  havgSP->SetBinContent(3,avgqAqC);
  havgSP->SetBinContent(4,avgqBqC);
  havgSP->Write();
  TH1D* hv2SP = new TH1D("hv2SP","v2SP",1,0,1);
  Double_t v2SP = avgqqAUps/sqrt(avgqAqB*avgqAqC/avgqBqC);
  cout << "v2SP = " << v2SP << endl;
  hv2SP->SetBinContent(1,v2SP);
  hv2SP->Write();
  avgSPFile->Close();

  //include the averages used to calculate the event plane resolution correction.
  TFile* resCorFile = new TFile(Form("averages/resCorFile_n%i_%s.root",nevt,binByBinStr.Data()),"RECREATE");
  TH1D* hRint = new TH1D("hRint","EP Resolution factor",1,0,1);
  avgCosAB = avgCosAB/ALLPASS;
  avgCosAC = avgCosAC/ALLPASS;
  avgCosBC = avgCosBC/ALLPASS;
  Double_t rmsCosAB = sqrt(sumsqrsCosAB/ALLPASS - pow(avgCosAB,2));
  Double_t rmsCosAC = sqrt(sumsqrsCosAC/ALLPASS - pow(avgCosAC,2));
  Double_t rmsCosBC = sqrt(sumsqrsCosBC/ALLPASS - pow(avgCosBC,2));
  cout << "avg(Cos(2*(deltaAB))) = " << avgCosAB << " +/- " << rmsCosAB << endl;
  cout << "avg(Cos(2*(deltaAC))) = " << avgCosAC << " +/- " << rmsCosAC << endl;
  cout << "avg(Cos(2*(deltaBC))) = " << avgCosBC << " +/- " << rmsCosBC << endl;
  rmsCosAB = rmsCosAB/sqrt(ALLPASS);
  rmsCosAC = rmsCosAC/sqrt(ALLPASS);
  rmsCosBC = rmsCosBC/sqrt(ALLPASS);
  Double_t RA = sqrt(avgCosAB*avgCosAC/avgCosBC);
  Double_t RAerr = 0.5*RA*sqrt(pow(rmsCosAB/avgCosAB,2) + pow(rmsCosAC/avgCosAC,2) + pow(rmsCosBC/avgCosBC,2));
  cout << "Event plane resolution factor = " << RA << " +/- " << RAerr << endl;
  hRint->SetBinContent(1,RA);
  hRint->SetBinError(1,RAerr);
  hRint->Write();

  cout << "filling pt histogram" << endl;
  for (int i=0; i<numptbins; i++) {
    avgCosABpt[i] = avgCosABpt[i]/ptPASS[i];
    avgCosACpt[i] = avgCosACpt[i]/ptPASS[i];
    avgCosBCpt[i] = avgCosBCpt[i]/ptPASS[i];
    Double_t rmsCosABpt = sqrt(sumsqrsCosABpt[i]/ptPASS[i] - pow(avgCosABpt[i],2));
    Double_t rmsCosACpt = sqrt(sumsqrsCosACpt[i]/ptPASS[i] - pow(avgCosACpt[i],2));
    Double_t rmsCosBCpt = sqrt(sumsqrsCosBCpt[i]/ptPASS[i] - pow(avgCosBCpt[i],2));
    rmsCosABpt = rmsCosABpt/sqrt(ptPASS[i]);
    rmsCosACpt = rmsCosACpt/sqrt(ptPASS[i]);
    rmsCosBCpt = rmsCosBCpt/sqrt(ptPASS[i]);
    Double_t RApt = sqrt(avgCosABpt[i]*avgCosACpt[i]/avgCosBCpt[i]);
    Double_t RApterr = 0.5*RApt*sqrt(pow(rmsCosABpt/avgCosABpt[i],2) + pow(rmsCosACpt/avgCosACpt[i],2) + pow(rmsCosBCpt/avgCosBCpt[i],2));
    hRpt->SetBinContent(i+1,RApt);
    hRpt->SetBinError(i+1,RApterr);
  }
  hRpt->Write();

  cout << "filling y histogram" << endl;
  for (int i=0; i<numybins; i++) {
    avgCosABy[i] = avgCosABy[i]/yPASS[i];
    avgCosACy[i] = avgCosACy[i]/yPASS[i];
    avgCosBCy[i] = avgCosBCy[i]/yPASS[i];
    Double_t rmsCosABy = sqrt(sumsqrsCosABy[i]/yPASS[i] - pow(avgCosABy[i],2));
    Double_t rmsCosACy = sqrt(sumsqrsCosACy[i]/yPASS[i] - pow(avgCosACy[i],2));
    Double_t rmsCosBCy = sqrt(sumsqrsCosBCy[i]/yPASS[i] - pow(avgCosBCy[i],2));
    rmsCosABy = rmsCosABy/sqrt(yPASS[i]);
    rmsCosACy = rmsCosACy/sqrt(yPASS[i]);
    rmsCosBCy = rmsCosBCy/sqrt(yPASS[i]);
    Double_t RAy = sqrt(avgCosABy[i]*avgCosACy[i]/avgCosBCy[i]);
    Double_t RAyerr = 0.5*RAy*sqrt(pow(rmsCosABy/avgCosABy[i],2) + pow(rmsCosACy/avgCosACy[i],2) + pow(rmsCosBCy/avgCosBCy[i],2));
    hRy->SetBinContent(i+1,RAy);
    hRy->SetBinError(i+1,RAyerr);
  }
  hRy->Write();

  cout << "filling cent histogram" << endl;
  for (int i=0; i<numcbins; i++) {
    avgCosABcent[i] = avgCosABcent[i]/centPASS[i];
    avgCosACcent[i] = avgCosACcent[i]/centPASS[i];
    avgCosBCcent[i] = avgCosBCcent[i]/centPASS[i];
    Double_t rmsCosABcent = sqrt(sumsqrsCosABcent[i]/centPASS[i] - pow(avgCosABcent[i],2));
    Double_t rmsCosACcent = sqrt(sumsqrsCosACcent[i]/centPASS[i] - pow(avgCosACcent[i],2));
    Double_t rmsCosBCcent = sqrt(sumsqrsCosBCcent[i]/centPASS[i] - pow(avgCosBCcent[i],2));
    rmsCosABcent = rmsCosABcent/sqrt(centPASS[i]);
    rmsCosACcent = rmsCosACcent/sqrt(centPASS[i]);
    rmsCosBCcent = rmsCosBCcent/sqrt(centPASS[i]);
    Double_t RAcent = sqrt(avgCosABcent[i]*avgCosACcent[i]/avgCosBCcent[i]);
    Double_t RAcenterr = 0.5*RAcent*sqrt(pow(rmsCosABcent/avgCosABcent[i],2) + pow(rmsCosACcent/avgCosACcent[i],2) + pow(rmsCosBCcent/avgCosBCcent[i],2));
    hRc->SetBinContent(i+1,RAcent);
    hRc->SetBinError(i+1,RAcenterr);
  }
  hRc->Write();

  hRpt->Sumw2();
  hRy->Sumw2();
  hRc->Sumw2();

  TCanvas* c1 = new TCanvas("c1","c1",0,0,900,300);
  c1->Divide(3);
  c1->cd(1);
  hRpt->Draw();
  c1->cd(2);
  hRy->Draw();
  c1->cd(3);
  hRc->Draw();

  hEpHF2->Write();
  hEpHFm2->Write();
  hEpHFp2->Write();
  hEptrackmid2->Write();

  TCanvas* c2 = new TCanvas("c2","c2",0,0,800,800);
  c2->Divide(2,2);
  c2->cd(1);
  hEpHF2->Draw();
  c2->cd(2);
  hEpHFm2->Draw();
  c2->cd(3);
  hEpHFp2->Draw();
  c2->cd(4);
  hEptrackmid2->Draw();

  c2->SaveAs(Form("plots/epHistos%s_n%i_%i.pdf",binByBinStr.Data(),nevt,dateStr));
  c2->SaveAs(Form("plots/epHistos%s_n%i_%i.png",binByBinStr.Data(),nevt,dateStr));

  TCanvas* c3 = new TCanvas("c3","c3",0,0,400,400);
  c3->cd();
  hEpHF2old->SetTitle("Event plane HF2");
  hEpHF2old->GetXaxis()->SetTitle("#psi");
  hEpHF2old->SetLineWidth(2);
  hEpHF2old->SetLineStyle(1);
  hEpHF2old->Draw();
  hEpHF2new->SetLineColor(2);
  hEpHF2new->SetLineWidth(2);
  hEpHF2new->SetLineStyle(7);
  hEpHF2new->Draw("same");
  hEpHF2->SetLineColor(3);
  hEpHF2->SetLineWidth(1);
  hEpHF2->SetLineStyle(1);
  hEpHF2->Draw("same");
  TLegend* legc3 = new TLegend(0.4,0.15,0.6,0.3); legc3->SetTextSize(12);
  legc3->SetTextFont(43);
  legc3->SetBorderSize(0);
  legc3->AddEntry(hEpHF2old,"ep raw","l");
  legc3->AddEntry(hEpHF2new,"ep recentered","l");
  legc3->AddEntry(hEpHF2,"ep rec+flattened","l");
  legc3->Draw("same");
  c3->SaveAs(Form("plots/EventPlaneHF2Change_n%i.pdf",nevt));
  c3->SaveAs(Form("plots/EventPlaneHF2Change_n%i.png",nevt));
  c3->SaveAs(Form("plots/EventPlaneHF2Change_n%i.C",nevt));

  TCanvas* c4 = new TCanvas("c4","c4",0,0,400,400);
  c4->cd();
  hEpHFm2old->SetTitle("Event plane HFm2");
  hEpHFm2old->SetLineWidth(2);
  hEpHFm2old->SetLineStyle(1);
  hEpHFm2old->Draw();
  hEpHFm2new->SetLineColor(2);
  hEpHFm2new->SetLineWidth(2);
  hEpHFm2new->SetLineStyle(7);
  hEpHFm2new->Draw("same");
  hEpHFm2->SetLineColor(3);
  hEpHFm2->SetLineWidth(1);
  hEpHFm2->SetLineStyle(1);
  hEpHFm2->Draw("same");
  TLegend* legc4 = new TLegend(0.4,0.15,0.6,0.3); legc4->SetTextSize(12);
  legc4->SetTextFont(43);
  legc4->SetBorderSize(0);
  legc4->AddEntry(hEpHFm2old,"ep raw","l");
  legc4->AddEntry(hEpHFm2new,"ep recentered","l");
  legc4->AddEntry(hEpHFm2,"ep rec+flattened","l");
  legc4->Draw("same");
  c4->SaveAs(Form("plots/EventPlaneHFm2Change_n%i.pdf",nevt));
  c4->SaveAs(Form("plots/EventPlaneHFm2Change_n%i.png",nevt));

  TCanvas* c5 = new TCanvas("c5","c5",0,0,400,400);
  c5->cd();
  hEpHFp2old->SetTitle("Event plane HFp2");
  hEpHFp2old->GetXaxis()->SetTitle("#psi");
  hEpHFp2old->SetLineWidth(2);
  hEpHFp2old->SetLineStyle(1);
  hEpHFp2old->Draw();
  hEpHFp2new->SetLineColor(2);
  hEpHFp2new->SetLineWidth(2);
  hEpHFp2new->SetLineStyle(7);
  hEpHFp2new->Draw("same");
  hEpHFp2->SetLineColor(3);
  hEpHFp2->SetLineWidth(1);
  hEpHFp2->SetLineStyle(1);
  hEpHFp2->Draw("same");
  TLegend* legc5 = new TLegend(0.4,0.15,0.6,0.3); legc5->SetTextSize(12);
  legc5->SetTextFont(43);
  legc5->SetBorderSize(0);
  legc5->AddEntry(hEpHFp2old,"ep raw","l");
  legc5->AddEntry(hEpHFp2new,"ep recentered","l");
  legc5->AddEntry(hEpHFp2,"ep rec+flattened","l");
  legc5->Draw("same");
  c5->SaveAs(Form("plots/EventPlaneHFp2Change_n%i.pdf",nevt));
  c5->SaveAs(Form("plots/EventPlaneHFp2Change_n%i.png",nevt));

  TCanvas* c6 = new TCanvas("c6","c6",0,0,400,400);
  c6->cd();
  hEptrackmid2old->SetTitle("Event plane trackmid2");
  hEptrackmid2old->GetXaxis()->SetTitle("#psi");
  hEptrackmid2old->SetLineWidth(2);
  hEptrackmid2old->SetLineStyle(1);
  hEptrackmid2old->Draw();
  hEptrackmid2new->SetLineColor(2);
  hEptrackmid2new->SetLineWidth(2);
  hEptrackmid2new->SetLineStyle(7);
  hEptrackmid2new->Draw("same");
  hEptrackmid2->SetLineColor(3);
  hEptrackmid2->SetLineWidth(1);
  hEptrackmid2->SetLineStyle(1);
  hEptrackmid2->Draw("same");
  TLegend* legc6 = new TLegend(0.4,0.15,0.6,0.3); legc6->SetTextSize(12);
  legc6->SetTextFont(43);
  legc6->SetBorderSize(0);
  legc6->AddEntry(hEptrackmid2old,"ep raw","l");
  legc6->AddEntry(hEptrackmid2new,"ep recentered","l");
  legc6->AddEntry(hEptrackmid2,"ep rec+flattened","l");
  legc6->Draw("same");
  c6->SaveAs(Form("plots/EventPlanetrackmid2Change_n%i.pdf",nevt));
  c6->SaveAs(Form("plots/EventPlanetrackmid2Change_n%i.png",nevt));

  resCorFile->Close();

  TCanvas* c7 = new TCanvas("c7","c7",0,0,400,400);
  c7->cd();
  //hdphiw->Draw();
  //hdphiw->SetLineColor(2);
  hdphi->Draw();

  TCanvas* c8 = new TCanvas("c8","c8",400,0,400,400);
  c8->cd();
  hphi->Draw();

  TCanvas* c9 = new TCanvas("c9","c9",0,0,400,400);
  c9->cd();
  //massHistoUpsilonCandidates->Scale(1,"width");
  massHistoUpsilonCandidates->SetFillColor(kYellow);
  massHistoUpsilonCandidates->SetTitle("Upsilon candidates");
  massHistoUpsilonCandidates->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  massHistoUpsilonCandidates->GetXaxis()->SetTitleOffset(1.05) ;
  massHistoUpsilonCandidates->GetXaxis()->CenterTitle();
  massHistoUpsilonCandidates->SetMinimum(0);
  massHistoUpsilonCandidates->Draw("hist");
  TString plotLabel = "PromptAOD_n1n2_pt0_SoftId_Vtx_HLTevt_HLTdm_NOTscaledByWidth";
  c9->SaveAs(Form("AllDimuons_%s.pdf",plotLabel.Data()));
  c9->SaveAs(Form("AllDimuons_%s.png",plotLabel.Data()));

  cout << "ptPASS = {";
  for (int i=0; i<numptbins; i++) {
    cout << ptPASS[i] << ",";
  }
  cout << "}" << endl;
  cout << "yPASS = {";
  for (int i=0; i<numybins; i++) {
    cout << yPASS[i] << ",";
  }
  cout << "}" << endl;
  cout << "centPASS = {";
  for (int i=0; i<numcbins; i++) {
    cout << centPASS[i] << ",";
  }
  cout << "}" << endl;
  cout << "ALLPASS = " << ALLPASS << endl;

  cout << endl << "All these averages should be zero after recentering:" << endl;
  avgqxrec = avgqxrec/ALLPASS;
  avgqyrec = avgqyrec/ALLPASS;
  cout << "avgqx-recentered = " << avgqxrec << endl;
  cout << "avgqy-recentered = " << avgqyrec << endl;
  avgqxHFm2rec = avgqxHFm2rec/ALLPASS;
  avgqyHFm2rec = avgqyHFm2rec/ALLPASS;
  cout << "avgqxHFm2-recentered = " << avgqxHFm2rec << endl;
  cout << "avgqyHFm2-recentered = " << avgqyHFm2rec << endl;
  avgqxHFp2rec = avgqxHFp2rec/ALLPASS;
  avgqyHFp2rec = avgqyHFp2rec/ALLPASS;
  cout << "avgqxHFp2-recentered = " << avgqxHFp2rec << endl;
  cout << "avgqyHFp2-recentered = " << avgqyHFp2rec << endl;
  avgqxtrackmid2rec = avgqxtrackmid2rec/ALLPASS;
  avgqytrackmid2rec = avgqytrackmid2rec/ALLPASS;
  cout << "avgqxtrackmid2-recentered = " << avgqxtrackmid2rec << endl;
  cout << "avgqytrackmid2-recentered = " << avgqytrackmid2rec << endl;
} 

