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


static const long MAXTREESIZE = 1000000000000;

void SkimMMTree_recentering_GetAverageEP_nords(int nevt=-1, int dateStr=20201118) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  float ptbins[5] = {0,3,6,10,50};
  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  Double_t cbins[4] = {10,30,50,90};
  const int numcbins = sizeof(cbins)/sizeof(double)-1;

  TH1D* hpt = new TH1D("hpt","hist vs pt",numptbins,ptbins);
  TH1D* hcent = new TH1D("hcent","hist vs cent",numcbins,cbins);

  gStyle->SetOptStat(0);
  //TH1D* hCosAC = new TH1D("hCosAC","cos(2*(psiA-psiC))",50,-1.2,1.2);
  TH1D* hRpt = new TH1D("hRpt","EP Resolution factor vs pt",numptbins,ptbins);
  TH1D* hRcent = new TH1D("hRcent","EP Resolution factor vs cent",numcbins,cbins);

  //Histograms to contain the distributions of qx and qy:
  TH1D* hqxold = new TH1D("hqxold","hqxold",100,-100,100);
  TH1D* hqxnew = new TH1D("hqxnew","hqxnew",100,-100,100);
  TH1D* hqyold = new TH1D("hqyold","hqyold",100,-100,100);
  TH1D* hqynew = new TH1D("hqynew","hqynew",100,-100,100);

  //Event planes
  TH1D* hepHF2old = new TH1D("hepHF2old","hepHF2old",50,-2,2);
  TH1D* hepHF2new = new TH1D("hepHF2new","hepHF2new",50,-2,2);
  TH1D* hepHFm2old = new TH1D("hepHFm2old","hepHFm2old",50,-2,2);
  TH1D* hepHFm2new = new TH1D("hepHFm2new","hepHFm2new",50,-2,2);
  TH1D* hepHFp2old = new TH1D("hepHFp2old","hepHFp2old",50,-2,2);
  TH1D* hepHFp2new = new TH1D("hepHFp2new","hepHFp2new",50,-2,2);
  TH1D* heptrackmid2old = new TH1D("heptrackmid2old","heptrackmid2old",50,-2,2);
  TH1D* heptrackmid2new = new TH1D("heptrackmid2new","heptrackmid2new",50,-2,2);

  TString inFileName = Form("skims/newOniaTree_Skim_UpsTrig_RD_RAW_n%i_20201118.root",nevt);
  if (nevt==-1) inFileName = "skims/newOniaTree_Skim_UpsTrig_RD_RAW_n-1_20201112.root";
  TFile* inFile = TFile::Open(inFileName,"READ");

  TTree* mytree = (TTree*)inFile->Get("mmep");
  TBranch* mm = (TBranch*)mytree->GetBranch("mm");

  TFile* newfile;
  newfile = new TFile(Form("skims/newOniaTree_Skim_UpsTrig_MM_RECENTERED_nords_n%i_%i.root",nevt,dateStr),"recreate");

  DiMuon dm;
  TTree* mmtree = new TTree("mmep","dimuonAndEventPlanes");
  mmtree->SetMaxTreeSize(MAXTREESIZE);
  mmtree->Branch("mm",&dm.run,branchString.Data());

  int ALLPASS = 0;

  if(nevt == -1) nevt = 56114317;

  int nevtReal = mytree->GetEntries();

  //averages for flattening
  const int flatOrder = 21;
  Double_t avgCosEp[flatOrder] = {0};
  Double_t avgSinEp[flatOrder] = {0};
  Double_t avgCosEpHFm2[flatOrder] = {0};
  Double_t avgSinEpHFm2[flatOrder] = {0};
  Double_t avgCosEpHFp2[flatOrder] = {0};
  Double_t avgSinEpHFp2[flatOrder] = {0};
  Double_t avgCosEptrackmid2[flatOrder] = {0};
  Double_t avgSinEptrackmid2[flatOrder] = {0};

  Double_t avgCosEppt[numptbins][flatOrder] = {0};
  Double_t avgSinEppt[numptbins][flatOrder] = {0};
  Double_t avgCosEpHFm2pt[numptbins][flatOrder] = {0};
  Double_t avgSinEpHFm2pt[numptbins][flatOrder] = {0};
  Double_t avgCosEpHFp2pt[numptbins][flatOrder] = {0};
  Double_t avgSinEpHFp2pt[numptbins][flatOrder] = {0};
  Double_t avgCosEptrackmid2pt[numptbins][flatOrder] = {0};
  Double_t avgSinEptrackmid2pt[numptbins][flatOrder] = {0};

  Double_t avgCosEpcent[numcbins][flatOrder] = {0};
  Double_t avgSinEpcent[numcbins][flatOrder] = {0};
  Double_t avgCosEpHFm2cent[numcbins][flatOrder] = {0};
  Double_t avgSinEpHFm2cent[numcbins][flatOrder] = {0};
  Double_t avgCosEpHFp2cent[numcbins][flatOrder] = {0};
  Double_t avgSinEpHFp2cent[numcbins][flatOrder] = {0};
  Double_t avgCosEptrackmid2cent[numcbins][flatOrder] = {0};
  Double_t avgSinEptrackmid2cent[numcbins][flatOrder] = {0};

  //averages for (unflattened) resolution correction
  Double_t avgCosAB = 0;
  Double_t avgCosAC = 0;
  Double_t avgCosBC = 0;
  Double_t avgCosABpt[numptbins] = {0};
  Double_t avgCosACpt[numptbins] = {0};
  Double_t avgCosBCpt[numptbins] = {0};
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
  Double_t sumsqrsCosABcent[numcbins] = {0};
  Double_t sumsqrsCosACcent[numcbins] = {0};
  Double_t sumsqrsCosBCcent[numcbins] = {0};

  //get values for re-centering
  cout << "getting values for re-centering..." << endl;
  //TFile* avgQFile = TFile::Open("avgQFile_fulldataset_2020_05_20.root","READ");
  TFile* avgQFile = TFile::Open(Form("averages/avgQFile%i.root",nevt),"READ");
  TH1D* havgqx = (TH1D*)avgQFile->Get("hqx;1");
  TH1D* havgqy = (TH1D*)avgQFile->Get("hqy;1");
  TH1D* havgqxHFm2 = (TH1D*)avgQFile->Get("hqxHFm2;1");
  TH1D* havgqyHFm2 = (TH1D*)avgQFile->Get("hqyHFm2;1");
  TH1D* havgqxHFp2 = (TH1D*)avgQFile->Get("hqxHFp2;1");
  TH1D* havgqyHFp2 = (TH1D*)avgQFile->Get("hqyHFp2;1");
  TH1D* havgqxtrackmid2 = (TH1D*)avgQFile->Get("hqxtrackmid2;1");
  TH1D* havgqytrackmid2 = (TH1D*)avgQFile->Get("hqytrackmid2;1");
  Double_t avgqx = havgqx->GetBinContent(1);
  Double_t avgqy = havgqy->GetBinContent(1);
  Double_t avgqxHFm2 = havgqxHFm2->GetBinContent(1);
  Double_t avgqyHFm2 = havgqyHFm2->GetBinContent(1);
  Double_t avgqxHFp2 = havgqxHFp2->GetBinContent(1);
  Double_t avgqyHFp2 = havgqyHFp2->GetBinContent(1);
  Double_t avgqxtrackmid2 = havgqxtrackmid2->GetBinContent(1);
  Double_t avgqytrackmid2 = havgqytrackmid2->GetBinContent(1);
  cout << "avg(qx) = " << avgqx << endl;
  cout << "avg(qy) = " << avgqy << endl;

  TH1D* havgqxpt = (TH1D*)avgQFile->Get("hqxpt;1");
  TH1D* havgqypt = (TH1D*)avgQFile->Get("hqypt;1");
  TH1D* havgqxHFm2pt = (TH1D*)avgQFile->Get("hqxHFm2pt;1");
  TH1D* havgqyHFm2pt = (TH1D*)avgQFile->Get("hqyHFm2pt;1");
  TH1D* havgqxHFp2pt = (TH1D*)avgQFile->Get("hqxHFp2pt;1");
  TH1D* havgqyHFp2pt = (TH1D*)avgQFile->Get("hqyHFp2pt;1");
  TH1D* havgqxtrackmid2pt = (TH1D*)avgQFile->Get("hqxtrackmid2pt;1");
  TH1D* havgqytrackmid2pt = (TH1D*)avgQFile->Get("hqytrackmid2pt;1");
  Double_t avgqxpt[numptbins];
  Double_t avgqypt[numptbins];
  Double_t avgqxHFm2pt[numptbins];
  Double_t avgqyHFm2pt[numptbins];
  Double_t avgqxHFp2pt[numptbins];
  Double_t avgqyHFp2pt[numptbins];
  Double_t avgqxtrackmid2pt[numptbins];
  Double_t avgqytrackmid2pt[numptbins];
  for (int i=0; i<numptbins; i++) {
    avgqxpt[i] = havgqxpt->GetBinContent(i+1);
    avgqypt[i] = havgqypt->GetBinContent(i+1);
    avgqxHFm2pt[i] = havgqxHFm2pt->GetBinContent(i+1);
    avgqyHFm2pt[i] = havgqyHFm2pt->GetBinContent(i+1);
    avgqxHFp2pt[i] = havgqxHFp2pt->GetBinContent(i+1);
    avgqyHFp2pt[i] = havgqyHFp2pt->GetBinContent(i+1);
    avgqxtrackmid2pt[i] = havgqxtrackmid2pt->GetBinContent(i+1);
    avgqytrackmid2pt[i] = havgqytrackmid2pt->GetBinContent(i+1);
  }

  TH1D* havgqxcent = (TH1D*)avgQFile->Get("hqxcent;1");
  TH1D* havgqycent = (TH1D*)avgQFile->Get("hqycent;1");
  TH1D* havgqxHFm2cent = (TH1D*)avgQFile->Get("hqxHFm2cent;1");
  TH1D* havgqyHFm2cent = (TH1D*)avgQFile->Get("hqyHFm2cent;1");
  TH1D* havgqxHFp2cent = (TH1D*)avgQFile->Get("hqxHFp2cent;1");
  TH1D* havgqyHFp2cent = (TH1D*)avgQFile->Get("hqyHFp2cent;1");
  TH1D* havgqxtrackmid2cent = (TH1D*)avgQFile->Get("hqxtrackmid2cent;1");
  TH1D* havgqytrackmid2cent = (TH1D*)avgQFile->Get("hqytrackmid2cent;1");
  Double_t avgqxcent[numcbins];
  Double_t avgqycent[numcbins];
  Double_t avgqxHFm2cent[numcbins];
  Double_t avgqyHFm2cent[numcbins];
  Double_t avgqxHFp2cent[numcbins];
  Double_t avgqyHFp2cent[numcbins];
  Double_t avgqxtrackmid2cent[numcbins];
  Double_t avgqytrackmid2cent[numcbins];
  for (int i=0; i<numcbins; i++) {
    avgqxcent[i] = havgqxcent->GetBinContent(i+1);
    avgqycent[i] = havgqycent->GetBinContent(i+1);
    avgqxHFm2cent[i] = havgqxHFm2cent->GetBinContent(i+1);
    avgqyHFm2cent[i] = havgqyHFm2cent->GetBinContent(i+1);
    avgqxHFp2cent[i] = havgqxHFp2cent->GetBinContent(i+1);
    avgqyHFp2cent[i] = havgqyHFp2cent->GetBinContent(i+1);
    avgqxtrackmid2cent[i] = havgqxtrackmid2cent->GetBinContent(i+1);
    avgqytrackmid2cent[i] = havgqytrackmid2cent->GetBinContent(i+1);
  }
  avgQFile->Close();
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

  newfile->cd();

  int ptPASS[numptbins] = {0};
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

      Double_t qxrec = qx[8]-avgqx;
      Double_t qyrec = qy[8]-avgqy;
      Double_t qxHFm2rec = qx[6]-avgqxHFm2;
      Double_t qyHFm2rec = qy[6]-avgqyHFm2;
      Double_t qxHFp2rec = qx[7]-avgqxHFp2;
      Double_t qyHFp2rec = qy[7]-avgqyHFp2;
      Double_t qxtrackmid2rec = qx[9]-avgqxtrackmid2;
      Double_t qytrackmid2rec = qy[9]-avgqytrackmid2;

      Double_t epHFm2 = atan2(qyHFm2rec,qxHFm2rec)/2;
      Double_t epHFp2 = atan2(qyHFp2rec,qxHFp2rec)/2;
      Double_t epHF2 = atan2(qyrec,qxrec)/2;
      Double_t eptrackmid2 = atan2(qytrackmid2rec,qxtrackmid2rec)/2;

      if (abs(epHF2)>5)
        continue;
      if (abs(epHFm2)>5)
        continue;
      if (abs(epHFp2)>5)
        continue;
      if (abs(eptrackmid2)>5)
        continue;

      ALLPASS++;

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

      hqxold->Fill(qx[8]);
      hqyold->Fill(qy[8]);
      hqxnew->Fill(qxrec);
      hqynew->Fill(qyrec);

      Double_t epHF2old = atan2(qy[8],qx[8])/2;
      Double_t epHFm2old = atan2(qy[6],qx[6])/2;
      Double_t epHFp2old = atan2(qy[7],qx[7])/2;
      Double_t eptrackmid2old = atan2(qy[9],qx[9])/2;
      hepHF2old->Fill(epHF2old);
      hepHF2new->Fill(epHF2);
      hepHFm2old->Fill(epHFm2old);
      hepHFm2new->Fill(epHFm2);
      hepHFp2old->Fill(epHFp2old);
      hepHFp2new->Fill(epHFp2);
      heptrackmid2old->Fill(eptrackmid2old);
      heptrackmid2new->Fill(eptrackmid2);

      avgqxrec += qxrec;
      avgqyrec += qyrec;
      avgqxHFm2rec += qxHFm2rec;
      avgqyHFm2rec += qyHFm2rec;
      avgqxHFp2rec += qxHFp2rec;
      avgqyHFp2rec += qyHFp2rec;
      avgqxtrackmid2rec += qxtrackmid2rec;
      avgqytrackmid2rec += qytrackmid2rec;

      //get average values for flattening later on.
      for (int n=1; n<=flatOrder; n++) {
        avgCosEp[n-1] += cos(2*n*dm.ep2);
        avgSinEp[n-1] += sin(2*n*dm.ep2);
        avgCosEpHFm2[n-1] += cos(2*n*epHFm2);
        avgSinEpHFm2[n-1] += sin(2*n*epHFm2);
        avgCosEpHFp2[n-1] += cos(2*n*epHFp2);
        avgSinEpHFp2[n-1] += sin(2*n*epHFp2);
        avgCosEptrackmid2[n-1] += cos(2*n*eptrackmid2);
        avgSinEptrackmid2[n-1] += sin(2*n*eptrackmid2);
      }

      //get values to use for event plane resolution.
      //epang[0] = HFm1
      //epang[1] = HFp1
      //epang[2] = HF1
      //epang[3] = trackm1
      //epang[4] = trackp1
      //epang[5] = Castor1
      //epHFm2 = HFm2*  -5<eta<-3
      //epHFp2 = HFp2*   3<eta<5
      //epHF2 = HF2
      //eptrackmid2 = trackmid2*   -0.8<eta<0.8
      //epang[10] = trackm2
      //epang[11] = trackp2
      //epang[12] = Castor2
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

      avgCosAB += cos(2*(deltaAB));
      avgCosAC += cos(2*(deltaAC));
      avgCosBC += cos(2*(deltaBC));

      //cout << "CosAB = " << cos(2*(deltaAB)) << endl;
      //hCosAC->Fill(cos(2*(deltaAC)));

      sumsqrsCosAB += pow(cos(2*(deltaAB)),2);
      sumsqrsCosAC += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBC += pow(cos(2*(deltaBC)),2);

      //binned averages:(without binned re-centering.)
      int whichptbin = hpt->FindBin(dm.pt)-1;
      int whichcbin = hcent->FindBin(dm.cBin/2)-1;

      for (int n=1; n<=flatOrder; n++) {
        avgCosEppt[whichptbin][n-1] += cos(2*n*dm.ep2);
        avgSinEppt[whichptbin][n-1] += sin(2*n*dm.ep2);
        avgCosEpHFm2pt[whichptbin][n-1] += cos(2*n*epHFm2);
        avgSinEpHFm2pt[whichptbin][n-1] += sin(2*n*epHFm2);
        avgCosEpHFp2pt[whichptbin][n-1] += cos(2*n*epHFp2);
        avgSinEpHFp2pt[whichptbin][n-1] += sin(2*n*epHFp2);
        avgCosEptrackmid2pt[whichptbin][n-1] += cos(2*n*eptrackmid2);
        avgSinEptrackmid2pt[whichptbin][n-1] += sin(2*n*eptrackmid2);
      }
      for (int n=1; n<=flatOrder; n++) {
        avgCosEpcent[whichcbin][n-1] += cos(2*n*dm.ep2);
        avgSinEpcent[whichcbin][n-1] += sin(2*n*dm.ep2);
        avgCosEpHFm2cent[whichcbin][n-1] += cos(2*n*epHFm2);
        avgSinEpHFm2cent[whichcbin][n-1] += sin(2*n*epHFm2);
        avgCosEpHFp2cent[whichcbin][n-1] += cos(2*n*epHFp2);
        avgSinEpHFp2cent[whichcbin][n-1] += sin(2*n*epHFp2);
        avgCosEptrackmid2cent[whichcbin][n-1] += cos(2*n*eptrackmid2);
        avgSinEptrackmid2cent[whichcbin][n-1] += sin(2*n*eptrackmid2);
      }


      avgCosABpt[whichptbin] += cos(2*(deltaAB));
      avgCosACpt[whichptbin] += cos(2*(deltaAC));
      avgCosBCpt[whichptbin] += cos(2*(deltaBC));
      sumsqrsCosABpt[whichptbin] += pow(cos(2*(deltaAB)),2);
      sumsqrsCosACpt[whichptbin] += pow(cos(2*(deltaAC)),2);
      sumsqrsCosBCpt[whichptbin] += pow(cos(2*(deltaBC)),2);
      ptPASS[whichptbin] += 1;

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

  //TCanvas* c2 = new TCanvas("c2","c2",0,0,400,400);
  //hCosAC->Draw();
  //hCosAC->SaveAs("hCosAC%i.pdf",nevt));

  //Save file with averages for flattening.
  TFile* avgFile = new TFile(Form("averages/avgEpFile%i.root",nevt),"RECREATE");
  const int avgEpArraySize = flatOrder+1;
  Double_t avgEpArray[avgEpArraySize] = {0};
  cout << "avgEpArray = {";
  for (int m=0; m<avgEpArraySize; m++) {
    if (m>0) cout << ",";
    avgEpArray[m] = m;
    cout << m;
  }
  cout << "}" << endl;

  TH1D* havgCosEp = new TH1D("havgCosEp","Average of Cos(nPsi)",flatOrder,avgEpArray);
  TH1D* havgSinEp = new TH1D("havgSinEp","Average of Sin(nPsi)",flatOrder,avgEpArray);
  TH1D* havgCosEpHFm2 = new TH1D("havgCosEpHFm2","Average of Cos(nPsi)",flatOrder,avgEpArray);
  TH1D* havgSinEpHFm2 = new TH1D("havgSinEpHFm2","Average of Sin(nPsi)",flatOrder,avgEpArray);
  TH1D* havgCosEpHFp2 = new TH1D("havgCosEpHFp2","Average of Cos(nPsi)",flatOrder,avgEpArray);
  TH1D* havgSinEpHFp2 = new TH1D("havgSinEpHFp2","Average of Sin(nPsi)",flatOrder,avgEpArray);
  TH1D* havgCosEptrackmid2 = new TH1D("havgCosEptrackmid2","Average of Cos(nPsi)",flatOrder,avgEpArray);
  TH1D* havgSinEptrackmid2 = new TH1D("havgSinEptrackmid2","Average of Sin(nPsi)",flatOrder,avgEpArray);

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

  for (int n=1; n<=flatOrder; n++) {
    cout << "sum(Cos(2*" << n << "*Psi)) = " << avgCosEp[n-1] << endl;
    cout << "sum(Sin(2*" << n << "*Psi)) = " << avgSinEp[n-1] << endl;
    avgCosEp[n-1] = avgCosEp[n-1]/ALLPASS;
    avgSinEp[n-1] = avgSinEp[n-1]/ALLPASS;
    avgCosEpHFm2[n-1] = avgCosEpHFm2[n-1]/ALLPASS;
    avgSinEpHFm2[n-1] = avgSinEpHFm2[n-1]/ALLPASS;
    avgCosEpHFp2[n-1] = avgCosEpHFp2[n-1]/ALLPASS;
    avgSinEpHFp2[n-1] = avgSinEpHFp2[n-1]/ALLPASS;
    avgCosEptrackmid2[n-1] = avgCosEptrackmid2[n-1]/ALLPASS;
    avgSinEptrackmid2[n-1] = avgSinEptrackmid2[n-1]/ALLPASS;
    cout << "avg(Cos(2*" << n << "*Psi)) = " << avgCosEp[n-1] << endl;
    cout << "avg(Sin(2*" << n << "*Psi)) = " << avgSinEp[n-1] << endl;
    havgCosEp->SetBinContent(n,avgCosEp[n-1]);
    havgSinEp->SetBinContent(n,avgSinEp[n-1]);
    havgCosEpHFm2->SetBinContent(n,avgCosEpHFm2[n-1]);
    havgSinEpHFm2->SetBinContent(n,avgSinEpHFm2[n-1]);
    havgCosEpHFp2->SetBinContent(n,avgCosEpHFp2[n-1]);
    havgSinEpHFp2->SetBinContent(n,avgSinEpHFp2[n-1]);
    havgCosEptrackmid2->SetBinContent(n,avgCosEptrackmid2[n-1]);
    havgSinEptrackmid2->SetBinContent(n,avgSinEptrackmid2[n-1]);
  }
  havgCosEp->Write();
  havgSinEp->Write();
  havgCosEpHFm2->Write();
  havgSinEpHFm2->Write();
  havgCosEpHFp2->Write();
  havgSinEpHFp2->Write();
  havgCosEptrackmid2->Write();
  havgSinEptrackmid2->Write();

  for (int i=0; i<numptbins; i++) {
    havgCosEppt[i] = new TH1D(Form("havgCosEppt[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEppt[i] = new TH1D(Form("havgSinEppt[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    havgCosEpHFm2pt[i] = new TH1D(Form("havgCosEpHFm2pt[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEpHFm2pt[i] = new TH1D(Form("havgSinEpHFm2pt[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    havgCosEpHFp2pt[i] = new TH1D(Form("havgCosEpHFp2pt[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEpHFp2pt[i] = new TH1D(Form("havgSinEpHFp2pt[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    havgCosEptrackmid2pt[i] = new TH1D(Form("havgCosEptrackmid2pt[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEptrackmid2pt[i] = new TH1D(Form("havgSinEptrackmid2pt[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    for (int n=1; n<=flatOrder; n++) {
      avgCosEppt[i][n-1] = avgCosEppt[i][n-1]/ptPASS[i];
      avgSinEppt[i][n-1] = avgSinEppt[i][n-1]/ptPASS[i];
      avgCosEpHFm2pt[i][n-1] = avgCosEpHFm2pt[i][n-1]/ptPASS[i];
      avgSinEpHFm2pt[i][n-1] = avgSinEpHFm2pt[i][n-1]/ptPASS[i];
      avgCosEpHFp2pt[i][n-1] = avgCosEpHFp2pt[i][n-1]/ptPASS[i];
      avgSinEpHFp2pt[i][n-1] = avgSinEpHFp2pt[i][n-1]/ptPASS[i];
      avgCosEptrackmid2pt[i][n-1] = avgCosEptrackmid2pt[i][n-1]/ptPASS[i];
      avgSinEptrackmid2pt[i][n-1] = avgSinEptrackmid2pt[i][n-1]/ptPASS[i];
      havgCosEppt[i]->SetBinContent(n,avgCosEppt[i][n-1]);
      havgSinEppt[i]->SetBinContent(n,avgSinEppt[i][n-1]);
      havgCosEpHFm2pt[i]->SetBinContent(n,avgCosEpHFm2pt[i][n-1]);
      havgSinEpHFm2pt[i]->SetBinContent(n,avgSinEpHFm2pt[i][n-1]);
      havgCosEpHFp2pt[i]->SetBinContent(n,avgCosEpHFp2pt[i][n-1]);
      havgSinEpHFp2pt[i]->SetBinContent(n,avgSinEpHFp2pt[i][n-1]);
      havgCosEptrackmid2pt[i]->SetBinContent(n,avgCosEptrackmid2pt[i][n-1]);
      havgSinEptrackmid2pt[i]->SetBinContent(n,avgSinEptrackmid2pt[i][n-1]);
    }
    havgCosEppt[i]->Write();
    havgSinEppt[i]->Write();
    havgCosEpHFm2pt[i]->Write();
    havgSinEpHFm2pt[i]->Write();
    havgCosEpHFp2pt[i]->Write();
    havgSinEpHFp2pt[i]->Write();
    havgCosEptrackmid2pt[i]->Write();
    havgSinEptrackmid2pt[i]->Write();
  }

  for (int i=0; i<numcbins; i++) {
    havgCosEpcent[i] = new TH1D(Form("havgCosEpcent[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEpcent[i] = new TH1D(Form("havgSinEpcent[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    havgCosEpHFm2cent[i] = new TH1D(Form("havgCosEpHFm2cent[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEpHFm2cent[i] = new TH1D(Form("havgSinEpHFm2cent[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    havgCosEpHFp2cent[i] = new TH1D(Form("havgCosEpHFp2cent[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEpHFp2cent[i] = new TH1D(Form("havgSinEpHFp2cent[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    havgCosEptrackmid2cent[i] = new TH1D(Form("havgCosEptrackmid2cent[%i]",i),"Average of Cos(nPsi)",flatOrder,avgEpArray);
    havgSinEptrackmid2cent[i] = new TH1D(Form("havgSinEptrackmid2cent[%i]",i),"Average of Sin(nPsi)",flatOrder,avgEpArray);
    for (int n=1; n<=flatOrder; n++) {
      avgCosEpcent[i][n-1] = avgCosEpcent[i][n-1]/centPASS[i];
      avgSinEpcent[i][n-1] = avgSinEpcent[i][n-1]/centPASS[i];
      avgCosEpHFm2cent[i][n-1] = avgCosEpHFm2cent[i][n-1]/centPASS[i];
      avgSinEpHFm2cent[i][n-1] = avgSinEpHFm2cent[i][n-1]/centPASS[i];
      avgCosEpHFp2cent[i][n-1] = avgCosEpHFp2cent[i][n-1]/centPASS[i];
      avgSinEpHFp2cent[i][n-1] = avgSinEpHFp2cent[i][n-1]/centPASS[i];
      avgCosEptrackmid2cent[i][n-1] = avgCosEptrackmid2cent[i][n-1]/centPASS[i];
      avgSinEptrackmid2cent[i][n-1] = avgSinEptrackmid2cent[i][n-1]/centPASS[i];
      havgCosEpcent[i]->SetBinContent(n,avgCosEpcent[i][n-1]);
      havgSinEpcent[i]->SetBinContent(n,avgSinEpcent[i][n-1]);
      havgCosEpHFm2cent[i]->SetBinContent(n,avgCosEpHFm2cent[i][n-1]);
      havgSinEpHFm2cent[i]->SetBinContent(n,avgSinEpHFm2cent[i][n-1]);
      havgCosEpHFp2cent[i]->SetBinContent(n,avgCosEpHFp2cent[i][n-1]);
      havgSinEpHFp2cent[i]->SetBinContent(n,avgSinEpHFp2cent[i][n-1]);
      havgCosEptrackmid2cent[i]->SetBinContent(n,avgCosEptrackmid2cent[i][n-1]);
      havgSinEptrackmid2cent[i]->SetBinContent(n,avgSinEptrackmid2cent[i][n-1]);
    }
    havgCosEpcent[i]->Write();
    havgSinEpcent[i]->Write();
    havgCosEpHFm2cent[i]->Write();
    havgSinEpHFm2cent[i]->Write();
    havgCosEpHFp2cent[i]->Write();
    havgSinEpHFp2cent[i]->Write();
    havgCosEptrackmid2cent[i]->Write();
    havgSinEptrackmid2cent[i]->Write();
  }

  //include the averages used to calculate the event plane resolution correction.
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
    hRcent->SetBinContent(i+1,RAcent);
    hRcent->SetBinError(i+1,RAcenterr);
  }
  hRcent->Write();

  hRpt->Sumw2();
  hRcent->Sumw2();

  TCanvas* c1 = new TCanvas("c1","c1",0,0,600,300);
  c1->Divide(2);
  c1->cd(1);
  hRpt->Draw();
  c1->cd(2);
  hRcent->Draw();

  hqxold->Write();
  hqxnew->Write();
  hqyold->Write();
  hqynew->Write();

  TCanvas* c2 = new TCanvas("c2","c2",0,0,800,400);
  c2->Divide(2);
  c2->cd(1);
  hqxold->SetTitle("qx recentering");
  hqxold->Draw();
  hqxnew->SetLineColor(2);
  hqxnew->Draw("same");

  TLegend* legqx = new TLegend(0.11,0.7,0.3,0.8); legqx->SetTextSize(12);
  legqx->SetTextFont(43);
  legqx->SetBorderSize(0);
  legqx->AddEntry(hqxold,Form("qx raw: avg=%.2f",avgqx),"l");
  legqx->AddEntry(hqxnew,Form("qx rec.: avg=%.2f",avgqxrec),"l");
  legqx->Draw("same");

  float xp = 0.15;
  float yp = 0.65;
  int textColor = kBlack;
  int textSize = 12;
  drawText(Form("RMS = %.2f",hqxold->GetRMS()), xp, yp, textColor, textSize);

  c2->cd(2);
  hqyold->SetTitle("qy recentering");
  hqyold->Draw();
  hqynew->SetLineColor(2);
  hqynew->Draw("same");

  TLegend* legqy = new TLegend(0.11,0.7,0.3,0.8); legqy->SetTextSize(12);
  legqy->SetTextFont(43);
  legqy->SetBorderSize(0);
  legqy->AddEntry(hqyold,Form("qy raw: avg=%.2f",avgqy),"l");
  legqy->AddEntry(hqynew,Form("qy rec.: avg=%.2f",avgqyrec),"l");
  legqy->Draw("same");

  drawText(Form("RMS = %.2f",hqyold->GetRMS()), xp, yp, textColor, textSize);

  c2->SaveAs(Form("plots/recenteringHistos_n%i.pdf",nevt));
  c2->SaveAs(Form("plots/recenteringHistos_n%i.png",nevt));

  hepHF2old->Write();
  hepHF2new->Write();
  hepHFm2old->Write();
  hepHFm2new->Write();
  hepHFp2old->Write();
  hepHFp2new->Write();
  heptrackmid2old->Write();
  heptrackmid2new->Write();

  TCanvas* c3 = new TCanvas("c3","c3",0,0,800,400);
  c3->Divide(2);
  c3->cd(1);
  hepHF2old->SetTitle("Event plane raw");
  hepHF2old->Draw();
  c3->cd(2);
  hepHF2new->SetTitle("Event plane recentered");
  hepHF2new->Draw();

  c3->SaveAs(Form("plots/recenteringEventPlane_n%i.pdf",nevt));
  c3->SaveAs(Form("plots/recenteringEventPlane_n%i.png",nevt));

  TCanvas* c3HFm2 = new TCanvas("c3HFm2","c3HFm2",0,0,800,400);
  c3HFm2->Divide(2);
  c3HFm2->cd(1);
  hepHFm2old->SetTitle("Event plane (HFm2) raw ");
  hepHFm2old->Draw();
  c3HFm2->cd(2);
  hepHFm2new->SetTitle("Event plane (HFm2) recentered");
  hepHFm2new->Draw();

  c3HFm2->SaveAs(Form("plots/recenteringEventPlaneHFm2_n%i.pdf",nevt));
  c3HFm2->SaveAs(Form("plots/recenteringEventPlaneHFm2_n%i.png",nevt));

  TCanvas* c3HFp2 = new TCanvas("c3HFp2","c3HFp2",0,0,800,400);
  c3HFp2->Divide(2);
  c3HFp2->cd(1);
  hepHFp2old->SetTitle("Event plane (HFp2) raw");
  hepHFp2old->Draw();
  c3HFp2->cd(2);
  hepHFp2new->SetTitle("Event plane (HFp2) recentered");
  hepHFp2new->Draw();

  c3HFp2->SaveAs(Form("plots/recenteringEventPlaneHFp2_n%i.pdf",nevt));
  c3HFp2->SaveAs(Form("plots/recenteringEventPlaneHFp2_n%i.png",nevt));

  TCanvas* c3trackmid2 = new TCanvas("c3trackmid2","c3trackmid2",0,0,800,400);
  c3trackmid2->Divide(2);
  c3trackmid2->cd(1);
  heptrackmid2old->SetTitle("Event plane (trackmid2) raw");
  heptrackmid2old->Draw();
  c3trackmid2->cd(2);
  heptrackmid2new->SetTitle("Event plane (trackmid2) recentered");
  heptrackmid2new->Draw();

  c3trackmid2->SaveAs(Form("plots/recenteringEventPlanetrackmid2_n%i.pdf",nevt));
  c3trackmid2->SaveAs(Form("plots/recenteringEventPlanetrackmid2_n%i.png",nevt));

  avgFile->Close();

  cout << "ptPASS = {";
  for (int i=0; i<numptbins; i++) {
    cout << ptPASS[i] << ",";
  }
  cout << "}" << endl;
  cout << "centPASS = {";
  for (int i=0; i<numcbins; i++) {
    cout << centPASS[i] << ",";
  }
  cout << "}" << endl;
  cout << "ALLPASS = " << ALLPASS << endl;

  cout << endl;
  cout << "avgqx = " << avgqx << endl;
  cout << "avgqy = " << avgqy << endl;
  cout << "avgqxHFm2 = " << avgqxHFm2 << endl;
  cout << "avgqyHFm2 = " << avgqyHFm2 << endl;
  cout << "avgqxHFp2 = " << avgqxHFp2 << endl;
  cout << "avgqyHFp2 = " << avgqyHFp2 << endl;
  cout << "avgqxtrackmid2 = " << avgqxtrackmid2 << endl;
  cout << "avgqytrackmid2 = " << avgqytrackmid2 << endl;

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

