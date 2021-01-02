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

void getAcceptanceFromSkim(int nevt=10000, int dateStr=20210102) 
{

  using namespace std;
  using namespace hi;

  Double_t ptbins[4] = {0,3,6,30};
  const int numptbins = sizeof(ptbins)/sizeof(Double_t)-1;
  Double_t ptbinsInt[2] = {0,30};
  const int numptbinsInt = 1;

  // Integrated bin:
  TH1D* hAccIntNoW = new TH1D("hAccIntNoW",";p_{T} (GeV/c);Acceptance of #Upsilon()",numptbinsInt,ptbinsInt);
  TH1D* hAccInt = new TH1D("hAccInt",";p_{T} (GeV/c);Acceptance of #Upsilon()",numptbinsInt,ptbinsInt);
  hAccIntNoW->Sumw2();
  hAccInt->Sumw2();
  hAccIntNoW->SetMarkerStyle(20);hAccIntNoW->SetMarkerColor(kBlue+2);
  hAccInt->SetMarkerStyle(20);hAccInt->SetMarkerColor(kBlue+2);
  
  // pT dependence:
  TH1D* hAccPtNoW = new TH1D("hAccPtNoW",";p_{T} (GeV/c);Acceptance of #Upsilon",numptbins,ptbins);
  TH1D* hAccPt = new TH1D("hAccPt",";p_{T} (GeV/c);Acceptance of #Upsilon",numptbins,ptbins);
  hAccPtNoW->Sumw2();
  hAccPt->Sumw2();
  hAccPtNoW->SetMarkerStyle(20);hAccPtNoW->SetMarkerColor(kBlue+2);
  hAccPt->SetMarkerStyle(20);hAccPt->SetMarkerColor(kBlue+2);

  TString inFileName = "/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_PbPb_502TeV_ClosureTest/skims/newOniaTree_Skim_UpsTrig_MC_RAW_20200615.root";
  TFile* inFile = TFile::Open(inFileName,"READ");

  TTree* mytree = (TTree*)inFile->Get("mmep");
  TBranch* mm = (TBranch*)mytree->GetBranch("mm");

  int DIMUIDPASS = 0;
  int ACCPASS = 0;

  int ptPASS[3] = {0};

  if(nevt == -1) nevt = 4023601;

  int nevtReal = mytree->GetEntries();

  cout << "Total events = " << nevtReal << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevtReal ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/nevtReal) << "%)" << endl;

    mytree->GetEntry(iev);
    TLeaf *ptLeaf = mm->GetLeaf("pt");
    TLeaf *yLeaf = mm->GetLeaf("y");
    TLeaf *pt1Leaf = mm->GetLeaf("pt1");
    TLeaf *eta1Leaf = mm->GetLeaf("eta1");
    TLeaf *pt2Leaf = mm->GetLeaf("pt2");
    TLeaf *eta2Leaf = mm->GetLeaf("eta2");

    Double_t pt = ptLeaf->GetValue();
    Double_t y = yLeaf->GetValue();
    Double_t pt1 = pt1Leaf->GetValue();
    Double_t eta1 = eta1Leaf->GetValue();
    Double_t pt2 = pt2Leaf->GetValue();
    Double_t eta2 = eta2Leaf->GetValue();

    //Count how many dimuons there are.
    if (pt>50 || abs(y)>2.4) continue;
    DIMUIDPASS++;

    //Count how many of them have both muons in the single-muon acceptance.

    if ( pt1>3.5 && pt2>3.5 && abs(eta1)<2.4 && abs(eta2)<2.4 ){
      ACCPASS++;
    }

  } //end of loop

  Double_t acceptance = (double)ACCPASS/(double)DIMUIDPASS;
  cout << "acceptance = " << acceptance << endl;
  hAccInt->Fill( acceptance );

}

