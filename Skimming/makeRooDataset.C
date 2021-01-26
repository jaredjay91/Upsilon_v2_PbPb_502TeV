#include <ctime>

#include <TLorentzVector.h>
#include "../HeaderFiles/commonUtility.h"
#include "../HeaderFiles/HiEvtPlaneList.h"
#include "../HeaderFiles/cutsAndBinUpsilonV2.h"
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


void makeRooDataset(int dateStr=20210125, bool NomAccTrue=kTRUE, bool NomEffTrue=kFALSE) 
{

  using namespace std;
  using namespace hi;
  using namespace RooFit;

  TFile* accFile = TFile::Open("../Corrections/Acceptance/acceptance_20210122.root");
  TH1D* hAccPt = (TH1D*)accFile->Get("hAccPt");
  TH1D* hAccPtNoW = (TH1D*)accFile->Get("hAccPtNoW");

  TFile* effFile = TFile::Open("../Corrections/Efficiency/efficiency_20210123.root");
  TH1D* hEffPtLowC = (TH1D*)effFile->Get("hEffPtLowC");
  TH1D* hEffPtMidC = (TH1D*)effFile->Get("hEffPtMidC");
  TH1D* hEffPtHighC = (TH1D*)effFile->Get("hEffPtHighC");
  TH1D* hEffPtLowCNoW = (TH1D*)effFile->Get("hEffPtLowCNoW");
  TH1D* hEffPtMidCNoW = (TH1D*)effFile->Get("hEffPtMidCNoW");
  TH1D* hEffPtHighCNoW = (TH1D*)effFile->Get("hEffPtHighCNoW");

  //TFile* inFile = TFile::Open("skims/newOniaTree_Skim_UpsTrig_RD_RAW_n-1_20201107.root","READ");
  //TFile* inFile = TFile::Open("condor/skims/newOniaTree_Skim_UpsTrig_RD_flattenedBinByBin_order21_n-1_20201114.root","READ");
  TFile* inFile = TFile::Open("condor/skims/newOniaTree_Skim_UpsTrig_MM_flattenedBinByBin_order21_n-1_20201119.root","READ");
  TTree* mytree = (TTree*)inFile->Get("mmep");
  TBranch* mm = (TBranch*)mytree->GetBranch("mm");

  int nevt = mytree->GetEntries();

  TString accstr = "";
  TString effstr = "";
  if (!NomAccTrue) accstr = "altAcc_";
  if (!NomEffTrue) effstr = "altEff_";
  TFile* newfile;
  //newfile = new TFile(Form("skims/newOniaTree_Skim_UpsTrig_RD_withDataset_%i.root",dateStr),"recreate");
  newfile = new TFile(Form("skims/newOniaTree_Skim_UpsTrig_MM_flattenedBinByBin_order21_n-1_withDataset_%s%s%i.root", accstr.Data(), effstr.Data(), dateStr),"recreate");
  
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","pt weight", 0, 10000,"");
  RooRealVar* recoQQsign = new RooRealVar("recoQQsign","qq sign",-1,3,"");
  RooRealVar* dphiEp2Var   = new RooRealVar("dphiEp2","Delta Phi from 2nd order event plane", -100,100,"");
  RooRealVar* qxVar   = new RooRealVar("qx","x-comp of 2nd order q vector", -100,100,"");
  RooRealVar* qyVar   = new RooRealVar("qy","y-comp of 2nd order q vector", -100,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var);
  argSet->add(*cBinVar); argSet->add(*ep2Var); argSet->add(*recoQQsign); argSet->add(*dphiEp2Var); argSet->add(*qxVar); argSet->add(*qyVar); argSet->add(*evtWeight);

  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);
  //RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet, WeightVar(*evtWeight));

  newfile->cd();

  // event loop start

  cout << "Total events = " << nevt << ", : " << mytree->GetEntries() << endl;

  double weight, pt, cBin;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    //"run/I:lumi:event:cBin:ep2/F:dphiEp2:vz:mass:pt:y:phi:eta:pt1:eta1:phi1:pt2:eta2:phi2:weight0:weight:oniaIndex/I:softFlag:highPtFlag:qxa/F:qya:qxb:qyb:qxc:qyc:qxdimu:qydimu"
    mytree->GetEntry(iev);
    TLeaf *massLeaf = mm->GetLeaf("mass");
    TLeaf *ptLeaf = mm->GetLeaf("pt");
    TLeaf *yLeaf = mm->GetLeaf("y");
    TLeaf *pt1Leaf = mm->GetLeaf("pt1");
    TLeaf *eta1Leaf = mm->GetLeaf("eta1");
    TLeaf *pt2Leaf = mm->GetLeaf("pt2");
    TLeaf *eta2Leaf = mm->GetLeaf("eta2");
    TLeaf *ep2Leaf = mm->GetLeaf("ep2");
    TLeaf *dphiEp2Leaf = mm->GetLeaf("dphiEp2");
    TLeaf *cBinLeaf = mm->GetLeaf("cBin");

    pt = (double)ptLeaf->GetValue();
    cBin = (double)cBinLeaf->GetValue();

    if (pt>50) continue;

    recoQQsign->setVal(0);
    ptVar->setVal(pt); 
    cBinVar->setVal(cBin);  
    massVar->setVal((double)massLeaf->GetValue());
    yVar->setVal((double)yLeaf->GetValue());
    pt1Var->setVal((double)pt1Leaf->GetValue());
    eta1Var->setVal((double)eta1Leaf->GetValue());
    pt2Var->setVal((double)pt2Leaf->GetValue());
    eta2Var->setVal((double)eta2Leaf->GetValue());
    ep2Var->setVal((double)ep2Leaf->GetValue());
    dphiEp2Var->setVal((double)dphiEp2Leaf->GetValue());

    //Apply acceptance and efficiency weights.
    weight = 1.0;
    int whichptbin = hAccPt->FindBin((double)ptLeaf->GetValue());
    if (NomAccTrue) weight = weight/(hAccPt->GetBinContent(whichptbin));
    else weight = weight/(hAccPtNoW->GetBinContent(whichptbin));

    whichptbin = hEffPtLowC->FindBin((double)ptLeaf->GetValue());
    if (cBin>20 && cBin<60) {//Centrality 10-30%
      if (NomEffTrue) weight = weight/(hEffPtLowC->GetBinContent(whichptbin));
      else weight = weight/(hEffPtLowCNoW->GetBinContent(whichptbin));
    }
    else if (cBin>=60 && cBin<100) {//Centrality 30-50%
      if (NomEffTrue) weight = weight/(hEffPtMidC->GetBinContent(whichptbin));
      else weight = weight/(hEffPtMidCNoW->GetBinContent(whichptbin));
    }
    else if (cBin>=100 && cBin<180) {//Centrality 50-90%
      if (NomEffTrue) weight = weight/(hEffPtHighC->GetBinContent(whichptbin));
      else weight = weight/(hEffPtHighCNoW->GetBinContent(whichptbin));
    }
    evtWeight->setVal( (double)weight ) ;

    if (weight>100) cout << "weight = " << weight << ", pt = " << pt << ", ptbin = " << whichptbin << endl;

    dataSet->add( *argSet);

  }

  cout << "Adding weights..." << endl;

  //Make sure the dataset is weighted. See https://root.cern.ch/doc/v610/rf403__weightedevts_8C.html
  RooDataSet* dataSetWeighted = new RooDataSet(dataSet->GetName(),dataSet->GetTitle(),dataSet,*dataSet->get(),0,evtWeight->GetName());

  //dataSet->Write();
  dataSetWeighted->Write();
  newfile->Close();

  cout << "Done." << endl;

} 

