#include <ctime>

#include <TLorentzVector.h>
#include "../../HeaderFiles/commonUtility.h"
#include "../../HeaderFiles/HiEvtPlaneList.h"
#include "../../HeaderFiles/cutsAndBinUpsilonV2.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"

#include "tnp_weight_lowptPbPb.h"

const Double_t pi = 3.141592653589;

//TF1* fdNdpTWgt = new TF1("fdNdpTWgt","([0]+[1]*x+[2]*x*x)/((x-[3])*(x-[3])*(x-[3]))");

double getSFs(double pt1, double eta1, double pt2, double eta2){
  double weight = 1.0;     
  weight = weight*tnp_weight_trg_pbpb(pt1,eta1);
  weight = weight*tnp_weight_muid_pbpb(pt1,eta1);
  weight = weight*tnp_weight_trk_pbpb(eta1);
  weight = weight*tnp_weight_trg_pbpb(pt2,eta2);
  weight = weight*tnp_weight_muid_pbpb(pt2,eta2);
  weight = weight*tnp_weight_trk_pbpb(eta2);
  return weight;
}

void getEfficiency(int nevt=-1, int dateStr=20210123, int weighted=0){

  gStyle->SetOptStat(0);

  using namespace std;
  using namespace hi;

  Double_t ptbins[33] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23, 26,29,32,35,38,41,44,47,50};
  int numptbins = sizeof(ptbins)/sizeof(Double_t)-1;

  // Integrated bin:
  TH1D* hEffInt = new TH1D("hEffInt",";p_{T} (GeV/c);Efficiency of Dimuons",1,0,50);
  hEffInt->Sumw2(); hEffInt->SetMinimum(0); hEffInt->SetMaximum(1.0);
  hEffInt->SetMarkerStyle(20);hEffInt->SetMarkerColor(kBlue+2);
  TH1D* hDenInt = new TH1D("hDenInt",";p_{T} (GeV/c);Efficiency of Dimuons",1,0,50);
  hDenInt->Sumw2(); hDenInt->SetMinimum(0); //hDenInt->SetMaximum(1.0);
  hDenInt->SetMarkerStyle(20);hDenInt->SetMarkerColor(kBlue+2);

  TH1D* hEffIntNoW = (TH1D*)hEffInt->Clone("hEffIntNoW");
  TH1D* hDenIntNoW = (TH1D*)hDenInt->Clone("hDenIntNoW");

  // pT dependence:
  TH1D* hEffPt = new TH1D("hEffPt",";p_{T} (GeV/c);Efficiency of Dimuons",numptbins,ptbins);
  hEffPt->Sumw2(); hEffPt->SetMinimum(0); hEffPt->SetMaximum(1.0);
  hEffPt->SetMarkerStyle(20);hEffPt->SetMarkerColor(kBlue+2);
  TH1D* hDenPt = new TH1D("hDenPt",";p_{T} (GeV/c);Efficiency of Dimuons",numptbins,ptbins);
  hDenPt->Sumw2(); hDenPt->SetMinimum(0); //hDenPt->SetMaximum(1.0);
  hDenPt->SetMarkerStyle(20);hDenPt->SetMarkerColor(kBlue+2);

  TH1D* hEffPtLowC = (TH1D*)hEffPt->Clone("hEffPtLowC");
  TH1D* hEffPtMidC = (TH1D*)hEffPt->Clone("hEffPtMidC");
  TH1D* hEffPtHighC = (TH1D*)hEffPt->Clone("hEffPtHighC");

  TH1D* hDenPtLowC = (TH1D*)hDenPt->Clone("hDenPtLowC");
  TH1D* hDenPtMidC = (TH1D*)hDenPt->Clone("hDenPtMidC");
  TH1D* hDenPtHighC = (TH1D*)hDenPt->Clone("hDenPtHighC");

  TH1D* hEffPtNoW = (TH1D*)hEffPt->Clone("hEffPtNoW");
  TH1D* hEffPtLowCNoW = (TH1D*)hEffPt->Clone("hEffPtLowCNoW");
  TH1D* hEffPtMidCNoW = (TH1D*)hEffPt->Clone("hEffPtMidCNoW");
  TH1D* hEffPtHighCNoW = (TH1D*)hEffPt->Clone("hEffPtHighCNoW");

  TH1D* hDenPtNoW = (TH1D*)hDenPt->Clone("hDenPtNoW");
  TH1D* hDenPtLowCNoW = (TH1D*)hDenPt->Clone("hDenPtLowCNoW");
  TH1D* hDenPtMidCNoW = (TH1D*)hDenPt->Clone("hDenPtMidCNoW");
  TH1D* hDenPtHighCNoW = (TH1D*)hDenPt->Clone("hDenPtHighCNoW");

  //pt weighting function
  TFile *fweight = TFile::Open("../Acceptance/Func_dNdpT_1S.root","READ");
  TF1* fitRatio = (TF1*)fweight->Get("fitRatio");

  TString inputMC1 = "/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_v1.root";
  TString inputMC2 = "/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_ext-v1.root";
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC1.Data());
  mytree->Add(inputMC2.Data());

  TChain* tree = new TChain("tree"); 
  tree->Add(inputMC1.Data());
  tree->Add(inputMC2.Data());

  //mytree->AddFriend(tree);

  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  Int_t           Reco_mu_size;
  Int_t           Reco_mu_whichGen[maxBranchSize];
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  Int_t           Gen_QQ_size;
  Int_t           Gen_mu_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_mu_whichGen;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_mu_4mom;   //!

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;

  Gen_QQ_4mom = 0;
  Gen_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);

  Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);


  Int_t           Gen_QQ_mupl_idx[maxBranchSize];
  Int_t           Gen_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx,&b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx,&b_Gen_QQ_mumi_idx);


  const int nEP = 29;  // number of event planes in the tree
  Double_t qx[nEP]; 
  Double_t qy[nEP]; 
  TBranch *b_qx;
  TBranch *b_qy;
  tree->SetBranchAddress("qx", qx, &b_qx);
  tree->SetBranchAddress("qy", qy, &b_qy);

  Double_t    	  epang[29];
  TBranch         *b_epang;
  tree->SetBranchAddress("epang", epang, &b_epang);
  
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  TLorentzVector* JP_Gen = new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;

  Double_t kTrigSel = 13;

  Double_t pt, rap, pt1, pt2, eta1, eta2, weight, ptweight;

  int nevtReal = mytree->GetEntries();

  if(nevt == -1) nevt = nevtReal;

  cout << "Total events = " << nevtReal << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/nevt) << "%)" << endl;

    mytree->GetEntry(iev);

    //Centrality cut 10-90%
    if (Centrality<20 || Centrality>180) continue;

    //Gen dimuon loop start
    for (Int_t irqq=0; irqq<Gen_QQ_size; ++irqq) 
    {

      //if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]] * Gen_mu_charge[Gen_QQ_mumi_idx[igen]] > 0) continue;

      //Get 4mom
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(irqq);
      pt = JP_Gen->Pt();
      rap = JP_Gen->Rapidity();

      if (pt>50 || fabs(rap)>2.4) continue;

      mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[irqq]);
      mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[irqq]);

      pt1 = mupl_Gen->Pt();
      eta1 = mupl_Gen->Eta();
      pt2 = mumi_Gen->Pt();
      eta2 = mumi_Gen->Eta();

      //cout << " Gen pt1=" << pt1 << ", eta1=" << eta1 << ", pt2=" << pt2 << ", eta2=" << eta2;

      //Apply acceptance cuts
      if ( pt1<3.5 || pt2<3.5 || fabs(eta1)>2.4 || fabs(eta2)>2.4 ) continue;

      //Get scale factors and weights
      weight = findNcoll(Centrality)*Gen_weight;

      //Apply pt reweighting function;
      ptweight = fitRatio->Eval(pt);

      //Increment DENOMINATOR
      //weight = weight*getSFs(pt1,eta1,pt2,eta2);
      hDenInt->Fill(pt,weight*ptweight);
      hDenPt->Fill(pt,weight*ptweight);
      hDenIntNoW->Fill(pt);
      hDenPtNoW->Fill(pt);

      if (Centrality<60) {
        hDenPtLowC->Fill(pt,weight*ptweight);//10-30%
        hDenPtLowCNoW->Fill(pt);
      }
      else if (Centrality<100) {
        hDenPtMidC->Fill(pt,weight*ptweight);//30-50%
        hDenPtMidCNoW->Fill(pt);
      }
      else {
        hDenPtHighC->Fill(pt,weight*ptweight);//50-90%
        hDenPtHighCNoW->Fill(pt);
      }
    }//end of gen dimuon loop

    //Apply event trigger
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    //Reco dimuon loop start
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {

      //Require opposite sign muons
      if (!(Reco_QQ_sign[irqq]==0)) continue;

      //Apply vertex cut
      if ( Reco_QQ_VtxProb[irqq] < 0.01 ) continue;

      //Apply dimuon trigger
      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
      
      //Apply muid
      bool passMuonTypePl = true;
      //passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      //passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      //passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      //passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

      if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
      if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      //cout << "Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]]=" << Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] << endl;
      //cout << "Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]]=" << Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] << endl;

      bool muplSoft = ( passMuonTypePl && //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) 
          ) ; 

      bool mumiSoft = ( passMuonTypeMi && //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  
          ) ; 

      if ( !(muplSoft && mumiSoft) ) 
        continue;   

      //Get 4mom
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      pt = JP_Reco->Pt();
      rap = JP_Reco->Rapidity();

      //cout << "Reco pt=" << pt << ", eta=" << eta;

      if (pt>50 || fabs(rap)>2.4) continue;

      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      pt1 = mupl_Reco->Pt();
      eta1 = mupl_Reco->Eta();
      pt2 = mumi_Reco->Pt();
      eta2 = mumi_Reco->Eta();

      //Apply acceptance cuts
      if ( pt1<3.5 || pt2<3.5 || fabs(eta1)>2.4 || fabs(eta2)>2.4 ) continue;

      //Apply pt reweighting function;
      weight = findNcoll(Centrality)*Gen_weight;
      ptweight = fitRatio->Eval(pt);

      //Get scale factors
      weight = weight*getSFs(pt1,eta1,pt2,eta2);

      //Increment NUMERATOR with scale factors
      hEffInt->Fill(pt,weight*ptweight);
      hEffPt->Fill(pt,weight*ptweight);
      hEffIntNoW->Fill(pt);
      hEffPtNoW->Fill(pt);

      if (Centrality<60) {
        hEffPtLowC->Fill(pt,weight*ptweight);//10-30%
        hEffPtLowCNoW->Fill(pt);
      }
      else if (Centrality<100) {
        hEffPtMidC->Fill(pt,weight*ptweight);//30-50%
        hEffPtMidCNoW->Fill(pt);
      }
      else {
        hEffPtHighC->Fill(pt,weight*ptweight);//50-90%
        hEffPtHighCNoW->Fill(pt);
      }

    }//end of reco dimuon loop

  } //end of event loop

  hEffInt->Divide(hDenInt);
  hEffPt->Divide(hDenPt);
  hEffPtLowC->Divide(hDenPtLowC);
  hEffPtMidC->Divide(hDenPtMidC);
  hEffPtHighC->Divide(hDenPtHighC);

  hEffIntNoW->Divide(hDenIntNoW);
  hEffPtNoW->Divide(hDenPtNoW);
  hEffPtLowCNoW->Divide(hDenPtLowCNoW);
  hEffPtMidCNoW->Divide(hDenPtMidCNoW);
  hEffPtHighCNoW->Divide(hDenPtHighCNoW);

  TCanvas* c1 = new TCanvas("c1","c1",0,0,1200,400);
  c1->Divide(3);
  c1->cd(1);
  hEffPtLowC->Draw("PE");
  hEffPtLowCNoW->SetMarkerStyle(24);
  hEffPtLowCNoW->Draw("same PE");
  c1->cd(2);
  hEffPtMidC->Draw("PE");
  hEffPtMidCNoW->SetMarkerStyle(24);
  hEffPtMidCNoW->Draw("same PE");
  c1->cd(3);
  hEffPtHighC->Draw("PE");
  hEffPtHighCNoW->SetMarkerStyle(24);
  hEffPtHighCNoW->Draw("same PE");

  TCanvas* c2 = new TCanvas("c2","c2",400,0,500,500);
  hDenPt->SetTitle("Gen: Centrality 50-90%");
  hDenPt->Draw("PE");

  c1->SaveAs(Form("EfficiencyPlot%s_%i.png",dateStr));
  c1->SaveAs(Form("EfficiencyPlot%s_%i.pdf",dateStr));

  TFile* effFile = new TFile(Form("efficiency%s_%i.root",dateStr),"RECREATE");
  hEffInt->Write();
  hEffPt->Write();
  hEffPtLowC->Write();
  hEffPtMidC->Write();
  hEffPtHighC->Write();
  hEffIntNoW->Write();
  hEffPtNoW->Write();
  hEffPtLowCNoW->Write();
  hEffPtMidCNoW->Write();
  hEffPtHighCNoW->Write();
  effFile->Close();
}

