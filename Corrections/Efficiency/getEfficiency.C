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

void getEfficiency(int nevt=-1, int dateStr=20210112, int weighted=0){

  gStyle->SetOptStat(0);

  using namespace std;
  using namespace hi;

  const int numptbins = 100;
  const int arraysize = numptbins+1;
  Double_t ptbins[arraysize];
  Double_t ptmax = 50.0;
  Double_t ptmin = 0.0;
  cout << "ptbins = {";
  for (int i=0; i<arraysize; i++){
    ptbins[i] = ptmax*i/((Double_t)numptbins);
    cout << ptbins[i] << ", ";
  }
  cout << "}" << endl;
  Double_t ptbinsInt[2] = {ptmin,ptmax};
  const int numptbinsInt = 1;

  // Integrated bin:
  TH1D* hEffInt = new TH1D("hEffInt",";p_{T} (GeV/c);Efficiency of Dimuons",numptbinsInt,ptbinsInt);
  hEffInt->Sumw2(); hEffInt->SetMinimum(0); hEffInt->SetMaximum(1.0);
  hEffInt->SetMarkerStyle(20);hEffInt->SetMarkerColor(kBlue+2);

  // pT dependence:
  TH1D* hEffPt = new TH1D("hEffPt",";p_{T} (GeV/c);Efficiency of Dimuons",numptbins,ptbins);
  hEffPt->Sumw2(); hEffPt->SetMinimum(0); hEffPt->SetMaximum(1.0);
  hEffPt->SetMarkerStyle(20);hEffPt->SetMarkerColor(kBlue+2);

  TString inputMC1 = "/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_v1.root";
  TString inputMC2 = "/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_ext-v1.root";
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC1.Data());
  mytree->Add(inputMC2.Data());

  TChain* tree = new TChain("tree"); 
  tree->Add(inputMC1.Data());
  tree->Add(inputMC2.Data());

  mytree->AddFriend(tree);

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
  Double_t NUMERATOR = 0;//# of RECO dimuons in acceptance with muId+Trig
  Double_t DENOMINATOR = 0;//# of GEN dimuons in acceptance

  Double_t NUMERATORPT[numptbins] = {0};//# of RECO dimuons in bin with muId+Trig
  Double_t DENOMINATORPT[numptbins] = {0};//# of GEN dimuons in bin

  Double_t pt, eta, pt1, pt2, eta1, eta2;

  int nevtReal = mytree->GetEntries();

  if(nevt == -1) nevt = nevtReal;

  cout << "Total events = " << nevtReal << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/nevt) << "%)" << endl;

    mytree->GetEntry(iev);

    //Apply event trigger
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    int Gen_QQ_pass = 0;
    int Reco_QQ_pass = 0;

      //cout << "Gen_QQ_size=" << Gen_QQ_size << ", Reco_QQ_size=" << Reco_QQ_size << ", Gen_mu_size=" << Gen_mu_size << ", Reco_mu_size=" << Reco_mu_size << endl;

    //Gen dimuon loop start
    for (Int_t irqq=0; irqq<Gen_QQ_size; ++irqq) 
    {

      //Get 4mom
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(irqq);
      pt = JP_Gen->Pt();
      eta = JP_Gen->Eta();

      if (pt>50 || abs(eta)>2.4) continue;

      mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[irqq]);
      mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[irqq]);

      pt1 = mupl_Gen->Pt();
      eta1 = mupl_Gen->Eta();
      pt2 = mumi_Gen->Pt();
      eta2 = mumi_Gen->Eta();

      //cout << " Gen pt1=" << pt1 << ", eta1=" << eta1 << ", pt2=" << pt2 << ", eta2=" << eta2;

      //Apply acceptance cuts
      if ( pt1<3.5 || pt2<3.5 || abs(eta1)>2.4 || abs(eta2)>2.4 ) continue;

      //Get scale factors

      //Increment DENOMINATOR with scale factors
      DENOMINATOR = DENOMINATOR + getSFs(pt1,eta1,pt2,eta2);

      //versus pt
      //cout << "Gen pt=" << pt << ", eta=" << eta << endl;
      int whichptbin = hEffPt->FindBin(pt)-1;
      //cout << "whichptbin=" << whichptbin << endl;
      DENOMINATORPT[whichptbin] = DENOMINATORPT[whichptbin] + getSFs(pt1,eta1,pt2,eta2);

      //cout << " passed!" << endl;
      Gen_QQ_pass++;
    }//end of gen dimuon loop

    //Reco dimuon loop start
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {

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
      eta = JP_Reco->Eta();

      //cout << "Reco pt=" << pt << ", eta=" << eta;

      if (pt>50 || abs(eta)>2.4) continue;

      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);

      pt1 = mupl_Reco->Pt();
      eta1 = mupl_Reco->Eta();
      pt2 = mumi_Reco->Pt();
      eta2 = mumi_Reco->Eta();

      //Apply acceptance cuts
      if ( pt1<3.5 || pt2<3.5 || abs(eta1)>2.4 || abs(eta2)>2.4 ) continue;

      //Get scale factors

      //Increment NUMERATOR with scale factors
      NUMERATOR = NUMERATOR + getSFs(pt1,eta1,pt2,eta2);

      //versus pt
      //cout << "Reco pt=" << pt << ", eta=" << eta;
      int whichptbin = hEffPt->FindBin(pt)-1;
      //cout << "whichptbin=" << whichptbin << endl;
      NUMERATORPT[whichptbin] = NUMERATORPT[whichptbin] + getSFs(pt1,eta1,pt2,eta2);

      //cout << " passed!" << endl;
      Reco_QQ_pass++;

    }//end of reco dimuon loop

    //cout << "Gen_QQ_pass=" << Gen_QQ_pass << ", Reco_QQ_pass=" << Reco_QQ_pass << endl;
    //if (Reco_QQ_pass>Gen_QQ_pass) {
    //  cout << "^^^^^^^ BAD ONE!!! ^^^^^^^" << endl << endl;
    //}
  } //end of event loop

  Double_t efficiency = (double)NUMERATOR/(double)DENOMINATOR;
  cout << "efficiency = " << NUMERATOR << "/" << DENOMINATOR << " = " << efficiency << endl;
  hEffInt->Fill( efficiency );

  for (int i=0; i<numptbins; i++) {
    efficiency = 0.0;
    if (DENOMINATORPT[i]>0) efficiency = (double)NUMERATORPT[i]/(double)DENOMINATORPT[i];
    cout << "efficiency = " << NUMERATORPT[i] << "/" << DENOMINATORPT[i] << " = " << efficiency << endl;
    hEffPt->SetBinContent(i+1,efficiency);
  }

  TCanvas* c1 = new TCanvas("c1","c1",0,0,500,500);
  hEffPt->Draw();

  TString weightstr = "NoW";

  c1->SaveAs(Form("EfficiencyPlot%s.png",weightstr.Data()));
  c1->SaveAs(Form("EfficiencyPlot%s.pdf",weightstr.Data()));

  TFile* effFile = new TFile(Form("efficiency%s.root",weightstr.Data()),"RECREATE");
  hEffInt->Write();
  hEffPt->Write();
  //feff->Write();
  effFile->Close();
}

