#include <ctime>
#include <TLorentzVector.h>
#include "../../HeaderFiles/commonUtility.h"
#include "../../HeaderFiles/HiEvtPlaneList.h"
#include "../../HeaderFiles/cutsAndBinUpsilonV2.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

const Double_t pi = 3.141592653589;

void plotAcceptance(int nevt=-1, int dateStr=20210702, int weighted=0){

  gStyle->SetOptStat(0);

  using namespace std;
  using namespace hi;

  Int_t nbinsy = 48;
  Double_t ylow = 0.0;
  Double_t yup = 12.0;

  Int_t nbinsx = 48;
  Double_t xlow = 0.0;
  Double_t xup = 2.4;

  TH2F* hReco = new TH2F("hReco","hReco",nbinsx, xlow, xup, nbinsy, ylow, yup);
  hReco->Sumw2();
  hReco->GetXaxis()->SetTitle("|#eta|");
  hReco->GetXaxis()->SetLabelSize(0.04);
  hReco->GetXaxis()->SetTitleSize(0.05);
  hReco->GetYaxis()->SetTitle("p_{T}");
  hReco->GetYaxis()->SetLabelSize(0.04);
  hReco->GetYaxis()->SetTitleSize(0.04);
  hReco->GetZaxis()->SetLabelSize(0.04);
  hReco->SetTitle("(Reco+Id)/Generated");

  TH2F* hGen = (TH2F*)hReco->Clone();

  TString inputMC1 = "/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_v1.root";
  TString inputMC2 = "/home/jared/Documents/Ubuntu_Overflow/DataTrees/2018PbPbMCOfficial/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_ext-v1.root";
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC1.Data());
  mytree->Add(inputMC2.Data());

  //TChain* tree = new TChain("tree"); 
  //tree->Add(inputMC1.Data());
  //tree->Add(inputMC2.Data());
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
  //Int_t           Reco_mu_idx[maxBranchSize];
  //TBranch        *b_Reco_mu_idx;
  //mytree->SetBranchAddress("Reco_mu_idx",Reco_mu_idx,&b_Reco_mu_idx);

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

  
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* mu_Reco = new TLorentzVector;
  //TLorentzVector* mu_Gen = new TLorentzVector;

  Double_t kTrigSel = 13;

  Double_t pt, eta, weight;

  int nevtReal = mytree->GetEntries();

  if(nevt == -1) nevt = nevtReal;

  cout << "Total events = " << nevtReal << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/nevt) << "%)" << endl;

    mytree->GetEntry(iev);

    //Apply event trigger
    //if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    //Gen muon loop start
    for (Int_t irmu=0; irmu<Reco_mu_size; ++irmu) 
    {

      //Get 4mom
      mu_Reco = (TLorentzVector*) Reco_mu_4mom->At(irmu);

      pt = mu_Reco->Pt();
      eta = mu_Reco->Eta();

      //Apply pt reweighting function;
      weight = 1.0;//findNcoll(Centrality)*Gen_weight;
      //ptweight = fitRatio->Eval(pt);

      //Increment DENOMINATOR
      hGen->Fill(eta,pt,weight);

    }//end of gen muon loop

    //Reco muon loop start
    for (Int_t irmu=0; irmu<Reco_mu_size; ++irmu) 
    {
      
      //Apply muid
      bool passMuonTypePl = true;
      bool passMuonTypeMi = true;

      if(Reco_mu_whichGen[irmu] == -1) continue;

      bool muSoft = ( passMuonTypePl && //(Reco_mu_TMOneStaTight[irmu]==true) &&
          (Reco_mu_nTrkWMea[irmu] > 5) &&
          (Reco_mu_nPixWMea[irmu] > 0) &&
          (fabs(Reco_mu_dxy[irmu])<0.3) &&
          (fabs(Reco_mu_dz[irmu])<20.) 
          ) ; 

      if ( ! muSoft ) continue;   

      //Get 4mom
      mu_Reco = (TLorentzVector*) Reco_mu_4mom->At(irmu);

      pt = mu_Reco->Pt();
      eta = mu_Reco->Eta();

      //Apply acceptance cuts
      //if ( pt1<3.5 || pt2<3.5 || fabs(eta1)>2.4 || fabs(eta2)>2.4 ) continue;

      //Apply pt reweighting function;
      weight = 1.0;//findNcoll(Centrality);//*Gen_weight;
      //ptweight = fitRatio->Eval(pt);

      //Increment NUMERATOR with scale factors
      hReco->Fill(eta,pt,weight);

    }//end of reco muon loop

  } //end of event loop

  TCanvas* cData = new TCanvas("cData","cData",50,50,400,400);
  cData->cd();
  cout << "Plotting data distribution" << endl;
  hReco->Divide(hGen);
  hReco->Draw("colz");
  gPad->SetRightMargin(0.12);

  //Acceptance cut line
  TLine* l1 = new TLine(0.0,3.5,2.4,3.5);
  l1->SetLineColor(kRed);
  l1->SetLineWidth(2);
  l1->Draw("same");

  cData->SaveAs("PbPbAcceptance.png");
  cData->SaveAs("PbPbAcceptance.pdf");

  //Save histograms.
  TFile* histfile2d = new TFile("histfile2d.root","RECREATE");
  hReco->Write();
  //hMC->Write();
  histfile2d->Close();
  delete histfile2d;

  //inputMC->Close();
  //delete inputMC;

}

