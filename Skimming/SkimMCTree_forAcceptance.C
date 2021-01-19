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


static const long MAXTREESIZE = 1000000000000;

void SkimMCTree_forAcceptance(int nevt=-1, int dateStr=20210111) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  Double_t ptbins[4] = {0,3,6,30};
  const int numptbins = sizeof(ptbins)/sizeof(Double_t)-1;

  Double_t centbins[5] = {0,10,30,50,100};
  const int numcentbins = sizeof(centbins)/sizeof(Double_t)-1;

  //TH1D* hPseudoQx = new TH1D("hPseudoQx","cos(2*psi)",50,-1.2,1.2);
  //TH1D* hPseudoQy = new TH1D("hPseudoQy","sin(2*psi)",50,-1.2,1.2);

  //TString fnameData1 = "../Oniatree_Ups1SMM_5p02TeV_TuneCP5_Embd_Gen_MC_190610.root";
  //TFile* MCfile = TFile::Open(fnameData1,"READ");
  //TTree *mytree = (TTree*)MCfile->Get("myTree");
  //TTree* tree = (TTree*)MCfile->Get("tree");
  //mytree->AddFriend(tree);

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
  Int_t           Gen_QQ_size;
  Int_t           Gen_mu_size;
  Int_t           Gen_mu_whichGen[maxBranchSize];
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  ULong64_t       Gen_QQ_trig[maxBranchSize];   //[Gen_QQ_size]
  Float_t         Gen_QQ_VtxProb[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_mu_whichGen;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_mu_4mom;   //!
  TBranch        *b_Gen_QQ_trig;   //!
  TBranch        *b_Gen_QQ_VtxProb;   //!

  Bool_t          Gen_mu_highPurity[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Gen_mu_highPurity;   //!
  mytree->SetBranchAddress("Gen_mu_highPurity", Gen_mu_highPurity, &b_Gen_mu_highPurity);

  Gen_QQ_4mom = 0;
  Gen_mu_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_mu_whichGen", Gen_mu_whichGen, &b_Gen_mu_whichGen);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
  mytree->SetBranchAddress("Gen_QQ_trig", Gen_QQ_trig, &b_Gen_QQ_trig);
  mytree->SetBranchAddress("Gen_QQ_VtxProb", Gen_QQ_VtxProb, &b_Gen_QQ_VtxProb);

  //  muon id 
  Int_t           Gen_QQ_mupl_idx[maxBranchSize];
  Int_t           Gen_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx,&b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx,&b_Gen_QQ_mumi_idx);

  Int_t           Gen_mu_nTrkHits[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Gen_mu_nTrkHits", Gen_mu_nTrkHits, &b_Gen_mu_nTrkHits);
  Float_t         Gen_mu_normChi2_global[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Gen_mu_normChi2_global", Gen_mu_normChi2_global, &b_Gen_mu_normChi2_global);
  Int_t           Gen_mu_nMuValHits[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Gen_mu_nMuValHits", Gen_mu_nMuValHits, &b_Gen_mu_nMuValHits);
  Int_t           Gen_mu_StationsMatched[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Gen_mu_StationsMatched", Gen_mu_StationsMatched, &b_Gen_mu_StationsMatched);
  Float_t         Gen_mu_dxy[maxBranchSize];   //[Gen_mu_size]
  Float_t         Gen_mu_dxyErr[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_dxy;   //!
  TBranch        *b_Gen_mu_dxyErr;   //!
  mytree->SetBranchAddress("Gen_mu_dxy", Gen_mu_dxy, &b_Gen_mu_dxy);
  mytree->SetBranchAddress("Gen_mu_dxyErr", Gen_mu_dxyErr, &b_Gen_mu_dxyErr);
  Float_t         Gen_mu_dz[maxBranchSize];   //[Gen_mu_size]
  Float_t         Gen_mu_dzErr[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_dz;   //!
  TBranch        *b_Gen_mu_dzErr;   //!
  mytree->SetBranchAddress("Gen_mu_dz", Gen_mu_dz, &b_Gen_mu_dz);
  mytree->SetBranchAddress("Gen_mu_dzErr", Gen_mu_dzErr, &b_Gen_mu_dzErr);
  Int_t           Gen_mu_nTrkWMea[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Gen_mu_nTrkWMea", Gen_mu_nTrkWMea, &b_Gen_mu_nTrkWMea);
  Bool_t          Gen_mu_TMOneStaTight[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_TMOneStaTight;   //!

  mytree->SetBranchAddress("Gen_mu_TMOneStaTight", Gen_mu_TMOneStaTight, &b_Gen_mu_TMOneStaTight);
  Int_t           Gen_mu_nPixWMea[maxBranchSize];   //[Gen_mu_size]
  TBranch        *b_Gen_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Gen_mu_nPixWMea", Gen_mu_nPixWMea, &b_Gen_mu_nPixWMea);
  Int_t           Gen_QQ_sign[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Gen_QQ_sign;   //!
  mytree->SetBranchAddress("Gen_QQ_sign", Gen_QQ_sign, &b_Gen_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Gen_mu_nPixValHits[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Gen_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Gen_mu_nPixValHits", Gen_mu_nPixValHits, &b_Gen_mu_nPixValHits);
  Float_t         Gen_mu_ptErr_global[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Gen_mu_ptErr_global;   //!
  mytree->SetBranchAddress("Gen_mu_ptErr_global", Gen_mu_ptErr_global, &b_Gen_mu_ptErr_global);

  Int_t           Gen_mu_SelectionType[maxBranchSize];
  TBranch        *b_Gen_mu_SelectionType;
  mytree->SetBranchAddress("Gen_mu_SelectionType", Gen_mu_SelectionType, &b_Gen_mu_SelectionType);


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
  
  TFile* newfile;
  newfile = new TFile(Form("skims/newOniaTree_Skim_NoTrig_MC_forAcceptance_%i.root",dateStr),"recreate");

  DiMuon dm;
  TTree* mmtree = new TTree("mmep","dimuonAndEventPlanes");
  mmtree->SetMaxTreeSize(MAXTREESIZE);
  mmtree->Branch("mm",&dm.run,branchString.Data());
  

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;


  int kTrigSel = 13;
  int DIMUIDPASS = 0;
  int ALLPASS = 0;

  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << ", : " << mytree->GetEntries() << endl;

  // event loop start
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    //if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    for (Int_t irqq=0; irqq<Gen_QQ_size; ++irqq) 
    {

      dm.clear();      // clear the output tree: 
      dm.run = runNb;
      dm.lumi = LS ;
      dm.event = eventNb ;
      dm.vz = zVtx;
      dm.cBin = Centrality ;

      JP_Reco = (TLorentzVector*) Gen_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[irqq]);

      dm.phi    = JP_Reco->Phi();
      dm.ep2 = epang[8];
      dm.dphiEp2 = getDPHI_Jared( dm.phi, dm.ep2);

      dm.mass   = JP_Reco->M();
      dm.pt     = JP_Reco->Pt();

      dm.y      = JP_Reco->Rapidity();
      dm.eta      = JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();
      dm.weight = 1.;


      mmtree->Fill();

    } // end of dimuon loop
  } //end of event loop

  mmtree->Write();  // Don't need to call Write() for trees
  newfile->Close();

}

