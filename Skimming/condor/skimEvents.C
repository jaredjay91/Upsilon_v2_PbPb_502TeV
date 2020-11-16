 
void skimEvents(int filenum=80)
{

  TString fnameDataReReco = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part*.root";

  // Get old file, old tree and set top branch address
  TChain* theChain = new TChain("myTree");
  theChain->Add(fnameDataReReco.Data());

  TChain* tree = new TChain("tree"); 
  tree->Add(fnameDataReReco.Data());

  theChain->AddFriend(tree);

  int kTrigSel = 13;
  const int maxBranchSize = 1000;

  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       HLTriggers;
  Int_t           Reco_mu_size;

  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_mu_size;

  theChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  theChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  theChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);

  const int nentries = theChain->GetEntries()/100;
  //const int nentries = 10000;
  const int firstevt = nentries*filenum;
  const int lastevt = nentries*(filenum+1);
  cout << "nentries = " << nentries << endl;
  const int stepSize = nentries/1000;

  //Event *event = nullptr;
  //theChain->SetBranchAddress("event", &event);
 
  // Create a new file + a clone of old tree in new file
  TFile newfile(Form("skims/skimmed_HLT13_DoubleMuon_ReReco_Oniatree_addvn%i.root",filenum), "recreate");
  //TFile newfile("/eos/cms/store/group/phys_heavyions/dileptons/jjay/DoubleQuarkonia/skimmedHLT13_4muons_PromptAOD_v1v2_Oniatree_addvn.root", "recreate");
  auto newmytree = theChain->CloneTree(0);
  auto newtree = tree->CloneTree(0);
 
  cout << "Now entering loop." << endl;
  for(int iev=firstevt; iev<lastevt; ++iev)
  {
    if(iev%stepSize==0) cout << ">>>>> EVENT " << iev << " / " << theChain->GetEntries() <<  " ("<<(int)(100.*(iev-firstevt)/nentries) << "%)" << endl;
    //cout << "Cloning entry " << iev << endl;
    theChain->GetEntry(iev);

    //if(Reco_mu_size<4) continue;
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    //if(!( (Reco_QQ_trig[iev]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;//this is a dimuon trigger, not an event trigger.

    newmytree->Fill();
    newtree->Fill();
    //event->Clear();
  }

  const int newentries = newtree->GetEntries();
  cout << "new tree has " << newentries << " events" << endl;
  cout << "That's a " << 100.0*(nentries-newentries)/(double)nentries  << "% reduction!" << endl;
  //newmytree->Print();
  newfile.Write();
}
