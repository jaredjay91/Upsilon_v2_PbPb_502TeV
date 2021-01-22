#include <ctime>

#include <TLorentzVector.h>
#include "../../HeaderFiles/commonUtility.h"
#include "../../HeaderFiles/HiEvtPlaneList.h"
#include "../../HeaderFiles/cutsAndBinUpsilonV2.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"

const Double_t pi = 3.141592653589;

TF1* fWgtAA1 = new TF1("fWgtAA1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fdNdpTWgt = new TF1("fdNdpTWgt","([0]+[1]*x+[2]*x*x)/((x-[3])*(x-[3])*(x-[3]))");

void getAcceptanceFrom19001Skim(int nevt=-1, int dateStr=20210122, int weighted=0) 
{

  gStyle->SetOptStat(0);

  using namespace std;
  using namespace hi;

  const int numptbins = 50;
  const int arraysize = numptbins+1;
  Double_t ptbins[arraysize];
  Double_t ptmax = 50.0;
  Double_t ptmin = 0.0;
  for (int i=0; i<arraysize; i++){
    ptbins[i] = ptmax*i/((Double_t)arraysize);
  }
  Double_t ptbinsInt[2] = {ptmin,ptmax};
  const int numptbinsInt = 1;

  // Integrated bin:
  TH1D* hDenInt = new TH1D("hDenInt",";p_{T} (GeV/c);Acceptance of Dimuons",1,0,50);
  TH1D* hAccInt = new TH1D("hAccInt",";p_{T} (GeV/c);Acceptance of Dimuons",1,0,50);
  hDenInt->Sumw2(); hDenInt->SetMinimum(0); hDenInt->SetMaximum(1.0);
  hAccInt->Sumw2(); hAccInt->SetMinimum(0); hAccInt->SetMaximum(1.0);
  hDenInt->SetMarkerStyle(20);hDenInt->SetMarkerColor(kBlue+2);
  hAccInt->SetMarkerStyle(20);hAccInt->SetMarkerColor(kBlue+2);

  TH1D* hDenIntNoW = new TH1D("hDenIntNoW",";p_{T} (GeV/c);Acceptance of Dimuons",1,0,50);
  TH1D* hAccIntNoW = new TH1D("hAccIntNoW",";p_{T} (GeV/c);Acceptance of Dimuons",1,0,50);
  hDenIntNoW->Sumw2(); hDenIntNoW->SetMinimum(0); hDenIntNoW->SetMaximum(1.0);
  hAccIntNoW->Sumw2(); hAccIntNoW->SetMinimum(0); hAccIntNoW->SetMaximum(1.0);
  hDenIntNoW->SetMarkerStyle(24);hDenIntNoW->SetMarkerColor(kBlue+2);
  hAccIntNoW->SetMarkerStyle(24);hAccIntNoW->SetMarkerColor(kBlue+2);

  // pT dependence:
  TH1D* hDenPt = new TH1D("hDenPt",";p_{T} (GeV/c);Acceptance of Dimuons",50,0,50);
  TH1D* hAccPt = new TH1D("hAccPt",";p_{T} (GeV/c);Acceptance of Dimuons",50,0,50);
  hDenPt->Sumw2(); hDenPt->SetMinimum(0); hDenPt->SetMaximum(1.0);
  hAccPt->Sumw2(); hAccPt->SetMinimum(0); hAccPt->SetMaximum(1.0);
  hDenPt->SetMarkerStyle(20);hDenPt->SetMarkerColor(kBlue+2);
  hAccPt->SetMarkerStyle(20);hAccPt->SetMarkerColor(kBlue+2);

  TH1D* hDenPtNoW = new TH1D("hDenPtNoW",";p_{T} (GeV/c);Acceptance of Dimuons",50,0,50);
  TH1D* hAccPtNoW = new TH1D("hAccPtNoW",";p_{T} (GeV/c);Acceptance of Dimuons",50,0,50);
  hDenPtNoW->Sumw2(); hDenPtNoW->SetMinimum(0); hDenPtNoW->SetMaximum(1.0);
  hAccPtNoW->Sumw2(); hAccPtNoW->SetMinimum(0); hAccPtNoW->SetMaximum(1.0);
  hDenPtNoW->SetMarkerStyle(24);hDenPtNoW->SetMarkerColor(kBlue+2);
  hAccPtNoW->SetMarkerStyle(24);hAccPtNoW->SetMarkerColor(kBlue+2);

  TFile *fweight = TFile::Open("Func_dNdpT_1S.root","READ");
  TF1* fitRatio = (TF1*)fweight->Get("fitRatio");
  double fdNdpTWgt_p0 = fitRatio->GetParameter(0);
  double fdNdpTWgt_p1 = fitRatio->GetParameter(1);
  double fdNdpTWgt_p2 = fitRatio->GetParameter(2);
  double fdNdpTWgt_p3 = fitRatio->GetParameter(3);
  fWgtAA1->SetParameters( 1.0001, 5.1, 2.0024, 12.4243);
  fdNdpTWgt->SetParameters( fdNdpTWgt_p0, fdNdpTWgt_p1, fdNdpTWgt_p2, fdNdpTWgt_p3);

  TString inFileName = "skimedForAcc_MC_Ups1S_20170808.root";
  TFile* inFile = TFile::Open(inFileName,"READ");
  TTree* mm = (TTree*)inFile->Get("mmGen");

  int nevtReal = mm->GetEntries();

  cout << "Total events = " << nevtReal << ", : " << mm->GetEntries() << endl;

  float mass = 0.0, pt = 0.0, phi = 0.0, y = 0.0, eta = 0.0;
  float pt1 = 0.0, phi1 = 0.0, eta1 = 0.0;
  float pt2 = 0.0, phi2 = 0.0, eta2 = 0.0;

  TBranch *b_mass;
  TBranch *b_y;
  TBranch *b_pt;
  TBranch *b_phi;
  TBranch *b_eta;
  TBranch *b_pt1;
  TBranch *b_phi1;
  TBranch *b_eta1;
  TBranch *b_pt2;
  TBranch *b_phi2;
  TBranch *b_eta2;

  mm->SetBranchAddress("mass", &mass, &b_mass);
  mm->SetBranchAddress("pt", &pt, &b_pt);
  mm->SetBranchAddress("phi", &phi, &b_phi);
  mm->SetBranchAddress("y", &y, &b_y);
  mm->SetBranchAddress("eta", &eta, &b_eta);
  mm->SetBranchAddress("pt1", &pt1, &b_pt1);
  mm->SetBranchAddress("eta1", &eta1, &b_eta1);
  mm->SetBranchAddress("phi1", &phi1, &b_phi1);
  mm->SetBranchAddress("pt2", &pt2, &b_pt2);
  mm->SetBranchAddress("eta2", &eta2, &b_eta2);
  mm->SetBranchAddress("phi2", &phi2, &b_phi2);

  // event loop start
  for(int iev=0; iev<nevtReal ; ++iev)
  {
    if(iev%10000==0) cout << ">>>>> EVENT " << iev << " / " << mm->GetEntries() <<  " ("<<(int)(100.*iev/nevtReal) << "%)" << endl;

    mm->GetEntry(iev);

    //Count how many dimuons there are.
    if (pt>=50 || fabs(y)>=2.4) continue;

    Double_t aawgt = fWgtAA1->Eval(pt); // apply weighting factor from functions defined above for AA
    aawgt=aawgt*fdNdpTWgt->Eval(pt);

    hDenInt->Fill(pt);
    hDenPt->Fill(pt,aawgt);
    hDenIntNoW->Fill(pt);
    hDenPtNoW->Fill(pt);

    //Count how many of them have both muons in the single-muon acceptance.
    if ( pt1>3.5 && pt2>3.5 && fabs(eta1)<2.4 && fabs(eta2)<2.4 ){
      //cout << "(pt,aawgt) = (" << pt << ", " << aawgt << ")" << endl;
      hAccInt->Fill(pt,aawgt);
      hAccPt->Fill(pt,aawgt);
      hAccIntNoW->Fill(pt);
      hAccPtNoW->Fill(pt);
    }


  } //end of loop

  hAccInt->Divide(hDenInt);
  hAccPt->Divide(hDenPt);
  hAccIntNoW->Divide(hDenIntNoW);
  hAccPtNoW->Divide(hDenPtNoW);

  TCanvas* c1 = new TCanvas("c1","c1",0,0,500,500);
  hAccPtNoW->Draw("PE");
  hAccPt->Draw("same PE");

  TLegend* leg = new TLegend(0.2,0.7,0.5,0.8);
  //leg->SetTextSize(19);
  //leg->SetTextFont(43);
  leg->SetBorderSize(0);
  leg->AddEntry(hAccPtNoW,"unweighted","pe");
  leg->AddEntry(hAccPt,"weighted","pe");
  leg->Draw("same");

  c1->SaveAs(Form("AcceptancePlot_%i.png",dateStr));
  c1->SaveAs(Form("AcceptancePlot_%i.pdf",dateStr));

  TFile* accFile = new TFile(Form("acceptance_%i.root",dateStr),"RECREATE");
  hAccInt->Write();
  hAccPt->Write();
  hAccIntNoW->Write();
  hAccPtNoW->Write();
  accFile->Close();
}

