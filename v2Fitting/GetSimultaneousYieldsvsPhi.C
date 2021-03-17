#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "TCanvas.h"


const int numphibins = 4;
float phibins[5] = {0.0, 0.125, 0.25, 0.375, 0.5};


void GetSimultaneousYieldsvsPhi(

  int collId = kAADATA,
  float ptLow = 0, float ptHigh = 3,
  float yLow = 2.1, float yHigh = 2.4,
  int cLow = 10, int cHigh = 90,
  float muPtCut = 3.5,
  int whichUpsilon = 1,
  int whichSyst = 0 )

{


  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  float dphiEp2Low; float dphiEp2High;

  TString fileLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);

  //Create histogram
  TH1D* yieldsVsPhi = new TH1D("yieldsVsPhi", "yieldsVsPhi",numphibins,phibins);

  dphiEp2Low = 0.0;
  dphiEp2High = 0.5;
  TString directory = "../SignalFitting/Simultaneous_dPhiFits/";
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString inFileName = Form("%ssim_dphi%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());
  cout << "Opening file: " << inFileName << endl;
  TFile* inFile = TFile::Open(inFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)inFile->Get("workspace");

  //Read in yields from fit files
  for (int iphi=0; iphi<numphibins; iphi++) {
    float yield = ws->var(Form("nSig%is[%i]",whichUpsilon,iphi))->getVal();
    float yielderr = ws->var(Form("nSig%is[%i]",whichUpsilon,iphi))->getError();
    //delete ws;
    cout << yield << " +/- " << yielderr << endl;
    yieldsVsPhi->SetBinContent(iphi+1, yield);
    yieldsVsPhi->SetBinError(iphi+1, yielderr);
  }

  inFile->Close("R");
  delete inFile;

  TCanvas* c1 = new TCanvas("c1","c1",50,50,550,550);
  yieldsVsPhi->Draw();
  yieldsVsPhi->SetMinimum(0);

  TFile* outFile = new TFile(Form("ExtractedYields/%syieldsVsPhi_%iS_%s.root", systStr.Data(), whichUpsilon, fileLabel.Data()),"RECREATE");
  yieldsVsPhi->Write();
  outFile->Close();

  delete yieldsVsPhi;
  delete c1;
  delete outFile;

  //cout << endl << "Here's what's in memory right now" << endl;
  //gDirectory->ls("-m");
  //cout << "that's all." << endl << endl;

}
