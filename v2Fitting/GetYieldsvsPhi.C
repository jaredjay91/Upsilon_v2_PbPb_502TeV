#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "TCanvas.h"


const int numphibins = 4;
float phibins[5] = {0.0, 0.125, 0.25, 0.375, 0.5};


void GetYieldsvsPhi(

  int collId = kAADATA,
  float ptLow = 0, float ptHigh = 6,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 100,
  float muPtCut = 4.0,
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

  //Read in yields from fit files
  TString directory = "../SignalFitting/dphiFits_R4a/";
  for (int iphi=0; iphi<numphibins; iphi++) {
    dphiEp2Low = phibins[iphi];
    dphiEp2High = phibins[iphi+1];
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString inFileName = Form("%s%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());
    cout << "Opening file: " << inFileName << endl;
    TFile* inFile = TFile::Open(inFileName,"READ");
    RooWorkspace *ws = (RooWorkspace*)inFile->Get("workspace");
    inFile->Close("R");
    float yield = ws->var(Form("nSig%is",whichUpsilon))->getVal();
    float yielderr = ws->var(Form("nSig%is",whichUpsilon))->getError();
    //delete ws;
    delete inFile;
    cout << yield << " +/- " << yielderr << endl;
    yieldsVsPhi->SetBinContent(iphi+1, yield);
    yieldsVsPhi->SetBinError(iphi+1, yielderr);
  }

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
