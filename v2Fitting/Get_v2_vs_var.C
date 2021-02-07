#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "Fitv2.C"


void Get_v2_vs_var(int whichUpsilon=1, int whichSyst=0) {

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";

  gStyle->SetOptStat(0);

  const int numybins = 3;
  float ybins[7] = {0.0, 0.8, 1.6, 2.4};

  const int numptbins = 4;
  float ptbins[5] = {0,3,6,10,50};

  const int numcbins = 3;
  float cbins[4] = {10,30,50,90};

  int collId = kAADATA;
  float muPtCut = 3.5;

  double v2Val; double v2Err;
  double* v2ValPtr; double* v2ErrPtr;
  v2ValPtr = &v2Val; v2ErrPtr = &v2Err;

  float ptLow; float ptHigh;
  float yLow; float yHigh;
  int cLow; int cHigh;

  TH1D* hv2pt = new TH1D("hv2pt","hv2pt",numptbins,ptbins);
  TH1D* hv2y = new TH1D("hv2y","hv2y",numybins,ybins);
  TH1D* hv2c = new TH1D("hv2c","hv2c",numcbins,cbins);

  //PT BINS
  for (int ipt=0; ipt<numptbins; ipt++) {
    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
    cout << "[" << ptLow << "," << ptHigh << "]" << endl;
    yLow = 0.0; yHigh = 2.4;
    cLow = 10; cHigh = 90;

    Fitv2(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, v2ValPtr, v2ErrPtr, whichUpsilon, whichSyst);
    hv2pt->SetBinContent(ipt+1,v2Val);
    hv2pt->SetBinError(ipt+1,v2Err);
  }
  
  TCanvas* cv2 = new TCanvas("cv2","cv2",1200,400);
  cv2->Divide(3,1);
  cv2->cd(1);
  hv2pt->Draw();

  //Y BINS
  for (int iy=0; iy<numybins; iy++) {
    yLow = ybins[iy];
    yHigh = ybins[iy+1];
    cout << "[" << yLow << "," << yHigh << "]" << endl;
    ptLow = 0; ptHigh = 50;
    cLow =10; cHigh = 90;

    Fitv2(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, v2ValPtr, v2ErrPtr, whichUpsilon, whichSyst);
    hv2y->SetBinContent(iy+1,v2Val);
    hv2y->SetBinError(iy+1,v2Err);
  }
  cv2->cd(2);
  hv2y->Draw();

  //CENTRALITY BINS
  for (int ic=0; ic<numcbins; ic++) {
    cLow = cbins[ic];
    cHigh = cbins[ic+1];
    cout << "[" << cLow << "," << cHigh << "]" << endl;
    ptLow = 0; ptHigh = 50;
    yLow =0.0; yHigh = 2.4;

    Fitv2(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, v2ValPtr, v2ErrPtr, whichUpsilon, whichSyst);
    hv2c->SetBinContent(ic+1,v2Val);
    hv2c->SetBinError(ic+1,v2Err);
  }
  cv2->cd(3);
  hv2c->Draw();

  TString outFileName = Form("Ups_%i_v2_cent10-90_%s.root", whichUpsilon, systStr.Data());
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hv2pt->Write();
  hv2y->Write();
  hv2c->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;

}
