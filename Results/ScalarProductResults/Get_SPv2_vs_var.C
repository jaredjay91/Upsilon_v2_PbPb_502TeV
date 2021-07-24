#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"

//{v2 point, top of stat. unc., top of syst. unc.}
double ptvals[4][3] = {
{0.00421687, 0.0232932, 0.0157631},
{-0.0178715, -0.0369478, -0.0223896},
{0.0579317, 0.0805221, 0.0619478},
{0.0182731, 0.0388554, 0.0247992}};
double centvals[3][3] = {
{0.00544933, 0.0215105, 0.0106119},
{0.0249522, 0.0427342, 0.03413},
{-0.016348, -0.0478967, -0.0272467}};

void get_v2_val_pt(int whichBin, double* v2ValPtr, double* v2ErrPtr, double* v2SystPtr) {
  *v2ValPtr = ptvals[whichBin][0];
  //The histogram draws the error bar with a size equal to whatever you set it to. It is not a percent error, but an absolute error.
  *v2ErrPtr = ptvals[whichBin][1]-ptvals[whichBin][0];
  *v2SystPtr = (ptvals[whichBin][2]-ptvals[whichBin][0])/ptvals[whichBin][0];
}

void get_v2_val_cent(int whichBin, double* v2ValPtr, double* v2ErrPtr, double* v2SystPtr) {
  *v2ValPtr = centvals[whichBin][0];
  *v2ErrPtr = centvals[whichBin][1]-centvals[whichBin][0];
  *v2SystPtr = (centvals[whichBin][2]-centvals[whichBin][0])/centvals[whichBin][0];
}

void Get_SPv2_vs_var(int whichUpsilon=1) {

  gStyle->SetOptStat(0);

  const int numptbins = 4;
  float ptbins[5] = {0,3,6,10,50};

  const int numcbins = 3;
  float cbins[4] = {10,30,50,90};

  double v2Val; double v2Err; double v2Syst;
  double* v2ValPtr; double* v2ErrPtr; double* v2SystPtr;
  v2ValPtr = &v2Val; v2ErrPtr = &v2Err; v2SystPtr = &v2Syst;

  float ptLow; float ptHigh;
  int cLow; int cHigh;

  TH1D* hv2pt = new TH1D("hv2pt","hv2pt",numptbins,ptbins);
  TH1D* hv2ptSyst = new TH1D("hv2ptSyst","hv2ptSyst",numptbins,ptbins);
  TH1D* hv2c = new TH1D("hv2c","hv2c",numcbins,cbins);
  TH1D* hv2cSyst = new TH1D("hv2cSyst","hv2cSyst",numcbins,cbins);

  //PT BINS
  for (int ipt=0; ipt<numptbins; ipt++) {
    get_v2_val_pt(ipt, v2ValPtr, v2ErrPtr, v2SystPtr);
    hv2pt->SetBinContent(ipt+1,v2Val);
    hv2pt->SetBinError(ipt+1,v2Err);
    hv2ptSyst->SetBinContent(ipt+1,v2Syst);
  }
  
  TCanvas* cv2 = new TCanvas("cv2","cv2",800,400);
  cv2->Divide(2,1);
  cv2->cd(1);
  hv2pt->SetMinimum(-0.05);
  hv2pt->SetMaximum(0.2);
  //hv2pt->Sumw2();
  hv2pt->Draw();

  //CENTRALITY BINS
  for (int ic=0; ic<numcbins; ic++) {
    get_v2_val_cent(ic, v2ValPtr, v2ErrPtr, v2SystPtr);
    hv2c->SetBinContent(ic+1,v2Val);
    hv2c->SetBinError(ic+1,v2Err);
    hv2cSyst->SetBinContent(ic+1,v2Syst);
  }
  cv2->cd(2);
  hv2c->SetMinimum(-0.15);
  hv2c->SetMaximum(0.15);
  hv2c->Draw();

  TString outFileName = Form("Ups_%i_SPv2_cent10-90.root", whichUpsilon);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hv2pt->Write();
  hv2ptSyst->Write();
  hv2c->Write();
  hv2cSyst->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;

}
