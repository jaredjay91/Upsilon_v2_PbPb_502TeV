#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "../HeaderFiles/cutsAndBin.h"


void Fitv2(

  int collId = kAADATA,
  float ptLow = 0, float ptHigh = 50,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 10, int cHigh = 90,
  float muPtCut = 3.5,
  double* v2ValPtr=0, double* v2ErrPtr=0,
  int whichUpsilon=1, int whichSyst=0 )

{

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";

  gStyle->SetOptStat(0);

  TString fileLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);

  //Load the yields vs phi
  TFile* yieldsFile = TFile::Open(Form("ExtractedYields/%syieldsVsPhi_%iS_%s.root", systStr.Data(), whichUpsilon, fileLabel.Data()), "READ");
  TH1D* yieldsVsPhi = (TH1D*)yieldsFile->Get("yieldsVsPhi;1");

  //Plot it
  TCanvas* c1 = new TCanvas("c1","c1",50,50,550,550);
  yieldsVsPhi->Draw();
  yieldsVsPhi->GetXaxis()->SetTitle("|#Delta#phi|/#pi");
  double totalintegral = yieldsVsPhi->Integral(1,4);
  yieldsVsPhi->Scale(1.0/totalintegral);

  //Fit it
  /*TF1* fitfunc = new TF1("fitfunc","[0]*( 1 + 2*[1]*cos(2*x*3.14159265) + 2*[2]*cos(3*x*3.14159265) + 2*[3]*cos(4*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("Amp","v2","v3","v4");*/
  TF1* fitfunc = new TF1("fitfunc","[0]*( 1 + 2*[1]*cos(2*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("Amp","v2");
  /*TF1* fitfunc = new TF1("fitfunc"," 1 + 2*[1]*cos(2*x*3.14159265)",0,0.5);
  fitfunc->SetParNames("v2");*/
  yieldsVsPhi->Fit("fitfunc");

  double v2Val = fitfunc->GetParameter(1);
  double v2Err = fitfunc->GetParError(1);
  cout << "v2 = " << v2Val << " +/- " << v2Err << endl;
  TLatex latex;
  latex.SetTextSize(0.05);
  //latex.SetTextAlign(13);
  latex.DrawLatex(0.25,50,Form("v_{2} = %.3f #pm %.3f",v2Val,v2Err));

  c1->Update();

  c1->SaveAs(Form("fitPlots/%s_v2_fit_pt%.1f-%.1f_y%.2f-%.2f_cent%i-%i.png", systStr.Data(), ptLow,ptHigh, yLow,yHigh, cLow, cHigh));

  *v2ValPtr = v2Val;
  *v2ErrPtr = v2Err;
  
  yieldsFile->Close();
  delete fitfunc;
  delete c1;
}
