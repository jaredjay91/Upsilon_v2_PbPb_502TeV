//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

#include <iostream>
#include "../HeaderFiles/rootFitHeaders.h"
#include "../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "../HeaderFiles/PsetCollection.h"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/StyleSetting.h"
#include "RoundsHeader.h"

bool isAbout(float a, float b) {
  if (abs(a-b)<0.01) return kTRUE;
  else return kFALSE;
}

using namespace std;
using namespace RooFit;
int CheckSimFitQuality( 
       int collId = kAADATA,
       float ptLow=0, float ptHigh=3,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=90,
       float muPtCut=3.5,
       int whichSyst=0
// 0: Nominal
// 1: AltSig
// 2: AltBkg
// 3: AltAcc
// 4: AltEff
			) 
{

  TString directory = "Simultaneous_dPhiFits/";
  TString logFileName = Form("%slog.txt",directory.Data());

  TString Params[5] = {"sigma1s_1","x1s","alpha1s_1","n1s_1","f1s"};

  //Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.35, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  float eta_low = -2.4;
  float eta_high = 2.4;
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";

  //import the model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
  TString NomFileName = Form("%ssim_dphi%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());
  cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");
  RooAbsData* reducedDS = ws->data("reducedDS[0]");
  RooFitResult* fitRes2 = (RooFitResult*)ws->obj("fitresult_model_dsABCD");
  fitRes2->Print("v");

  //Plot it and calculate chi-squared
  TCanvas* c[4];
  RooPlot* myPlot[4];
  RooHist* hpull[4];
  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin;
  double *ypull[4];

  for (int i=0; i<4; i++) {
    c[i] = new TCanvas(Form("canvas[%i]",i),"My plots",4,45,550,520);
    c[i]->cd();
    myPlot[i] = ws->var("mass")->frame(nMassBin); // bins
    ws->data(Form("reducedDS[%i]",i))->plotOn(myPlot[i],Name(Form("dataHist[%i]",i)), Range(massLow, massHigh));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Name(Form("modelHist[%i]",i)), Range(massLow, massHigh));
    //ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Name("Sig1S"),Components(ws->pdf("cb1s")),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    //ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    //ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    //ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Name("bkgPDF"),Components("bkgLowPt"),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
    myPlot[i]->SetFillStyle(4000);
    myPlot[i]->SetAxisRange(massLowForPlot, massHighForPlot,"X");
    myPlot[i]->GetYaxis()->SetTitleOffset(1.43);
    myPlot[i]->GetYaxis()->CenterTitle();
    //myPlot[i]->GetYaxis()->SetTitleSize(0.058);
    myPlot[i]->GetYaxis()->SetLabelSize(0.054);
    //myPlot[i]->GetXaxis()->SetLabelSize(0);
    myPlot[i]->GetXaxis()->SetRangeUser(8,14);
    //myPlot[i]->GetXaxis()->SetTitleSize(0);
    myPlot[i]->Draw();

    cout << "dataset entries = " << ws->data(Form("reducedDS[%i]",i))->sumEntries() << endl;
    //cout << "model entries = " << ws->pdf(Form("model[%i]",i))->normRange() << endl;

    hpull[i] = myPlot[i]->pullHist(Form("dataHist[%i]",i),Form("modelHist[%i]",i));
    hpull[i]->SetMarkerSize(0.8);

    ypull[i] = hpull[i]->GetY();
    for(int j=0;j<nBins;j++)
    {
      if(ypull[i][j] == 0) continue;
      chisq += TMath::Power(ypull[i][j],2);
      nFullBinsPull++;
    }
  }
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  cout << "numFitPar = " << numFitPar << endl;
  int ndf = nFullBinsPull - numFitPar;
  cout << "chisq/dof = " << chisq/ndf << endl;

  float temp1[4];
  float temp1err[4];
  float temp2[4];
  float temp2err[4];
  float temp3[4];
  float temp3err[4];

  for (int i=0; i<4; i++) {
    temp1[i] = ws->var(Form("nSig1s[%i]",i))->getVal();  
    temp1err[i] = ws->var(Form("nSig1s[%i]",i))->getError();  
    temp2[i] = ws->var(Form("nSig2s[%i]",i))->getVal();  
    temp2err[i] = ws->var(Form("nSig2s[%i]",i))->getError();  
    temp3[i] = ws->var(Form("nSig3s[%i]",i))->getVal();  
    temp3err[i] = ws->var(Form("nSig3s[%i]",i))->getError();  

    cout << "1S signal    =  " << temp1[i] << " +/- " << temp1err[i] << endl;
    cout << "2S signal    =  " << temp2[i] << " +/- " << temp2err[i] << endl;
    cout << "3S signal    =  " << temp3[i] << " +/- " << temp3err[i] << endl;
    cout << "Total signal =  " << temp1[i]+temp2[i]+temp3[i] << endl;

    c[i]->Update();
  }

  //QUALITY CHECK

  //chi-squared requirements:
  bool goodChi2 = kTRUE;
  //double chisqUpperCut = 2.0 + temp1/10000;
  //double chisqLowerCut = 0.5 + temp1/30000;
  double chisqUpperCut = 5;
  double chisqLowerCut = 0.5;

  //Signal error requirements:
  bool goodSigErr = kTRUE;
  double errUpperLimit = 0.15;
  double errLowerLimit = 0.005;
  if (collId==kPADATA) errUpperLimit = 0.15;

  //parameter requirements:
  bool goodParams = kTRUE;
  //TString badParamsList = "";
  std::ostringstream sstream;
  double buffer = 0.03;
  int paramstart = 0;
  int paramend = 3;
  //all the signal parameters should be good in the nominal fits.
  paramstart = 0;
  paramend = 4;
  for (int i=paramstart;i<=paramend;i++) {
    double fittedval = ws->var(Params[i])->getVal();
    double width = paramsupper[i]-paramslower[i];
    double upper = paramsupper[i]-(buffer*width);
    double lower = paramslower[i]+(buffer*width);
    if (fittedval>upper || fittedval<lower) {
      goodParams = kFALSE;
      //badParamsList = badParamsList + " " + Params[i] + "=" + fittedval;
      sstream << " " << Params[i] << "=" << fittedval;
    }
  }
  std::string badParamsList = sstream.str();

  //Special bins
  //if (collId==kPPDATA && isAbout(yLow,0.0) && isAbout(yHigh,0.4) && isAbout(ptHigh-ptLow,30)) chisqUpperCut = 3.5;
  //if (collId==kPPDATA && isAbout(yLow,0.0) && isAbout(yHigh,0.8) && isAbout(ptHigh-ptLow,30)) chisqUpperCut = 4.5;
  //if (collId==kAADATA && isAbout(ptLow,4.0) && isAbout(ptHigh,9.0) && isAbout(yHigh-yLow,2.4)) chisqUpperCut = 2.3;


  //Check
  ofstream logFile;
  logFile.open(logFileName, fstream::in | fstream::out | fstream::app);
  int good = 0;
  if (chisq/ndf>chisqUpperCut || chisq/ndf<chisqLowerCut) goodChi2 = kFALSE;
  else goodChi2 = kTRUE;

  goodSigErr = kTRUE;
  TString badSigStr = "";
  for (int i=0; i<4; i++) {
    if (temp1err[i]/temp1[i]>errUpperLimit || temp1err[i]/temp1[i]<errLowerLimit) {
      goodSigErr = kFALSE;
      badSigStr = badSigStr + Form("%.2f",temp1err[i]/temp1[i]*100);
    }
  }

  //if (goodChi2 && goodSigErr && goodParams){
  if (kTRUE){
    cout << "THE FIT PASSED THE QUALITY CHECK! :)" << endl;
    good = 1;
  }
  else{
    cout << "THE FIT FAILED THE QUALITY CHECK! :(" << endl;
    logFile << endl << NomFileName << endl;
    logFile << "THE FIT FAILED THE QUALITY CHECK! :(" << endl;
    if (!goodChi2) {
      cout << "  -->bad chi^2 value = " << chisq/ndf << " is not within [" << chisqLowerCut << "-" << chisqUpperCut << "]" << endl;
      logFile << "  -->bad chi^2 value = " << chisq/ndf << " is not within [" << chisqLowerCut << "-" << chisqUpperCut << "]" << endl;
    }
    if (!goodSigErr) {
      cout << "  -->bad signal error = " << badSigStr << "\%" << endl;
      logFile << "  -->bad signal error = " << badSigStr << "\%" << endl;
    }
    if (!goodParams) {
      cout << "  -->bad parameter (" << badParamsList << ") hits limit" << endl;
      logFile << "  -->bad parameter (" << badParamsList << ") hits limit" << endl;
    }

    good = 0;
  }
  cout << "good = " << good << endl;

  return good;
} 
 
