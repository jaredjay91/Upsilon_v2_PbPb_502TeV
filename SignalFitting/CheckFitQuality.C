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
int CheckFitQuality( 
       int collId = kAADATA,
       float ptLow=16, float ptHigh=30,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=60,
       float muPtCut=4.0,
       float dphiEp2Low = 0,//multiplied by PI
       float dphiEp2High = 0.5,
       int whichSyst=0,
// 0: Nominal
// 1: AltSig
// 2: AltBkg
// 3: AltAcc
// 4: AltEff
       int whichRound=R4a
			) 
{

  TString directory = "AllParamFree/";
  if (dphiEp2High-dphiEp2Low < 0.5) {
    directory = Form("dphiFits_%s/",roundLabel[whichRound].Data());
  }
  else if (whichRound>0) directory = Form("RoundFits_%s/",roundLabel[whichRound].Data());
  TString logFileName = Form("%slog.txt",directory.Data());

  TString Params[5] = {"sigma1s_1","x1s","alpha1s_1","n1s_1","f1s"};

  //Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.25, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
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
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("%s%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());
  //TString NomFileName = Form("FitsWithSimICs_run1fnxalphafixed/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");
  RooAbsData* reducedDS = ws->data("reducedDS");
  RooFitResult* fitRes2 = (RooFitResult*)ws->obj("fitresult_model_reducedDS");

  //Plot it
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  ws->pdf("model")->plotOn(myPlot,Name("modelHist"));

  // PULL
  RooHist* hpull = myPlot->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);

  //calculate chi-squared
  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin; 
  double *ypull = hpull->GetY();
  for(int i=0;i<nBins;i++)
  {
    if(ypull[i] == 0) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  cout << "numFitPar = " << numFitPar << endl;
  int ndf = nFullBinsPull - numFitPar;
  cout << "chisq/dof = " << chisq/ndf << endl;

  float temp1 = ws->var("nSig1s")->getVal();  
  float temp1err = ws->var("nSig1s")->getError();  
  float temp2 = ws->var("nSig2s")->getVal();  
  float temp2err = ws->var("nSig2s")->getError();  
  float temp3 = ws->var("nSig3s")->getVal();  
  float temp3err = ws->var("nSig3s")->getError();  
  
  cout << "1S signal    =  " << temp1 << " +/- " << temp1err << endl;
  cout << "2S signal    =  " << temp2 << " +/- " << temp2err << endl;
  cout << "3S signal    =  " << temp3 << " +/- " << temp3err << endl;
  cout << "Total signal =  " << temp1+temp2+temp3 << endl;

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
  if (whichRound==R4a || whichRound==R4b) {//all the signal parameters should be good in the nominal fits.
    paramstart = 0;
    paramend = 4;
  }
  else if (whichRound==R3a || whichRound==R3b) {//x has to be good in round 3 so we can get a good average for x.
    paramstart = 1;
    paramend = 3;
  }
  else {//Only alpha and n need to be good in the AllParamFree fits, i.e. indices 2 and 3.
    paramstart = 2;
    paramend = 3;
  }
  //for (int i=0;i<5;i++) {
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
  if (temp1err/temp1>errUpperLimit || temp1err/temp1<errLowerLimit) goodSigErr = kFALSE;
  else goodSigErr = kTRUE;

  if (goodChi2 && goodSigErr && goodParams){
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
      cout << "  -->bad signal error = " << temp1err/temp1*100 << "\%" << endl;
      logFile << "  -->bad signal error = " << temp1err/temp1*100 << "\%" << endl;
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
 
