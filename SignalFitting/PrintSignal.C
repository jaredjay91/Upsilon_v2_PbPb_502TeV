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


using namespace std;
using namespace RooFit;


bool isAbout(float a, float b) {
  if (abs(a-b)<0.01) return kTRUE;
  else return kFALSE;
}


void PrintSignalOneBin(
       int collId = kAADATA,
       float ptLow=0, float ptHigh=50,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=90,
       float muPtCut=3.5
			) 
{
  //import the model
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
  TString NomFileName = Form("RoundFits_R4a/nomfitresults_upsilon_%s.root",kineLabel.Data());
  //cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");

  float nomDataSize = ws->data("reducedDS")->sumEntries();
  float signal1 = ws->var("nSig1s")->getVal();
  float signal1err = ws->var("nSig1s")->getError();
  float signal2 = ws->var("nSig2s")->getVal();
  float signal2err = ws->var("nSig2s")->getError();
  float signal3 = ws->var("nSig3s")->getVal();
  float signal3err = ws->var("nSig3s")->getError();

  cout << " & " << signal1 << " $\\pm$ " << signal1err;
  cout << " & " << signal2 << " $\\pm$ " << signal2err;
  cout << " & " << signal3 << " $\\pm$ " << signal3err;
  cout << " \\\\" << endl;
}


void PrintSignal() 
{

  int collId = kAADATA;
  float ptLow=0; float ptHigh=50;
  float yLow=0.0; float yHigh=2.4;
  int cLow=10; int cHigh=90;
  float muPtCut=3.5;

  float eta_low = -2.4;
  float eta_high = 2.4;

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  const int numybins = 3;
  float ybins[7] = {0.0, 0.8, 1.6, 2.4};

  const int numptbins = 4;
  float ptbins[5] = {0,3,6,10,50};

  const int numcbins = 3;
  float cbins[4] = {10,30,50,90};

  //PT BINS
  cout << "\\hline" << endl;
  for (int ipt=0; ipt<numptbins; ipt++) {
    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
    cout << "$" << ptLow << "\\le p_T < " << ptHigh << "$~GeV";
    yLow = 0.0; yHigh = 2.4;
    cLow = 10; cHigh = 90;

    PrintSignalOneBin(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut);
  }
  
  //Y BINS
  cout << "\\hline" << endl;
  for (int iy=0; iy<numybins; iy++) {
    yLow = ybins[iy];
    yHigh = ybins[iy+1];
    cout << "$" << yLow << "\\le y < " << yHigh << "$";
    ptLow = 0; ptHigh = 50;
    cLow =10; cHigh = 90;

    PrintSignalOneBin(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut);
  }

  //CENTRALITY BINS
  cout << "\\hline" << endl;
  for (int ic=0; ic<numcbins; ic++) {
    cLow = cbins[ic];
    cHigh = cbins[ic+1];
    cout << "Centrality " << cLow << "-" << cHigh << "\\%";
    ptLow = 0; ptHigh = 50;
    yLow =0.0; yHigh = 2.4;

    PrintSignalOneBin(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut);
  }

  cout << "\\hline" << endl;
  cout << "Centrality 10-90\\%";
  PrintSignalOneBin(collId, 0, 50, 0.0, 2.4, 10, 90, muPtCut);

}
 
