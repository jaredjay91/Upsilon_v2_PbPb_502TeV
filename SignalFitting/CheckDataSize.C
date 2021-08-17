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
void CheckDataSize( 
       int collId = kAADATA,
       float ptLow=0, float ptHigh=3,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=90,
       float muPtCut=3.5
			) 
{

  TString directory = "Simultaneous_dPhiFits/";

  float eta_low = -2.4;
  float eta_high = 2.4;

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;


  //import the simultaneous model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
  TString SimFileName = Form("%ssim_dphinomfitresults_upsilon_%s.root", directory.Data(), kineLabel.Data());
  cout << SimFileName << endl;
  if (gSystem->AccessPathName(SimFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* SimFile = TFile::Open(SimFileName,"READ");
  RooWorkspace *Simws = (RooWorkspace*)SimFile->Get("workspace");

  float simDataSize = 0.0;
  float simModelSize = 0.0;
  float simSignal = 0.0;
  for (int i=0; i<4; i++) {
    float thisSize = Simws->data(Form("reducedDS[%i]",i))->sumEntries();
    float thisSignal = Simws->var(Form("nSig1s[%i]",i))->getVal();
    float thisModelSize = Simws->var(Form("nSig1s[%i]",i))->getVal() + Simws->var(Form("nSig2s[%i]",i))->getVal() + Simws->var(Form("nSig3s[%i]",i))->getVal() + Simws->var(Form("nBkg[%i]",i))->getVal();
    cout << Form("reducedDS[%i] has ",i) << thisSize << " events and " << thisSignal << "Y(1S)" << endl;
    simDataSize = simDataSize + thisSize;
    simSignal = simSignal + thisSignal;
    simModelSize = simModelSize + thisModelSize;
  }
  cout << "simDataSize = " << simDataSize << " and simModelSize = " << simModelSize << " and simSignal = " << simSignal << endl;


  //import the regular model
  //cout << "Importing workspace" << endl;
  TString NomFileName = Form("RoundFits_R4a/nomfitresults_upsilon_%s.root",kineLabel.Data());
  //cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");

  float nomDataSize = ws->data("reducedDS")->sumEntries();
  float nomSignal = ws->var("nSig1s")->getVal();
  float nomModelSize = ws->var("nSig1s")->getVal() + ws->var("nSig2s")->getVal() + ws->var("nSig3s")->getVal() + ws->var("nBkg")->getVal();
  cout << "nomDataSize = " << nomDataSize << " and nomModelSize = " << nomModelSize << " and nomSignal = " << nomSignal << endl;

  cout << "nomSignal/simSignal = " << nomSignal/simSignal << endl;
  cout << "nomModelSize/simModelSize = " << nomModelSize/simModelSize << endl;
}
 
