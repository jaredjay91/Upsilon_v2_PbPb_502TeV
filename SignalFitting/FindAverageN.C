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


using namespace std;
using namespace RooFit;
void FindAverageN(int collId=kPADATA) {

  int cLow = 0; int cHigh = 100;
  float dphiEp2Low = 0; float dphiEp2High = 0.5;
  float muPtCut=4.0;

  /*float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[7] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[4] = {0.0,0.8,1.6,2.4};
  float ptbins3[3] = {0,6,30};
  float ybins3[3] = {0.0,1.2,2.4};*/

  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ptbins3[3] = {0,6,30};
  float ybins3[3] = {-1.93,0.0,1.93};

  int numptbins, numybins;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow = 0.0;
  float yHigh = 2.4;

  //start at y=0 if using PP data


  //Define pointers to sets of bins
  float *ptbinsptr;
  float *ybinsptr;
  int ystart;

  double nsumpt = 0;
  double nsumy = 0;
  double xsumpt = 0;
  double xsumy = 0;
  double alphasumpt = 0;
  double alphasumy = 0;
  int totalptbins = 0;
  int totalybins = 0;
  double nsumptweighted = 0;
  double nsumyweighted = 0;
  double totalptbinsize = 30;
  double totalybinsize = 2.4;

  double nsumsqry = 0;
  double xsumsqry = 0;
  double alphasumsqry = 0;

  //TString directory = "AllParamFree/";
  //TString directory = "Fits_with_n_fixed/";
  //TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithSimICs_navgnew/";
  //TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/OfficialNominalFitsFeb2018/";
  //TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithSimICs_Feb202018/";
  TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithSimICs_AllParamFree_March15/";

  //loop through all the bins and do the fit in each one.
  for (int whichUpsilon =1; whichUpsilon<=3; whichUpsilon++) {

    if (whichUpsilon==1) {
      numptbins = 6;
      numybins = sizeof(ybins1)/sizeof(float)-1;
      ptbinsptr = &ptbins1[0];
      ybinsptr = &ybins1[0];
    }
    else if (whichUpsilon==2) {
      numptbins = 3;
      numybins = sizeof(ybins2)/sizeof(float)-1;
      ptbinsptr = &ptbins2[0];
      ybinsptr = &ybins2[0];
    }
    else if (whichUpsilon==3) {
      numptbins = 2;
      numybins = sizeof(ybins3)/sizeof(float)-1;
      ptbinsptr = &ptbins3[0];
      ybinsptr = &ybins3[0];
    }

    ystart = 0;

    /*for (int ipt = 0; ipt<numptbins; ipt++) {
      ptLow = *(ptbinsptr+ipt);
      ptHigh = *(ptbinsptr+ipt+1);
      yLow = *(ybinsptr+ystart);//-1.93 for pPb, 0.0 for pp
      yHigh = 1.93;

      //import fit result
      TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
      //cout << NomFileName << endl;
      cout << Form("pt%.1f-%.1f: ",ptLow,ptHigh);
      TFile* NomFile = TFile::Open(NomFileName,"READ");
      RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
      NomFile->Close("R");

      double nval = Nomws->var("n1s_1")->getVal();
      cout << "n1s_1 = " << nval << endl;
      double xval = Nomws->var("x1s")->getVal();
      cout << "x1s = " << xval << endl;
      double alphaval = Nomws->var("alpha1s_1")->getVal();
      cout << "alpha1s_1 = " << alphaval << endl;
      nsumpt += nval;
      xsumpt += xval;
      alphasumpt += alphaval;
      totalptbins += 1;
      nsumptweighted += nval*(ptHigh-ptLow);

      delete NomFile;
      delete Nomws;
    }*/

    for (int iy = ystart; iy<numybins; iy++) {
      ptLow = 0;
      ptHigh = 30;
      yLow = *(ybinsptr+iy);
      yHigh = *(ybinsptr+iy+1);

      //import fit result
      TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
      //cout << NomFileName << endl;
      cout << Form("y%.2f-%.2f: ",yLow,yHigh);
      TFile* NomFile = TFile::Open(NomFileName,"READ");
      RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
      NomFile->Close("R");

      double nval = Nomws->var("n1s_1")->getVal();
      cout << "n1s_1 = " << nval << endl;
      double xval = Nomws->var("x1s")->getVal();
      cout << "x1s = " << xval << endl;
      double alphaval = Nomws->var("alpha1s_1")->getVal();
      cout << "alpha1s_1 = " << alphaval << endl;
      nsumy += nval;
      xsumy += xval;
      alphasumy += alphaval;
      nsumsqry += nval*nval;
      xsumsqry += xval*xval;
      alphasumsqry += alphaval*alphaval;
      totalybins += 1;
      nsumyweighted += nval*(yHigh-yLow);

      delete NomFile;
      delete Nomws;
    }
  }

  //Calculate averages
  double navgpt = nsumpt/totalptbins;
  double navgy = nsumy/totalybins;
  double nmsy = nsumsqry/totalybins;
  double nsd = sqrt(nmsy-navgy*navgy);
  //cout << "navgpt = " << nsumpt << "/" << totalptbins << " = " << navgpt << endl;
  cout << "navgy = " << nsumy << "/" << totalybins << " = " << navgy << " +/- " << nsd << endl;
  nsumptweighted = nsumptweighted*totalptbins/totalptbinsize;
  nsumyweighted = nsumyweighted*totalybins/totalybinsize;
  double navgptweighted = nsumptweighted/totalptbins;
  double navgyweighted = nsumyweighted/totalybins;
  //cout << "navgptweighted = " << nsumptweighted << "/" << totalptbins << " = " << navgptweighted << endl;
  //cout << "navgyweighted = " << nsumyweighted << "/" << totalybins << " = " << navgyweighted << endl;

  double xavgpt = xsumpt/totalptbins;
  double xavgy = xsumy/totalybins;
  double xmsy = xsumsqry/totalybins;
  double xsd = sqrt(xmsy-xavgy*xavgy);
  //cout << "xavgpt = " << xsumpt << "/" << totalptbins << " = " << xavgpt << endl;
  cout << "xavgy = " << xsumy << "/" << totalybins << " = " << xavgy << " +/- " << xsd << endl;

  double alphaavgpt = alphasumpt/totalptbins;
  double alphaavgy = alphasumy/totalybins;
  double alphamsy = alphasumsqry/totalybins;
  double alphasd = sqrt(alphamsy-alphaavgy*alphaavgy);
  //cout << "alphaavgpt = " << alphasumpt << "/" << totalptbins << " = " << alphaavgpt << endl;
  cout << "alphaavgy = " << alphasumy << "/" << totalybins << " = " << alphaavgy << " +/- " << alphasd << endl;
}
