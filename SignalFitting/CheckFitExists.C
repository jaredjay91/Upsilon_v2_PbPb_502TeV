#include <iostream>
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"


using namespace std;
int CheckFitExists( 
       int collId = kAADATA,
       float ptLow=16, float ptHigh=30,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=60,
       float muPtCut=4.0,
       float dphiEp2Low = 0,//multiplied by PI
       float dphiEp2High = 0.5,
       int whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  TString directory = "AllParamFree/";
  if (dphiEp2High-dphiEp2Low < 0.5) {
    directory = "dphiFits/";
  }

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
  cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  else {
    return 1;
  }
} 
 
