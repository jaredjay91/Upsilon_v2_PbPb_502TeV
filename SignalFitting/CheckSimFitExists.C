#include <iostream>
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "RoundsHeader.h"


using namespace std;
int CheckSimFitExists( 
       int collId = kAADATA,
       float ptLow=0, float ptHigh=3,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=90,
       float muPtCut=3.5,
       int whichSyst=3
// 0: Nominal
// 1: AltSig
// 2: AltBkg
// 3: AltAcc
// 4: AltEff
// 5: AltConst
			) 
{

  TString directory = "Simultaneous_dPhiFits/";
  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
  TString NomFileName = Form("%ssim_dphi%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());
  cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  else {
    cout << "THE FIT EXISTS! :)";
    return 1;
  }
} 
 
