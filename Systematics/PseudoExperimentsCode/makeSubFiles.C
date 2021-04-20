#include <iostream>
#include <stdio.h>
//#include "../../HeaderFiles/rootFitHeaders.h"
//#include "../../HeaderFiles/commonUtility.h"
#include "../../HeaderFiles/cutsAndBin.h"

bool fileExists(TString filename) {
  if (gSystem->AccessPathName(filename)) return false;
  else return true;
}

bool fileIsGood(TString filename) {
  TFile* f = new TFile(filename);
  if (f->GetListOfKeys()->Contains("ntuple")) return true;
  else return false;
}

void writeFile(TString filename,
	  float ptLow = 0, float ptHigh = 50,
	  float yLow = 0.0, float yHigh = 2.4,
	  int cLow = 10, int cHigh = 90,
	  int whichSyst = 5,
          int whichUpsilon = 1) {

  remove(filename.Data());

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  int collId = kAADATA;
  float muPtCut = 3.5;
  float dphiEp2Low = 0.0;
  float dphiEp2High = 0.5;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);

  ofstream subFile;
  subFile.open(filename.Data(), fstream::in | fstream::out | fstream::app);
  subFile << "executable            = condorExec.sh" << endl;
  subFile << "output                = output/condor.$(ClusterId).$(ProcId).out" << endl;
  subFile << "error                 = error/condor.$(ClusterId).$(ProcId).err" << endl;
  subFile << "log                   = log/condor.$(ClusterId).log" << endl;
//  subFile << "+JobFlavour           = \"longlunch\" " << endl;
  subFile << "+JobFlavour           = \"workday\" " << endl;
  for (int i=1; i<101; i++) {
    TString outFileName = Form("PseudoExpResults/%sPseudoExpResults_%i_%s_%i.root", systStr.Data(), whichUpsilon, kineLabel.Data(), i);
    if (fileExists(outFileName)) {
      cout << "File " << outFileName << " exists";
      if (fileIsGood(outFileName)) {
        cout << " and it's good. :)" << endl;
        continue;
      }
      else cout << " but it's bad. :(" << endl;
    }
    cout << "File " << outFileName << " does not exist. :(" << endl;
    subFile << "arguments = " << ptLow << " " << ptHigh << " " << yLow << " " << yHigh << " " << cLow << " " << cHigh << " " << whichSyst << " " << i << " " << whichUpsilon << endl;
    subFile << "queue" << endl << endl;
  }
}

void makeSubFiles(int whichSyst=5, int whichUpsilon=2) {

  writeFile("condorJobs_pt_0_3.sub", 0, 3, 0.0, 2.4, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_pt_3_6.sub", 3, 6, 0.0, 2.4, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_pt_6_10.sub", 6, 10, 0.0, 2.4, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_pt_10_50.sub", 10, 50, 0.0, 2.4, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_y_00_08.sub", 0, 50, 0.0, 0.8, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_y_08_16.sub", 0, 50, 0.8, 1.6, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_y_16_24.sub", 0, 50, 1.6, 2.4, 10, 90, whichSyst, whichUpsilon);
  writeFile("condorJobs_c_10_30.sub", 0, 50, 0.0, 2.4, 10, 30, whichSyst, whichUpsilon);
  writeFile("condorJobs_c_30_50.sub", 0, 50, 0.0, 2.4, 30, 50, whichSyst, whichUpsilon);
  writeFile("condorJobs_c_50_90.sub", 0, 50, 0.0, 2.4, 50, 90, whichSyst, whichUpsilon);
}
