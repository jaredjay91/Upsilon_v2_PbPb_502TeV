#include "GetYieldsvsPhi.C"

//const int numphibins = 4;
//float phibins[5] = {0.0, 0.125, 0.25, 0.375, 0.5};

const int numybins = 3;
float ybins[7] = {0.0, 0.8, 1.6, 2.4};

const int numptbins = 4;
float ptbins[5] = {0,3,6,10,50};

const int numcbins = 3;
float cbins[4] = {10,30,50,90};

void GetAllYieldsvsPhi(int whichUpsilon=1) {

  int collId = kAADATA;
  int whichModel = 0;
  float muPtCut=3.5;
  float ptLow, ptHigh, yLow, yHigh;
  int cLow, cHigh;

  bool INTBIN = kTRUE;
  bool PTBINS = kTRUE;
  bool YBINS = kTRUE;
  bool CBINS = kTRUE;

  //integrated bin
  if (INTBIN) GetYieldsvsPhi(collId, 0, 50, 0, 2.4, 10, 90, muPtCut, whichUpsilon);

  //Pt bins
  if (PTBINS) {
    yLow = 0.0;
    yHigh = 2.4;
    cLow = 10;
    cHigh = 90;
    for (int ipt=0; ipt<numptbins; ipt++) {
      ptLow = ptbins[ipt];
      ptHigh = ptbins[ipt+1];
      cout << "[" << ptLow << "," << ptHigh << "]" << endl;
      GetYieldsvsPhi(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, whichUpsilon);
    }
  }

  //Rapidity bins
  if (YBINS) {
    ptLow = 0;
    ptHigh = 50;
    cLow = 10;
    cHigh = 90;
    for (int iy=0; iy<numybins; iy++) {
      yLow = ybins[iy];
      yHigh = ybins[iy+1];
      cout << "[" << yLow << "," << yHigh << "]" << endl;
      GetYieldsvsPhi(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, whichUpsilon);
    }
  }
 
  //Centrality bins
  if (CBINS) {
    ptLow = 0;
    ptHigh = 50;
    yLow = 0.0;
    yHigh = 2.4;
    for (int ic=0; ic<numcbins; ic++) {
      cLow = cbins[ic];
      cHigh = cbins[ic+1];
      cout << "[" << cLow << "," << cHigh << "]" << endl;
      GetYieldsvsPhi(collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, whichUpsilon);
    }
  }
}



