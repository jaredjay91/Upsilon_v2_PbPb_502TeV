#include "../../HeaderFiles/cutsAndBin.h"

void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=18, float textFont=43){
   TLatex *tex = new TLatex(xp,yp,text);
   tex->SetTextFont(textFont);
   //   if(bold)tex->SetTextFont(43);
   tex->SetTextSize(textSize);
   tex->SetTextColor(textColor);
   tex->SetLineWidth(1);
   tex->SetNDC();
   tex->Draw();
}

double getNominalV2(int whichBin=0, float ptLow=0, float ptHigh=50, float yLow=0.0, float yHigh=2.4, int cLow=10, int cHigh=90){

  TString histlabel;
  if (ptHigh-ptLow<49) histlabel = "pt";
  else if (yHigh-yLow<2.3) histlabel = "y";
  else if (cHigh-cLow<79) histlabel = "c";

  TFile* inFile = new TFile("../../v2Fitting/Ups_1_v2_cent10-90_nom.root","READ");
  TH1D* hv2 = (TH1D*)inFile->Get(Form("hv2%s",histlabel.Data()));

  //Apply event plane resolution correction
  bool EpResCorrection = kTRUE;
  TString resCorFileName = "../../Skimming/condor/averages/resCorFile_n56114317_BinByBin.root";
  TFile* EPRFile = new TFile(resCorFileName,"READ");
  TH1D* hR = (TH1D*)EPRFile->Get(Form("hR%s",histlabel.Data()));
  hv2->Sumw2(); hR->Sumw2();
  if (EpResCorrection) {
    hv2->Divide(hR);
  }

  return hv2->GetBinContent(whichBin);
}

void getSystOneBin(TH1D* hSyst, int whichBin=0, float ptLow=0, float ptHigh=50, float yLow=0.0, float yHigh=2.4, int cLow=10, int cHigh=90, TString systStr="") {

  if (whichBin<1) {
    cout << "ERROR: invalid bin number" << endl;
    return;
  }

  int collId = kAADATA;
  float muPtCut = 3.5;
  float dphiEp2Low = 0.0;
  float dphiEp2High = 0.5;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString pseudoExpFileName = Form("PseudoExpResults/%sPseudoExpResults_%s.root", systStr.Data(), kineLabel.Data());
  TFile* pseudoExpFile = TFile::Open(pseudoExpFileName, "READ");
  TNtuple* ntuple = (TNtuple*)pseudoExpFile->Get("ntuple");

  TCanvas* cSyst = new TCanvas("cSyst","cSyst",1200,400);
  cSyst->Divide(3,1);

  int numv2bins = 50;
  double v2max = 0.1;
  TH1D* hNom = new TH1D("hNom","nominal v2",numv2bins,-v2max,v2max);
  TH1D* hAlt = new TH1D("hAlt","alternate v2",numv2bins,-v2max,v2max);
  TString perc = "%";
  TH1D* hDiff = new TH1D("hDiff",Form("%s difference",perc.Data()),numv2bins,-v2max,v2max);
  TH2D* hchisq = new TH2D("hchisq","chisq check",numv2bins,0,6000,100,0,20);

  cSyst->cd(1);
  ntuple->Draw("v2Nom>>hNom");

  cSyst->cd(2);
  ntuple->Draw("v2Alt>>hAlt");

  cSyst->cd(3);
  hDiff->Draw();

  TLeaf *v2NomLeaf = ntuple->GetLeaf("v2Nom");
  TLeaf *v2AltLeaf = ntuple->GetLeaf("v2Alt");
  TLeaf *chisqNomLeaf = ntuple->GetLeaf("chisqNom");
  TLeaf *chisqAltLeaf = ntuple->GetLeaf("chisqAlt");
  int numEntries = ntuple->GetEntries();
  double avgDiff = 0;
  double avgAbsDiff = 0;
  double sumSqrs = 0;
  int numBelow = 0;
  double numList[numEntries];
  for (int j=0; j<numEntries; j++) {
    ntuple->GetEntry(j);
    double v2Nom = v2NomLeaf->GetValue();
    double v2Alt = v2AltLeaf->GetValue();
    //double percDiff = fabs((v2Alt-v2Nom)/v2Nom)*100;
    //double percDiff = ((v2Alt-v2Nom)/v2Nom)*100;
    double percDiff = v2Alt-v2Nom;
    if (fabs(percDiff)>0.1) {
      cout << "Omitting this outlier:" << endl;
      cout << "v2Nom=" << v2Nom << endl;
      cout << "v2Alt=" << v2Alt << endl;
      cout << "percDiff=" << percDiff << endl;
      continue;
    }
    numList[j] = fabs(percDiff);
    avgDiff = avgDiff + percDiff;
    avgAbsDiff = avgAbsDiff + fabs(percDiff);
    sumSqrs = sumSqrs + pow(percDiff,2);
    hDiff->Fill(percDiff);
    double chisqNom = chisqNomLeaf->GetValue();
    double chisqAlt = chisqAltLeaf->GetValue();
    hchisq->Fill(percDiff, chisqNom+chisqAlt);
  }

  sort(numList, numList+numEntries);
  double median;
  if (numEntries%2==0) median = (numList[numEntries/2]+numList[numEntries/2-1])/2;
  else median = numList[(numEntries+1)/2];
  double onesigma = numList[67];
  avgDiff = avgDiff/numEntries;
  avgAbsDiff = avgAbsDiff/numEntries;
  double stdDev = sqrt((sumSqrs-pow(avgDiff,2))/numEntries);
  cout << endl << avgDiff << " = " << hDiff->GetMean() << endl << endl;
  //double systVal = avgAbsDiff/100;
  //double systVal = fabs(avgDiff)/100;
  //double systVal = median/100;
  //double systVal = onesigma/100;
  double systVal = stdDev;
  if (fabs(avgDiff)>systVal) systVal = fabs(avgDiff);

  //Convert to relative uncertainty:
  systVal = systVal/fabs(getNominalV2(whichBin, ptLow, ptHigh, yLow, yHigh, cLow, cHigh));
  cout << "systVal = " << systVal << endl;
  hSyst->SetBinContent(whichBin, systVal);
  cout << "hSyst->GetBinContent(whichBin) = " << hSyst->GetBinContent(whichBin) << endl;

  cSyst->cd(3);
  float pos_text_x = 0.15;
  float pos_text_y = 0.85;
  float pos_y_diff = 0.06;
  float text_size = 16;
  int text_color = 1;
  drawText(Form("mean = %.4f", fabs(avgDiff)), pos_text_x,pos_text_y,2,text_size);
  drawText(Form("st. dev. = %.4f", stdDev), pos_text_x,pos_text_y-pos_y_diff,2,text_size);
  drawText(Form("mean(abs) = %.4f", fabs(avgAbsDiff)), pos_text_x,pos_text_y-pos_y_diff*2,2,text_size);
  drawText(Form("median(abs) = %.4f", median), pos_text_x,pos_text_y-pos_y_diff*3,2,text_size);
  drawText(Form("68th(abs) = %.4f", onesigma), pos_text_x,pos_text_y-pos_y_diff*4,2,text_size);

  cSyst->Update();
  cSyst->SaveAs(Form("Plots/%sPseudoExpResults_%s.png", systStr.Data(), kineLabel.Data()));
  cSyst->SaveAs(Form("Plots/%sPseudoExpResults_%s.pdf", systStr.Data(), kineLabel.Data()));

  TCanvas* cchisq = new TCanvas("cchisq","cchisq",400,400);
  cchisq->cd();
  hchisq->SetMarkerStyle(8);
  hchisq->Draw();
  cchisq->SaveAs(Form("Plots/%sChisqCheck_%s.png", systStr.Data(), kineLabel.Data()));

  delete hNom;
  delete hAlt;
  delete hDiff;
  delete cSyst;
  delete ntuple;
  delete pseudoExpFile;
}


void getSystematicsFromPseudoExps(int whichUpsilon=1, int whichSyst=1) {

  if (! (whichSyst==1 || whichSyst==2 || whichSyst==3 || whichSyst==4 || whichSyst==5) ) {
    cout << "Error: Invalid value of 'whichSyst'" << endl;
    return;
  }

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  gStyle->SetOptStat(0);

  const int numybins = 3;
  float ybins[4] = {0.0, 0.8, 1.6, 2.4};

  const int numptbins = 4;
  float ptbins[5] = {0,3,6,10,50};

  const int numcbins = 3;
  float cbins[4] = {10,30,50,90};

  TH1D* hSystpt = new TH1D("hv2pt","hv2pt",numptbins,ptbins);
  TH1D* hSysty = new TH1D("hv2y","hv2y",numybins,ybins);
  TH1D* hSystc = new TH1D("hv2c","hv2c",numcbins,cbins);

  //Pt bins
  float ptLow = 0.0; float ptHigh = 50;
  float yLow = 0.0; float yHigh = 2.4;
  int cLow = 10; int cHigh = 90;
  for (int i=0; i<numptbins; i++) {
    ptLow = ptbins[i];
    ptHigh = ptbins[i+1];
    getSystOneBin(hSystpt, i+1, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, systStr);
  }

  //rapidity bins
  ptLow = 0.0; ptHigh = 50;
  yLow = 0.0; yHigh = 2.4;
  cLow = 10; cHigh = 90;
  for (int i=0; i<numybins; i++) {
    yLow = ybins[i];
    yHigh = ybins[i+1];
    getSystOneBin(hSysty, i+1, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, systStr);
  }

  //centrality bins
  ptLow = 0.0; ptHigh = 50;
  yLow = 0.0; yHigh = 2.4;
  cLow = 10; cHigh = 90;
  for (int i=0; i<numcbins; i++) {
    cLow = cbins[i];
    cHigh = cbins[i+1];
    getSystOneBin(hSystc, i+1, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, systStr);
  }

  cout << "Drawing final histograms..." << endl;
  TCanvas* cSyst = new TCanvas("cSyst","cSyst",1200,400);
  cSyst->Divide(3,1);
  cSyst->cd(1);
  hSystpt->Draw();
  cSyst->cd(2);
  hSysty->Draw();
  cSyst->cd(3);
  hSystc->Draw();

  cSyst->Update();

  cout << "Saving plots..." << endl;
  cSyst->SaveAs(Form("Plots/Ups_%i_v2_cent10-90_%sSyst.png", whichUpsilon, systStr.Data()));
  cSyst->SaveAs(Form("Plots/Ups_%i_v2_cent10-90_%sSyst.pdf", whichUpsilon, systStr.Data()));

  cout << "Creating output file..." << endl;
  TString outFileName = Form("Ups_%i_v2_cent10-90_%sSyst.root", whichUpsilon, systStr.Data());
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hSystpt->Write();
  hSysty->Write();
  hSystc->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;

}
