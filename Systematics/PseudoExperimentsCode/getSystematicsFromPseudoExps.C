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

double getNominalV2(int whichUpsilon=1, int whichBin=0, float ptLow=0, float ptHigh=50, float yLow=0.0, float yHigh=2.4, int cLow=10, int cHigh=90){

  TString histlabel;
  if (ptHigh-ptLow<49) histlabel = "pt";
  else if (yHigh-yLow<2.3) histlabel = "y";
  else if (cHigh-cLow<79) histlabel = "c";

  TFile* inFile = new TFile(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/v2Fitting/Ups_%i_v2_cent10-90_nom.root", whichUpsilon),"READ");
  TH1D* hv2 = (TH1D*)inFile->Get(Form("hv2%s",histlabel.Data()));

  //Apply event plane resolution correction
  bool EpResCorrection = kTRUE;
  TString resCorFileName = "/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Skimming/condor/averages/resCorFile_n56114317_BinByBin.root";
  TFile* EPRFile = new TFile(resCorFileName,"READ");
  TH1D* hR = (TH1D*)EPRFile->Get(Form("hR%s",histlabel.Data()));
  hv2->Sumw2(); hR->Sumw2();
  if (EpResCorrection) {
    hv2->Divide(hR);
  }

  return hv2->GetBinContent(whichBin);
}

void getSystOneBin(TH1D* hSyst, int whichBin=0, float ptLow=0, float ptHigh=50, float yLow=0.0, float yHigh=2.4, int cLow=10, int cHigh=90, TString systStr="", int whichUpsilon=1, TString nomFunc="", TString altFunc="") {

  if (whichBin<1) {
    cout << "ERROR: invalid bin number" << endl;
    return;
  }

  if (systStr=="altBkg") {
    if (ptLow >= 5) {
      nomFunc = "Exp";
      altFunc = "Power Law";
    }
  }

  TString whichUpsString = "1_";
  if (whichUpsilon==2) whichUpsString = "2_";
  int collId = kAADATA;
  float muPtCut = 3.5;
  float dphiEp2Low = 0.0;
  float dphiEp2High = 0.5;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString pseudoExpFileName = Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/PseudoExpResults/%sPseudoExpResults_%s%s.root", systStr.Data(), whichUpsString.Data(), kineLabel.Data());
  TFile* pseudoExpFile = TFile::Open(pseudoExpFileName, "READ");
  TNtuple* ntuple = (TNtuple*)pseudoExpFile->Get("ntuple");

  TCanvas* cSyst1 = new TCanvas("cSyst1","cSyst1",400,400);
  //TCanvas* cSyst2 = new TCanvas("cSyst2","cSyst2",400,400);
  TCanvas* cSyst3 = new TCanvas("cSyst3","cSyst3",400,400);

  int numv2bins = 100;
  double v2min = -0.045;
  double v2max = 0.062;
  double percmax = 0.026;
  TH1D* hNom = new TH1D("hNom","nominal v_{2}",numv2bins,v2min,v2max);
  TH1D* hAlt = new TH1D("hAlt","alternate v_{2}",numv2bins,v2min,v2max);
  TString perc = "%";
  TH1D* hDiff = new TH1D("hDiff",Form("%s difference","%"),numv2bins,-percmax,percmax);
  TH2D* hchisq = new TH2D("hchisq","chisq check",numv2bins,0,6000,100,0,20);

  float labelSize = 0.05;
  float legSize = 0.045;
  hNom->SetTitle("");
  hNom->GetXaxis()->SetTitle("v_{2}");
  hNom->GetXaxis()->SetTitleSize(labelSize);
  hNom->GetYaxis()->SetTitle("counts");
  hNom->GetYaxis()->SetTitleSize(labelSize);
  hNom->GetXaxis()->SetLabelSize(labelSize);
  hNom->GetYaxis()->SetLabelSize(labelSize);
  hNom->SetMaximum(119);
  hNom->SetFillStyle(1001);
  hNom->SetFillColor(kBlue);
  hNom->SetLineColor(kBlue);
  hNom->SetFillColorAlpha(kBlue,0.1);
  //hNom->GetXaxis()->ChangeLabel("-1","","-.6","","-.2","0",".2","",".6","","1");
  //hNom->GetXaxis()->SetNdivisions(-5,-5);
  hAlt->SetTitle("");
  hAlt->GetXaxis()->SetTitle("v_{2,alt}");
  hAlt->GetXaxis()->SetTitleSize(labelSize);
  hAlt->GetXaxis()->SetLabelSize(labelSize);
  hAlt->GetYaxis()->SetTitle("counts");
  hAlt->GetYaxis()->SetTitleSize(labelSize);
  hAlt->GetYaxis()->SetLabelSize(labelSize);
  hAlt->SetMaximum(119);
  hAlt->SetFillStyle(1001);
  hAlt->SetFillColor(kRed);
  hAlt->SetLineColor(kRed);
  hAlt->SetFillColorAlpha(kRed,0.1);
  hDiff->SetTitle("");
  hDiff->GetXaxis()->SetTitle("x = v_{2,alt} - v_{2,nom}");
  hDiff->GetXaxis()->SetTitleSize(labelSize);
  hDiff->GetXaxis()->SetLabelSize(labelSize);
  hDiff->GetYaxis()->SetTitle("counts");
  hDiff->GetYaxis()->SetTitleSize(labelSize);
  hDiff->GetYaxis()->SetLabelSize(labelSize);
  hDiff->SetMaximum(119);
  hDiff->SetFillStyle(1001);
  hDiff->SetFillColor(kGreen+3);
  hDiff->SetLineColor(kGreen+3);
  hDiff->SetFillColorAlpha(kGreen+3,0.25);

  cSyst1->cd();
  ntuple->Draw("v2Nom>>hNom");
  cSyst1->SetLeftMargin(0.14);
  cSyst1->SetBottomMargin(0.11);
  //gPad->SetTicks(1,1);

  //cSyst2->cd();
  ntuple->Draw("v2Alt>>hAlt");
  //cSyst2->SetLeftMargin(0.14);
  //cSyst2->SetBottomMargin(0.11);

  hNom->Draw();
  hAlt->Draw("same");

  gPad->RedrawAxis();

  cSyst3->cd();
  hDiff->Draw();
  cSyst3->SetLeftMargin(0.14);
  cSyst3->SetBottomMargin(0.11);

  gPad->RedrawAxis();

  TLeaf *v2NomLeaf = ntuple->GetLeaf("v2Nom");
  TLeaf *v2AltLeaf = ntuple->GetLeaf("v2Alt");
  TLeaf *chisqNomLeaf = ntuple->GetLeaf("chisqNom");
  TLeaf *chisqAltLeaf = ntuple->GetLeaf("chisqAlt");
  int numEntries = ntuple->GetEntries();
  double avgv2 = 0;
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
    avgv2 = avgv2 + v2Nom;
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
  avgv2 = avgv2/numEntries;
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
  double nominalv2 = getNominalV2(whichUpsilon, whichBin, ptLow, ptHigh, yLow, yHigh, cLow, cHigh);
  systVal = systVal/fabs(nominalv2);
  cout << "systVal = " << systVal << endl;
  hSyst->SetBinContent(whichBin, systVal);
  cout << "hSyst->GetBinContent(whichBin) = " << hSyst->GetBinContent(whichBin) << endl;

  cSyst3->cd();
  float pos_text_x = 0.18;
  float pos_text_y = 0.83;
  float pos_y_diff = 0.06;
  float text_size = 16;
  int text_color = 1;
  drawText(Form("#LT x #GT = %.4f", fabs(avgDiff)), pos_text_x,pos_text_y,text_color,text_size);
  drawText(Form("#sigma_{x} = %.4f", stdDev), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  //drawText(Form("#LT |x| #GT = %.4f", fabs(avgAbsDiff)), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  //drawText(Form("med(|x|) = %.4f", median), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  //drawText(Form("68th(|x|) = %.4f", onesigma), pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);
  cSyst1->cd();
  drawText(Form("gen. v_{2} = %.3f", avgv2), pos_text_x,pos_text_y,text_color,text_size);

  cSyst1->cd();
  TLegend* l1 = new TLegend(0.5,0.78,0.89,0.89);
  l1->SetBorderSize(0);
  l1->SetTextSize(legSize);
  l1->AddEntry(hNom,nomFunc.Data(),"l");
  l1->AddEntry(hAlt,altFunc.Data(),"l");
  l1->Draw("same");

  //cSyst2->cd();
  //TLegend* l2 = new TLegend(0.5,0.82,0.89,0.89);
  //l2->SetBorderSize(0);
  //l2->SetTextSize(legSize);
  //l2->AddEntry(hAlt,altFunc.Data(),"l");
  //l2->Draw("same");

  cSyst3->cd();
  TLegend* l3 = new TLegend(0.5,0.82,0.89,0.89);
  l3->SetBorderSize(0);
  l3->SetTextSize(legSize);
  l3->AddEntry(hDiff,"Difference","l");
  l3->Draw("same");

  cSyst1->Update();
  cSyst1->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/%sPseudoExpResults_%s%s_1.png", systStr.Data(), whichUpsString.Data(), kineLabel.Data()));
  cSyst1->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/%sPseudoExpResults_%s%s_1.pdf", systStr.Data(), whichUpsString.Data(), kineLabel.Data()));
  //cSyst2->Update();
  //cSyst2->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/%sPseudoExpResults_%s%s_2.png", systStr.Data(), whichUpsString.Data(), kineLabel.Data()));
  //cSyst2->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/%sPseudoExpResults_%s%s_2.pdf", systStr.Data(), whichUpsString.Data(), kineLabel.Data()));
  cSyst3->Update();
  cSyst3->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/%sPseudoExpResults_%s%s_3.png", systStr.Data(), whichUpsString.Data(), kineLabel.Data()));
  cSyst3->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/%sPseudoExpResults_%s%s_3.pdf", systStr.Data(), whichUpsString.Data(), kineLabel.Data()));

  TCanvas* cchisq = new TCanvas("cchisq","cchisq",400,400);
  cchisq->cd();
  hchisq->SetMarkerStyle(8);
  hchisq->Draw();
  cchisq->SaveAs(Form("Plots/%sChisqCheck_%s.png", systStr.Data(), kineLabel.Data()));

  delete hNom;
  delete hAlt;
  delete hDiff;
  delete cSyst1;
  //delete cSyst2;
  delete cSyst3;
  delete ntuple;
  delete pseudoExpFile;
}


void getSystematicsFromPseudoExps(int whichUpsilon=1, int whichSyst=5) {

  if (! (whichSyst==1 || whichSyst==2 || whichSyst==3 || whichSyst==4 || whichSyst==5) ) {
    cout << "Error: Invalid value of 'whichSyst'" << endl;
    return;
  }

  TString systStr = "";
  TString nomFunc = "";
  TString altFunc = "";
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) {
    systStr = "altSig";
    nomFunc = "DCB";
    altFunc = "CB+Gaus";
  }
  else if (whichSyst==2) {
    systStr = "altBkg";
    nomFunc = "Erf*Exp";
    altFunc = "Chebychev";
  }
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) {
    systStr = "altConst";
    nomFunc = "Nom. constraints";
    altFunc = "Alt. constraints";
  }

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
    getSystOneBin(hSystpt, i+1, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, systStr, whichUpsilon, nomFunc, altFunc);
  }

  //rapidity bins
  ptLow = 0.0; ptHigh = 50;
  yLow = 0.0; yHigh = 2.4;
  cLow = 10; cHigh = 90;
  for (int i=0; i<numybins; i++) {
    yLow = ybins[i];
    yHigh = ybins[i+1];
    getSystOneBin(hSysty, i+1, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, systStr, whichUpsilon, nomFunc, altFunc);
  }

  //centrality bins
  ptLow = 0.0; ptHigh = 50;
  yLow = 0.0; yHigh = 2.4;
  cLow = 10; cHigh = 90;
  for (int i=0; i<numcbins; i++) {
    cLow = cbins[i];
    cHigh = cbins[i+1];
    getSystOneBin(hSystc, i+1, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, systStr, whichUpsilon, nomFunc, altFunc);
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
  cSyst->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/Ups_%i_v2_cent10-90_%sSyst.png", whichUpsilon, systStr.Data()));
  cSyst->SaveAs(Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Plots/Ups_%i_v2_cent10-90_%sSyst.pdf", whichUpsilon, systStr.Data()));

  cout << "Creating output file..." << endl;
  TString outFileName = Form("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Systematics/PseudoExperimentsCode/Ups_%i_v2_cent10-90_%sSyst.root", whichUpsilon, systStr.Data());
  cout << outFileName << endl;
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hSystpt->Write();
  hSysty->Write();
  hSystc->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;

}
