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


void PlotFittedParams(int whichUpsilon=1, int collId=kAADATA) {

  //TString directory = "AllParamFree/";
  TString directory = "RoundFits_R4b/";

  float scale = whichUpsilon*40;
  if (collId==kPPDATA) scale = scale*8;
  float muPtCut=3.5;

  float dphiEp2Low = 0;
  float dphiEp2High = 0.5;

  //arrays of upper and lower limits.
  //{0:sigma1s, 1:x1s, 2:alpha, 3:n1s, 4:f1s, 5:errsigma, 6:errmu, 7:lambda}
  double paramsupper[8] = {0.3, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double plotlimupper[8] = {0.3, 1.2, 5.5, 5.5, 1.2, 15.0, 15.0, 25.0};
  double plotlimlower[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  //choose a set of bins
  float ptbins[5] = {0,3,6,10,50};
  float ybins[4] = {0.0,0.8,1.6,2.4};
  float cbins[4] = {10,30,50,90};

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numcbins = sizeof(cbins)/sizeof(float)-1;

  float cMin = 10; float cMax = 90;
  float ptMin = 0; float ptMax = 50;
  float yMin = 0.0; float yMax = 2.4;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");

  //set up canvases
  TCanvas *calpha = new TCanvas("calpha","calpha",4,45,1200,400);
  calpha->Divide(3,1);
  TCanvas *cf1s = new TCanvas("cf1s","cf1s",4,45,1200,400);
  cf1s->Divide(3,1);
  TCanvas *cmass = new TCanvas("cmass","cmass",4,45,1200,400);
  cmass->Divide(3,1);
  TCanvas *cn1s = new TCanvas("cn1s","cn1s",4,45,1200,400);
  cn1s->Divide(3,1);
  TCanvas *csigma1s = new TCanvas("csigma1s","csigma1s",4,45,1200,400);
  csigma1s->Divide(3,1);
  TCanvas *cx1s = new TCanvas("cx1s","cx1s",4,45,1200,400);
  cx1s->Divide(3,1);
  TCanvas *cerrmu = new TCanvas("cerrmu","cerrmu",4,45,1200,400);
  cerrmu->Divide(3,1);
  TCanvas *cerrsigma = new TCanvas("cerrsigma","cerrsigma",4,45,1200,400);
  cerrsigma->Divide(3,1);
  TCanvas *clambda = new TCanvas("clambda","clambda",4,45,1200,400);
  clambda->Divide(3,1);
  TCanvas *cnSig1s = new TCanvas("cnSig1s","cnSig1s",4,45,1200,400);
  cnSig1s->Divide(3,1);
  TCanvas *cnSig2s = new TCanvas("cnSig2s","cnSig2s",4,45,1200,400);
  cnSig2s->Divide(3,1);
  TCanvas *cnSig3s = new TCanvas("cnSig3s","cnSig3s",4,45,1200,400);
  cnSig3s->Divide(3,1);
  TCanvas *cnBkg = new TCanvas("cnBkg","cnBkg",4,45,1200,400);
  cnBkg->Divide(3,1);

  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs pt",numptbins,ptbins);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs pt",numptbins,ptbins);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs pt",numptbins,ptbins);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs pt",numptbins,ptbins);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs pt",numptbins,ptbins);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs pt",numptbins,ptbins);
  TH1F* herrmupt = new TH1F("herrmupt","errmu vs pt",numptbins,ptbins);
  TH1F* herrsigmapt = new TH1F("herrsigmapt","errsigma vs pt",numptbins,ptbins);
  TH1F* hlambdapt = new TH1F("hlambdapt","errlambda vs pt",numptbins,ptbins);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs pt",numptbins,ptbins);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs pt",numptbins,ptbins);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs pt",numptbins,ptbins);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs pt",numptbins,ptbins);

  TH1F* halphay = new TH1F("halphay","alpha vs y",numybins,ybins);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs y",numybins,ybins);
  TH1F* hmassy = new TH1F("hmassy","mass vs y",numybins,ybins);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs y",numybins,ybins);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs y",numybins,ybins);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs y",numybins,ybins);
  TH1F* herrmuy = new TH1F("herrmuy","errmu vs y",numybins,ybins);
  TH1F* herrsigmay = new TH1F("herrsigmay","errsigma vs y",numybins,ybins);
  TH1F* hlambday = new TH1F("hlambday","errlambda vs y",numybins,ybins);
  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs y",numybins,ybins);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs y",numybins,ybins);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs y",numybins,ybins);
  TH1F* hnBkgy = new TH1F("hnBkgy","nBkg vs y",numybins,ybins);

  TH1F* halphac = new TH1F("halphac","alpha vs c",numcbins,cbins);
  TH1F* hf1sc = new TH1F("hf1sc","f1s vs c",numcbins,cbins);
  TH1F* hmassc = new TH1F("hmassc","mass vs c",numcbins,cbins);
  TH1F* hn1sc = new TH1F("hn1sc","n1s vs c",numcbins,cbins);
  TH1F* hsigma1sc = new TH1F("hsigma1sc","sigma vs c",numcbins,cbins);
  TH1F* hx1sc = new TH1F("hx1sc","x1s vs c",numcbins,cbins);
  TH1F* herrmuc = new TH1F("herrmuc","errmu vs c",numcbins,cbins);
  TH1F* herrsigmac = new TH1F("herrsigmac","errsigma vs c",numcbins,cbins);
  TH1F* hlambdac = new TH1F("hlambdac","errlambda vs c",numcbins,cbins);
  TH1F* hnSig1sc = new TH1F("hnSig1sc","nSig1s vs c",numcbins,cbins);
  TH1F* hnSig2sc = new TH1F("hnSig2sc","nSig2s vs c",numcbins,cbins);
  TH1F* hnSig3sc = new TH1F("hnSig3sc","nSig3s vs c",numcbins,cbins);
  TH1F* hnBkgc = new TH1F("hnBkgc","nBkg vs c",numcbins,cbins);

  TString kineLabel, NomFileName;
  TFile* NomFile;
  float ptLow, ptHigh, yLow, yHigh, cLow, cHigh;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
    if (collId==kAADATA) {
      yLow = 0.0;
      yHigh = 2.4;
    }
    cLow = 10;
    cHigh = 90;

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float ywidth = yHigh-yLow;
    if (collId==kPPDATA) ywidth = 2*ywidth;
    tempalpha = ws->var("alpha1s_1")->getVal();  
    tempalphaerr = ws->var("alpha1s_1")->getError();
    tempf1s = ws->var("f1s")->getVal();  
    tempf1serr = ws->var("f1s")->getError();
    tempmass = ws->var("m_{#Upsilon(1S)}")->getVal();  
    tempmasserr = ws->var("m_{#Upsilon(1S)}")->getError();
    tempn1s = ws->var("n1s_1")->getVal();  
    tempn1serr = ws->var("n1s_1")->getError();
    tempsigma1s = ws->var("sigma1s_1")->getVal();  
    tempsigma1serr = ws->var("sigma1s_1")->getError();
    tempx1s = ws->var("x1s")->getVal();  
    tempx1serr = ws->var("x1s")->getError();
    if (ptLow<5) {
      temperrmu = ws->var("#mu")->getVal();  
      temperrmuerr = ws->var("#mu")->getError();
      temperrsigma = ws->var("#sigma")->getVal();  
      temperrsigmaerr = ws->var("#sigma")->getError();
    }
    else {
      temperrmu = 0;
      temperrmuerr = 0;
      temperrsigma = 0;  
      temperrsigmaerr = 0;
    }
    templambda = ws->var("#lambda")->getVal();  
    templambdaerr = ws->var("#lambda")->getError();
    tempnSig1s = ws->var("nSig1s")->getVal();//ywidth;  
    tempnSig1serr = ws->var("nSig1s")->getError();
    tempnSig2s = ws->var("nSig2s")->getVal();//ywidth;  
    tempnSig2serr = ws->var("nSig2s")->getError();
    tempnSig3s = ws->var("nSig3s")->getVal();//ywidth;  
    tempnSig3serr = ws->var("nSig3s")->getError();
    tempnBkg = ws->var("nBkg")->getVal();//ywidth;  
    tempnBkgerr = ws->var("nBkg")->getError();

    //fill histograms
    halphapt->SetBinContent(ipt+1, tempalpha);
    halphapt->SetBinError  (ipt+1, tempalphaerr);
    hf1spt->SetBinContent(ipt+1, tempf1s);
    hf1spt->SetBinError  (ipt+1, tempf1serr);
    hmasspt->SetBinContent(ipt+1, tempmass);
    hmasspt->SetBinError  (ipt+1, tempmasserr);
    hn1spt->SetBinContent(ipt+1, tempn1s);
    hn1spt->SetBinError  (ipt+1, tempn1serr);
    hsigma1spt->SetBinContent(ipt+1, tempsigma1s);
    hsigma1spt->SetBinError  (ipt+1, tempsigma1serr);
    hx1spt->SetBinContent(ipt+1, tempx1s);
    hx1spt->SetBinError  (ipt+1, tempx1serr);
    if (ptLow<5) {
      herrmupt->SetBinContent(ipt+1, temperrmu);
      herrmupt->SetBinError  (ipt+1, temperrmuerr);
      herrsigmapt->SetBinContent(ipt+1, temperrsigma);
      herrsigmapt->SetBinError  (ipt+1, temperrsigmaerr);
    }
    hlambdapt->SetBinContent(ipt+1, templambda);
    hlambdapt->SetBinError  (ipt+1, templambdaerr);
    hnSig1spt->SetBinContent(ipt+1, tempnSig1s);
    hnSig1spt->SetBinError  (ipt+1, tempnSig1serr);
    hnSig2spt->SetBinContent(ipt+1, tempnSig2s);
    hnSig2spt->SetBinError  (ipt+1, tempnSig2serr);
    hnSig3spt->SetBinContent(ipt+1, tempnSig3s);
    hnSig3spt->SetBinError  (ipt+1, tempnSig3serr);
    hnBkgpt->SetBinContent(ipt+1, tempnBkg);
    hnBkgpt->SetBinError  (ipt+1, tempnBkgerr);
    delete ws;
    delete NomFile;

  }

  //draw pt plots
  calpha->cd(1);
  halphapt->SetXTitle("pT");
  halphapt->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphapt->Draw();
  //halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(0,paramslower[2],50,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,paramsupper[2],50,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd(1);
  hf1spt->SetXTitle("pT");
  hf1spt->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1spt->Draw();
  //hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(0,paramslower[4],50,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,paramsupper[4],50,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd(1);
  hmasspt->SetXTitle("pT");
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
  //hmasspt->Fit("pol1");
  TLine *massptlinepdg = new TLine(0,9.46,50,9.46);
  massptlinepdg->SetLineColor(3);
  massptlinepdg->Draw();
  TLine *massptlineLow = new TLine(0,9.36,50,9.36);
  massptlineLow->SetLineColor(kRed);
  massptlineLow->Draw();
  TLine *massptlineHigh = new TLine(0,9.56,50,9.56);
  massptlineHigh->SetLineColor(kRed);
  massptlineHigh->Draw();

  cn1s->cd(1);
  hn1spt->SetXTitle("pT");
  hn1spt->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1spt->Draw();
  //hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(0,paramslower[3],50,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,paramsupper[3],50,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd(1);
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1spt->Draw();
  //hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(0,paramslower[0],50,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],50,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd(1);
  hx1spt->SetXTitle("pT");
  hx1spt->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1spt->Draw();
  //hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(0,paramslower[1],50,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,paramsupper[1],50,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

  cerrmu->cd(1);
  herrmupt->SetXTitle("pT");
  herrmupt->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmupt->Draw();
  //TLine *errmuptlineLow = new TLine(0,paramslower[6],50,paramslower[6]);
  //errmuptlineLow->SetLineColor(kRed);
  //errmuptlineLow->Draw();
  //TLine *errmuptlineHigh = new TLine(0,paramsupper[6],50,paramsupper[6]);
  //errmuptlineHigh->SetLineColor(kRed);
  //errmuptlineHigh->Draw();

  cerrsigma->cd(1);
  herrsigmapt->SetXTitle("pT");
  herrsigmapt->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmapt->Draw();
  //TLine *errsigmaptlineLow = new TLine(0,paramslower[5],50,paramslower[5]);
  //errsigmaptlineLow->SetLineColor(kRed);
  //errsigmaptlineLow->Draw();
  //TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],50,paramsupper[5]);
  //errsigmaptlineHigh->SetLineColor(kRed);
  //errsigmaptlineHigh->Draw();

  clambda->cd(1);
  hlambdapt->SetXTitle("pT");
  hlambdapt->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambdapt->Draw();
  //TLine *lambdaptlineLow = new TLine(0,paramslower[7],50,paramslower[7]);
  //lambdaptlineLow->SetLineColor(kRed);
  //lambdaptlineLow->Draw();
  //TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],50,paramsupper[7]);
  //lambdaptlineHigh->SetLineColor(kRed);
  //lambdaptlineHigh->Draw();

  cnSig1s->cd(1);
  //gPad->SetLogy();
  hnSig1spt->SetXTitle("pT");
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->Draw();
  //hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(0,0,50,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd(1);
  //gPad->SetLogy();
  hnSig2spt->SetXTitle("pT");
  hnSig2spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig2spt->Draw();
  //hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(0,0,50,0);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd(1);
  //gPad->SetLogy();
  hnSig3spt->SetXTitle("pT");
  hnSig3spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig3spt->Draw();
  //hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(0,0,50,0);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  cnBkg->cd(1);
  //gPad->SetLogy();
  hnBkgpt->SetXTitle("pT");
  hnBkgpt->GetYaxis()->SetRangeUser(10,40000*scale);
  hnBkgpt->Draw();
  //hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(0,0,50,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    ptLow = 0;
    ptHigh = 50;
    yLow = ybins[iy];
    yHigh = ybins[iy+1];
    if (collId==kPPDATA && yLow<0) {
      yLow = TMath::Abs(ybins[iy+1]);
      yHigh = TMath::Abs(ybins[iy]);
    }
    cLow = 10;
    cHigh = 90;

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *yws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float ywidth = yHigh-yLow;
    if (collId==kPPDATA) ywidth = 2*ywidth;
    tempalpha = yws->var("alpha1s_1")->getVal();  
    tempalphaerr = yws->var("alpha1s_1")->getError();
    tempf1s = yws->var("f1s")->getVal();  
    tempf1serr = yws->var("f1s")->getError();
    tempmass = yws->var("m_{#Upsilon(1S)}")->getVal();  
    tempmasserr = yws->var("m_{#Upsilon(1S)}")->getError();
    tempn1s = yws->var("n1s_1")->getVal();  
    tempn1serr = yws->var("n1s_1")->getError();
    tempsigma1s = yws->var("sigma1s_1")->getVal();  
    tempsigma1serr = yws->var("sigma1s_1")->getError();
    tempx1s = yws->var("x1s")->getVal();  
    tempx1serr = yws->var("x1s")->getError();
    if (ptLow<5) {
      temperrmu = yws->var("#mu")->getVal();  
      temperrmuerr = yws->var("#mu")->getError();
      temperrsigma = yws->var("#sigma")->getVal();  
      temperrsigmaerr = yws->var("#sigma")->getError();
    }
    templambda = yws->var("#lambda")->getVal();  
    templambdaerr = yws->var("#lambda")->getError();
    tempnSig1s = yws->var("nSig1s")->getVal();//ywidth;  
    tempnSig1serr = yws->var("nSig1s")->getError();
    tempnSig2s = yws->var("nSig2s")->getVal();//ywidth;  
    tempnSig2serr = yws->var("nSig2s")->getError();
    tempnSig3s = yws->var("nSig3s")->getVal();//ywidth;  
    tempnSig3serr = yws->var("nSig3s")->getError();
    tempnBkg = yws->var("nBkg")->getVal();//ywidth;  
    tempnBkgerr = yws->var("nBkg")->getError();

    //fill histograms
    halphay->SetBinContent(iy+1, tempalpha);
    halphay->SetBinError  (iy+1, tempalphaerr);
    hf1sy->SetBinContent(iy+1, tempf1s);
    hf1sy->SetBinError  (iy+1, tempf1serr);
    hmassy->SetBinContent(iy+1, tempmass);
    hmassy->SetBinError  (iy+1, tempmasserr);
    hn1sy->SetBinContent(iy+1, tempn1s);
    hn1sy->SetBinError  (iy+1, tempn1serr);
    hsigma1sy->SetBinContent(iy+1, tempsigma1s);
    hsigma1sy->SetBinError  (iy+1, tempsigma1serr);
    hx1sy->SetBinContent(iy+1, tempx1s);
    hx1sy->SetBinError  (iy+1, tempx1serr);
    if (ptLow<5) {
      herrmuy->SetBinContent(iy+1, temperrmu);
      herrmuy->SetBinError  (iy+1, temperrmuerr);
      herrsigmay->SetBinContent(iy+1, temperrsigma);
      herrsigmay->SetBinError  (iy+1, temperrsigmaerr);
    }
    hlambday->SetBinContent(iy+1, templambda);
    hlambday->SetBinError  (iy+1, templambdaerr);
    hnSig1sy->SetBinContent(iy+1, tempnSig1s);
    hnSig1sy->SetBinError  (iy+1, tempnSig1serr);
    hnSig2sy->SetBinContent(iy+1, tempnSig2s);
    hnSig2sy->SetBinError  (iy+1, tempnSig2serr);
    hnSig3sy->SetBinContent(iy+1, tempnSig3s);
    hnSig3sy->SetBinError  (iy+1, tempnSig3serr);
    hnBkgy->SetBinContent(iy+1, tempnBkg);
    hnBkgy->SetBinError  (iy+1, tempnBkgerr);
    delete yws;
    delete NomFile;
  }

  //draw rapidity plots
  calpha->cd(2);
  halphay->SetXTitle("y");
  halphay->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphay->Draw();
  //halphay->Fit("pol1");
  TLine *alphaylineLow = new TLine(0.0,paramslower[2],2.4,paramslower[2]);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(0.0,paramsupper[2],2.4,paramsupper[2]);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();

  cf1s->cd(2);
  hf1sy->SetXTitle("y");
  hf1sy->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1sy->Draw();
  //hf1sy->Fit("pol1");
  TLine *f1sylineLow = new TLine(0.0,paramslower[4],2.4,paramslower[4]);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(0.0,paramsupper[4],2.4,paramsupper[4]);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();

  cmass->cd(2);
  hmassy->SetXTitle("y");
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->Draw();
  //hmassy->Fit("pol1");
  TLine *massylinepdg = new TLine(0.0,9.46,2.4,9.46);
  massylinepdg->SetLineColor(3);
  massylinepdg->Draw();
  TLine *massylineLow = new TLine(0.0,9.36,2.4,9.36);
  massylineLow->SetLineColor(kRed);
  massylineLow->Draw();
  TLine *massylineHigh = new TLine(0.0,9.56,2.4,9.56);
  massylineHigh->SetLineColor(kRed);
  massylineHigh->Draw();

  cn1s->cd(2);
  hn1sy->SetXTitle("y");
  hn1sy->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1sy->Draw();
  //hn1sy->Fit("pol1");
  TLine *n1sylineLow = new TLine(0.0,paramslower[3],2.4,paramslower[3]);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(0.0,paramsupper[3],2.4,paramsupper[3]);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();

  csigma1s->cd(2);
  hsigma1sy->SetXTitle("y");
  hsigma1sy->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1sy->Draw();
  //hsigma1sy->Fit("pol1");
  TLine *sigmaylineLow = new TLine(0.0,paramslower[0],2.4,paramslower[0]);
  sigmaylineLow->SetLineColor(kRed);
  sigmaylineLow->Draw();
  TLine *sigmaylineHigh = new TLine(0.0,paramsupper[0],2.4,paramsupper[0]);
  sigmaylineHigh->SetLineColor(kRed);
  sigmaylineHigh->Draw();

  cx1s->cd(2);
  hx1sy->SetXTitle("y");
  hx1sy->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1sy->Draw();
  //hx1sy->Fit("pol1");
  TLine *x1sylineLow = new TLine(0.0,paramslower[1],2.4,paramslower[1]);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(0.0,paramsupper[1],2.4,paramsupper[1]);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();

  cerrmu->cd(2);
  herrmuy->SetXTitle("y");
  herrmuy->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmuy->Draw();
  //TLine *errmuylineLow = new TLine(0,paramslower[6],50,paramslower[6]);
  //errmuylineLow->SetLineColor(kRed);
  //errmuylineLow->Draw();
  //TLine *errmuylineHigh = new TLine(0,paramsupper[6],50,paramsupper[6]);
  //errmuylineHigh->SetLineColor(kRed);
  //errmuylineHigh->Draw();

  cerrsigma->cd(2);
  herrsigmay->SetXTitle("y");
  herrsigmay->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmay->Draw();
  //TLine *errsigmaylineLow = new TLine(0,paramslower[5],50,paramslower[5]);
  //errsigmaylineLow->SetLineColor(kRed);
  //errsigmaylineLow->Draw();
  //TLine *errsigmaylineHigh = new TLine(0,paramsupper[5],50,paramsupper[5]);
  //errsigmaylineHigh->SetLineColor(kRed);
  //errsigmaylineHigh->Draw();

  clambda->cd(2);
  hlambday->SetXTitle("y");
  hlambday->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambday->Draw();
  //TLine *lambdaylineLow = new TLine(0,paramslower[7],50,paramslower[7]);
  //lambdaylineLow->SetLineColor(kRed);
  //lambdaylineLow->Draw();
  //TLine *lambdaylineHigh = new TLine(0,paramsupper[7],50,paramsupper[7]);
  //lambdaylineHigh->SetLineColor(kRed);
  //lambdaylineHigh->Draw();

  cnSig1s->cd(2);
  //gPad->SetLogy();
  hnSig1sy->SetXTitle("y");
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1sy->Draw();
  //hnSig1sy->Fit("pol1");
  TLine *nSig1sylineLow = new TLine(0.0,0,2.4,0);
  nSig1sylineLow->SetLineColor(kRed);
  nSig1sylineLow->Draw();

  cnSig2s->cd(2);
  //gPad->SetLogy();
  hnSig2sy->SetXTitle("y");
  hnSig2sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig2sy->Draw();
  //hnSig2sy->Fit("pol1");
  TLine *nSig2sylineLow = new TLine(0.0,0,2.4,0);
  nSig2sylineLow->SetLineColor(kRed);
  nSig2sylineLow->Draw();

  cnSig3s->cd(2);
  //gPad->SetLogy();
  hnSig3sy->SetXTitle("y");
  hnSig3sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig3sy->Draw();
  //hnSig3sy->Fit("pol1");
  TLine *nSig3sylineLow = new TLine(0.0,0,2.4,0);
  nSig3sylineLow->SetLineColor(kRed);
  nSig3sylineLow->Draw();

  cnBkg->cd(2);
  //gPad->SetLogy();
  hnBkgy->SetXTitle("y");
  hnBkgy->GetYaxis()->SetRangeUser(10,40000*scale);
  hnBkgy->Draw();
  //hnBkgy->Fit("pol1");
  TLine *nBkgylineLow = new TLine(0.0,0,2.4,0);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();


  //1S centrality loop
  for (int ic = 0; ic<numcbins; ic++) {
    ptLow = 0;
    ptHigh = 50;
    cLow = cbins[ic];
    cHigh = cbins[ic+1];
    yLow = 0.0;
    yHigh = 2.4;

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *cws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float cwidth = cHigh-cLow;
    if (collId==kPPDATA) cwidth = 2*cwidth;
    tempalpha = cws->var("alpha1s_1")->getVal();  
    tempalphaerr = cws->var("alpha1s_1")->getError();
    tempf1s = cws->var("f1s")->getVal();  
    tempf1serr = cws->var("f1s")->getError();
    tempmass = cws->var("m_{#Upsilon(1S)}")->getVal();  
    tempmasserr = cws->var("m_{#Upsilon(1S)}")->getError();
    tempn1s = cws->var("n1s_1")->getVal();  
    tempn1serr = cws->var("n1s_1")->getError();
    tempsigma1s = cws->var("sigma1s_1")->getVal();  
    tempsigma1serr = cws->var("sigma1s_1")->getError();
    tempx1s = cws->var("x1s")->getVal();  
    tempx1serr = cws->var("x1s")->getError();
    if (ptLow<5) {
      temperrmu = cws->var("#mu")->getVal();  
      temperrmuerr = cws->var("#mu")->getError();
      temperrsigma = cws->var("#sigma")->getVal();  
      temperrsigmaerr = cws->var("#sigma")->getError();
    }
    templambda = cws->var("#lambda")->getVal();  
    templambdaerr = cws->var("#lambda")->getError();
    tempnSig1s = cws->var("nSig1s")->getVal();//cwidth;  
    tempnSig1serr = cws->var("nSig1s")->getError();
    tempnSig2s = cws->var("nSig2s")->getVal();//cwidth;  
    tempnSig2serr = cws->var("nSig2s")->getError();
    tempnSig3s = cws->var("nSig3s")->getVal();//cwidth;  
    tempnSig3serr = cws->var("nSig3s")->getError();
    tempnBkg = cws->var("nBkg")->getVal();//cwidth;  
    tempnBkgerr = cws->var("nBkg")->getError();

    //fill histograms
    halphac->SetBinContent(ic+1, tempalpha);
    halphac->SetBinError  (ic+1, tempalphaerr);
    hf1sc->SetBinContent(ic+1, tempf1s);
    hf1sc->SetBinError  (ic+1, tempf1serr);
    hmassc->SetBinContent(ic+1, tempmass);
    hmassc->SetBinError  (ic+1, tempmasserr);
    hn1sc->SetBinContent(ic+1, tempn1s);
    hn1sc->SetBinError  (ic+1, tempn1serr);
    hsigma1sc->SetBinContent(ic+1, tempsigma1s);
    hsigma1sc->SetBinError  (ic+1, tempsigma1serr);
    hx1sc->SetBinContent(ic+1, tempx1s);
    hx1sc->SetBinError  (ic+1, tempx1serr);
    if (ptLow<5) {
      herrmuc->SetBinContent(ic+1, temperrmu);
      herrmuc->SetBinError  (ic+1, temperrmuerr);
      herrsigmac->SetBinContent(ic+1, temperrsigma);
      herrsigmac->SetBinError  (ic+1, temperrsigmaerr);
    }
    hlambdac->SetBinContent(ic+1, templambda);
    hlambdac->SetBinError  (ic+1, templambdaerr);
    hnSig1sc->SetBinContent(ic+1, tempnSig1s);
    hnSig1sc->SetBinError  (ic+1, tempnSig1serr);
    hnSig2sc->SetBinContent(ic+1, tempnSig2s);
    hnSig2sc->SetBinError  (ic+1, tempnSig2serr);
    hnSig3sc->SetBinContent(ic+1, tempnSig3s);
    hnSig3sc->SetBinError  (ic+1, tempnSig3serr);
    hnBkgc->SetBinContent(ic+1, tempnBkg);
    hnBkgc->SetBinError  (ic+1, tempnBkgerr);
    delete cws;
    delete NomFile;
  }

  //draw centrality plots
  calpha->cd(3);
  halphac->SetXTitle("c");
  halphac->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphac->Draw();
  //halphac->Fit("pol1");
  TLine *alphaclineLow = new TLine(cMin,paramslower[2],cMax,paramslower[2]);
  alphaclineLow->SetLineColor(kRed);
  alphaclineLow->Draw();
  TLine *alphaclineHigh = new TLine(cMin,paramsupper[2],cMax,paramsupper[2]);
  alphaclineHigh->SetLineColor(kRed);
  alphaclineHigh->Draw();

  cf1s->cd(3);
  hf1sc->SetXTitle("c");
  hf1sc->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1sc->Draw();
  //hf1sc->Fit("pol1");
  TLine *f1sclineLow = new TLine(cMin,paramslower[4],cMax,paramslower[4]);
  f1sclineLow->SetLineColor(kRed);
  f1sclineLow->Draw();
  TLine *f1sclineHigh = new TLine(cMin,paramsupper[4],cMax,paramsupper[4]);
  f1sclineHigh->SetLineColor(kRed);
  f1sclineHigh->Draw();

  cmass->cd(3);
  hmassc->SetXTitle("c");
  hmassc->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassc->Draw();
  //hmassc->Fit("pol1");
  TLine *massclinepdg = new TLine(cMin,9.46,cMax,9.46);
  massclinepdg->SetLineColor(3);
  massclinepdg->Draw();
  TLine *massclineLow = new TLine(cMin,9.36,cMax,9.36);
  massclineLow->SetLineColor(kRed);
  massclineLow->Draw();
  TLine *massclineHigh = new TLine(cMin,9.56,cMax,9.56);
  massclineHigh->SetLineColor(kRed);
  massclineHigh->Draw();

  cn1s->cd(3);
  hn1sc->SetXTitle("c");
  hn1sc->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1sc->Draw();
  //hn1sc->Fit("pol1");
  TLine *n1sclineLow = new TLine(cMin,paramslower[3],cMax,paramslower[3]);
  n1sclineLow->SetLineColor(kRed);
  n1sclineLow->Draw();
  TLine *n1sclineHigh = new TLine(cMin,paramsupper[3],cMax,paramsupper[3]);
  n1sclineHigh->SetLineColor(kRed);
  n1sclineHigh->Draw();

  csigma1s->cd(3);
  hsigma1sc->SetXTitle("c");
  hsigma1sc->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1sc->Draw();
  //hsigma1sc->Fit("pol1");
  TLine *sigmaclineLow = new TLine(cMin,paramslower[0],cMax,paramslower[0]);
  sigmaclineLow->SetLineColor(kRed);
  sigmaclineLow->Draw();
  TLine *sigmaclineHigh = new TLine(cMin,paramsupper[0],cMax,paramsupper[0]);
  sigmaclineHigh->SetLineColor(kRed);
  sigmaclineHigh->Draw();

  cx1s->cd(3);
  hx1sc->SetXTitle("c");
  hx1sc->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1sc->Draw();
  //hx1sc->Fit("pol1");
  TLine *x1sclineLow = new TLine(cMin,paramslower[1],cMax,paramslower[1]);
  x1sclineLow->SetLineColor(kRed);
  x1sclineLow->Draw();
  TLine *x1sclineHigh = new TLine(cMin,paramsupper[1],cMax,paramsupper[1]);
  x1sclineHigh->SetLineColor(kRed);
  x1sclineHigh->Draw();

  cerrmu->cd(3);
  herrmuc->SetXTitle("c");
  herrmuc->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmuc->Draw();
  //TLine *errmuclineLow = new TLine(0,paramslower[6],50,paramslower[6]);
  //errmuclineLow->SetLineColor(kRed);
  //errmuclineLow->Draw();
  //TLine *errmuclineHigh = new TLine(0,paramsupper[6],50,paramsupper[6]);
  //errmuclineHigh->SetLineColor(kRed);
  //errmuclineHigh->Draw();

  cerrsigma->cd(3);
  herrsigmac->SetXTitle("c");
  herrsigmac->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmac->Draw();
  //TLine *errsigmaclineLow = new TLine(0,paramslower[5],50,paramslower[5]);
  //errsigmaclineLow->SetLineColor(kRed);
  //errsigmaclineLow->Draw();
  //TLine *errsigmaclineHigh = new TLine(0,paramsupper[5],50,paramsupper[5]);
  //errsigmaclineHigh->SetLineColor(kRed);
  //errsigmaclineHigh->Draw();

  clambda->cd(3);
  hlambdac->SetXTitle("c");
  hlambdac->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambdac->Draw();
  //TLine *lambdaclineLow = new TLine(0,paramslower[7],50,paramslower[7]);
  //lambdaclineLow->SetLineColor(kRed);
  //lambdaclineLow->Draw();
  //TLine *lambdaclineHigh = new TLine(0,paramsupper[7],50,paramsupper[7]);
  //lambdaclineHigh->SetLineColor(kRed);
  //lambdaclineHigh->Draw();

  cnSig1s->cd(3);
  //gPad->SetLogc();
  hnSig1sc->SetXTitle("c");
  hnSig1sc->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1sc->Draw();
  //hnSig1sc->Fit("pol1");
  TLine *nSig1sclineLow = new TLine(cMin,0,cMax,0);
  nSig1sclineLow->SetLineColor(kRed);
  nSig1sclineLow->Draw();

  cnSig2s->cd(3);
  //gPad->SetLogc();
  hnSig2sc->SetXTitle("c");
  hnSig2sc->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig2sc->Draw();
  //hnSig2sc->Fit("pol1");
  TLine *nSig2sclineLow = new TLine(cMin,0,cMax,0);
  nSig2sclineLow->SetLineColor(kRed);
  nSig2sclineLow->Draw();

  cnSig3s->cd(3);
  //gPad->SetLogc();
  hnSig3sc->SetXTitle("c");
  hnSig3sc->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig3sc->Draw();
  //hnSig3sc->Fit("pol1");
  TLine *nSig3sclineLow = new TLine(cMin,0,cMax,0);
  nSig3sclineLow->SetLineColor(kRed);
  nSig3sclineLow->Draw();

  cnBkg->cd(3);
  //gPad->SetLogc();
  hnBkgc->SetXTitle("c");
  hnBkgc->GetYaxis()->SetRangeUser(10,40000*scale);
  hnBkgc->Draw();
  //hnBkgc->Fit("pol1");
  TLine *nBkgclineLow = new TLine(cMin,0,cMax,0);
  nBkgclineLow->SetLineColor(kRed);
  nBkgclineLow->Draw();


 TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";

  //save plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins.png",strId.Data(),whichUpsilon));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins.png",strId.Data(),whichUpsilon));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins.png",strId.Data(),whichUpsilon));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins.png",strId.Data(),whichUpsilon));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins.png",strId.Data(),whichUpsilon));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins.png",strId.Data(),whichUpsilon));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins.png",strId.Data(),whichUpsilon));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins.png",strId.Data(),whichUpsilon));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins.png",strId.Data(),whichUpsilon));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins.png",strId.Data(),whichUpsilon));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins.png",strId.Data(),whichUpsilon));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins.png",strId.Data(),whichUpsilon));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins.png",strId.Data(),whichUpsilon));

/*  calpha->Print("h1.pdf(","pdf");
  cf1s->Print("h1.pdf","pdf");
  cmass->Print("h1.pdf","pdf");
  cn1s->Print("h1.pdf","pdf");
  csigma1s->Print("h1.pdf","pdf");
  cx1s->Print("h1.pdf","pdf");
  cerrmu->Print("h1.pdf","pdf");
  cerrsigma->Print("h1.pdf","pdf");
  clambda->Print("h1.pdf","pdf");
  cnSig1s->Print("h1.pdf","pdf");
  cnSig2s->Print("h1.pdf","pdf");
  cnSig3s->Print("h1.pdf","pdf");
  cnBkg->Print("h1.pdf)","pdf");
*/
}
