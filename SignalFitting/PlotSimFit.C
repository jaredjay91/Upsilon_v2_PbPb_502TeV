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
void PlotSimFit( 
       int collId = kAADATA,
       float ptLow=10, float ptHigh=50,
       float yLow=0.0, float yHigh=2.4,
       int cLow=10, int cHigh=90,
       float muPtCut=3.5,
       int whichSyst=0
// 0: Nominal
// 1: AltSig
// 2: AltBkg
// 3: AltAcc
// 4: AltEff
// 5: AltConst
			) 
{

  TString directory = "Simultaneous_dPhiFits/";

  TString Params[5] = {"sigma1s_1","x1s","alpha1s_1","n1s_1","f1s"};

  TString paramNames[13] = {"n1s_1","alpha1s_1","x1s","f1s","sigma1s_1","m_{#Upsilon(1S)}","#lambda","#mu","#sigma"};

  //Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.35, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  float eta_low = -2.4;
  float eta_high = 2.4;

  float dphiEp2bins[5] = {0.0, 0.125, 0.25, 0.375, 0.5};
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  //import the model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
  TString SimFileName = Form("%ssim_dphi%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());
  cout << SimFileName << endl;
  if (gSystem->AccessPathName(SimFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* SimFile = TFile::Open(SimFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)SimFile->Get("workspace");
  RooAbsData* reducedDS = ws->data("reducedDS[0]");
  RooFitResult* fitRes2 = (RooFitResult*)ws->obj("fitresult_model_dsABCD");
  fitRes2->Print("v");

  //Plot it and calculate chi-squared
  TCanvas* c[4];
  TPad *pad1[4];
  TPad *pad2[4];
  RooPlot* myPlot[4];
  RooHist* hpull[4];
  RooPlot* pullFrame[4];
  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin;
  double *ypull[4];
  TLine l1[4];

  //Scale the yield parameters to match the true value (no change in relative yields, so no change in v2. This just ensures that the proper yields are printed on the plots.)
  float simSignal = 0.0;
  for (int i=0; i<4; i++) {
    simSignal = simSignal + ws->var(Form("nSig1s[%i]",i))->getVal();
  }
  TString NomFileName = Form("RoundFits_R4a/nomfitresults_upsilon_%s.root",kineLabel.Data());
  //cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  float nomSignal = Nomws->var("nSig1s")->getVal();
  double scaleFactor = nomSignal/simSignal;

  RooAbsPdf* cb1s = ws->pdf("cb1s");
  RooAbsPdf* cb2s = ws->pdf("cb2s");
  RooAbsPdf* cb3s = ws->pdf("cb3s");
  RooAbsPdf* bkg;
  if (ptLow<5) bkg = ws->pdf("bkgLowPt");
  else bkg = ws->pdf("bkgHighPt");

    //setTDRStyle();
    writeExtraText = false;
    extraText = "Preliminary";

  TString perc = "%";

  for (int i=0; i<4; i++) {
    c[i] = new TCanvas(Form("canvas[%i]",i),"My plots",4,45,550,520);
    c[i]->cd();
    pad1[i] = new TPad(Form("pad1[%i]",i), Form("pad1[%i]",i), 0, 0.25, 0.98, 1.0);
    pad1[i]->SetTicks(1,1);
    pad1[i]->Draw(); pad1[i]->cd();
    myPlot[i] = ws->var("mass")->frame(nMassBin); // bins
    ws->data(Form("reducedDS[%i]",i))->plotOn(myPlot[i],Name(Form("dataHist[%i]",i)), Range(massLow, massHigh));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Name(Form("modelHist[%i]",i)), Range(massLow, massHigh));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Name("Sig1S"),Components(RooArgSet(*cb1s)), Range(massLow, massHigh),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Components(RooArgSet(*cb2s)), Range(massLow, massHigh),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Components(RooArgSet(*cb3s)), Range(massLow, massHigh),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot[i],Name("bkgPDF"),Components(RooArgSet(*bkg)), Range(massLow, massHigh),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
    myPlot[i]->SetFillStyle(4000);
    myPlot[i]->SetTitle("");
    myPlot[i]->SetAxisRange(massLowForPlot, massHighForPlot,"X");
    myPlot[i]->GetYaxis()->SetTitleOffset(1.43);
    myPlot[i]->GetYaxis()->CenterTitle();
    myPlot[i]->GetYaxis()->SetTitleSize(0.058);
    myPlot[i]->GetYaxis()->SetLabelSize(0.054);
    myPlot[i]->GetXaxis()->SetLabelSize(0);
    myPlot[i]->GetXaxis()->SetRangeUser(8,14);
    myPlot[i]->GetXaxis()->SetTitleSize(0);
    myPlot[i]->Draw();

    //float datasetSize = ws->data(Form("reducedDS[%i]",i))->sumEntries();
    //float modelSize = ws->var(Form("nSig1s[%i]",i))->getVal() + ws->var(Form("nSig2s[%i]",i))->getVal() + ws->var(Form("nSig3s[%i]",i))->getVal() + ws->var(Form("nBkg[%i]",i))->getVal();
    //cout << "dataset entries = " << datasetSize << endl;
    //cout << "model entries = " << ws->pdf(Form("model[%i]",i))->normRange() << endl;
    //cout << "Normalization = " << modelSize << endl;
    //scaleFactor[i] = datasetSize/modelSize;
    //cout << "Scale factor = " << scaleFactor[i] << endl;

    float dphiEp2Low = dphiEp2bins[i];
    float dphiEp2High = dphiEp2bins[i+1];

    float pos_text_x = 0.43;
    float pos_text_y = 0.816;
    float pos_y_diff = 0.05;
    float text_size = 12;
    int text_color = 1;
    if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
    else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
    if (collId==kPPDATA || collId==kAADATA) {
      if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
      else drawText(Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
      }
    else if (collId==kPADATA) {
      if(yLow==-yHigh) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
      else drawText(Form("%.2f < y_{CM}^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
      }
    drawText(Form("p_{T}^{#mu} > %.1f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
    drawText(Form("|#eta^{#mu}| < %.1f",eta_high), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
    if (collId==kAADATA) {
      drawText(Form("Centrality %d-%d%s",cLow,cHigh,perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);
      if (dphiEp2Low==0) drawText(Form("|#Delta#phi| < %.2f#pi", dphiEp2High), pos_text_x,pos_text_y-pos_y_diff*5,text_color,text_size);
      else drawText(Form("%.2f#pi < |#Delta#phi| < %.2f#pi", dphiEp2Low, dphiEp2High), pos_text_x,pos_text_y-pos_y_diff*5,text_color,text_size);
    }

    //Print fitted parameter values
    float pos_text_x_params = 0.7;
    float pos_text_y_params = 0.82;
    float pos_y_diff_params = 0.04;
    float text_size_params = 12;
    int text_color_params = 1;
    int numParams = 9;
    TString paramText;
    for (int iparam = 0; iparam<numParams; iparam++) {
      if (ptLow>5) {
        if (iparam>6) break;
      }
      if (iparam <= -1) {
        paramText = Form("%s = %.3f", paramNames[iparam].Data(), ws->var(paramNames[iparam].Data())->getVal());
      }
      else {
        paramText = Form("%s = %.3f #pm %.3f", paramNames[iparam].Data(), ws->var(paramNames[iparam].Data())->getVal(), ws->var(paramNames[iparam].Data())->getError());
      }
      drawText(paramText, pos_text_x_params,pos_text_y_params-pos_y_diff_params*iparam,text_color_params,text_size_params);
    }
    paramText = Form("#Upsilon_{1S} signal = %.0f #pm %.0f", ws->var(Form("nSig1s[%i]",i))->getVal()*scaleFactor, ws->var(Form("nSig1s[%i]",i))->getError()*scaleFactor);
    drawText(paramText, pos_text_x_params, pos_text_y_params-pos_y_diff_params*9, text_color_params, text_size_params);

    // PULL
    pad2[i] = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.30);
    pad2[i]->SetTopMargin(0); // Upper and lower plot are joined
    pad2[i]->SetBottomMargin(0.67);
    pad1[i]->SetLeftMargin(0.18);
    pad1[i]->SetRightMargin(0.02);
    pad2[i]->SetRightMargin(0.02);
    pad2[i]->SetLeftMargin(0.18);
    pad2[i]->SetTicks(1,1);
    pad2[i]->cd();
    hpull[i] = myPlot[i]->pullHist(Form("dataHist[%i]",i),Form("modelHist[%i]",i));
    hpull[i]->SetMarkerSize(0.8);
    pullFrame[i] = ws->var("mass")->frame(Title("Pull Distribution")) ;
    pullFrame[i]->addPlotable(hpull[i],"P") ;
    pullFrame[i]->SetTitleSize(0);
    pullFrame[i]->GetYaxis()->SetTitleOffset(0.43) ;
    pullFrame[i]->GetYaxis()->SetTitle("Pull") ;
    pullFrame[i]->GetYaxis()->SetTitleSize(0.18) ; //19
    pullFrame[i]->GetYaxis()->SetLabelSize(0.113) ; // 113
    pullFrame[i]->GetYaxis()->SetRangeUser(-3.8,3.8) ;
    pullFrame[i]->GetYaxis()->CenterTitle();

    pullFrame[i]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
    pullFrame[i]->GetXaxis()->SetTitleOffset(1.05) ;
    pullFrame[i]->GetXaxis()->SetLabelOffset(0.04) ;
    pullFrame[i]->GetXaxis()->SetLabelSize(0.20) ; //23
    pullFrame[i]->GetXaxis()->SetTitleSize(0.25) ;  //28
    pullFrame[i]->GetXaxis()->CenterTitle();

    pullFrame[i]->GetYaxis()->SetTickSize(0.04);
    pullFrame[i]->GetYaxis()->SetNdivisions(404);
    pullFrame[i]->GetXaxis()->SetTickSize(0.03);
    pullFrame[i]->Draw();

    ypull[i] = hpull[i]->GetY();
    for(int j=0;j<nBins;j++)
    {
      if(ypull[i][j] == 0) continue;
      chisq += TMath::Power(ypull[i][j],2);
      nFullBinsPull++;
    }

    l1[i] = TLine(massLow,0,massHigh,0);
    l1[i].SetLineStyle(9);
    l1[i].Draw("same");
    pad1[i]->Update();
              
    TString label;
    label="";
    if(collId == kPPDATA) CMS_lumi(pad1[i], 1 ,33);
    else if(collId == kAADATA && cLow < 60) CMS_lumi(pad1[i], 2 ,10);
    else if(collId == kPADATA) CMS_lumi(pad1[i], 3 ,33);
    else if(collId == kAADATA && cLow>=60) CMS_lumi(pad1[i], 21 ,33);

    pad1[i]->Update();
    pad2[i]->Update();

    c[i]->cd();
    pad1[i]->Draw();
    pad2[i]->Draw();

    pad1[i]->Update();
    pad2[i]->Update();
  }
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  cout << "numFitPar = " << numFitPar << endl;
  int ndf = nFullBinsPull - numFitPar;
  cout << "chisq/dof = " << chisq/ndf << endl;

  float temp1[4];
  float temp1err[4];
  float temp2[4];
  float temp2err[4];
  float temp3[4];
  float temp3err[4];

  for (int i=0; i<4; i++) {
    temp1[i] = ws->var(Form("nSig1s[%i]",i))->getVal()*scaleFactor;  
    temp1err[i] = ws->var(Form("nSig1s[%i]",i))->getError()*scaleFactor;  
    temp2[i] = ws->var(Form("nSig2s[%i]",i))->getVal()*scaleFactor;  
    temp2err[i] = ws->var(Form("nSig2s[%i]",i))->getError()*scaleFactor;  
    temp3[i] = ws->var(Form("nSig3s[%i]",i))->getVal()*scaleFactor;  
    temp3err[i] = ws->var(Form("nSig3s[%i]",i))->getError()*scaleFactor;  

    cout << "1S signal    =  " << temp1[i] << " +/- " << temp1err[i] << endl;
    cout << "2S signal    =  " << temp2[i] << " +/- " << temp2err[i] << endl;
    cout << "3S signal    =  " << temp3[i] << " +/- " << temp3err[i] << endl;
    cout << "Total signal =  " << temp1[i]+temp2[i]+temp3[i] << endl;

    c[i]->Update();

    c[i]->SaveAs(Form("%ssim_dphi%sfitresults_upsilon_%s_%i.png", directory.Data(), systStr.Data(), kineLabel.Data(), i));
    c[i]->SaveAs(Form("%ssim_dphi%sfitresults_upsilon_%s_%i.pdf", directory.Data(), systStr.Data(), kineLabel.Data(), i ));

  }

} 
 
