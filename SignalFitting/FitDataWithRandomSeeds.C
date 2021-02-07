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
#include "TCanvas.h"
#include "TPad.h"
#include "TRandom3.h"

const double pi = 3.14159265;

//Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

using namespace std;
using namespace RooFit;
double FitDataWithRandomSeeds( 
       int collId = kAADATA,
       float ptLow=0, float ptHigh=50,
       float yLow=2.1, float yHigh=2.4,//Run 1 has p going in -z direction
       int cLow=10, int cHigh=90,//%centrality
       float muPtCut=3.5,
       float dphiEp2Low = 0.0,//units of PI
       float dphiEp2High = 0.12,//range is [0,0.5]
       int whichSyst=0,   
// 0: Nominal
// 1: AltSig
// 2: AltBkg
// 3: AltAcc
// 4: AltEff
       bool randomSeeds=kTRUE
			) 
{

  TString directory = "AllParamFree/";
  //TString directory = "Fits_with_n_fixed/";

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";

  float eta_low = -2.4;
  float eta_high = 2.4;
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  TFile* f1;
  TFile* f2;
  float yLowLab;
  float yHighLab;

  TString kineCut;

  //Select Data Set
  if (whichSyst==3){//altAcc
    f1 = new TFile("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Skimming/skims/newOniaTree_Skim_UpsTrig_MM_flattenedBinByBin_order21_n-1_withDataset_altAcc_20210125.root");
  }
  else if (whichSyst==4){//AltEff
    f1 = new TFile("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Skimming/skims/newOniaTree_Skim_UpsTrig_MM_flattenedBinByBin_order21_n-1_withDataset_altEff_20210125.root");
  }
  else {//nominal dataset
    f1 = new TFile("/home/jared/Documents/Ubuntu_Overflow/Upsilon_v2_502TeV_ThesisNew/Skimming/skims/newOniaTree_Skim_UpsTrig_MM_flattenedBinByBin_order21_n-1_withDataset_20210125.root");
  }
  yLowLab = yLow;
  yHighLab = yHigh;
  kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i && ((abs(dphiEp2)>%.2f && abs(dphiEp2)<%.2f) || (abs(dphiEp2)>%.2f && abs(dphiEp2)<%.2f))",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low, cLow*2,cHigh*2, dphiEp2Low*pi,dphiEp2High*pi, (1-dphiEp2High)*pi,(1-dphiEp2Low)*pi);


  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  if (collId==kPADATA) {
    RooDataSet *dataset2 = (RooDataSet*)f2->Get("dataset");
    dataset->append(*dataset2);
    delete dataset2;
  }
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  delete dataset;
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );


  //**************SIGNAL****************:
  double sigma1s_1_init = 0.07;
  double x1s_init = 0.5;
  double alpha1s_1_init = 2.5;
  double n1s_1_init = 2.5;
  double f1s_init = 0.6;

  //randomize ICs:
  if (randomSeeds) {
    double rangeScale = 0.5;
    TRandom3 rnd3(0);
    sigma1s_1_init = rangeScale*(paramsupper[0]-paramslower[0])*(rnd3.Rndm()-0.5) + (paramsupper[0]+paramslower[0])/2;
    x1s_init = rangeScale*(paramsupper[1]-paramslower[1])*(rnd3.Rndm()-0.5) + (paramsupper[1]+paramslower[1])/2;
    alpha1s_1_init = rangeScale*(paramsupper[2]-paramslower[2])*(rnd3.Rndm()-0.5) + (paramsupper[2]+paramslower[2])/2;
    n1s_1_init = rangeScale*(paramsupper[3]-paramslower[3])*(rnd3.Rndm()-0.5) + (paramsupper[3]+paramslower[3])/2;
    f1s_init = rangeScale*(paramsupper[4]-paramslower[4])*(rnd3.Rndm()-0.5) + (paramsupper[4]+paramslower[4])/2;
  }

  //Print ICs:
  cout << endl;
  cout << "SEEDS:" << endl;
  cout << "sigma1s_1_init = " << sigma1s_1_init << endl;
  cout << "x1s_init = " << x1s_init << endl;
  cout << "alpha1s_1_init = " << alpha1s_1_init << endl;
  cout << "n1s_1_init = " << n1s_1_init << endl;
  cout << "f1s_init = " << f1s_init << endl << endl;

  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, paramslower[0], paramsupper[0]);
  RooFormulaVar sigma2s_1("sigma2s_1","@0*@1",RooArgList(sigma1s_1,mRatio21) );
  RooFormulaVar sigma3s_1("sigma3s_1","@0*@1",RooArgList(sigma1s_1,mRatio31) );

  RooRealVar x1s = RooRealVar("x1s","sigma ratio ", x1s_init, paramslower[1], paramsupper[1]);

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, x1s) );
  RooFormulaVar sigma2s_2("sigma2s_2","@0*@1",RooArgList(sigma1s_2,mRatio21) );
  RooFormulaVar sigma3s_2("sigma3s_2","@0*@1",RooArgList(sigma1s_2,mRatio31) );

  RooRealVar alpha1s_1("alpha1s_1","tail shift", alpha1s_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar alpha2s_1("alpha2s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_1("alpha3s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha2s_2("alpha2s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_2("alpha3s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", n1s_1_init, paramslower[3], paramsupper[3]);
  RooFormulaVar n2s_1("n2s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_1("n3s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n2s_2("n2s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_2("n3s_2","1.0*@0",RooArgList(n1s_1) );

  RooRealVar f1s = RooRealVar("f1s","1S CB fraction", f1s_init, paramslower[4], paramsupper[4]);
  RooFormulaVar f2s("f2s","1.0*@0",RooArgList(f1s) );
  RooFormulaVar f3s("f3s","1.0*@0",RooArgList(f1s) );

  // Set up crystal ball shapes
  RooCBShape cb1s_1 = RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape cb2s_1 = RooCBShape("cball2s_1", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape cb3s_1 = RooCBShape("cball3s_1", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_1, alpha3s_1, n3s_1);

  RooGaussian gauss1s = RooGaussian("gauss1s","gaussian PDF",*(ws->var("mass")),mean1s,sigma1s_2);
  RooGaussian gauss2s = RooGaussian("gauss2s","gaussian PDF",*(ws->var("mass")),mean2s,sigma2s_2);
  RooGaussian gauss3s = RooGaussian("gauss3s","gaussian PDF",*(ws->var("mass")),mean3s,sigma3s_2);

  RooCBShape cb1s_2 = RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape cb2s_2 = RooCBShape("cball2s_2", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape cb3s_2 = RooCBShape("cball3s_2", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_2, alpha3s_2, n3s_2);

  RooAddPdf* cb1s;
  RooAddPdf* cb2s;
  RooAddPdf* cb3s;
if (whichSyst==1) {
  //AltSig: CB+GAUSSIAN
  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(cb1s_1,gauss1s), RooArgList(f1s) );
  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(cb2s_1,gauss2s), RooArgList(f1s) );
  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(cb3s_1,gauss3s), RooArgList(f1s) );
}
else {
  //Nominal: DOUBLE CRYSTAL BALL
  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(cb1s_1,cb1s_2), RooArgList(f1s) );
  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(cb2s_1,cb2s_2), RooArgList(f1s) );
  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(cb3s_1,cb3s_2), RooArgList(f1s) );
}

  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",50000,0,1000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",1000,0,360000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",10,0,260000);


  //**********BACKGROUND****************

  //CHEBYCHEV
  RooRealVar ach1("Ach1","Acheb1",0,-1,1);
  RooRealVar ach2("Ach2","Acheb2",-0.1,-1,1);
  RooRealVar ach3("Ach3","Acheb3",0,-1,1);
  RooRealVar ach4("Ach4","Acheb4",0,-1,1);
  RooChebychev* bkgCheb = new RooChebychev("bkgLowPt","Background",*(ws->var("mass")),RooArgList(ach1,ach2,ach3,ach4));
  
  //POWER LAW
  RooRealVar m0("m0","m0",1,0,100);
  RooRealVar pow("pow","pow",10,0,100);
  RooRealVar mpow("mpow","mpow",0,0,100);
  RooGenericPdf* bkgPow = new RooGenericPdf("bkgHighPt","Background","TMath::Power(@0,@3)/TMath::Power(1+@0/@1,@2)",RooArgList(*(ws->var("mass")),m0,pow,mpow));

  //NOMINAL BACKGROUND
  double m_lambda_init = 8;
  double err_mu_init = 8;
  double err_sigma_init = 1;
  if (randomSeeds) {
    double rangeScale = 0.5;
    TRandom3 rnd3(0);
    m_lambda_init = rangeScale*(paramsupper[7]-paramslower[7])*(rnd3.Rndm()-0.5) + (paramsupper[7]+paramslower[7])/2;
    err_sigma_init = rangeScale*(paramsupper[6]-paramslower[6])*(rnd3.Rndm()-0.5) + (paramsupper[6]+paramslower[6])/2;
    err_mu_init = rangeScale*(paramsupper[5]-paramslower[5])*(rnd3.Rndm()-0.5) + (paramsupper[5]+paramslower[5])/2;
  }
  RooRealVar err_mu("#mu","err_mu", err_mu_init, paramslower[5], paramsupper[5]) ;
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);
  
  RooGenericPdf *bkgErfExp = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );
  RooGenericPdf *bkgExp = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));

  RooAbsPdf *bkg;
  RooAbsPdf *bkgLowPt;
  RooAbsPdf *bkgHighPt;

  if (whichSyst==2) {//alternative background
    //low pt: chebychev
    bkgLowPt = bkgCheb;
    //high pt: power law
    bkgHighPt = bkgPow;
  }
  else {//nominal background:
    //low pt: erf*exponential
    bkgLowPt = bkgErfExp;
    //high pt: exponential
    bkgHighPt = bkgExp;
  }

  if (ptLow >= 5) bkg = bkgHighPt;
  else bkg = bkgLowPt;

  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",25000,0,5000000);

  //fix parameters
  float nfix;
  float xfix;
  float alphafix;
  float ffix;
  float sigmafix;
  float lambdafix;
  float errmufix;
  float errsigmafix;
  if (collId==kPPDATA) float nfix = 2.14973;
  else if (collId==kAADATA) float nfix = 3.59489;
  //n1s_1.setVal(nfix);
  //x1s->setVal(0.6);
  //alpha1s_1.setVal(1.5);
  if (dphiEp2High-dphiEp2Low < 0.5) {
    directory = "dphiFits/";
    //import ICs
    cout << "Importing workspace" << endl;
    TString kineLabelICs = getKineLabel(collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
    TString NomFileName = Form("AllParamFree/%sfitresults_upsilon_%s.root", systStr.Data(), kineLabelICs.Data());
    cout << NomFileName << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    nfix = Nomws->var("n1s_1")->getVal();
    xfix = Nomws->var("x1s")->getVal();
    alphafix = Nomws->var("alpha1s_1")->getVal();
    ffix = Nomws->var("f1s")->getVal();
    sigmafix = Nomws->var("sigma1s_1")->getVal();
    //lambdafix = Nomws->var("#lambda")->getVal();
    //errmufix = Nomws->var("#mu")->getVal();
    //errsigmafix = Nomws->var("#sigma")->getVal();
    //cout << "lambdafix = " << lambdafix << endl;
    n1s_1.setVal(nfix);
    n1s_1.setConstant();
    x1s.setVal(xfix);
    x1s.setConstant();
    alpha1s_1.setVal(alphafix);
    alpha1s_1.setConstant();
    f1s.setVal(ffix);
    f1s.setConstant();
    sigma1s_1.setVal(sigmafix);
    sigma1s_1.setConstant();
    //err_mu.setVal(errmufix);
    //err_mu.setConstant();
    //err_sigma.setVal(errsigmafix);
    //err_sigma.setConstant();
    //m_lambda.setVal(lambdafix);
    //m_lambda.setConstant();
    delete Nomws;
    NomFile->Close("R");
    delete NomFile;
    cout << endl;
    cout << "IMPORTED SEEDS:" << endl;
    cout << "sigma1s_1_init = " << sigmafix << endl;
    cout << "x1s_init = " << xfix << endl;
    cout << "alpha1s_1_init = " << alphafix << endl;
    cout << "n1s_1_init = " << nfix << endl;
    cout << "f1s_init = " << ffix << endl << endl;
  }

  //Build the model
  RooAddPdf* model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));

  ws->import(*model);

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  //Fit the model to the data
  cout << "Now fitting..." << endl;
  RooFitResult* fitRes2 = (RooFitResult*)ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  ws->import(*fitRes2);

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.43);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.054);
  myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();
  fitRes2->Print("v");
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

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

  TLegend fitleg = TLegend(0.76,0.4,0.91,0.7); fitleg.SetTextSize(19);
  fitleg.SetTextFont(43);
  fitleg.SetBorderSize(0);
  fitleg.AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
  fitleg.AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  fitleg.AddEntry(myPlot2->findObject("Sig1S"),"Signal","l");
  fitleg.AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
  fitleg.Draw("same");

  // PULL
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.30);
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad1->SetLeftMargin(0.18);
  pad1->SetRightMargin(0.02);
  pad2->SetRightMargin(0.02);
  pad2->SetLeftMargin(0.18);
  pad2->SetTicks(1,1);
  pad2->cd();
  
  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.43) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.18) ; //19
  pullFrame->GetYaxis()->SetLabelSize(0.113) ; // 113
  pullFrame->GetYaxis()->SetRangeUser(-3.8,3.8) ;
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetLabelSize(0.20) ; //23
  pullFrame->GetXaxis()->SetTitleSize(0.25) ;  //28
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw() ;

  //calculate chi-squared
  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin; 
  double *ypull = hpull->GetY();
  for(int i=0;i<nBins;i++)
  {
    if(ypull[i] == 0) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;
  cout << "chisq/dof = " << chisq << "/" << ndf << " = " << chisq/ndf << endl;

  //continue beautifying the plot and print out results
  TLine l1 = TLine(massLow,0,massHigh,0);
  l1.SetLineStyle(9);
  l1.Draw("same");
  pad1->Update();
              
  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";

  TString label;
  label="";
  if(collId == kPPDATA) CMS_lumi(pad1, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(pad1, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(pad1, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(pad1, 21 ,33);

  pad1->Update();
  pad2->Update();

  c1->cd();
  pad1->Draw();
  pad2->Draw();

  pad1->Update();
  pad2->Update();

  c1->SaveAs(Form("%s%sfitresults_upsilon_%s.png", directory.Data(), systStr.Data(), kineLabel.Data() ));
  c1->SaveAs(Form("%s%sfitresults_upsilon_%s.pdf", directory.Data(), systStr.Data(), kineLabel.Data() ));

  float temp1 = ws->var("nSig1s")->getVal();  
  float temp1err = ws->var("nSig1s")->getError();  
  float temp2 = ws->var("nSig2s")->getVal();  
  float temp2err = ws->var("nSig2s")->getError();  
  float temp3 = ws->var("nSig3s")->getVal();  
  float temp3err = ws->var("nSig3s")->getError();  

  cout << "1S signal    =  " << temp1 << " +/- " << temp1err << endl;
  cout << "2S signal    =  " << temp2 << " +/- " << temp2err << endl;
  cout << "3S signal    =  " << temp3 << " +/- " << temp3err << endl;

	cout << "if ( binMatched( "<<muPtCut<<",  " << ptLow <<", "<< ptHigh<<", "<< yLow<<", "<< yHigh << ", " << cLow << ", " << cHigh << ") ) " ; 
  cout << "  { setSignalParMC( " ;
  cout <<  ws->var("n1s_1")->getVal() << ", " <<  ws->var("alpha1s_1")->getVal() << ", "<<  ws->var("sigma1s_1")->getVal() << ", " ;
  cout <<  ws->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  ws->var("f1s")->getVal() << ", "<<  ws->var("x1s")->getVal() << " );} " << endl;

TString outFileName = Form("%s%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data() );
  TFile* outf = new TFile(outFileName,"recreate");
  c1->Write();
  ws->Write();
  outf->Close();

  cout << "here1" << endl;
  delete f1;
  cout << "here2" << endl;
  if (collId==kPADATA) delete f2;
  cout << "here3" << endl;
  //delete ws;
  cout << "here4" << endl;
  delete cb1s;
  cout << "here5" << endl;
  delete cb2s;
  cout << "here6" << endl;
  delete cb3s;
  cout << "here7" << endl;
  delete nSig1s;
  cout << "here8" << endl;
  delete nSig2s;
  cout << "here9" << endl;
  delete nSig3s;
  cout << "here10" << endl;
  delete nBkg;
  cout << "here11" << endl;
  delete bkgLowPt;
  cout << "here12" << endl;
  delete bkgHighPt;
  cout << "here13" << endl;
  delete model;
  cout << "here14" << endl;
  delete fitRes2;
  cout << "here15" << endl;
  delete reducedDS;
  cout << "here16" << endl;
  //delete myPlot;
  cout << "here17" << endl;
  //delete myPlot2;
  cout << "here18" << endl;
  delete pad2;
  cout << "here19" << endl;
  delete pad1;
  cout << "here20" << endl;
  delete c1;
  cout << "here21" << endl;

  //Before closing root, you can type: gDirectory->ls("-m") to see what objects are still in memory.
  cout << endl << "Here's what's in memory right now" << endl;
  gDirectory->ls("-m");
  cout << "that's all." << endl << endl;

  return chisq/ndf;
  //return 0.0;
} 
 
