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
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RoundsHeader.h"

const double pi = 3.14159265;

//Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
double paramsupper[8] = {0.35, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
double paramslower[8] = {0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

double randomIC(double lo=0, double up=1) {
  double rangeScale = 0.5;
  TRandom3 rnd3(0);
  return rangeScale*(up-lo)*(rnd3.Rndm()-0.5) + (up+lo)/2;
}

using namespace std;
using namespace RooFit;
void Fit4dphiBinsSimultaneously( 
       int collId = kAADATA,  
       float ptLow=0, float ptHigh=3, 
       float yLow=2.1, float yHigh=2.4,//Run 1 has p going in -z direction
       int cLow=10, int cHigh=90,
       float muPtCut=3.5,
       int whichSyst=4,   
// 0: Nominal
// 1: AltSig
// 2: AltBkg
// 3: AltAcc
// 4: AltEff
// 5: AltConst
       bool randomSeeds=kTRUE,
       int whichRound=R4a
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

  if (whichSyst==5) whichRound = R4b;

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

  TString kineCut[4];
  float dphibins[5] = {0.0,0.125,0.25,0.375,0.5};

  cout << endl << "PERFORMING SIMULTANEOUS FIT IN DPHI BINS [" << dphibins[0] << "," << dphibins[1] << "], [" << dphibins[1] << "," << dphibins[2] << "], [" << dphibins[2] << "," << dphibins[3] << "], [" << dphibins[3] << "," << dphibins[4] << "]" << endl << endl;

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
  for (int i=0; i<4; i++) {
    kineCut[i] = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && cBin>%i && cBin<%i && ((abs(dphiEp2)>=%.2f && abs(dphiEp2)<%.2f) || (abs(dphiEp2)>=%.2f && abs(dphiEp2)<%.2f))",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low, cLow*2,cHigh*2, dphibins[i]*pi, dphibins[i+1]*pi, (1-dphibins[i+1])*pi,(1-dphibins[i])*pi);
    if (muPtCut>0) kineCut[i] = kineCut[i] + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  }

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  dataset->Print();
  cout << "####################################" << endl;
  RooDataSet *reducedDS[4];
  for (int i=0; i<4; i++) {
    reducedDS[i] = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut[i].Data() );
    reducedDS[i]->SetName(Form("reducedDS[%i]",i));
    ws->import(*reducedDS[i]);
    reducedDS[i]->Print();
  }
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  RooCategory tp("tp","tp");
  tp.defineType("A");
  tp.defineType("B");
  tp.defineType("C");
  tp.defineType("D");

  // Create a dataset that imports contents of all the above datasets mapped by index category tp
  RooDataSet* dsABCD = new RooDataSet("dsABCD","dsABCD",RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS[0]),Import("B",*reducedDS[1]),Import("C",*reducedDS[2]),Import("D",*reducedDS[3]));
  cout << "******** New Combined Dataset ***********" << endl;
  dsABCD->Print();
  ws->import(*dsABCD);

  //import info from full bin
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
  TString NomFileName = Form("RoundFits_R4a/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");
  double totSig1 = Nomws->var("nSig1s")->getVal();
  double totSig2 = Nomws->var("nSig2s")->getVal();
  double totSig3 = Nomws->var("nSig3s")->getVal();
  double totBkg = Nomws->var("nBkg")->getVal();

  //SIGNAL:
  double sigma1s_1_init = 0.1;
  double x1s_init = 0.6;
  double alpha1s_1_init = 1.5;
  double n1s_1_init = 2.0;
  double f1s_init = 0.5;

  //randomize ICs:
  if (randomSeeds) {
    sigma1s_1_init = randomIC(paramslower[0],paramsupper[0]);
    x1s_init = randomIC(paramslower[1],paramsupper[1]);
    alpha1s_1_init = randomIC(paramslower[2],paramsupper[2]);
    n1s_1_init = randomIC(paramslower[3],paramsupper[3]);
    f1s_init = randomIC(paramslower[4],paramsupper[4]);
  }

  //Print ICs:
  cout << endl;
  cout << "SEEDS:" << endl;
  cout << "sigma1s_1_init = " << sigma1s_1_init << endl;
  cout << "x1s_init = " << x1s_init << endl;
  cout << "alpha1s_1_init = " << alpha1s_1_init << endl;
  cout << "n1s_1_init = " << n1s_1_init << endl;
  cout << "f1s_init = " << f1s_init << endl << endl;

  TCanvas* c[4];
  RooPlot* myPlot[4];

  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );


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

  RooRealVar *nSig1s[4];
  RooRealVar *nSig2s[4];
  RooRealVar *nSig3s[4];


  //**********BACKGROUND****************

  //CHEBYCHEV
  double ach1_init = 0;
  double ach2_init = -0.1;
  double ach3_init = 0;
  double ach4_init = 0;
  if (randomSeeds) {
    ach1_init = randomIC(-1,1);
    ach2_init = randomIC(-1,1);
    ach3_init = randomIC(-1,1);
    ach4_init = randomIC(-1,1);
  }
  RooRealVar ach1("Ach1","Acheb1",ach1_init,-1,1);
  RooRealVar ach2("Ach2","Acheb2",ach2_init,-1,1);
  RooRealVar ach3("Ach3","Acheb3",ach3_init,-1,1);
  RooRealVar ach4("Ach4","Acheb4",ach4_init,-1,1);
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
    m_lambda_init = randomIC(paramslower[7],paramsupper[7]);
    err_sigma_init = randomIC(paramslower[6],paramsupper[6]);
    err_mu_init = randomIC(paramslower[5],paramsupper[5]);
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

  RooRealVar *nBkg[4];

  for (int i=0; i<4; i++) {
    c[i] =  new TCanvas(Form("canvas[%i]",i),"My plots",4,45,550,520);
    c[i]->cd();
    myPlot[i] = ws->var("mass")->frame(nMassBin); // bins
    ws->data(Form("reducedDS[%i]",i))->plotOn(myPlot[i],Name(Form("dataHist[%i]",i)));

    nSig1s[i]= new RooRealVar(Form("nSig1s[%i]",i)," 1S signals",totSig1/4,0,totSig1);
    nSig2s[i]= new RooRealVar(Form("nSig2s[%i]",i)," 2S signals",totSig2/4,-20,totSig2);
    nSig3s[i]= new RooRealVar(Form("nSig3s[%i]",i)," 3S signals",totSig3/4,-50,totSig3);

    nBkg[i] = new RooRealVar(Form("nBkg[%i]",i),"fraction of component 1 in bkg",totBkg/4,0,totBkg);
  }

  // Construct Gaussian constraints on parameters for the final round of fitting
  //bool badbin = kFALSE;
  //if (abs(ptLow-6)<0.1 && abs(ptHigh-9)<0.1 && abs(yLow+2.87)<0.1 && abs(yHigh-1.93)<0.1) badbin = kTRUE;
  float nmu, xmu, alphamu, fmu;
  float ndev, xdev, alphadev, fdev;
  //nominal constraints derived from path a
  if (whichRound==R4a) {
    alphamu = 1.225790;
    nmu = 3.828210;
    xmu = 0.426153;
    alphadev = 0.230156;
    ndev = 0.720624;
    xdev = 0.110304;
    TString kineLabelICs = getKineLabel(collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
    TString NomFileName = Form("RoundFits_R3a/nomfitresults_upsilon_%s.root", kineLabelICs.Data());
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    fmu = Nomws->var("f1s")->getVal();
    fdev = Nomws->var("f1s")->getError();
    delete Nomws;
    NomFile->Close("R");
    delete NomFile;
  }
  else if (whichRound==R4b) {
    //alternate constraints derived from path b
    alphamu = 1.798211;
    nmu = 1.784743;
    xmu = 0.451104;
    alphadev = 0.942203;
    ndev = 1.218936;
    xdev = 0.088950;
    TString kineLabelICs = getKineLabel(collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);
    TString NomFileName = Form("RoundFits_R3b/nomfitresults_upsilon_%s.root", kineLabelICs.Data());
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    fmu = Nomws->var("f1s")->getVal();
    fdev = Nomws->var("f1s")->getError();
    delete Nomws;
    NomFile->Close("R");
    delete NomFile;
  }

  RooGaussian nconstraint("nconstraint","nconstraint", n1s_1,RooConst(nmu),RooConst(ndev));
  RooGaussian alphaconstraint("alphaconstraint","alphaconstraint", alpha1s_1,RooConst(alphamu),RooConst(alphadev));
  RooGaussian xconstraint("xconstraint","xconstraint", x1s,RooConst(xmu),RooConst(xdev));
  RooGaussian fconstraint("fconstraint","fconstraint", f1s,RooConst(fmu),RooConst(fdev));

  RooArgSet *allConstraints;
  RooArgSet *constParams;
  allConstraints = new RooArgSet(nconstraint, alphaconstraint, xconstraint, fconstraint);
  constParams = new RooArgSet(n1s_1,alpha1s_1,x1s,f1s);
  //allConstraints = new RooArgSet(nconstraint, alphaconstraint, xconstraint);
  //constParams = new RooArgSet(n1s_1,alpha1s_1,x1s);
  //ws->import(*allConstraints);

  //Build the model
  RooAddPdf* model[4];
  RooPlot* myPlot2[4];
  for (int i=0; i<4; i++) {
    model[i] = new RooAddPdf(Form("model[%i]",i),"1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s[i],*nSig2s[i],*nSig3s[i],*nBkg[i]));
    //ws->import(*model[i]);
    c[i]->cd();
    myPlot2[i] = (RooPlot*)myPlot[i]->Clone();
    ws->data(Form("reducedDS[%i]",i))->plotOn(myPlot2[i],Name(Form("dataOS_FIT[%i]",i)),MarkerSize(.8));
  }

  // Construct simultaneous PDF
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",tp);
  simPdf->addPdf(*model[0],"A") ;
  simPdf->addPdf(*model[1],"B") ;
  simPdf->addPdf(*model[2],"C") ;
  simPdf->addPdf(*model[3],"D") ;
  RooProdPdf modelC("model","model with constraints",RooArgSet(*simPdf,*allConstraints));
  //ws->import(*simPdf);
  ws->import(modelC);

  cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("model")->fitTo(*dsABCD, Constrain(*constParams), Save(), Hesse(kTRUE), Range(massLow,massHigh), Timer(kTRUE), Extended(kTRUE));
  cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;
  cout << endl << "Importing fit result..." << endl;
  ws->import(*fitResSim);

  cout << endl << "Plotting fit functions..." << endl;
  for (int i=0; i<4; i++) {
    c[i]->cd();
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Name(Form("modelHist[%i]",i)), Range(massLow, massHigh));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Name("Sig1S"),Components(RooArgSet(*cb1s)), Range(massLow, massHigh),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Components(RooArgSet(*cb2s)), Range(massLow, massHigh),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Components(RooArgSet(*cb3s)), Range(massLow, massHigh),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
    ws->pdf(Form("model[%i]",i))->plotOn(myPlot2[i],Name("bkgPDF"),Components(RooArgSet(*bkg)), Range(massLow, massHigh),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

    //make a pretty plot
    myPlot2[i]->SetFillStyle(4000);
    myPlot2[i]->SetAxisRange(massLowForPlot, massHighForPlot,"X");
    myPlot2[i]->GetYaxis()->SetTitleOffset(1.43);
    myPlot2[i]->GetYaxis()->CenterTitle();
    //myPlot2[i]->GetYaxis()->SetTitleSize(0.058);
    myPlot2[i]->GetYaxis()->SetLabelSize(0.054);
    //myPlot2[i]->GetXaxis()->SetLabelSize(0);
    myPlot2[i]->GetXaxis()->SetRangeUser(8,14);
    //myPlot2[i]->GetXaxis()->SetTitleSize(0);
    myPlot2[i]->Draw();

    c[i]->Update();
    c[i]->SaveAs(Form("%ssim_dphi%sfitresults_upsilon_%s_%i.png", directory.Data(), systStr.Data(), kineLabel.Data(), i));
    c[i]->SaveAs(Form("%ssim_dphi%sfitresults_upsilon_%s_%i.pdf", directory.Data(), systStr.Data(), kineLabel.Data(), i ));
  }

  fitResSim->Print("v");

  cout << endl << "Creating output file..." << endl;
  TString outFileName = Form("%ssim_dphi%sfitresults_upsilon_%s.root", directory.Data(), systStr.Data(), kineLabel.Data());

  TFile* outf = new TFile(outFileName,"recreate");
  ws->Write();
  outf->Close();

  cout << endl << "Finished everything." << endl;

//    delete cb1s;
//    delete cb2s;
//    delete cb3s;
//    delete bkgLowPt;
//    delete bkgHighPt;
//  for (int i=0; i<4; i++) {
//    delete c[i];
//    delete nSig1s[i];
//    delete nSig2s[i];
//    delete nSig3s[i];
//    delete nBkg[i];
//    delete model[i];
//    delete reducedDS[i];
//  }

  delete f1;
//  delete simPdf;
//  delete dsABCD;
//  delete outf;

  cout << endl << "What's left in memory:" << endl;
  gDirectory->ls("-m");
} 
 
