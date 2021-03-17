#include "../HeaderFiles/SONGKYO.h"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/cutsAndBin.h"



void draw_v2(int whichUpsilon=1) {

  float ptbins[5] = {0,3,6,10,50};
  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  float ybins[4] = {0.0,0.8,1.6,2.4};
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  float cbins[4] = {10,30,50,90};
  const int numcbins = sizeof(cbins)/sizeof(float)-1;

  setTDRStyle();
  writeExtraText = true;

  //Get v2 histograms
  TFile* inFile = new TFile(Form("../v2Fitting/Ups_%i_v2_cent10-90_nom.root",whichUpsilon),"READ");
  TH1D* hv2pt = (TH1D*)inFile->Get("hv2pt");
  TH1D* hv2y = (TH1D*)inFile->Get("hv2y");
  TH1D* hv2c = (TH1D*)inFile->Get("hv2c");

  //Apply event plane resolution correction
  bool EpResCorrection = kTRUE;
  TString resCorFileName = "../Skimming/condor/averages/resCorFile_n56114317_BinByBin.root";
  TFile* EPRFile = new TFile(resCorFileName,"READ");
  TH1D* hRpt = (TH1D*)EPRFile->Get("hRpt");
  TH1D* hRy = (TH1D*)EPRFile->Get("hRy");
  TH1D* hRc = (TH1D*)EPRFile->Get("hRc");
  hv2pt->Sumw2();
  hRpt->Sumw2();
  hRy->Sumw2();
  hRc->Sumw2();
  if (EpResCorrection) {
    cout << " Applying resoulution correction" << endl;
    hv2pt->Divide(hRpt);
    hv2y->Divide(hRy);
    hv2c->Divide(hRc);
  }

  //make histograms of systematics
  TString systFileName = "../Systematics/Ups_1_v2_cent10-90_SystCombined.root";
  TFile* systFile = new TFile(systFileName,"READ");
  TH1D* hsyspt = (TH1D*)systFile->Get("hv2pt");
  TH1D* hsysy = (TH1D*)systFile->Get("hv2y");
  TH1D* hsysc = (TH1D*)systFile->Get("hv2c");

  //make TGraphs
  TGraphErrors* gv2pt = new TGraphErrors(hv2pt);
  TGraphErrors* gv2pt_sys = new TGraphErrors(hv2pt);
  for (int ipt=0; ipt<numptbins; ipt++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2pt->GetPoint(ipt, pxtmp, pytmp);
    extmp=gv2pt->GetErrorX(ipt);
    eytmp=gv2pt->GetErrorY(ipt);
    relsys=hsyspt->GetBinContent(ipt+1);
    //remove ex
    gv2pt->SetPointError(ipt, 0, eytmp);
    //set ey for systematic error
    gv2pt_sys->SetPointError(ipt, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2pt, 1, 1); 
  SetGraphStyleSys(gv2pt_sys, 1);

  float ptmin=0; float ptmax=50;
  gv2pt_sys->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  gv2pt_sys->GetXaxis()->CenterTitle();
  gv2pt_sys->GetXaxis()->SetTitleOffset(1.);
  gv2pt_sys->GetXaxis()->SetLimits(ptmin,ptmax);
  gv2pt_sys->GetXaxis()->SetTitleSize(0.06);
  gv2pt_sys->GetYaxis()->SetTitle("v_{2}");
  gv2pt_sys->GetYaxis()->CenterTitle();
  gv2pt_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2pt_sys->GetYaxis()->SetTitleSize(0.06);
  gv2pt_sys->SetMinimum(-0.1);
  gv2pt_sys->SetMaximum(0.2);

  TGraphErrors* gv2y = new TGraphErrors(hv2y);
  TGraphErrors* gv2y_sys = new TGraphErrors(hv2y);
  for (int iy=0; iy<numybins; iy++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2y->GetPoint(iy, pxtmp, pytmp);
    extmp=gv2y->GetErrorX(iy);
    eytmp=gv2y->GetErrorY(iy);
    relsys=hsysy->GetBinContent(iy+1);
    //remove ex
    gv2y->SetPointError(iy, 0, eytmp);
    //set ey for systematic error
    gv2y_sys->SetPointError(iy, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2y, 1, 1); 
  SetGraphStyleSys(gv2y_sys, 1);

  float ymin=0; float ymax=2.4;
  gv2y_sys->GetXaxis()->SetTitle("y^{#varUpsilon}");
  gv2y_sys->GetXaxis()->CenterTitle();
  gv2y_sys->GetXaxis()->SetTitleOffset(1.);
  gv2y_sys->GetXaxis()->SetLimits(ymin,ymax);
  gv2y_sys->GetXaxis()->SetTitleSize(0.06);
  gv2y_sys->GetYaxis()->SetTitle("v_{2}");
  gv2y_sys->GetYaxis()->CenterTitle();
  gv2y_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2y_sys->GetYaxis()->SetTitleSize(0.06);
  gv2y_sys->SetMinimum(-0.1);
  gv2y_sys->SetMaximum(0.2);

  TGraphErrors* gv2c = new TGraphErrors(hv2c);
  TGraphErrors* gv2c_sys = new TGraphErrors(hv2c);
  for (int ic=0; ic<numcbins; ic++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2c->GetPoint(ic, pxtmp, pytmp);
    extmp=gv2c->GetErrorX(ic);
    eytmp=gv2c->GetErrorY(ic);
    relsys=hsysc->GetBinContent(ic+1);
    //remove ex
    gv2c->SetPointError(ic, 0, eytmp);
    //set ey for systematic error
    gv2c_sys->SetPointError(ic, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2c, 1, 1); 
  SetGraphStyleSys(gv2c_sys, 1);

  float cmin=10; float cmax=90;
  gv2c_sys->GetXaxis()->SetTitle("Centrality (%)");
  gv2c_sys->GetXaxis()->CenterTitle();
  gv2c_sys->GetXaxis()->SetTitleOffset(1.);
  gv2c_sys->GetXaxis()->SetLimits(cmin,cmax);
  gv2c_sys->GetXaxis()->SetTitleSize(0.06);
  gv2c_sys->GetYaxis()->SetTitle("v_{2}");
  gv2c_sys->GetYaxis()->CenterTitle();
  gv2c_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2c_sys->GetYaxis()->SetTitleSize(0.06);
  gv2c_sys->SetMinimum(-0.1);
  gv2c_sys->SetMaximum(0.2);

  TCanvas* cpt = new TCanvas("cpt","cpt",40,40,600,600);

  gv2pt_sys->Draw("A5");
  gv2pt->Draw("P");

  TLegend *legpt= new TLegend(0.7, 0.6, 0.9, 0.7);
  SetLegendStyle(legpt);
  legpt->AddEntry(gv2pt,Form(" #Upsilon(%iS)",whichUpsilon),"lp");

  legpt->Draw("same");
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  TLine *lpt = new TLine(0,0,50,0);
  lpt->SetLineStyle(9);
  lpt->Draw();

  int collId = kAADATA;
  int cLow = 0;
  if(collId == kPPDATA) CMS_lumi(cpt, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(cpt, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(cpt, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(cpt, 21 ,33);

  cpt->SaveAs(Form("Plots/v2_vs_pt_%is.png",whichUpsilon));
  cpt->SaveAs(Form("Plots/v2_vs_pt_%is.pdf",whichUpsilon));

  TCanvas* cy = new TCanvas("cy","cy",40,40,600,600);

  gv2y_sys->Draw("A5");
  gv2y->Draw("P");

  TLegend *legy= new TLegend(0.7, 0.6, 0.9, 0.7);
  SetLegendStyle(legy);
  legy->AddEntry(gv2y,Form(" #Upsilon(%iS)",whichUpsilon),"lp");

  legy->Draw("same");
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  TLine *ly = new TLine(0,0,2.4,0);
  ly->SetLineStyle(9);
  ly->Draw();

  if(collId == kPPDATA) CMS_lumi(cy, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(cy, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(cy, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(cy, 21 ,33);

  cy->SaveAs(Form("Plots/v2_vs_y_%is.png",whichUpsilon));
  cy->SaveAs(Form("Plots/v2_vs_y_%is.pdf",whichUpsilon));

  TCanvas* cc = new TCanvas("cc","cc",40,40,600,600);

  gv2c_sys->Draw("A5");
  gv2c->Draw("P");

  TLegend *legc= new TLegend(0.7, 0.6, 0.9, 0.7);
  SetLegendStyle(legc);
  legc->AddEntry(gv2c,Form(" #Upsilon(%iS)",whichUpsilon),"lp");

  legc->Draw("same");
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  TLine *lc = new TLine(10,0,90,0);
  lc->SetLineStyle(9);
  lc->Draw();

  if(collId == kPPDATA) CMS_lumi(cc, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(cc, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(cc, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(cc, 21 ,33);

  cc->SaveAs(Form("Plots/v2_vs_c_%is.png",whichUpsilon));
  cc->SaveAs(Form("Plots/v2_vs_c_%is.pdf",whichUpsilon));

}
