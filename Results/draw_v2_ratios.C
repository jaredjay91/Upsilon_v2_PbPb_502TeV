#include "../HeaderFiles/SONGKYO.h"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/cutsAndBin.h"


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

void draw_v2_ratios() {

  float ptbins[5] = {0,3,6,10,50};
  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  float ybins[4] = {0.0,0.8,1.6,2.4};
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  float cbins[4] = {10,30,50,90};
  const int numcbins = sizeof(cbins)/sizeof(float)-1;

  float plotmin = -6;
  float plotmax = 10;

  setTDRStyle();
  writeExtraText = false;

  float legx1=0.6;
  float legy1=0.25;
  float legx2=0.8;
  float legy2=0.35;

  int styleParam = 3;

  //Get v2 histograms
  cout << "getting 1s histos..." << endl;
  TFile* inFile1s = new TFile(Form("../v2Fitting/Ups_%i_v2_cent10-90_nom.root",1),"READ");
  TH1D* hv2pt1s = (TH1D*)inFile1s->Get("hv2pt");
  TH1D* hv2y1s = (TH1D*)inFile1s->Get("hv2y");
  TH1D* hv2c1s = (TH1D*)inFile1s->Get("hv2c");
  //inFile1s->Close();
  //delete inFile1s;

  cout << "getting 2s histos..." << endl;
  TFile* inFile2s = new TFile(Form("../v2Fitting/Ups_%i_v2_cent10-90_nom.root",2),"READ");
  TH1D* hv2pt2s = (TH1D*)inFile2s->Get("hv2pt");
  TH1D* hv2y2s = (TH1D*)inFile2s->Get("hv2y");
  TH1D* hv2c2s = (TH1D*)inFile2s->Get("hv2c");
  //inFile2s->Close();
  //delete inFile2s;

  //The resolution correction is unnecessary because it cancels out in the ratio.

  //make ratio histos
  cout << "dividing histos..." << endl;
  TH1D* hratiopt = (TH1D*)hv2pt2s->Clone("hratiopt");
  TH1D* hratioy = (TH1D*)hv2y2s->Clone("hratioy");
  TH1D* hratioc = (TH1D*)hv2c2s->Clone("hratioc");
  hratiopt->Divide(hv2pt1s);
  hratioy->Divide(hv2y1s);
  hratioc->Divide(hv2c1s);

  //make histograms of systematics
  cout << "Getting systematics..." << endl;
  TString systFileName = "../Systematics/Ups_1_v2_cent10-90_SystCombined.root";
  TFile* systFile = new TFile(systFileName,"READ");
  TH1D* hsyspt1 = (TH1D*)systFile->Get("hv2pt");
  TH1D* hsysy1 = (TH1D*)systFile->Get("hv2y");
  TH1D* hsysc1 = (TH1D*)systFile->Get("hv2c");
  TString systFile2Name = "../Systematics/Ups_2_v2_cent10-90_SystCombined.root";
  TFile* systFile2 = new TFile(systFile2Name,"READ");
  TH1D* hsyspt2 = (TH1D*)systFile2->Get("hv2pt");
  TH1D* hsysy2 = (TH1D*)systFile2->Get("hv2y");
  TH1D* hsysc2 = (TH1D*)systFile2->Get("hv2c");

  //make TGraphs
  cout << "making pt TGraphs..." << endl;
  TGraphErrors* gv2pt = new TGraphErrors(hratiopt);
  TGraphErrors* gv2pt_sys = new TGraphErrors(hratiopt);
  for (int ipt=0; ipt<numptbins; ipt++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2pt->GetPoint(ipt, pxtmp, pytmp);
    extmp=gv2pt->GetErrorX(ipt);
    eytmp=gv2pt->GetErrorY(ipt);
    relsys=sqrt(pow(hsyspt1->GetBinContent(ipt+1),2) + pow(hsyspt2->GetBinContent(ipt+1),2)); //Uncertainty is found as the sum of squares of the two sources of uncertainty.
    //relsys=sqrt(2)*hsyspt->GetBinContent(ipt+1);//Error is multiplied by the square root of 2 because we are taking the ratio of two data points with the same relative error.
    //remove ex
    gv2pt->SetPointError(ipt, 0, eytmp);
    //set ey for systematic error
    gv2pt_sys->SetPointError(ipt, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2pt, styleParam, styleParam); 
  SetGraphStyleSys(gv2pt_sys, styleParam);

  float ptmin=0; float ptmax=50;
  gv2pt_sys->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  gv2pt_sys->GetXaxis()->CenterTitle();
  gv2pt_sys->GetXaxis()->SetTitleOffset(1.);
  gv2pt_sys->GetXaxis()->SetLimits(ptmin,ptmax);
  gv2pt_sys->GetXaxis()->SetTitleSize(0.06);
  gv2pt_sys->GetYaxis()->SetTitle("v_{2}^{#varUpsilon(2S)}/v_{2}^{#varUpsilon(1S)}");
  gv2pt_sys->GetYaxis()->CenterTitle();
  gv2pt_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2pt_sys->GetYaxis()->SetTitleSize(0.06);
  gv2pt_sys->SetMinimum(plotmin);
  gv2pt_sys->SetMaximum(plotmax);

  cout << "making y TGraphs..." << endl;
  TGraphErrors* gv2y = new TGraphErrors(hratioy);
  TGraphErrors* gv2y_sys = new TGraphErrors(hratioy);
  for (int iy=0; iy<numybins; iy++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2y->GetPoint(iy, pxtmp, pytmp);
    extmp=gv2y->GetErrorX(iy);
    eytmp=gv2y->GetErrorY(iy);
    relsys=sqrt(pow(hsysy1->GetBinContent(iy+1),2) + pow(hsysy2->GetBinContent(iy+1),2)); //Uncertainty is found as the sum of squares of the two sources of uncertainty.
    //remove ex
    gv2y->SetPointError(iy, 0, eytmp);
    //set ey for systematic error
    gv2y_sys->SetPointError(iy, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2y, styleParam, styleParam); 
  SetGraphStyleSys(gv2y_sys, styleParam);

  float ymin=0; float ymax=2.4;
  gv2y_sys->GetXaxis()->SetTitle("|y^{#varUpsilon}|");
  gv2y_sys->GetXaxis()->CenterTitle();
  gv2y_sys->GetXaxis()->SetTitleOffset(1.);
  gv2y_sys->GetXaxis()->SetLimits(ymin,ymax);
  gv2y_sys->GetXaxis()->SetTitleSize(0.06);
  gv2y_sys->GetYaxis()->SetTitle("v_{2}^{#varUpsilon(2S)}/v_{2}^{#varUpsilon(1S)}");
  gv2y_sys->GetYaxis()->CenterTitle();
  gv2y_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2y_sys->GetYaxis()->SetTitleSize(0.06);
  gv2y_sys->SetMinimum(plotmin);
  gv2y_sys->SetMaximum(plotmax);

  cout << "making c TGraphs..." << endl;
  TGraphErrors* gv2c = new TGraphErrors(hratioc);
  TGraphErrors* gv2c_sys = new TGraphErrors(hratioc);
  for (int ic=0; ic<numcbins; ic++) {
    double pxtmp=0; double pytmp=0; double extmp=0; double eytmp=0; double relsys=0;
    gv2c->GetPoint(ic, pxtmp, pytmp);
    extmp=gv2c->GetErrorX(ic);
    eytmp=gv2c->GetErrorY(ic);
    relsys=sqrt(pow(hsysc1->GetBinContent(ic+1),2) + pow(hsysc2->GetBinContent(ic+1),2)); //Uncertainty is found as the sum of squares of the two sources of uncertainty.
    //remove ex
    gv2c->SetPointError(ic, 0, eytmp);
    //set ey for systematic error
    gv2c_sys->SetPointError(ic, extmp, pytmp*relsys);
  }

  SetGraphStyle(gv2c, styleParam, styleParam); 
  SetGraphStyleSys(gv2c_sys, styleParam);

  float cmin=10; float cmax=90;
  gv2c_sys->GetXaxis()->SetTitle("Centrality (%)");
  gv2c_sys->GetXaxis()->CenterTitle();
  gv2c_sys->GetXaxis()->SetTitleOffset(1.);
  gv2c_sys->GetXaxis()->SetLimits(cmin,cmax);
  gv2c_sys->GetXaxis()->SetTitleSize(0.06);
  gv2c_sys->GetYaxis()->SetTitle("v_{2}^{#varUpsilon(2S)}/v_{2}^{#varUpsilon(1S)}");
  gv2c_sys->GetYaxis()->CenterTitle();
  gv2c_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2c_sys->GetYaxis()->SetTitleSize(0.06);
  gv2c_sys->SetMinimum(plotmin);
  gv2c_sys->SetMaximum(plotmax);

  TCanvas* cpt = new TCanvas("cpt","cpt",40,40,600,600);

  gv2pt_sys->Draw("A5");
  gv2pt->Draw("P");

  TString perc = "%";
  float pos_text_x = 0.25;
  float pos_text_y = 0.82;
  float pos_y_diff = 0.05;
  float text_size = 18;
  int text_color = 1;
  float muPtCut = 3.5;
  float yCut = 2.4;
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", muPtCut ), pos_text_x,pos_text_y,text_color,text_size);
  drawText(Form("|#eta^{#mu}| < %.2f", yCut ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  drawText(Form("|y^{#varUpsilon}| < %.2f", yCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("Centrality %i-%i%s", 10,90, perc.Data() ), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);

  TLegend *legpt= new TLegend(legx1,legy1,legx2,legy2);
  SetLegendStyle(legpt);
  legpt->AddEntry(gv2pt,"v_{2}^{#varUpsilon(2S)}/v_{2}^{#varUpsilon(1S)}","lp");

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

  cpt->SaveAs("Plots/ratio_vs_pt.png");
  cpt->SaveAs("Plots/ratio_vs_pt.pdf");

  TCanvas* cy = new TCanvas("cy","cy",40,40,600,600);

  gv2y_sys->Draw("A5");
  gv2y->Draw("P");

  TLegend *legy= new TLegend(legx1,legy1,legx2,legy2);
  SetLegendStyle(legy);
  legy->AddEntry(gv2y,"v_{2}^{#varUpsilon(2S)}/v_{2}^{#varUpsilon(1S)}","lp");

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

  cy->SaveAs("Plots/ratio_vs_y.png");
  cy->SaveAs("Plots/ratio_vs_y.pdf");

  TCanvas* cc = new TCanvas("cc","cc",40,40,600,600);

  gv2c_sys->Draw("A5");
  gv2c->Draw("P");

  TLegend *legc= new TLegend(legx1,legy1,legx2,legy2);
  SetLegendStyle(legc);
  legc->AddEntry(gv2c,"v_{2}^{#varUpsilon(2S)}/v_{2}^{#varUpsilon(1S)}","lp");

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

  cc->SaveAs("Plots/ratio_vs_c.png");
  cc->SaveAs("Plots/ratio_vs_c.pdf");

}
