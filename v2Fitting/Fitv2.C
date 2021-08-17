#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
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

void Fitv2(

  int collId = kAADATA,
  float ptLow = 0, float ptHigh = 50,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 10, int cHigh = 90,
  float muPtCut = 3.5,
  double* v2ValPtr=0, double* v2ErrPtr=0,
  int whichUpsilon=1, int whichSyst=0 )

{

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TString fileLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, 0.0, 0.5);

  //Load the yields vs phi
  TFile* yieldsFile = TFile::Open(Form("ExtractedYields/%syieldsVsPhi_%iS_%s.root", systStr.Data(), whichUpsilon, fileLabel.Data()), "READ");
  TH1D* yieldsVsPhi = (TH1D*)yieldsFile->Get("yieldsVsPhi;1");

  //Print some information about the yields
  float tot_yield = 0.0;
  float tot_yielderr = 0.0;
  for (int i=1; i<5; i++) {
    float yield = yieldsVsPhi->GetBinContent(i);
    float yielderr = yieldsVsPhi->GetBinError(i);
    cout << "Yield = " << yield << " +/- " << yielderr << endl;
    tot_yield += yield;
    tot_yielderr += pow(yielderr,2);
  }
  cout << "Total yield = " << tot_yield << " $\\pm$ " << sqrt(tot_yielderr) << endl;

  //Plot it
  TCanvas* c1 = new TCanvas("c1","c1",50,50,550,550);
  yieldsVsPhi->Draw();
  yieldsVsPhi->SetTitle("");
  yieldsVsPhi->GetXaxis()->SetTitle("|#Delta#phi|/#pi");
  yieldsVsPhi->GetXaxis()->SetTitleSize(0.045);
  yieldsVsPhi->GetXaxis()->SetLabelSize(0.045);
  yieldsVsPhi->GetYaxis()->SetTitle(Form("Normalized yield of #Upsilon(%iS)",whichUpsilon));
  yieldsVsPhi->GetYaxis()->SetTitleOffset(1.5);
  yieldsVsPhi->GetYaxis()->SetTitleSize(0.045);
  yieldsVsPhi->GetYaxis()->SetLabelSize(0.045);
  yieldsVsPhi->SetMarkerStyle(8);
  yieldsVsPhi->SetMarkerColor(1);
  yieldsVsPhi->SetLineColor(1);
  double totalintegral = yieldsVsPhi->Integral(1,4);
  yieldsVsPhi->Scale(1.0/totalintegral);
  yieldsVsPhi->SetMinimum(0);
  yieldsVsPhi->SetMaximum(0.5);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);

  //Fit it
  /*TF1* fitfunc = new TF1("fitfunc","[0]*( 1 + 2*[1]*cos(2*x*3.14159265) + 2*[2]*cos(3*x*3.14159265) + 2*[3]*cos(4*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("Amp","v2","v3","v4");*/
  /*TF1* fitfunc = new TF1("fitfunc","[0]*( 1 + 2*[1]*cos(2*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("Amp","v2");*/
  TF1* fitfunc = new TF1("fitfunc"," 0.25*(1 + 2*[0]*cos(2*x*3.14159265))",0,0.5);
  fitfunc->SetParNames("v2");
  yieldsVsPhi->Fit("fitfunc");

  TLegend* leg1 = new TLegend(0.55,0.7,0.87,0.85);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->AddEntry(yieldsVsPhi,"Normalized yield","pel");
  leg1->AddEntry(fitfunc,"Fit","l");
  leg1->Draw("same");

  double v2Val = fitfunc->GetParameter(0);
  double v2Err = fitfunc->GetParError(0);
  cout << "v2 = " << v2Val << " +/- " << v2Err << endl;
  TLatex latex;
  latex.SetTextSize(0.05);
  //latex.SetTextAlign(13);
  latex.DrawLatex(0.25,50,Form("v_{2} = %.3f #pm %.3f",v2Val,v2Err));

  TString perc = "%";
  float pos_text_x = 0.2;
  float pos_text_y = 0.8;
  float pos_y_diff = 0.06;
  float text_size = 20;
  int text_color = 1;
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", muPtCut ), pos_text_x,pos_text_y,text_color,text_size);
  drawText(Form("|#eta^{#mu}| < %.1f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  pos_text_y = 0.27;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else drawText(Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  drawText(Form("Centrality %i-%i%s", cLow,cHigh, perc.Data() ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  pos_text_x = 0.6;
  drawText(Form("v_{2} = %.3f #pm %.3f", v2Val,v2Err), pos_text_x,pos_text_y-pos_y_diff*2,2,text_size);

  c1->Update();

  c1->SaveAs(Form("fitPlots/%s%is_v2_fit_pt%.1f-%.1f_y%.2f-%.2f_cent%i-%i.png", systStr.Data(), whichUpsilon, ptLow,ptHigh, yLow,yHigh, cLow, cHigh));

  *v2ValPtr = v2Val;
  *v2ErrPtr = v2Err;
  
  yieldsFile->Close();
  delete fitfunc;
  delete c1;
}
