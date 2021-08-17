#include <ctime>

#include "../../HeaderFiles/commonUtility.h"
#include "../../HeaderFiles/HiEvtPlaneList.h"
#include "../../HeaderFiles/cutsAndBinUpsilonV2.h"

const Double_t pi = 3.141592653589;


void PlotEfficiency(int dateStr=20210811){

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  float textSize = 0.055;

  TFile* effFile = TFile::Open(Form("efficiency_%i.root",dateStr),"READ");
  TH1D* hEffInt = (TH1D*)effFile->Get("hEffInt");
  TH1D* hEffPt = (TH1D*)effFile->Get("hEffPt");
  TH1D* hEffPtLowC = (TH1D*)effFile->Get("hEffPtLowC");
  TH1D* hEffPtMidC = (TH1D*)effFile->Get("hEffPtMidC");
  TH1D* hEffPtHighC = (TH1D*)effFile->Get("hEffPtHighC");
  TH1D* hEffIntNoW = (TH1D*)effFile->Get("hEffIntNoW");
  TH1D* hEffPtNoW = (TH1D*)effFile->Get("hEffPtNoW");
  TH1D* hEffPtLowCNoW = (TH1D*)effFile->Get("hEffPtLowCNoW");
  TH1D* hEffPtMidCNoW = (TH1D*)effFile->Get("hEffPtMidCNoW");
  TH1D* hEffPtHighCNoW = (TH1D*)effFile->Get("hEffPtHighCNoW");

  TCanvas* c1 = new TCanvas("c1","c1",0,0,400,400);
  hEffPtLowC->GetXaxis()->SetTitleSize(textSize);
  hEffPtLowC->GetXaxis()->SetLabelSize(textSize);
  hEffPtLowC->GetYaxis()->SetTitleSize(textSize);
  hEffPtLowC->GetYaxis()->SetLabelSize(textSize);
  hEffPtLowC->Draw("PE");
  hEffPtLowCNoW->SetMarkerStyle(24);
  hEffPtLowCNoW->Draw("same PE");

  TCanvas* c2 = new TCanvas("c2","c2",0,0,400,400);
  hEffPtMidC->GetXaxis()->SetTitleSize(textSize);
  hEffPtMidC->GetXaxis()->SetLabelSize(textSize);
  hEffPtMidC->GetYaxis()->SetTitleSize(textSize);
  hEffPtMidC->GetYaxis()->SetLabelSize(textSize);
  hEffPtMidC->Draw("PE");
  hEffPtMidCNoW->SetMarkerStyle(24);
  hEffPtMidCNoW->Draw("same PE");

  TCanvas* c3 = new TCanvas("c3","c3",0,0,400,400);
  hEffPtHighC->GetXaxis()->SetTitleSize(textSize);
  hEffPtHighC->GetXaxis()->SetLabelSize(textSize);
  hEffPtHighC->GetYaxis()->SetTitleSize(textSize);
  hEffPtHighC->GetYaxis()->SetLabelSize(textSize);
  hEffPtHighC->Draw("PE");
  hEffPtHighCNoW->SetMarkerStyle(24);
  hEffPtHighCNoW->Draw("same PE");

  float pos_text_x = 0.18;
  float pos_text_y = 0.18;
  float pos_y_diff = 0.06;
  float text_size = 20;
  int text_color = 1;

  c1->cd();
  TLegend* leg1 = new TLegend(0.5,0.3,0.89,0.4);
  leg1->SetTextSize(textSize);
  //leg1->SetTextFont(43);
  leg1->SetBorderSize(0);
  leg1->AddEntry(hEffPtLowCNoW,"unweighted","pe");
  leg1->AddEntry(hEffPtLowC,"weighted","pe");
  leg1->Draw("same");
  drawText("Centrality 10-30\%", pos_text_x,pos_text_y,text_color,text_size);
  c1->SetBottomMargin(0.13);
  c1->SetLeftMargin(0.15);

  c2->cd();
  TLegend* leg2 = new TLegend(0.5,0.3,0.89,0.4);
  leg2->SetTextSize(textSize);
  //leg2->SetTextFont(43);
  leg2->SetBorderSize(0);
  leg2->AddEntry(hEffPtMidCNoW,"unweighted","pe");
  leg2->AddEntry(hEffPtMidC,"weighted","pe");
  leg2->Draw("same");
  drawText("Centrality 30-50\%", pos_text_x,pos_text_y,text_color,text_size);
  c2->SetBottomMargin(0.13);
  c2->SetLeftMargin(0.15);

  c3->cd();
  TLegend* leg3 = new TLegend(0.5,0.3,0.89,0.4);
  leg3->SetTextSize(textSize);
  //leg3->SetTextFont(43);
  leg3->SetBorderSize(0);
  leg3->AddEntry(hEffPtHighCNoW,"unweighted","pe");
  leg3->AddEntry(hEffPtHighC,"weighted","pe");
  leg3->Draw("same");
  drawText("Centrality 50-90\%", pos_text_x,pos_text_y,text_color,text_size);
  c3->SetBottomMargin(0.13);
  c3->SetLeftMargin(0.15);

  c1->SaveAs(Form("EfficiencyPlot_LowC_%i.png",dateStr));
  c1->SaveAs(Form("EfficiencyPlot_LowC_%i.pdf",dateStr));

  c2->SaveAs(Form("EfficiencyPlot_MidC_%i.png",dateStr));
  c2->SaveAs(Form("EfficiencyPlot_MidC_%i.pdf",dateStr));

  c3->SaveAs(Form("EfficiencyPlot_HighC_%i.png",dateStr));
  c3->SaveAs(Form("EfficiencyPlot_HighC_%i.pdf",dateStr));

}

