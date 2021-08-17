
void drawResCor() {

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  float textSize = 0.055;

  TFile* inFile = TFile::Open("averages/resCorFile_n56114317_BinByBin.root","READ");
  TH1D* hRpt = (TH1D*)inFile->Get("hRpt");
  TH1D* hRy = (TH1D*)inFile->Get("hRy");
  TH1D* hRc = (TH1D*)inFile->Get("hRc");

  hRpt->SetTitle("");
  hRpt->SetTitleSize(textSize);
  hRpt->GetXaxis()->SetTitle("p_{T} (GeV)");
  hRpt->GetXaxis()->SetTitleSize(textSize);
  hRpt->GetXaxis()->SetLabelSize(textSize);
  hRpt->GetYaxis()->SetTitle("EP Resolution Correction");
  hRpt->GetYaxis()->SetTitleSize(textSize);
  hRpt->GetYaxis()->SetLabelSize(textSize);
  hRpt->SetMarkerStyle(8);
  hRpt->SetMarkerSize(1.5);
  hRpt->SetMarkerColor(kBlue);
  hRpt->SetLineColor(kBlue);
  hRpt->SetMinimum(0.35);
  hRpt->SetMaximum(1.0);
  hRpt->Sumw2();

  hRy->SetTitle("");
  hRy->SetTitleSize(textSize);
  hRy->GetXaxis()->SetTitle("y");
  hRy->GetXaxis()->SetTitleSize(textSize);
  hRy->GetXaxis()->SetLabelSize(textSize);
  hRy->GetYaxis()->SetTitle("EP Resolution Correction");
  hRy->GetYaxis()->SetTitleSize(textSize);
  hRy->GetYaxis()->SetLabelSize(textSize);
  hRy->SetMarkerStyle(8);
  hRy->SetMarkerSize(1.5);
  hRy->SetMarkerColor(kBlue);
  hRy->SetLineColor(kBlue);
  hRy->SetMinimum(0.35);
  hRy->SetMaximum(1.0);
  hRy->Sumw2();

  hRc->SetTitle("");
  hRc->SetTitleSize(textSize);
  hRc->GetXaxis()->SetTitle("centrality (%)");
  hRc->GetXaxis()->SetTitleSize(textSize);
  hRc->GetXaxis()->SetLabelSize(textSize);
  hRc->GetYaxis()->SetTitle("EP Resolution Correction");
  hRc->GetYaxis()->SetTitleSize(textSize);
  hRc->GetYaxis()->SetLabelSize(textSize);
  hRc->SetMarkerStyle(8);
  hRc->SetMarkerSize(1.5);
  hRc->SetMarkerColor(kBlue);
  hRc->SetLineColor(kBlue);
  hRc->SetMinimum(0.35);
  hRc->SetMaximum(1.0);
  hRc->Sumw2();



  TCanvas* cpt = new TCanvas("cpt","cpt",0,0,500,500);
  hRpt->Draw();
  TLegend* lpt = new TLegend(0.2,0.2,0.6,0.3);
  lpt->SetBorderSize(0);
  lpt->SetTextSize(textSize);
  lpt->AddEntry(hRpt,"Resolution correction","pel");
  lpt->AddEntry(hRpt,"factor","");
  //lpt->Draw("same");
  cpt->SetBottomMargin(0.15);
  cpt->SetLeftMargin(0.15);

  TCanvas* cy = new TCanvas("cy","cy",0,0,500,500);
  hRy->Draw();
  TLegend* ly = new TLegend(0.2,0.2,0.6,0.3);
  ly->SetBorderSize(0);
  ly->SetTextSize(textSize);
  ly->AddEntry(hRy,"Resolution correction","pel");
  ly->AddEntry(hRy,"factor","");
  //ly->Draw("same");
  cy->SetBottomMargin(0.15);
  cy->SetLeftMargin(0.15);

  TCanvas* cc = new TCanvas("cc","cc",0,0,500,500);
  hRc->Draw();
  TLegend* lc = new TLegend(0.2,0.2,0.6,0.3);
  lc->SetBorderSize(0);
  lc->SetTextSize(textSize);
  lc->AddEntry(hRc,"Resolution correction","pel");
  lc->AddEntry(hRc,"factor","");
  //lc->Draw("same");
  cc->SetBottomMargin(0.15);
  cc->SetLeftMargin(0.15);

  cpt->SaveAs("resCorpt.pdf");
  cy->SaveAs("resCory.pdf");
  cc->SaveAs("resCorc.pdf");
}
