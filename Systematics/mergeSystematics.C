

void mergeSystematics(int whichUpsilon=1){

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  double etaGapSensitivity = 0.025; //From https://arxiv.org/pdf/1204.1850.pdf
  double evtPlaneRes = 0.01;
  double evtPlaneUnc = evtPlaneRes+etaGapSensitivity;

  double singleMuTrkUnc = 0.006; //For PbPb (from Jpsi analysis)  
  double doubleMuTrkUnc = 2*singleMuTrkUnc;

  TString systStr[5] = {"Signal function","Background function","Acceptance weighting","Efficiency weighting","Parameter constraints"};

  TString systFileName[5] = {
	Form("PseudoExperimentsCode/Ups_%i_v2_cent10-90_altSigSyst.root", whichUpsilon),
	Form("PseudoExperimentsCode/Ups_%i_v2_cent10-90_altBkgSyst.root", whichUpsilon),
	Form("Ups_%i_v2_cent10-90_altAccSyst.root", whichUpsilon),
	Form("Ups_%i_v2_cent10-90_altEffSyst.root", whichUpsilon),
	Form("PseudoExperimentsCode/Ups_%i_v2_cent10-90_altConstSyst.root", whichUpsilon)};

  int color[5] = {2, 3, 4, 6, 7};

  TFile* altFile[5];
  TH1D* hSystpt[5];
  TH1D* hSysty[5];
  TH1D* hSystc[5];

  for (int i=0; i<5; i++) {
    altFile[i] = TFile::Open(systFileName[i],"READ");
    hSystpt[i] = (TH1D*)altFile[i]->Get("hv2pt");
    hSysty[i] = (TH1D*)altFile[i]->Get("hv2y");
    hSystc[i] = (TH1D*)altFile[i]->Get("hv2c");
  }

  //cout << "pow(3,2) = " << pow(3,2) << endl;

  //Combine systematics via sum of squares
  int nptbins = hSystpt[0]->GetNbinsX();
  int nybins = hSysty[0]->GetNbinsX();
  int ncbins = hSystc[0]->GetNbinsX();

  double sumsqrs;
  double systval;

  // PT BINS
  TH1D* hSystptCombined = (TH1D*)hSystpt[0]->Clone();
  for (int i=0; i<nptbins; i++) {
    sumsqrs = 0.0;
    for (int jSyst=0; jSyst<5; jSyst++) {
      systval = hSystpt[jSyst]->GetBinContent(i+1);
      sumsqrs = sumsqrs + pow(systval,2);
    }
    sumsqrs = sumsqrs + pow(evtPlaneUnc,2) + pow(doubleMuTrkUnc,2);
    sumsqrs = sqrt(sumsqrs);
    hSystptCombined->SetBinContent(i+1,sumsqrs);
  }

  // Y BINS
  TH1D* hSystyCombined = (TH1D*)hSysty[0]->Clone();
  for (int i=0; i<nybins; i++) {
    sumsqrs = 0.0;
    for (int jSyst=0; jSyst<5; jSyst++) {
      systval = hSysty[jSyst]->GetBinContent(i+1);
      sumsqrs = sumsqrs + pow(systval,2);
    }
    sumsqrs = sumsqrs + pow(evtPlaneUnc,2) + pow(doubleMuTrkUnc,2);
    sumsqrs = sqrt(sumsqrs);
    hSystyCombined->SetBinContent(i+1,sumsqrs);
  }

  // CENTRALITY BINS
  TH1D* hSystcCombined = (TH1D*)hSystc[0]->Clone();
  for (int i=0; i<ncbins; i++) {
    sumsqrs = 0.0;
    for (int jSyst=0; jSyst<5; jSyst++) {
      systval = hSystc[jSyst]->GetBinContent(i+1);
      sumsqrs = sumsqrs + pow(systval,2);
    }
    sumsqrs = sumsqrs + pow(evtPlaneUnc,2) + pow(doubleMuTrkUnc,2);
    sumsqrs = sqrt(sumsqrs);
    hSystcCombined->SetBinContent(i+1,sumsqrs);
  }

  float labelSize = 0.05;
  float legx1 = 0.35;
  float legx2 = 0.89;
  float legy1 = 0.55;
  float legy2 = 0.89;
  TCanvas* cpt = new TCanvas("cpt","cpt",0,0,500,500);
  cpt->cd();
  hSystptCombined->SetTitle("");
  hSystptCombined->GetXaxis()->SetTitle("p_{T}^{#varUpsilon}");
  hSystptCombined->GetXaxis()->SetTitleSize(labelSize);
  hSystptCombined->GetXaxis()->SetLabelSize(labelSize);
  hSystptCombined->GetYaxis()->SetTitle("Relative Systematic Uncertainty in v_{2}^{#varUpsilon}");
  hSystptCombined->GetYaxis()->SetTitleSize(labelSize);
  hSystptCombined->GetYaxis()->SetLabelSize(labelSize);
  hSystptCombined->SetMinimum(0);
  hSystptCombined->SetMaximum(2.7);
  hSystptCombined->SetLineColor(1);
  hSystptCombined->SetLineWidth(2);
  hSystptCombined->Draw();
  TLegend* legpt = new TLegend(legx1,legy1,legx2,legy2);
  legpt->SetTextSize(20);
  legpt->SetTextFont(43);
  legpt->SetBorderSize(0);
  legpt->AddEntry(hSystptCombined,"Total","l");
  for (int i=0; i<5; i++) {
    hSystpt[i]->Draw("hist same");
    hSystpt[i]->SetLineColor(color[i]);
    hSystpt[i]->SetLineWidth(2);
    legpt->AddEntry(hSystpt[i],systStr[i].Data(),"l");
  }
  legpt->Draw("same");
  hSystptCombined->Draw("same");
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->RedrawAxis();


  TCanvas* cy = new TCanvas("cy","cy",400,0,500,500);
  cy->cd();
  hSystyCombined->SetTitle("");
  hSystyCombined->GetXaxis()->SetTitle("y^{#varUpsilon}");
  hSystyCombined->GetXaxis()->SetTitleSize(labelSize);
  hSystyCombined->GetXaxis()->SetLabelSize(labelSize);
  hSystyCombined->GetYaxis()->SetTitle("Relative Systematic Uncertainty in v_{2}^{#varUpsilon}");
  hSystyCombined->GetYaxis()->SetTitleSize(labelSize);
  hSystyCombined->GetYaxis()->SetLabelSize(labelSize);
  hSystyCombined->SetMinimum(0);
  hSystyCombined->SetMaximum(2.7);
  hSystyCombined->SetLineColor(1);
  hSystyCombined->SetLineWidth(2);
  hSystyCombined->Draw();
  TLegend* legy = new TLegend(legx1,legy1,legx2,legy2);
  legy->SetTextSize(20);
  legy->SetTextFont(43);
  legy->SetBorderSize(0);
  legy->AddEntry(hSystyCombined,"Total","l");
  for (int i=0; i<5; i++) {
    hSysty[i]->Draw("hist same");
    hSysty[i]->SetLineColor(color[i]);
    hSysty[i]->SetLineWidth(2);
    legy->AddEntry(hSysty[i],systStr[i].Data(),"l");
  }
  legy->Draw("same");
  hSystyCombined->Draw("same");
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->RedrawAxis();

  TCanvas* cc = new TCanvas("cc","cc",800,0,500,500);
  cc->cd();
  hSystcCombined->SetTitle("");
  hSystcCombined->GetXaxis()->SetTitle("Centrality (%)");
  hSystcCombined->GetXaxis()->SetTitleSize(labelSize);
  hSystcCombined->GetXaxis()->SetLabelSize(labelSize);
  hSystcCombined->GetYaxis()->SetTitle("Relative Systematic Uncertainty in v_{2}^{#varUpsilon}");
  hSystcCombined->GetYaxis()->SetTitleSize(labelSize);
  hSystcCombined->GetYaxis()->SetLabelSize(labelSize);
  hSystcCombined->SetMinimum(0);
  hSystcCombined->SetMaximum(2.7);
  hSystcCombined->SetLineColor(1);
  hSystcCombined->SetLineWidth(2);
  hSystcCombined->Draw();
  TLegend* legc = new TLegend(legx1,legy1,legx2,legy2);
  legc->SetTextSize(20);
  legc->SetTextFont(43);
  legc->SetBorderSize(0);
  legc->AddEntry(hSystcCombined,"Total","l");
  for (int i=0; i<5; i++) {
    hSystc[i]->Draw("hist same");
    hSystc[i]->SetLineColor(color[i]);
    hSystc[i]->SetLineWidth(2);
    legc->AddEntry(hSystc[i],systStr[i].Data(),"l");
  }
  legc->Draw("same");
  hSystcCombined->Draw("same");
  hSystcCombined->GetYaxis()->SetTitleOffset(1.4);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  gPad->RedrawAxis();

  cpt->SaveAs(Form("Plots/Systpt_%is.pdf", whichUpsilon));
  cpt->SaveAs(Form("Plots/Systpt_%is.png", whichUpsilon));
  cy->SaveAs(Form("Plots/Systy_%is.pdf", whichUpsilon));
  cy->SaveAs(Form("Plots/Systy_%is.png", whichUpsilon));
  cc->SaveAs(Form("Plots/Systc_%is.pdf", whichUpsilon));
  cc->SaveAs(Form("Plots/Systc_%is.png", whichUpsilon));


  TString outFileName = Form("Ups_%i_v2_cent10-90_SystCombined.root", whichUpsilon);
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hSystptCombined->Write();
  hSystyCombined->Write();
  hSystcCombined->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;
}
