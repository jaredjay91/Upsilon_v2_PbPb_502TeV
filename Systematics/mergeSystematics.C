

void mergeSystematics(){

  double etaGapSensitivity = 0.025; //From https://arxiv.org/pdf/1204.1850.pdf
  double evtPlaneRes = 0.01;
  double evtPlaneUnc = evtPlaneRes+etaGapSensitivity;

  double singleMuTrkUnc = 0.006; //For PbPb (from Jpsi analysis)  
  double doubleMuTrkUnc = 2*singleMuTrkUnc;

  TString systStr[5] = {"altSig","altBkg","altAcc","altEff","altConst"};

  TString systFileName[5] = {
	"PseudoExperimentsCode/Ups_1_v2_cent10-90_altSigSyst.root",
	"PseudoExperimentsCode/Ups_1_v2_cent10-90_altBkgSyst.root",
	"Ups_1_v2_cent10-90_altAccSyst.root",
	"Ups_1_v2_cent10-90_altEffSyst.root",
	"Ups_1_v2_cent10-90_altConstSyst.root"};

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

  TCanvas* cpt = new TCanvas("cpt","cpt",0,0,500,500);
  cpt->cd();
  hSystptCombined->SetTitle("Systematics vs p_{T}");
  hSystptCombined->GetXaxis()->SetTitle("p_{T}^{#Upsilon}");
  hSystptCombined->SetMinimum(0);
  hSystptCombined->SetLineColor(1);
  hSystptCombined->Draw();
  TLegend* legpt = new TLegend(0.7,0.6,0.89,0.8);
  legpt->SetTextSize(19);
  legpt->SetTextFont(43);
  legpt->SetBorderSize(0);
  legpt->AddEntry(hSystptCombined,"Total","l");
  for (int i=0; i<5; i++) {
    hSystpt[i]->Draw("hist same");
    hSystpt[i]->SetLineColor(color[i]);
    legpt->AddEntry(hSystpt[i],systStr[i].Data(),"l");
  }
  legpt->Draw("same");
  hSystptCombined->Draw("same");
  gPad->RedrawAxis();


  TCanvas* cy = new TCanvas("cy","cy",400,0,500,500);
  cy->cd();
  hSystyCombined->SetTitle("Systematics vs Rapidity");
  hSystyCombined->GetXaxis()->SetTitle("y^{#Upsilon}");
  hSystyCombined->SetMinimum(0);
  hSystyCombined->SetLineColor(1);
  hSystyCombined->Draw();
  TLegend* legy = new TLegend(0.7,0.6,0.89,0.8);
  legy->SetTextSize(19);
  legy->SetTextFont(43);
  legy->SetBorderSize(0);
  legy->AddEntry(hSystyCombined,"Total","l");
  for (int i=0; i<5; i++) {
    hSysty[i]->Draw("hist same");
    hSysty[i]->SetLineColor(color[i]);
    legy->AddEntry(hSysty[i],systStr[i].Data(),"l");
  }
  legy->Draw("same");
  hSystyCombined->Draw("same");
  gPad->RedrawAxis();

  TCanvas* cc = new TCanvas("cc","cc",800,0,500,500);
  cc->cd();
  hSystcCombined->SetTitle("Systematics vs Centrality");
  hSystcCombined->GetXaxis()->SetTitle("Centrality (%)");
  hSystcCombined->SetMinimum(0);
  hSystcCombined->SetLineColor(1);
  hSystcCombined->Draw();
  TLegend* legc = new TLegend(0.7,0.6,0.89,0.8);
  legc->SetTextSize(19);
  legc->SetTextFont(43);
  legc->SetBorderSize(0);
  legc->AddEntry(hSystcCombined,"Total","l");
  for (int i=0; i<5; i++) {
    hSystc[i]->Draw("hist same");
    hSystc[i]->SetLineColor(color[i]);
    legc->AddEntry(hSystc[i],systStr[i].Data(),"l");
  }
  legc->Draw("same");
  hSystcCombined->Draw("same");
  gPad->RedrawAxis();

  cpt->SaveAs("Plots/Systpt.pdf");
  cpt->SaveAs("Plots/Systpt.png");
  cy->SaveAs("Plots/Systy.pdf");
  cy->SaveAs("Plots/Systy.png");
  cc->SaveAs("Plots/Systc.pdf");
  cc->SaveAs("Plots/Systc.png");


  TString outFileName = Form("Ups_1_v2_cent10-90_SystCombined.root");
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hSystptCombined->Write();
  hSystyCombined->Write();
  hSystcCombined->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;
}
