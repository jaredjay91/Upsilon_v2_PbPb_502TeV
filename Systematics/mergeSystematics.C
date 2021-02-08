

void mergeSystematics(){

  TString systStr[4] = {"altSig","altBkg","altAcc","altEff"};

  TFile* altFile[4];
  TH1D* hSystpt[4];
  TH1D* hSysty[4];
  TH1D* hSystc[4];

  for (int i=0; i<4; i++) {
    TString altFileName = Form("Ups_1_v2_cent10-90_%sSyst.root", systStr[i].Data());
    altFile[i] = TFile::Open(altFileName,"READ");
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
    for (int jSyst=0; jSyst<4; jSyst++) {
      systval = hSystpt[jSyst]->GetBinContent(i);
      sumsqrs = sumsqrs + pow(systval,2);
    }
    sumsqrs = sqrt(sumsqrs);
    hSystptCombined->SetBinContent(i,sumsqrs);
  }

  // Y BINS
  TH1D* hSystyCombined = (TH1D*)hSysty[0]->Clone();
  for (int i=0; i<nybins; i++) {
    sumsqrs = 0.0;
    for (int jSyst=0; jSyst<4; jSyst++) {
      systval = hSysty[jSyst]->GetBinContent(i);
      sumsqrs = sumsqrs + pow(systval,2);
    }
    sumsqrs = sqrt(sumsqrs);
    hSystyCombined->SetBinContent(i,sumsqrs);
  }

  // CENTRALITY BINS
  TH1D* hSystcCombined = (TH1D*)hSystc[0]->Clone();
  for (int i=0; i<ncbins; i++) {
    sumsqrs = 0.0;
    for (int jSyst=0; jSyst<4; jSyst++) {
      systval = hSystc[jSyst]->GetBinContent(i);
      sumsqrs = sumsqrs + pow(systval,2);
    }
    sumsqrs = sqrt(sumsqrs);
    hSystcCombined->SetBinContent(i,sumsqrs);
  }

  TString outFileName = Form("Ups_1_v2_cent10-90_SystCombined.root");
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hSystptCombined->Write();
  hSystyCombined->Write();
  hSystcCombined->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;
}
