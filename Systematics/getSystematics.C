

void subtractHists(TH1D* hNom, TH1D* hAlt, TH1D* hSyst) {

  int nbinsNom = hNom->GetNbinsX();
  int nbinsAlt = hAlt->GetNbinsX();
  int nbinsSyst = hSyst->GetNbinsX();
  if (!(nbinsAlt==nbinsNom && nbinsSyst==nbinsNom)) {
    cout << "Error: histograms have different number of bins!" << endl;
    return;
  }

  for (int i=1; i<=nbinsNom; i++) {
    double valNom = hNom->GetBinContent(i);
    double valAlt = hAlt->GetBinContent(i);
    double valSyst = fabs( (valAlt-valNom)/valNom );
    hSyst->SetBinContent(i, valSyst);
  }

}

void getSystematics(int whichUpsilon=1, int whichSyst=4) {

  if (! (whichSyst==1 || whichSyst==2 || whichSyst==3 || whichSyst==4 || whichSyst==5) ) {
    cout << "Error: Invalid value of 'whichSyst'" << endl;
    return;
  }

  TString systStr;
  if (whichSyst==0) systStr = "nom";
  else if (whichSyst==1) systStr = "altSig";
  else if (whichSyst==2) systStr = "altBkg";
  else if (whichSyst==3) systStr = "altAcc";
  else if (whichSyst==4) systStr = "altEff";
  else if (whichSyst==5) systStr = "altConst";

  TString nomFileName = Form("../v2Fitting/Ups_%i_v2_cent10-90_nom.root", whichUpsilon);
  TFile* nomFile = TFile::Open(nomFileName,"READ");
  TH1D* hNompt = (TH1D*)nomFile->Get("hv2pt");
  TH1D* hNomy = (TH1D*)nomFile->Get("hv2y");
  TH1D* hNomc = (TH1D*)nomFile->Get("hv2c");

  TString altFileName = Form("../v2Fitting/Ups_%i_v2_cent10-90_%s.root", whichUpsilon, systStr.Data());
  TFile* altFile = TFile::Open(altFileName,"READ");
  TH1D* hAltpt = (TH1D*)altFile->Get("hv2pt");
  TH1D* hAlty = (TH1D*)altFile->Get("hv2y");
  TH1D* hAltc = (TH1D*)altFile->Get("hv2c");

  TH1D* hSystpt = (TH1D*)hNompt->Clone();
  TH1D* hSysty = (TH1D*)hNomy->Clone();
  TH1D* hSystc = (TH1D*)hNomc->Clone();

  subtractHists(hNompt,hAltpt,hSystpt);
  subtractHists(hNomy,hAlty,hSysty);
  subtractHists(hNomc,hAltc,hSystc);

  TCanvas* cSyst = new TCanvas("cSyst","cSyst",1200,400);
  cSyst->Divide(3,1);
  cSyst->cd(1);
  hSystpt->Draw();
  cSyst->cd(2);
  hSysty->Draw();
  cSyst->cd(3);
  hSystc->Draw();

  cSyst->SaveAs(Form("Plots/Ups_%i_v2_cent10-90_%sSyst.png", whichUpsilon, systStr.Data()));

  TString outFileName = Form("Ups_%i_v2_cent10-90_%sSyst.root", whichUpsilon, systStr.Data());
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hSystpt->Write();
  hSysty->Write();
  hSystc->Write();
  outFile->Close();
  cout << "File created: " << outFileName << endl;

}
