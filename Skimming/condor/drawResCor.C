
void drawResCor() {

  gStyle->SetOptStat(0);

  TFile* inFile = TFile::Open("averages/resCorFile_n56114317_BinByBin.root","READ");
  TH1D* hRpt = (TH1D*)inFile->Get("hRpt");
  TH1D* hRy = (TH1D*)inFile->Get("hRy");
  TH1D* hRc = (TH1D*)inFile->Get("hRc");

  hRpt->GetXaxis()->SetTitle("p_T");
  hRpt->SetMinimum(0.4);
  hRpt->SetMaximum(1.0);
  //hRpt->GetXaxis()->SetTitleSize(0.06);
  hRpt->GetYaxis()->SetTitleSize(0.06);

  hRy->GetXaxis()->SetTitle("y");
  hRy->SetMinimum(0.4);
  hRy->SetMaximum(1.0);
  //hRpt->GetXaxis()->SetTitleSize(0.06);
  hRy->GetYaxis()->SetTitleSize(0.06);

  hRc->GetXaxis()->SetTitle("centrality");
  hRc->SetMinimum(0.4);
  hRc->SetMaximum(1.0);
  //hRc->GetXaxis()->SetTitleSize(0.06);
  hRc->GetYaxis()->SetTitleSize(0.06);

  hRpt->Sumw2();
  hRc->Sumw2();

  TCanvas* cpt = new TCanvas("cpt","cpt",0,0,500,500);
  hRpt->Draw();

  TCanvas* cy = new TCanvas("cy","cy",0,0,500,500);
  hRy->Draw();

  TCanvas* cc = new TCanvas("cc","cc",0,0,500,500);
  hRc->Draw();

  cpt->SaveAs("resCorpt.pdf");
  cy->SaveAs("resCory.pdf");
  cc->SaveAs("resCorc.pdf");
}
