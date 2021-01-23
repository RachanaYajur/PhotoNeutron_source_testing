{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  //TFile f("Sb124_BERT_BeO_1E9.root");      
  TFile f("gamma_EM_10MeV_5mm_1e9.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  //c1->SetFillColor(30);
  //c1.Setlogy();
  //TImageDump *imgdump = new TImageDump("test.png");
  //c1->Paint();
/*  TH1D* hist4 = (TH1D*)f.Get("4");
  hist4->Draw("HIST");
  hist4->Setlogy();
  hist4->GetXaxis()->SetRange(0,3);
  hist4->SetLineColor(2);
 // hist4->Draw("HIST");
  
  */
  //imgdump->Close();
  TH1D* hist6 = (TH1D*)f.Get("6");
  //hist6->Setlogy();
  //hist6->GetXaxis()->SetRangeUser(0,0.02);
  hist6->Draw("HIST");
  c1->SaveAs("gamma_EM_10MeV_5mm_1E9.jpg");
}
