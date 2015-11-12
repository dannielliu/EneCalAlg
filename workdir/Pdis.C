{
  TTree *treekpi;
  TTree *treekpipi;
  TTree *treepsip;
  TTree *treekk;
  TTree *tree6pi;
  TTree *treekkpipi;

  TFile *file1 = new TFile("plot_kpi_D0D0bar.root");
  file1->GetObject("vars",treekpi);
  TFile *file2 = new TFile("plot_kpipi_DmDp.root");
  file2->GetObject("vars",treekpipi);
  TFile *file3 = new TFile("plot_pipill.root");
  file3->GetObject("vars",treepsip);
  TFile *file4 = new TFile("plot_kk.root");
  file4->GetObject("vars",treekk);
  TFile *file5 = new TFile("plot_6pi.root");
  file5->GetObject("vars",tree6pi);
  TFile *file6 = new TFile("plot_kkpipi.root");
  file6->GetObject("vars",treekkpipi);

  TH1D *hppi_kpi   = new TH1D("hppi_kpi"   ,"momentum distribution",200,0,2);
  hppi_kpi->SetLineColor(1);
  TH1D *hppi_kpipi = new TH1D("hppi_kpipi" ,"momentum distribution",200,0,2);
  hppi_kpipi->SetLineColor(2);
  TH1D *hppi_pipill= new TH1D("hppi_pipill","momentum distribution",200,0,2);
  hppi_pipill->SetLineColor(3);
  TH1D *hppi_6pi   = new TH1D("hppi_6pi"   ,"momentum distribution",200,0,2);
  hppi_6pi->SetLineColor(4);
  TH1D *hppi_kkpipi= new TH1D("hppi_kkpipi","momentum distribution",200,0,2);
  hppi_kkpipi->SetLineColor(5);

  TH1D *hpk_kpi   = new TH1D("hpk_kpi"   ,"momentum distribution",200,0,2.2);
  hpk_kpi->SetLineColor(1);
  TH1D *hpk_kpipi = new TH1D("hpk_kpipi" ,"momentum distribution",200,0,2.2);
  hpk_kpipi->SetLineColor(2);
  TH1D *hpk_kk    = new TH1D("hpk_kk"    ,"momentum distribution",200,0,2.2);
  hpk_kk->SetLineColor(3);
  TH1D *hpk_kkpipi= new TH1D("hpk_kkpipi","momentum distribution",200,0,2.2);
  hpk_kkpipi->SetLineColor(4);
  
  double ppi1,ppi2,ppi3,ppi6pi,ppikkpipi;
  double pk1, pk2, pk3, pkkkpipi;
  treekpi->SetBranchAddress("p2",&ppi1);
  treekpi->SetBranchAddress("p1",&pk1);
  
  treekpipi->SetBranchAddress("ppi1",&ppi2);
  treekpipi->SetBranchAddress("pk1",&pk2);
  
  treepsip->SetBranchAddress("p1",&ppi3);
  treekk->SetBranchAddress("p1",&pk3);
  double mass;
  treekk->SetBranchAddress("mass",&mass);

  tree6pi->SetBranchAddress("p1",&ppi6pi);

  treekkpipi->SetBranchAddress("ppi1",&ppikkpipi);
  treekkpipi->SetBranchAddress("pk1",&pkkkpipi);

  for (int i=0; i<treekpi->GetEntries(); i++){
    treekpi->GetEntry(i);
	hppi_kpi->Fill(ppi1);
	hpk_kpi->Fill(pk1);
  }
  for (int i=0; i<treekpipi->GetEntries(); i++){
    treekpipi->GetEntry(i);
	hppi_kpipi->Fill(ppi2);
	hpk_kpipi->Fill(pk2);
  }
  for (int i=0; i<treepsip->GetEntries(); i++){
    treepsip->GetEntry(i);
	hppi_pipill->Fill(ppi3);
  }
  for (int i=0; i<treekk->GetEntries(); i++){
    treekk->GetEntry(i);
	if (mass > 1.2) hpk_kk->Fill(pk3);
  }
  for (int i=0; i<tree6pi->GetEntries(); i++){
    tree6pi->GetEntry(i);
	hppi_6pi->Fill(ppi6pi);
  }
  for (int i=0; i<treekkpipi->GetEntries(); i++){
    treekkpipi->GetEntry(i);
	hppi_kkpipi->Fill(ppikkpipi);
	hpk_kkpipi->Fill(pkkkpipi);
  }
  // p pi distribution
  new TCanvas();
  gStyle->SetOptStat(0);
  hppi_kpi->Draw();
  hppi_kpipi->Scale(hppi_kpi->GetEntries()/hppi_kpipi->GetEntries());
  hppi_kpipi->Draw("same");
  hppi_pipill->Scale(hppi_kpi->GetEntries()/hppi_pipill->GetEntries()/2);
  hppi_pipill->Draw("same");
  hppi_6pi->Scale(hppi_kpi->GetEntries()/hppi_6pi->GetEntries());
  hppi_6pi->Draw("same");
    hppi_kkpipi->Scale(hppi_kpi->GetEntries()/hppi_kkpipi->GetEntries());
    hppi_kkpipi->Draw("same");
  
  TLegend *legend = new TLegend(0.65,0.65,0.85,0.85);
  legend->AddEntry(hppi_kpi, "one #pi in K #pi");
  legend->AddEntry(hppi_kpipi, "one #pi in K #pi #pi");
  legend->AddEntry(hppi_pipill, "one #pi in #pi #pi J/#psi");
  legend->AddEntry(hppi_6pi, "one #pi in 6 #pi");
  legend->AddEntry(hppi_kkpipi, "one #pi in K K #pi #pi");
  legend->Draw();

  // p k distribution
  new TCanvas();
  //gStyle->SetOptStat(0);
  hpk_kpi->Draw();
  hpk_kpipi->Scale(hpk_kpi->GetEntries()/hpk_kpipi->GetEntries());
  hpk_kpipi->Draw("same");
  //hpk_kk->Scale(hpk_kpi->GetEntries()/hpk_kk->GetEntries()/1.5);
  //hpk_kk->Draw("same");
  hpk_kkpipi->Scale(hpk_kpi->GetEntries()/hpk_kkpipi->GetEntries());
  hpk_kkpipi->Draw("same");
  
  legend = new TLegend(0.65,0.65,0.85,0.85);
  legend->AddEntry(hpk_kpi, "one K in K #pi");
  legend->AddEntry(hpk_kpipi, "one K in K #pi #pi");
  //legend->AddEntry(hpk_kk, "one K in #phi #rightarrow K K");
  legend->AddEntry(hpk_kkpipi, "one K in K K #pi #pi");
  legend->Draw();

}
