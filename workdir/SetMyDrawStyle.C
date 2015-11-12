void SetMyDrawStyleGra(TCanvas *c1)
{
  //gStyle->SetFitFormat(".6f");
  //gStyle->SetOptFit(1111);
 	int n = c1->GetListOfPrimitives()->GetSize();
	std::string histname;
	for (int i=0; i<n; i++){
		if (!strncmp(c1->GetListOfPrimitives()->At(i)->ClassName(), "TGraph",6)){
			histname = (c1->GetListOfPrimitives()->At(i))->GetName();
			break;
		}
	} 
  //((TGraph*)c1->GetPrimitive("Graph"))->SetTitle("factor at each energy point");
  //((TPaveText*)c1->GetPrimitive("title"))->SetTextSize(0.05);
  // void SetMargin(Float_t left, Float_t right, Float_t bottom, Float_t top)
  c1->SetMargin(0.12,0.12,0.12,0.1);

  //((TPaveText*)c1->GetPrimitive("title"))->SetTextSize(0.04);

  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetNdivisions(505);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetLabelSize(0.04);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTickLength(0.04);
  //((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitle("CM energy (GeV)");
  //((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitle("Factor");
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleSize(0.04);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleOffset(1.2);

  //((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetRangeUser(0.996, 1.006);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetNdivisions(505);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetLabelSize(0.04);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTickLength(0.04);
  //((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitle("Factor");
  //((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitle("#delta E");
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleSize(0.04);
  ((TGraph*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleOffset(1.2);
}

void SetMyDrawStyleRoo(TCanvas *c1)
{
  //gStyle->SetFitFormat(".6f");
  //gStyle->SetOptFit(1111);
  
  //((TGraph*)c1->GetPrimitive("Graph"))->SetTitle("factor at each energy point");
  //((TPaveText*)c1->GetPrimitive("title"))->SetTextSize(0.05);
  // void SetMargin(Float_t left, Float_t right, Float_t bottom, Float_t top)
  //c1->SetMargin(0.15,0.15,0.15,0.15);
  c1->SetMargin(0.12,0.12,0.12,0.1);
  	
	int n = c1->GetListOfPrimitives()->GetSize();
	std::string histname;
	for (int i=0; i<n; i++){
		if (!strncmp(c1->GetListOfPrimitives()->At(i)->ClassName(), "TH1",3)){
			histname = (c1->GetListOfPrimitives()->At(i))->GetName();
			break;
		}
	}

  //((TPaveText*)c1->GetPrimitive("title"))->SetTextSize(0.04);
  
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetNdivisions(505);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetLabelSize(0.04);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTickLength(0.04);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleSize(0.04);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleOffset(1.2);

  //((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetRangeUser(0.996, 1.006);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetNdivisions(505);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetLabelSize(0.04);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTickLength(0.04);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleSize(0.04);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleOffset(1.2);

  int Nbin = ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetNbinsX();
  double Xmin = ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->GetXmin();
  double Xmax = ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->GetXmax();
  double xgap = (Xmax-Xmin)/Nbin*1000;
  char title[100];
  sprintf(title,"Events / %2.1f MeV", xgap);
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->SetTitle("Invariant mass spectrum");
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitle("Mass (GeV)");
  ((TH1D*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitle(title);
}


void ShowErrorRange(TGraphErrors* graph)
{
  // de = 0
  TF1 *fcon = new TF1("fcon","0",0.99,1.01);
  fcon->SetLineColor(3);
  
  // fit function
  TF1 *fit = graph->GetFunction("facfit");
  double fac = fit->GetParameter(0);
  double slo = fit->GetParameter(1);
  double de = graph->GetErrorY(graph->GetN()/2);
  
  // construct edge functions
  TF1 *f2 = new TF1("f2","[1]*(x-[0])",0.99,1.01);
  TF1 *f3 = new TF1("f3","[1]*(x-[0])",0.99,1.01);
  f2->SetParameters(fac-de/slo,slo);
  f3->SetParameters(fac+de/slo,slo);
  
  f2->SetLineStyle(2);
  f3->SetLineStyle(2);

  double fmin = f2->GetX(0);
  double fmax = f3->GetX(0);
  TArrow *arr0 = new TArrow(fac ,-0.002,fac ,0);
  TArrow *arr1 = new TArrow(fmin,-0.002,fmin,0);
  TArrow *arr2 = new TArrow(fmax,-0.002,fmax,0);
  arr0->SetArrowSize(0.035);
  arr1->SetArrowSize(0.035);
  arr2->SetArrowSize(0.035);
  arr0->SetLineWidth(2);
  arr1->SetLineWidth(2);
  arr2->SetLineWidth(2);

  graph->Draw();
  f2->Draw("same");
  f3->Draw("same");
  fcon->Draw("same");
  
  arr0->Draw();
  arr1->Draw();
  arr2->Draw();
}

void SetMyDrawStyleTF1(TCanvas *c1)
{
  //gStyle->SetFitFormat(".6f");
  //gStyle->SetOptFit(1111);
 	int n = c1->GetListOfPrimitives()->GetSize();
	std::string histname;
 	for (int i=0; i<n; i++){
 		if (!strncmp(c1->GetListOfPrimitives()->At(i)->ClassName(), "TF1",3)){
			histname = (c1->GetListOfPrimitives()->At(i))->GetName();
			break;
		}
	} 
  //((TGraph*)c1->GetPrimitive("Graph"))->SetTitle("factor at each energy point");
  //((TPaveText*)c1->GetPrimitive("title"))->SetTextSize(0.05);
  // void SetMargin(Float_t left, Float_t right, Float_t bottom, Float_t top)
  //c1->SetMargin(0.15,0.15,0.15,0.15);
  c1->SetMargin(0.12,0.12,0.12,0.1);

  //((TPaveText*)c1->GetPrimitive("title"))->SetTextSize(0.04);

  ((TF1*)c1->GetPrimitive(histname.c_str()))->SetTitle("Get Factor");
  //((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetNdivisions(505);
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetLabelSize(0.045);
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTickLength(0.045);
  //((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitle("CM energy (GeV)");
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitle("f_{K}");
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleSize(0.045);
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetXaxis()->SetTitleOffset(1.3);

  //((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetRangeUser(0.996, 1.006);
  //((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetNdivisions(505);
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetLabelSize(0.045);
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTickLength(0.045);
  //((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitle("Factor");
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitle("f_{#pi}");
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleSize(0.045);
  ((TF1*)c1->GetPrimitive(histname.c_str()))->GetYaxis()->SetTitleOffset(1.3);
}


