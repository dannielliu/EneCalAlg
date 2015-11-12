{
  double x[5]  = {4230,  4260,  4360,  4420,  4600};
  double xe[5] = {0};
//// pipill
//double y[5]  = {1.00090,1.00079, 1.00078, 1.00124, 1.00132};
//double ye[5] = {7.8e-5, 9.12e-5, 0.00012, 9.28e-5, 0.00016};
  // kpi
  double y[5]  = {1.00054, 1.00053, 1.00055, 1.00069, 1.00066};
  double ye[5] = {2.16e-5, 3.07e-5, 2.81e-5, 1.79e-5, 3.05e-5};

  TCanvas *c1 = new TCanvas();
  gStyle->SetOptFit(1111);
  gStyle->SetFitFormat("5.6g");
  TGraphErrors *graph = new TGraphErrors(5,x,y,xe,ye);
  graph->Draw("AP");
  graph->Fit("pol0");

}
