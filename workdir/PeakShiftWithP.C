
int PeakShiftWithP()
{

  // D -> K pi, peak shift with pK
  double x[5] = {0.5, 0.7, 0.9, 1.1, 1.3};
  double xe[5]= {0.1, 0.1, 0.1, 0.1, 0.1};
  double y[5] = {1.863733, 1.863818, 1.863911, 1.864223, 1.864516};
  double ye[5]= {0.000450, 0.000109, 0.000082, 0.000118, 0.000315};
  
  TGraphErrors* ps = new TGraphErrors(5,x,y,xe,ye);
  TF1* f1 = new TF1("f1","1.863990",0.4,1.4);
  
  TCanvas *c1 = new TCanvas();
  ps->Draw("APL");
  f1->SetLineColor(kRed);
  f1->Draw("same");
}
