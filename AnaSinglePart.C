#include <fstream>
#include <iostream>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
using namespace std;

int main(int argc, char** argv)
{
  std::string filename;
  std::string outdir=".";
  if (argc>1) filename = argv[1];
  if (argc>2) outdir = argv[2];
  ifstream inpar(filename.c_str());
  if (!inpar){
    std::cout<<"can not open:"<<filename<<std::endl;
    return -1;
  }

  double p[100], fac[100], face[100];
  double pe[100]={0};
  int count=0;
  char tmpchr[100];
  
  while (!inpar.eof()){
    count++;
    inpar >>p[count-1]>>fac[count-1]>>face[count-1];
    inpar.read(tmpchr,2);
    inpar.seekg(-2,std::ios::cur);
  }
  TCanvas *c1 = new TCanvas();
  TGraphErrors *graph = new TGraphErrors(count,p,fac,pe,face);
  graph->SetTitle("factors");
  graph->GetYaxis()->SetRangeUser(0.99,1.01);
  graph->SetMarkerStyle(5);
  gStyle->SetOptFit(1111);
  graph->Draw("AP");
  graph->Fit("pol0","","",p[0],p[count-1]);
  TF1 *f1 = graph->GetFunction("pol0");
  TPaveText *pt1 = new TPaveText(0.12,0.70,0.5,0.80,"BRNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(4000);
  pt1->SetTextAlign(12);
  pt1->SetTextFont(42);
  pt1->SetTextSize(0.035);
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",f1->GetParameter(0),f1->GetParError(0));
  pt1->AddText(tmpchr);
  pt1->Draw();
  
  sprintf(tmpchr,"%s/Ks_factors.pdf",outdir.c_str());
  c1->Print(tmpchr);
  return 0;
}
