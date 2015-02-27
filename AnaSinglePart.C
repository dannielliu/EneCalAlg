#include <fstream>
#include <iostream>
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TFile.h"
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

  double p[1000], fac[1000], face[1000];
  double cos[1000];
  double pe[1000]={0};
  double cose[1000]={0};
  int count=0;
  char tmpchr[1000];
  
  std::cout<<"start reading file."<<std::endl;
  while (!inpar.eof()){
    count++;
    //std::cout<<"reading ..."<<std::endl;
    inpar >>p[count-1]>>cos[count-1]>>fac[count-1]>>face[count-1];
    //std::cout<<"checking ..."<<std::endl;
    inpar.read(tmpchr,2);
    inpar.seekg(-2,std::ios::cur);
  }
  std::cout<<"reading file end."<<std::endl;
  std::cout<<"there are "<<count<<" points."<<std::endl;

  TFile *f = new TFile("aaaa.root","recreate");
  TCanvas *c1 = new TCanvas();
  TGraph2DErrors *graph = new TGraph2DErrors(count,p,cos,fac,pe,cose,face);
  graph->SetTitle("factors");
  graph->GetZaxis()->SetRangeUser(0.95,1.05);
  graph->SetMarkerStyle(5);
  //gStyle->SetOptFit(1111);
  graph->Draw("p0");
  //graph->Fit("pol0","","",p[0],p[count-1]);
  //TF1 *f1 = graph->GetFunction("pol0");
  //TPaveText *pt1 = new TPaveText(0.12,0.70,0.5,0.80,"BRNDC");
  //pt1->SetBorderSize(0);
  //pt1->SetFillStyle(4000);
  //pt1->SetTextAlign(12);
  //pt1->SetTextFont(42);
  //pt1->SetTextSize(0.035);
  //sprintf(tmpchr,"factor = %1.6f #pm %1.6f",f1->GetParameter(0),f1->GetParError(0));
  //pt1->AddText(tmpchr);
  //pt1->Draw();
  c1->Write();

  sprintf(tmpchr,"%s/Ks_factors.pdf",outdir.c_str());
  c1->Print(tmpchr);
  return 0;
}
