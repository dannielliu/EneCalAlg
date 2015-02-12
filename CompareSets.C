#include <iostream>
#include <fstream>
#include <string>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"

using namespace std;

int main(int argc, char** argv)
{

  char filenames[200][100];
  int fileNo = 0;
  if (argc<2){
    cout<<"no input file"<<endl;
    return -1;
  }
  if (argc>200){
    cout<<"too many input files."<<endl;
    return 1;
  }
  cout<<"There are "<<argc-1<<" input files."<<endl;;
  for (int i=1;i<argc;i++){
    sprintf(filenames[i-1],"%s",argv[i]);
    fileNo++;
  }

  //int ppart;
  int validNo;
  int np[200];
  double p[200][20];
  double pe[200][20];
  double fac[200][20];
  double face[200][20];

  char tmpchr[100];
  //double p1,fac1,face1;
  //int count[20]={0}; // count vaild part
  //double facv[20][200]; // buffer for pars
  //double facev[20][200];
  double pcut[21];
  pcut[0] =0.0;
  pcut[1] =0.05;
  pcut[2] =0.10;
  pcut[3] =0.15;
  pcut[4] =0.20;
  pcut[5] =0.25;
  pcut[6] =0.30;
  pcut[7] =0.35;
  pcut[8] =0.40;
  pcut[9] =0.45;
  pcut[10]=0.50;
  pcut[11]=0.60;
  pcut[12]=0.70;
  pcut[13]=0.80;
  pcut[14]=0.90;
  pcut[15]=1.00;
  pcut[16]=1.20;
  pcut[17]=1.40;
  pcut[18]=1.60;
  pcut[19]=1.80;
  pcut[20]=2.00;

  validNo = fileNo;
  int fid=0;
  for (int parti=0;parti<fileNo;parti++)
  {
    ifstream inpar(filenames[parti]);
    if (!inpar){
      cout<<"Warning: can not open file: "<<filenames[parti]<<endl;
      validNo--;
      continue;
    }
    cout<<"reading file: "<<filenames[parti]<<endl;
    np[fid]=0;
    while(!inpar.eof()){
      inpar>>p[fid][np[fid]]>>fac[fid][np[fid]]>>face[fid][np[fid]];
      pe[fid][np[fid]]=0;
      np[fid]++;
      
      inpar.read(tmpchr,2);
      inpar.seekg(-2,std::ios::cur);
    }
    fid++;
  }

  TCanvas *c1 = new TCanvas();
  TFile *f = new TFile("aaa.root","RECREATE");
  char lname[200][100]={"Rvalue","xyz4230","xyz4260","xyz4360","xyz4420","xyz4600"};
  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
  for (fid=0;fid<validNo;fid++){
    std::cout<<"file "<<fid<<" have "<<np[fid]<<" parts."<<std::endl;
    TGraphErrors *graph = new TGraphErrors(np[fid],p[fid],fac[fid],pe[fid],face[fid]);
    graph->SetFillColor(0);
    graph->SetMarkerStyle(21+fid);
    graph->SetLineColor(2+fid);
    graph->SetMarkerColor(2+fid);
    if (fid == 0){
      graph->SetTitle("factors");
      graph->GetYaxis()->SetRangeUser(0.995,1.005);
      graph->GetXaxis()->SetRangeUser(0.,1.5);
      graph->SetMarkerStyle(5);
      gStyle->SetOptFit(1111);
      graph->Draw("APC");
    }
    else graph->Draw("CP");
    //sprintf(tmpchr,"sample %d",fid);
    legend->AddEntry(graph,lname[fid]);

  } 
  legend->Draw(); 
  sprintf(tmpchr,"./Ks_compare_factors_%d.pdf",fid);
  c1->Print(tmpchr);
  c1->Write();
  f->Close();
/*
  TGraphErrors *graph = new TGraphErrors(cp,fp,ffac,fpe,fface);
  graph->SetTitle("factors");
  graph->GetYaxis()->SetRangeUser(0.995,1.005);
  graph->SetMarkerStyle(5);
  gStyle->SetOptFit(1111);
  graph->Draw("AP");
  graph->Fit("pol0","","",fp[0],fp[cp-1]);
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
*/  
  return 0;

}
