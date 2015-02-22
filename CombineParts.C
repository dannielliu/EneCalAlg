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

  int ppart;
  double p[20];
  double pe[20];
  double fac[20];
  double face[20];

  char tmpchr[100];
  double p1,fac1,face1;
  int count[20]={0}; // count vaild part
  double facv[20][200]; // buffer for pars
  double facev[20][200];
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

  for (int parti=0;parti<fileNo;parti++)
  {
    ifstream inpar(filenames[parti]);
    if (!inpar){
      cout<<"Warning: can not open file: "<<filenames[parti]<<endl;
      continue;
    }
    cout<<"reading file: "<<filenames[parti]<<endl;
    while(!inpar.eof()){
      inpar>>p1>>fac1>>face1;
      int id=-1;
      for (int pi=0;pi<20;pi++){
        if (p1>=pcut[pi] && p1<pcut[pi+1]){
          id = pi;
          break;
        }
      }
      facv[id][count[id]] = fac1;
      facev[id][count[id]]= face1;
      count[id]++;
      
      inpar.read(tmpchr,2);
      inpar.seekg(-2,std::ios::cur);
    }
  }

  for (int pi=0;pi<20;pi++){
    double sumerr=0;// sum of 1/error^2
    double sumfac=0;
    cout<<"part  "<<pi<<" have "<<count[pi]<<" samples!"<<endl;
    if (count[pi]==0) continue;
    for (int i=0;i<count[pi];i++){
      sumerr += 1./(facev[pi][i]*facev[pi][i]);
      sumfac += 1./(facev[pi][i]*facev[pi][i])*facv[pi][i];
    }
    p[pi] = pcut[pi] + (pcut[pi+1]-pcut[pi])/2;
    pe[pi] = 0;
    if (sumerr!=0){
      fac[pi] = sumfac/sumerr;
      face[pi]= TMath::Sqrt(1./sumerr);
    }
    else fac[pi]=0;
    cout<<"part  "<<pi<<" error "<<face[pi]<<" factor "<<fac[pi]<<endl;
  }
  
  double fp[20],ffac[20],fface[20];
  double fpe[20]={0};
  int cp=0;
  sprintf(tmpchr,"combinedpar.txt");
  ofstream of(tmpchr);
  for (int i=0;i<20;i++){
    if (count[i]==0) continue;
    fp[cp] = p[i];
    ffac[cp]=fac[i];
    fface[cp]=face[i];

    of<<fp[cp]<<"\t"<<ffac[cp]<<"\t"<<fface[cp]<<endl;
    cout<<fp[cp]<<"\t"<<ffac[cp]<<"\t"<<fface[cp]<<endl;

    cp++;
  }
 
  TCanvas *c1 = new TCanvas();
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
  
  sprintf(tmpchr,"./Ks_factors_sum.pdf");
  c1->Print(tmpchr);
  return 0;

}
