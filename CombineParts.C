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
  char filenames[200][1000];
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

  const int Npart = 25;
  //const int Ncos = 10;
  char tmpchr[100];
  double p1,fac1,face1;
  int count[Npart]={0}; // count vaild part
  double facv[Npart][200]; // buffer for pars
  double facev[Npart][200];
  double p[Npart];
  double pe[Npart];
//double cos[Ncos];
//double cose[Ncos];
  double fac[Npart];
  double face[Npart];

  double pcut[Npart+1];
//double coscut[Ncos+1];
  pcut[0] =0.0  ;//coscut[0] = -1.0;
  pcut[1] =0.05 ;//coscut[2] = -0.6;
  pcut[2] =0.10 ;//coscut[4] = -0.2;
  pcut[3] =0.15 ;//coscut[6] =  0.2;
  pcut[4] =0.20 ;//coscut[8] =  0.6;
  pcut[5] =0.25 ;//coscut[10]=  1.0;
  pcut[6] =0.30 ;
  pcut[7] =0.35 ;
  pcut[8] =0.40 ;
  pcut[9] =0.45 ;
  pcut[10]=0.50 ;
  pcut[11]=0.60 ;
  pcut[12]=0.70 ;
  pcut[13]=0.80 ;
  pcut[14]=0.90 ;
  pcut[15]=1.00 ;
  pcut[16]=1.10 ;
  pcut[17]=1.20 ;
  pcut[18]=1.30 ;
  pcut[19]=1.40 ;
  pcut[20]=1.50 ;
  pcut[21]=1.60 ;
  pcut[22]=1.70 ;
  pcut[23]=1.80 ;
  pcut[24]=1.90 ;
  pcut[25]=2.00 ;

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
	  int jd=-1;
      for (int pi=0;pi<Npart;pi++){
        if (p1>=pcut[pi] && p1<pcut[pi+1]){
          id = pi;
          break;
        }
      }
	//for (int ci=0;ci<Ncos;ci++){
	//  if (cos1>=coscut[ci] && cos1<coscut[ci+1]){
	//    jd = ci;
	//    break;
	//  }
	//}
      facv[id][count[id]] = fac1;
      facev[id][count[id]]= face1;
      count[id]++;
      
      inpar.read(tmpchr,2);
      inpar.seekg(-2,std::ios::cur);
    }
  }

  for (int pi=0;pi<Npart;pi++){
  //for (int ci=0;ci<Ncos;ci++){
    double sumerr=0;// sum of 1/error^2
    double sumfac=0;
    cout<<"part "<<pi<<" have "<<count[pi]<<" samples!"<<endl;
    if (count[pi]==0) continue;
    for (int i=0;i<count[pi];i++){
      sumerr += 1./(facev[pi][i]*facev[pi][i]);
      sumfac += 1./(facev[pi][i]*facev[pi][i])*facv[pi][i];
    }
    p[pi] = pcut[pi] + (pcut[pi+1]-pcut[pi])/2;
    pe[pi] = 0;
////cos[ci] = coscut[ci] + (coscut[ci+1]-coscut[ci])/2;
////cose[ci] = 0;
    if (sumerr!=0){
      fac[pi] = sumfac/sumerr;
      face[pi]= TMath::Sqrt(1./sumerr);
    }
    else fac[pi]=0;
    cout<<"part "<<pi<<" error "<<face[pi]<<" factor "<<fac[pi]<<endl;
  //}// ci end
  }// pi end
  
  //double fp[Npart],ffac[Npart],fface[Npart];
  //double fpe[Npart]={0};
  //int cp=0;
  sprintf(tmpchr,"combinedpar.txt");
  ofstream of(tmpchr);
  for (int i=0;i<Npart;i++){
  //for (int j=0;j<Ncos;j++){
    if (count[i]==0) continue;
    //fp[cp] = p[i];
    //ffac[cp]=fac[i];
    //fface[cp]=face[i];

    of<<p[i]<<"\t"<<fac[i]<<"\t"<<face[i]<<endl;
    //cout<<fp[cp]<<"\t"<<ffac[cp]<<"\t"<<fface[cp]<<endl;

    //cp++;
  //}
  }
 /*
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
  */
  return 0;

}
