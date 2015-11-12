#include <iostream>
#include <fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include <sstream>

int PeakCompare2()
{
  //char filename[1000]={"f6picutp"};
  char filename[1000]={"fkkpipi"};
  gStyle->SetOptFit(1111);
  gStyle->SetFitFormat(".6f");

  double ene,peak,peakerr;
  //double x[200],xe[200];
  double x[10][200],xe[10][200];
  double y[10][200],ye[10][200];
//double x2[200],x2e[200];
//double x3[200],x3e[200];
//double x4[200],x4e[200];
//double y1[200],y1e[200];
//double y2[200],y2e[200];
//double y3[200],y3e[200];
//double y4[200],y4e[200];
  
  ifstream in(filename);
  if(!in){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  int id=0;
  int setsize1 = 0;
  int setsize2 = 0;
  char tmp[1000];
  bool valid=false;
  istringstream iss;
  int pos_cut=0;
  int pos_nocut=0;
  // find position of "#cut"
  in.clear();
  while (!in.eof())
  {
    in.getline(tmp,1000);
    iss.clear();
    iss.str(tmp);
    if (strncmp(tmp,"#cut",4)==0) pos_cut = in.tellg();
    if (strncmp(tmp,"#nocut",6)==0) pos_nocut = in.tellg();
  }
  if (pos_cut==0){
    std::cout<<"Can not find #cut"<<std::endl;
    return -1;
  }
  if (pos_nocut==0){
    std::cout<<"Can not find #nocut"<<std::endl;
    return -1;
  }

  in.clear();
  in.seekg(pos_cut,std::ios::beg);
  id=0;
  int enecount=0;
  bool loop = true;
  while (loop){
     in.getline(tmp,1000);
     iss.clear(); // clear error status
     iss.str(tmp);
     if (in.eof()) break; // file end
     if (strncmp(tmp,"#EOS",4)==0) break;// data set end
     if (iss.peek()==-1) continue;// empty line
     if (iss.peek()=='#') continue;// invalid line
     
     iss>>ene>>peak>>peakerr;
     if (ene == x[enecount][id]) enecount ++;
     else {enecount = 0; id++;}
     x[enecount][id] = ene;///1000;
     xe[enecount][id]=0;
     y[enecount][id] = peak - ene;
     ye[enecount][id] = peakerr;
  }
  setsize1 = id;
  
  in.clear();
  in.seekg(pos_nocut,std::ios::beg);
  id=0;
  enecount=5;
  loop = true;
  while (loop){
     in.getline(tmp,1000);
     iss.clear(); // clear error status
     iss.str(tmp);
     if (in.eof()) break; // file end
     if (strncmp(tmp,"#EOS",4)==0) break;// data set end
     if (iss.peek()==-1) continue;// empty line
     if (iss.peek()=='#') continue;// invalid line
     
     iss>>ene>>peak>>peakerr;
     if (ene == x[enecount][id]) enecount ++;
     else {enecount = 5; id++;}
     x[enecount][id] = ene;///1000;
     xe[enecount][id]=0;
     y[enecount][id] = peak - ene;
     ye[enecount][id] = peakerr;
  }
  setsize2 = id;

  //id --;
  std::cout<<"finish reading data"<<std::endl;
  TCanvas *c1  = new TCanvas("c1","c1",800,600);
  gStyle->SetOptFit(0);

  double mean1,mean2,mean3,mean4;
  TGraphErrors *graph1 = new TGraphErrors(setsize1,&x[0][1],&y[0][1],&xe[0][1],&ye[0][1]);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph1->GetYaxis()->SetTitle("#delta E (GeV)");
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(4);
  graph1->SetLineColor(1);
  graph1->SetFillColor(0);
  graph1->Fit("pol0");
  graph1->GetFunction("pol0")->SetLineColor(1);
  mean1 = graph1->GetFunction("pol0")->GetParameter(0);
 
  TGraphErrors *graph2 = new TGraphErrors(setsize1,&x[1][1],&y[1][1],&xe[1][1],&ye[1][1]);
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(5);
  graph2->SetLineColor(2);
  graph2->SetFillColor(0);
  graph2->Fit("pol0");
  graph2->GetFunction("pol0")->SetLineColor(2);
  mean2 = graph2->GetFunction("pol0")->GetParameter(0);
  
  TGraphErrors *graph3 = new TGraphErrors(setsize2,&x[5][1],&y[5][1],&xe[5][1],&ye[5][1]);
  graph3->SetMarkerColor(3);
  graph3->SetMarkerStyle(6);
  graph3->SetLineColor(3);
  graph3->SetFillColor(0);
  graph3->Fit("pol0");
  graph3->GetFunction("pol0")->SetLineColor(3);
  mean3 = graph3->GetFunction("pol0")->GetParameter(0);
 
  TGraphErrors *graph4 = new TGraphErrors(setsize2,&x[6][1],&y[6][1],&xe[6][1],&ye[6][1]);
  graph4->SetMarkerColor(4);
  graph4->SetMarkerStyle(7);
  graph4->SetLineColor(4);
  graph4->SetFillColor(0);
  graph4->Fit("pol0");
  graph4->GetFunction("pol0")->SetLineColor(4);
  mean4 = graph4->GetFunction("pol0")->GetParameter(0);
   
  graph1->Draw("AP");
  graph2->Draw("P");
  graph3->Draw("P");
  graph4->Draw("P");
  //std::cout<<"finish raw graph"<<std::endl;

  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
  char label[100];
  sprintf(label,"p < 1.2 GeV raw data, #mu = %.5f",mean1);
  legend->AddEntry(graph1,label);
  sprintf(label,"p < 1.2 GeV corrected data, #mu = %.5f",mean2);
  legend->AddEntry(graph2,label);
  sprintf(label,"all p raw data, #mu = %.5f",mean3);
  legend->AddEntry(graph3,label);
  sprintf(label,"all p corrected data, #mu = %.5f",mean4);
  legend->AddEntry(graph4,label);
  legend->Draw();
  //c1->Write();
  //

  // calculate chi2 and error
  // chi2 = dEi^2*Wi/n, 
  // dEi = E1 - E2
  // Wi   = 1/sqrt(sigma_1^2 + sigma_2^2),
  // 
  // error = sum( dEi^2 * Wi )/(n*Wi)
  if (setsize1 != setsize2){
    std::cout<<"data set with different size: "<<setsize1 <<", "<<setsize2<<std::endl;
    std::cout<<"will not calculate chi2 and error"<<std::endl;
    return -2;
  }
  double wi=0;
  double sumw=0;
  double dEi=0;
  double dEw=0;
  for (int i=0; i<setsize1;i++){
     dEi = y[1][i] - y[6][i];
     wi  = sqrt(ye[1][i]*ye[1][i]+ye[6][i]*ye[6][i]);
     sumw += wi;
     dEw  += dEi*dEi*wi;
  }
  double chi2 = dEw/setsize1;
  double syserr = sqrt(chi2/sumw);
  std::cout<<"Chi2 when comparing correct energy is "<< chi2<<std::endl;
  std::cout<<"sys err when comparing correct energy is "<< syserr<<std::endl;
 
  return 0;
}
