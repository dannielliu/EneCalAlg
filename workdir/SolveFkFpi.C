#include <iostream>
#include <fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"

int main(int argc, char** argv)
{
  char filename[1000];
  if (argc<3) return -1;
  sprintf(filename, "%s",argv[1]);
  char outroot[1000]={"bbb.root"};
  //if (argc>2) sprintf(outroot,"%s",argv[2]);
  bool checksize = true;
  bool checkenergy = true;

  TFile file(outroot,"recreate");
  gStyle->SetOptFit(0);
  //gStyle->SetFitFormat("5.7g");
  double ene,peak,peakerr;
  //double y3[200],y3e[200];
  //double y4[200],y4e[200];
  double e1[200], e2[200];
  double a1[200], b1[200], c11[200], c12[200];
  double a2[200], b2[200], c21[200], c22[200];
  double a1e[200], b1e[200], c11e[200], c12e[200];
  double a2e[200], b2e[200], c21e[200], c22e[200];

  ifstream in1(argv[1]);
  if(!in1){
    std::cout<<"Can not open file: "<<argv[1]<<std::endl;
  }
  ifstream in2(argv[2]);
  if(!in2){
    std::cout<<"Can not open file: "<<argv[2]<<std::endl;
  }
  int id=0; 
  int size1, size2;
  char tmp[1000];
  bool valid=false;
  while (!in1.eof() && id<104){
  	if (!valid){
  	  in1.getline(tmp,1000);
	  valid = !strncmp(tmp,"#new",4);
	  //std::cout<<valid<<'\t'<<tmp<<std::endl;
	  continue;
	}
    	double trashvar;
	in1>>ene;
	e1[id] = ene;///1000;
	in1>>trashvar;
	in1>>trashvar;
	in1>>trashvar;
	in1.read(tmp,6);
	in1>>a1[id];
	in1.read(tmp,9);
	in1>>c11[id];
	in1.read(tmp,7);
	in1>>a1e[id];
	in1.read(tmp,6);
	in1>>c11e[id];
	in1.read(tmp,9);
	in1>>c11e[id];
	
	in1>>ene;
	//std::cout<<"ene @"<<ene<<std::endl;
	if (ene != e1[id]){
	  std::cout<<"Error: energy not equal! "<<ene<<std::endl;
	  return 1;
	}
	in1>>trashvar;
	in1>>trashvar;
	in1>>trashvar;
    	in1.read(tmp,6);
	in1>>b1[id];
	in1.read(tmp,9);
	in1>>c12[id];
	in1.read(tmp,7);
	in1>>b1e[id];
	in1.read(tmp,6);
	in1>>c12e[id];
	in1.read(tmp,9);
	in1>>c12e[id];
	
	id++;
	//std::cout<<ene<<"\t"<<id<<std::endl;
	in1.getline(tmp,1000);
	in1.read(tmp,10);
	in1.seekg(-10,std::ios::cur);
	bool endreading = false;
	for (int i=0; i<10; i++){
	  if (tmp[i] == '#') {
	    endreading = true;
	    in1.seekg(i,std::ios::cur);
		in1.getline(tmp,1000);
		std::cout<<tmp<<std::endl;
		break;
	  }
	}
	if (endreading) break;
  }
  in1.close();
  in1.clear(std::ios::goodbit);
  size1 = id;
  id = 0;
  valid = false;
  while (!in2.eof() && id<104){
  	if (!valid){
  	  in2.getline(tmp,1000);
	  valid = !strncmp(tmp,"#new",4);
	  continue;
	}
	//if(strncmp(tmp,"#new",4)) continue; // if success, cmp return 0, else return other value;
    	double trashvar;
	in2>>ene;
	e2[id] = ene;///1000;
	in2>>trashvar;
	in2>>trashvar;
	in2>>trashvar;
	in2.read(tmp,6);
	in2>>a2[id];
	in2.read(tmp,9);
	in2>>c21[id];
	in2.read(tmp,7);
	in2>>a2e[id];
	in2.read(tmp,6);
	in2>>c21e[id];
	in2.read(tmp,9);
	in2>>c21e[id];
	
	in2>>ene;
	//std::cout<<"ene2 @"<<ene<<std::endl;
	if (ene != e2[id]){
	  std::cout<<"Error: energy not equal!"<<std::endl;
	  return 1;
	}
	in2>>trashvar;
	in2>>trashvar;
	in2>>trashvar;
    	in2.read(tmp,6);
	in2>>b2[id];
	in2.read(tmp,9);
	in2>>c22[id];
	in2.read(tmp,7);
	in2>>b2e[id];
	in2.read(tmp,6);
	in2>>c22e[id];
	in2.read(tmp,9);
	in2>>c22e[id];
	
	//std::cout<<"file 2 "<<ene<<std::endl;
	
	id++;
	std::cout<<ene<<"\t"<<id<<std::endl;
	in2.getline(tmp,1000);
	in2.read(tmp,10);
	in2.seekg(-10,std::ios::cur);
	bool endreading = false;
	for (int i=0; i<10; i++){
	  if (tmp[i] == '#') {
	    endreading = true;
	    in2.seekg(i,std::ios::cur);
		in2.getline(tmp,1000);
		std::cout<<tmp<<std::endl;
		break;
	  }
	}
	if (endreading) break;
  }
  in2.close();
  in2.clear();
  size2 = id;
  id = 0;

  double x[200],xe[200];
  double y1[200],y1e[200];
  double y2[200],y2e[200];
  if (checksize)
  if (size1 != size2) {
    std::cout<<"Please check data, size not equal, size 1 is "<<size1<<" size 2 is "<< size2<<std::endl;
    return 2;
  }
  for (id = 0; id < size1; id++){
    if (checkenergy)
    	if (e1[id]!=e2[id]){
		std::cout<<"energy not equal "<< e1[id]<<" "<<e2[id]<<" id is "<<id<<std::endl;
		return 3;
	}
    x[id] = e1[id];
    xe[id] = 0;
    double c1 = (c11[id]+c12[id])/2. ;
    double c2 = (c21[id]+c22[id])/2. ;
    double c1e = (c11e[id]+c12e[id])/2. ;
    double c2e = (c21e[id]+c22e[id])/2. ;
    y1[id] = (a1[id]*c2-a2[id]*c1)/(a2[id]*b1[id]-a1[id]*b2[id]);
    y2[id] = (b1[id]*c2-b2[id]*c1)/(b2[id]*a1[id]-b1[id]*a2[id]);
    y1e[id]= sqrt( 
    			TMath::Power(a2[id]*y2[id]*a1e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a1[id]*y2[id]*a2e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a2[id]*y1[id]*b1e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a1[id]*y1[id]*b2e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a2[id]*c1e/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a1[id]*c2e/(a2[id]*b1[id]-a1[id]*b2[id]),2)
			);
    //y2e[id]= ;
    y2e[id]= sqrt( 
    			TMath::Power(a1[id]*y1[id]*a2e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a2[id]*y1[id]*a1e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a1[id]*y2[id]*b2e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a2[id]*y2[id]*b1e[id]/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a1[id]*c2e/(a2[id]*b1[id]-a1[id]*b2[id]),2)
		+	TMath::Power(a2[id]*c1e/(a2[id]*b1[id]-a1[id]*b2[id]),2)
			);
    y1[id] += 1;
    y2[id] += 1;
    std::cout<<x[id]<<'\t'<<y1[id]<<'\t'<<y2[id]<<'\t'<<y1e[id]<<'\t'<<y2e[id]<<std::endl;
  }

  
  double mean1, chi2_1;
  double mean2, chi2_2;
  //id --;
  std::cout<<"Finish reading data, point No. is "<< size1<<std::endl;
  TCanvas *c1  = new TCanvas("c1","c1",800,600);
  TGraphErrors *graph1 = new TGraphErrors(id,&x[0],&y1[0],&xe[0],&y1e[0]);
  graph1->GetXaxis()->SetRangeUser(3.7,4.7);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(4);
  graph1->SetLineColor(1);
  graph1->SetFillColor(0);
  graph1->Draw("AP");
  graph1->Fit("pol0");
  graph1->GetFunction("pol0")->SetLineColor(1);
  mean1 = graph1->GetFunction("pol0")->GetParameter(0);
  chi2_1 = graph1->Chisquare(graph1->GetFunction("pol0"))/103;
  std::cout<<"finish raw graph "<< id <<std::endl;
  
//c1->Write();
//return 0;

  TGraphErrors *graph2 = new TGraphErrors(id,&x[0],&y2[0],&xe[0],&y2e[0]);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(5);
  graph2->SetLineColor(2);
  graph2->SetFillColor(0);
  graph2->Draw("P");
  graph2->Fit("pol0");
  graph2->GetFunction("pol0")->SetLineColor(2);
  mean2 = graph2->GetFunction("pol0")->GetParameter(0);
  chi2_2 = graph2->Chisquare(graph2->GetFunction("pol0"))/103;
  std::cout<<"finish corrected graph"<<std::endl;

  TLegend *legend = new TLegend(0.75,0.75,0.95,0.9);
  char label[100];
  sprintf(label,"f_{#pi}, #mu=%.6f, #chi^{2}=%.3f",mean1,chi2_1);
  legend->AddEntry(graph1,label);
  sprintf(label,"f_{K}, #mu=%.6f, #chi^{2}=%.3f",mean2,chi2_2);
  legend->AddEntry(graph2,label);
  legend->Draw();
  std::cout<<"finish legend"<<std::endl;

  c1->Write();
  std::cout<<"finish write canvas"<<std::endl;
  file.Close();

  return 0;

}
