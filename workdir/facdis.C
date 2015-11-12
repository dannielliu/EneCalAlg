#include <iostream>
#include <fstream>
#include <sstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
using std::istringstream;

int main(int argc, char** argv)
{
  char filename[1000];
  if (argc<2) return -1;
  sprintf(filename, "%s",argv[1]);
  char outroot[1000]={"bbb.root"};
  if (argc>2) sprintf(outroot,"%s",argv[2]);

  TFile file(outroot,"recreate");
  gStyle->SetOptFit(1111);
  gStyle->SetFitFormat("5.6g");
  double ene,peak,peakerr;
  double x[200],xe[200];
  double y1[200],y1e[200];
  double y2[200],y2e[200];
  //double y3[200],y3e[200];
  //double y4[200],y4e[200];
  
  ifstream in(filename);
  if(!in){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  int id=0;
  char tmp[1000];
  bool valid=false;
  while (!in.eof() && id<104){
    	double trashvar;
  	in.getline(tmp,1000);
  	if (!valid){
	  valid = !strncmp(tmp,"#new",4);
	  //std::cout<<valid<<'\t'<<tmp<<std::endl;
	  continue;
	}

	istringstream iss;
	iss.str(tmp);

	iss >> ene >> peak >> peakerr;

	x[id] = ene;///1000;
	xe[id]=0;
	y1[id] = peak ;
	y1e[id] = fabs(peakerr);
    
	id++;

	iss >> y2[id] >> y2e[id];
	
	

	//std::cout<<ene<<"\t"<<id<<std::endl;
	in.read(tmp,10);
	in.seekg(-10,std::ios::cur);
    
	bool endreading = false;
	for (int i=0; i<10; i++){
	  if (tmp[i] == '#') {
	    endreading = true;
	    in.seekg(i,std::ios::cur);
		in.getline(tmp,1000);
		std::cout<<tmp<<std::endl;
		break;
	  }
	}
	if (endreading) break;
  }
  //id --;
  std::cout<<"finish reading data, point no. is "<< id<<std::endl;
  TCanvas *c1  = new TCanvas("c1","c1",800,600);
  TGraphErrors *graph1 = new TGraphErrors(id,x,y1,xe,y1e);
  graph1->GetXaxis()->SetRangeUser(3.7,4.7);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(4);
  graph1->SetLineColor(1);
  graph1->SetFillColor(0);
  graph1->Draw("AP");
  graph1->Fit("pol0");
  graph1->GetFunction("pol0")->SetLineColor(1);
  std::cout<<"finish raw graph "<< id <<std::endl;
  
  c1->Write();
  
  TCanvas *c2 = new TCanvas();
  c2->SetName("c2");
  TGraphErrors *graph2 = new TGraphErrors(id,x,y2,xe,y2e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(5);
  graph2->SetLineColor(2);
  graph2->SetFillColor(0);
  graph2->Draw("AP");
  graph2->Fit("pol0");
  graph2->GetFunction("pol0")->SetLineColor(2);
  std::cout<<"finish corrected graph"<<std::endl;
 
  c2->Write();
  return 0;
/*
  id=0;
  while (!in.eof() && id<104){
    double trashvar;
	in>>ene;
	x[id] = ene;///1000;
	xe[id]=0;
	//in>>trashvar;
	in>>peak>>peakerr;
	y2[id] = peak ;
	y2e[id] = peakerr;
    
	id++;
	//std::cout<<ene<<"\t"<<id<<std::endl;
	in.getline(tmp,1000);
    //in.getline(tmp,1000);
	in.read(tmp,10);
	in.seekg(-10,std::ios::cur);
    
	bool endreading = false;
	for (int i=0; i<10; i++){
	  if (tmp[i] == '#') {
	    endreading = true;
		in.seekg(i,std::ios::cur);
		in.getline(tmp,1000);
	  }
	}
	if (endreading) break;
  }
  std::cout<<"finish raw graph "<< id <<std::endl;

  TGraphErrors *graph2 = new TGraphErrors(id,x,y2,xe,y2e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(5);
  graph2->SetLineColor(2);
  graph2->SetFillColor(0);
  graph2->Draw("P");
  graph2->Fit("pol0");
  graph2->GetFunction("pol0")->SetLineColor(2);
  std::cout<<"finish corrected graph"<<std::endl;

  TLegend *legend = new TLegend(0.75,0.75,0.95,0.9);
  legend->AddEntry(graph1,"R scan data");
  legend->AddEntry(graph2,"xyz data");
  legend->Draw();
  std::cout<<"finish legend"<<std::endl;
*/
  c1->Write();
  std::cout<<"finish write canvas"<<std::endl;
  file.Close();

  return 0;

}
