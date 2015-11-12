#include <iostream>
#include <fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
int main(int argc, char** argv)
{
  char filename[1000];
  if (argc<3) return -1;
  sprintf(filename, "%s",argv[1]);
  char outroot[1000]={"bbb.root"};
  //if (argc>2) sprintf(outroot,"%s",argv[2]);

  TFile file(outroot,"recreate");
  gStyle->SetOptFit(1111);
  gStyle->SetFitFormat("5.6g");

  double ene,peak,peakerr;
  double x[200],xe[200];
  double y1[200],y1e[200];
  double y2[200],y2e[200];
  double y3[200],y3e[200];
  double y4[200],y4e[200];
  double y5[200],y5e[200];
  
  ifstream in(filename);
  if(!in){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  ifstream in2(argv[2]);
  if(!in2){
    std::cout<<"Can not open file: "<<argv[2]<<std::endl;
  }
  int id=0;
  char tmp[1000];
  bool valid=false;
  while (!in.eof() && id<104){
  	if (!valid){
  	  in.getline(tmp,1000);
	  valid = !strncmp(tmp,"#new",4);
	  continue;
	}
    	in>>ene;
	x[id] = ene;///1000;
	xe[id]=0;
	in>>peak>>peakerr;
	y1[id] = peak - ene;
	y1e[id] = peakerr;
	in.getline(tmp,1000);

    	in>>ene;
	in>>peak>>peakerr;
	y2[id] = peak - ene;
	y2e[id] = peakerr;
	in.getline(tmp,1000);
	
	// data check
	
	if (y1e[id] < 2e-4 && y2e[id] < 2e-4) {
	  y1e[id] = 2e-3; 
	  y2e[id] = y1e[id];
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	if (y1e[id] < 5e-4) {
	  y1e[id] = y2e[id];
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	if (y2e[id] < 5e-4) {
	  y2e[id] = y1e[id];
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}

	
	id++;
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
  std::cout<<"finish reading data"<<std::endl;
  TCanvas *c1  = new TCanvas("c1","c1",800,600);
  TGraphErrors *graph1 = new TGraphErrors(id,x,y1,xe,y1e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(4);
  graph1->SetLineColor(1);
  graph1->SetFillColor(0);
  graph1->Draw("AP");
  graph1->Fit("pol0");
  graph1->GetFunction("pol0")->SetLineColor(1);
  std::cout<<"finish raw graph"<<std::endl;

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


  id=0;
  valid = false;
  while (!in2.eof() && id<104){
    double trashvar;
  	if (!valid){
  	  in2.getline(tmp,1000);
	  valid = !strncmp(tmp,"#new",4);
	  continue;
	}
	in2>>ene;
	x[id] = ene;///1000;
	xe[id]=0;
	in2>>peak>>peakerr;
	y3[id] = peak - ene;
	y3e[id] = peakerr;
	in2.getline(tmp,1000);

        in2>>ene;
	in2>>peak>>peakerr;
	y4[id] = peak -ene;
	y4e[id] = peakerr;
	in2.getline(tmp,1000);
    	
	// data check
	if (y3e[id] < 2e-4 && y4e[id] < 2e-4) {
	  y3e[id] = 2e-3; 
	  y4e[id] = y3e[id];
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	if (y3e[id] < 5e-4) {
	  y3e[id] = y4e[id];
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	if (y4e[id] < 5e-4) {
	  y4e[id] = y3e[id];
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}


	id++;
	in2.read(tmp,10);
	in2.seekg(-10,std::ios::cur);
    
	bool endreading = false;
	for (int i=0; i<10; i++){
	  if (tmp[i] == '#') {
	    endreading = true;
		in2.seekg(i,std::ios::cur);
		in2.getline(tmp,1000);
	  }
	}
	if (endreading) break;
  }
  std::cout<<"finish raw graph "<< id <<std::endl;
//TGraphErrors *graph3 = new TGraphErrors(id,x,y3,xe,y3e);
////graph1->GetXaxis()->SetRangeUser(3.5,5.0);
////graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
//graph3->SetMarkerColor(3);
//graph3->SetMarkerStyle(6);
//graph3->SetLineColor(3);
//graph3->SetFillColor(0);
//graph3->Draw("P");
//graph3->Fit("pol0");
//graph3->GetFunction("pol0")->SetLineColor(3);
//std::cout<<"finish corrected graph"<<std::endl;



  TLegend *legend = new TLegend(0.75,0.75,0.95,0.9);
  legend->AddEntry(graph1,"raw");
  legend->AddEntry(graph2,"cor");
 // legend->AddEntry(graph3,"cor2");
  legend->Draw();
  std::cout<<"finish legend"<<std::endl;
  
  for (int i=0; i<104; i++){
    double sum_mw = y2[i]/(y2e[i]*y2e[i]) + y4[i]/(y4e[i]*y4e[i]);
    double sum_w  = 1./(y2e[i]*y2e[i]) + 1./(y4e[i]*y4e[i]);
    y5[i] = sum_mw/sum_w;
    y5e[i]= 1./sqrt(sum_w);
    std::cout<< x[i] <<'\t'<<y5[i]<<'\t'<<y5e[i]<<std::endl;
  }
  TCanvas *c2  = new TCanvas("c2","c2",800,600);
  TGraphErrors *graph5 = new TGraphErrors(104,x,y5,xe,y5e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph5->SetMarkerColor(1);
  graph5->SetMarkerStyle(4);
  graph5->SetLineColor(1);
  graph5->SetFillColor(0);
  graph5->Draw("AP");
  graph5->Fit("pol0");
  graph5->GetFunction("pol0")->SetLineColor(1);


  c1->Write();
  c2->Write();
  std::cout<<"finish write canvas"<<std::endl;
  file.Close();

  return 0;

}
