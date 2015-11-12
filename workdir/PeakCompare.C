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
  gStyle->SetFitFormat(".6f");

  double ene,peak,peakerr;
  double x[200],xe[200];
  double x1[200],x1e[200];
  double x2[200],x2e[200];
  double x3[200],x3e[200];
  double x4[200],x4e[200];
  double y1[200],y1e[200];
  double y2[200],y2e[200];
  double y3[200],y3e[200];
  double y4[200],y4e[200];
  
  ifstream in(filename);
  if(!in){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  ifstream in2(argv[2]);
  if(!in2){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  int id=0;
  char tmp[1000];
  bool valid=false;
  while (!in.eof() && id<104){
  	if (!valid){
  	  in.getline(tmp,1000);
	  //valid = !strncmp(tmp,"#xyz",4);
	  valid = !strncmp(tmp,"#last",5);
	  //std::cout<<valid<<'\t'<<tmp<<std::endl;
	  continue;
	}
    in>>ene;
	x1[id] = ene;///1000;
	x1e[id]=0;
	in>>peak>>peakerr;
	y1[id] = peak - ene;
	y1e[id] = peakerr;
	in.getline(tmp,1000);
	
	in>>ene;
	in>>peak>>peakerr;
	x2[id] = ene;///1000;
	x2e[id]=0;
	y2[id] = peak-ene;
	y2e[id] = peakerr;
	in.getline(tmp,1000);

	//in.getline(tmp,1000);
	//in.getline(tmp,1000);

	if (y1e[id] < 1e-4 && y2e[id]<1e-4 ) {
	  y1e[id] = 1e-3; 
	  y2e[id] = 1e-3;
	  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	if (y1e[id] < 1e-4 ) {
	  y1e[id] = y2e[id]; 
	  //std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	if (y2e[id] < 1e-4 ) {
	  y2e[id] = y1e[id]; 
	  //std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
	}
	
	id++;
	//in.getline(tmp,1000);
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
  
  TGraphErrors *graph1 = new TGraphErrors(id,x1,y1,x1e,y1e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(4);
  graph1->SetLineColor(1);
  graph1->SetFillColor(0);
  graph1->Fit("pol0");
  graph1->GetFunction("pol0")->SetLineColor(1);
 
  TGraphErrors *graph2 = new TGraphErrors(id,x2,y2,x2e,y2e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph2->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph2->GetYaxis()->SetTitle("#delta E (GeV)");
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(5);
  graph2->SetLineColor(2);
  graph2->SetFillColor(0);
  graph2->Fit("pol0");
  graph2->GetFunction("pol0")->SetLineColor(2);
   
  graph1->Draw("AP");
  graph2->Draw("P");
  //std::cout<<"finish raw graph"<<std::endl;

  TLegend *legend = new TLegend(0.75,0.75,0.95,0.9);
  legend->AddEntry(graph1,"raw");
  legend->AddEntry(graph2,"cor");
  legend->Draw();
  c1->Write();
  
  id = 0;
  valid = false;
  while (!in2.eof() && id<104){
  	if (!valid){
  	  in2.getline(tmp,1000);
	  valid = !strncmp(tmp,"#last",5);
	  //std::cout<<valid<<'\t'<<tmp<<std::endl;
	  continue;
	}
    	in2>>ene;
	x3[id] = ene/1e3;///1000;
	x3e[id]=0;
	in2>>peak>>peakerr;
	y3[id] = peak/1e3 - ene/1e3;
	y3e[id] = peakerr/1e3;

	in2>>peak>>peakerr;
	x4[id] = ene/1e3;
	x4e[id]=0;
	y4[id] = peak/1e3 - ene/1e3;
	y4e[id] = peakerr/1e3;
	in2.getline(tmp,1000);

	//std::cout<<ene<<'\t'<<y1[id]<<'\t'<<y1e[id]<<std::endl;
	// data check
	
////////if (y3e[id] < 5e-4 ) {
////////  y3e[id] = 2e-3; 
////////  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
////////}
////////if (y3e[id] < 5e-4 ) {
////////  y3e[id] = 2e-3; 
////////  std::cout<<"error is abnormal at energy point: "<<ene<<std::endl;
////////}
	
	id++;
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

  TCanvas *c2 = new TCanvas("c2","c2");
  TGraphErrors *graph3 = new TGraphErrors(id,x3,y3,x3e,y3e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph3->SetMarkerColor(3);
  graph3->SetMarkerStyle(6);
  graph3->SetLineColor(3);
  graph3->SetFillColor(0);
  graph3->Fit("pol0");
  graph3->GetFunction("pol0")->SetLineColor(3);
  graph3->Draw("AP");
  //std::cout<<"finish corrected graph"<<std::endl;

  TGraphErrors *graph4 = new TGraphErrors(id,x4,y4,x4e,y4e);
  //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
  //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
  graph4->SetMarkerColor(4);
  graph4->SetMarkerStyle(7);
  graph4->SetLineColor(4);
  graph4->SetFillColor(0);
  graph4->Fit("pol0");
  graph4->GetFunction("pol0")->SetLineColor(4);
  graph4->Draw("P");

  TLegend *legend2 = new TLegend(0.75,0.75,0.95,0.9);
  legend2->AddEntry(graph3,"raw");
  legend2->AddEntry(graph4,"cor");
  legend2->Draw();
  c2->Write();

  TCanvas *c3 = new TCanvas("c3","c3");
  graph2->Draw("AP");
  graph4->Draw("P");

  TLegend *legend3 = new TLegend(0.75,0.75,0.95,0.9);
  legend3->AddEntry(graph2,"result with this method");
  legend3->AddEntry(graph4,"previous result");
  legend3->Draw();
  std::cout<<"finish legend"<<std::endl;
  c3->Write();

  graph1->Write();
  graph2->Write();
  graph3->Write();
  graph4->Write();
  
  std::cout<<"finish write canvas"<<std::endl;
  file.Close();

  return 0;

}
