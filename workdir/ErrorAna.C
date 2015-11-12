#include <iostream>
#include <fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iomanip>
using namespace std;

int main(int argc, char** argv)
{
  char filename[1000];
  if (argc<2) return -1;
  sprintf(filename, "%s",argv[1]);
  //char outroot[1000]={"bbb.root"};
  //if (argc>2) sprintf(outroot,"%s",argv[2]);

  //TFile file(outroot,"recreate");
  gStyle->SetOptFit(1111);
  gStyle->SetFitFormat("5.6g");

  double ene,peak,peakerr;
  double x[200],xe[200];
  double y1[200],y1e[200];
  double y2[200],y2e[200];
  double y3[200],y3e[200];
  double y4[200],y4e[200];
  
  ifstream in(filename);
  if(!in){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  ofstream out("table.tex");
  
  out<<"\\begin{longtable}{ccccc}%{|c|c|c|c|}"<<std::endl;
  out<<"\\caption{summary of energy calibration (unit: MeV)}\\\\"<<std::endl;
  out<<"\\toprule"<<std::endl;
  out<<"nominal energy & observed energy & corrected energy & statistic error & system error \\\\"<<std::endl;
  out<<"\\midrule"<<std::endl;
  out<<"\\endfirsthead"<<std::endl;
  out<<"\\midrule"<<std::endl;
  out<<"nominal energy & observed energy & corrected energy & statistic error & system error \\\\"<<std::endl;
  out<<"\\midrule"<<std::endl;
  out<<"\\endhead"<<std::endl;
  out<<"\\midrule"<<std::endl;
  out<<"\\multicolumn{5}{r}{continue \\dots} \\\\"<<std::endl;
  out<<"\\endfoot"<<std::endl;
  out<<"\\bottomrule"<<std::endl;
  out<<"\\endlastfoot"<<std::endl;
  //out<<"nominal energy & corrected energy & statistic error & system error \\\\"<<std::endl;
  
  out<<setiosflags(ios::fixed)<<setprecision(6);
  int id=0;
  char tmp[1000];
  bool valid=false;
  while (!in.eof() && id<104){
  	if (!valid){
  	  in.getline(tmp,1000);
	  valid = !strncmp(tmp,"#system_error",13);
	  //std::cout<<valid<<'\t'<<tmp<<std::endl;
	  continue;
	}
    	in>>ene;
	x[id] = ene*1000;///1000;
	xe[id]=0;
	in>>peak>>peakerr;
	y1[id] = peak*1000;
	y1e[id] = peakerr*1000;
	in.getline(tmp,1000);

    	in>>ene;
        in>>peak>>peakerr;
        y2[id] = peak*1000;
        y2e[id] = peakerr*1000;
        in.getline(tmp,1000);
	
    in>>ene;
    in>>peak>>peakerr;
    y3[id] = peak*1000;
    y3e[id] = peakerr*1000;
    in.getline(tmp,1000);
	
    in>>ene;
    in>>peak>>peakerr;
    y4[id] = peak*1000;
    y4e[id] = peakerr*1000;
    in.getline(tmp,1000);
	//std::cout<<ene<<'\t'<<y1[id]<<'\t'<<y1e[id]<<std::endl;
	// data check
	//y1e[id]<y2e[id]? y1e[id] = y2e[id] : y2e[id] = y1e[id];
    //out<<"\\hline"<<std::endl;
    out<<setiosflags(ios::fixed)<<setprecision(0)<<x[id]<<" & "<<setprecision(2)<<y1[id]<<" & "<<y2[id]<<" & "<<y1e[id]<<" & "<<(fabs(y2[id]-y3[id])+fabs(y4[id]-y2[id]))/2.<<"\\\\"<<std::endl;
	if (y1e[id] < 5e-4 && y2e[id] < 5e-4) {
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
  //out<<"\\hline"<<std::endl;
  out<<"\\end{longtable}"<<std::endl;

  return 0;

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

//c1->Write();
//return 0;

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
  
  c1->Write();
  std::cout<<"finish write canvas"<<std::endl;
  //file.Close();

  return 0;

}
