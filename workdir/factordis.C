#include <iostream>
#include <fstream>
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

int main(int argc, char** argv)
//int factordis()
{
  char filename[1000]={"factorcol"};
  //if (argc<2) return -1;
  //sprintf(filename, "%s",argv[1]);

  TFile file("bbb.root","recreate");
  double ene,peak,peakerr;
  double x[200],xe[200];
  double y[20][200],ye[20][200];
//double y2[200],y2e[200];
//double y3[200],y3e[200];
//double y4[200],y4e[200];
  
  ifstream in(filename);
  if(!in){
    std::cout<<"Can not open file: "<<filename<<std::endl;
  }
  int id=0;
  char tmp[1000];
  char flag[] = "cut";
  int iarray=0;
	for(int i=0; i<10; i++){
	  x[i] = 2./10*i-1+0.1;
	  xe[i] = 0.1;
	}
  std::string line;
  //while (!in.eof())
  for (int cnt =0;cnt<6; cnt++)
  {
    in.getline(&tmp[0],1000);
	bool notcut = strncmp(tmp,flag,3);
	std::cout<<tmp<<'\t'<<in.tellg()<<'\t'<<notcut<<std::endl;
	if (notcut) in.getline(&tmp[0],1000);
	notcut = strncmp(tmp,flag,3);
	if (notcut) in.getline(&tmp[0],1000);
	notcut = strncmp(tmp,flag,3);
	std::cout<<tmp<<'\t'<<in.tellg()<<'\t'<<notcut<<std::endl;
	for(int i=0; i<10; i++){
	  do {in.read(tmp,1);}
	  while(tmp[0]!='\t' && tmp[0]!=' ');
	  in>>y[iarray][i]>>ye[iarray][i];
	  in.getline(tmp,1000);
	  //std::cout<<"\t"<<y[iarray][i]<<'\t'<<ye[iarray][i]<<std::endl;
	}
	iarray++;
	in.read(tmp,2);
	in.seekg(-2,std::ios::cur);
  }

  TCanvas *c1  = new TCanvas();
  TGraphErrors *graph[20];
  for (int i=0; i<6; i++){
    graph[i] = new TGraphErrors(10,x,y[i],xe,ye[i]);
    //graph1->GetXaxis()->SetRangeUser(3.5,5.0);
    //graph1->GetYaxis()->SetRangeUser(-0.01,0.01);
    graph[i]->SetMarkerColor(i+1);
    graph[i]->SetMarkerStyle(i+4);
    graph[i]->SetLineColor(i+1);
    graph[i]->SetFillColor(0);
    if (i==0) graph[i]->Draw("AP");
	else graph[i]->Draw("P");
 // graph[i]->Fit("pol0");
 // graph[i]->GetFunction("pol0")->SetLineColor(1);
  }
  
  TLegend *legend = new TLegend(0.75,0.75,0.95,0.9);
  legend->AddEntry(graph[0],"#psi' #rightarrow #pi #pi l l cut #pi^{+}");
  legend->AddEntry(graph[1],"#psi' #rightarrow #pi #pi l l cut #pi^{-}");
  legend->AddEntry(graph[2],"D^{0} #rightarrow K K cut K^{+}");
  legend->AddEntry(graph[3],"D^{0} #rightarrow K K cut K^{-}");
  legend->AddEntry(graph[4],"D^{0} #rightarrow K^{-} #pi^{+} cut K^{-}");
  legend->AddEntry(graph[5],"D^{0} #rightarrow K^{-} #pi^{+} cut #pi^{+}");
  legend->Draw();
  
  c1->Write();
  file.Close();

  return 0;

}
