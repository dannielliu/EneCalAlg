#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TGaxis.h"
#include "TPad.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
using RooFit::Title;
using RooFit::Components;
using RooFit::LineStyle;
using RooFit::LineColor;
using RooFit::Range;
extern std::string outputdir;
extern std::vector<double> px1,py1,pz1,px2,py2,pz2;
extern std::vector<double> mass1;
extern double m0;
extern double mparticle,mparticle2,mparticle3,mparticle4;
extern double sigma;
extern double width;
extern double weight;
extern double slope;

std::string outputdir=".";
void FitAndSave(TH1D **hmass, 
	      int &part,
	      RooPlot *&xframe,
                RooDataHist* &data_k, 
	      RooAddPdf* &sum,
	      RooGaussian &gaus,
	      RooChebychev &bkg,
	      RooRealVar &x,
	      RooRealVar &mean,
	      RooRealVar &sigma1,
	      RooRealVar &co1,
	      RooRealVar &signal,
	      RooRealVar &background,
	      double &philow,
	      double &phiup,
	      double *Ps,
	      ofstream &ofpar,
	      TFile *&f);

int main(int argc,char** argv)
{

   ifstream sample("sample.txt");
   ofstream ofpar;
   ofpar.open("parkk.txt",std::ios::app);
   ofstream detail;
   detail.open("detailkk.txt",std::ios::app);
   detail<<"kk algrithm: will give factors for kaon"<<std::endl;
   TFile *f=new TFile("plot.root","RECREATE");
   
   double args[10];
  //double mpeak=1.019455;
  double &mpeak=args[1]=1.019455;
  double &msigma=args[2]=0.0025;
  double &pcut=args[3]=1.0;
  double &runmode=args[4]=0;
  //double mpeak=3.0969;
  //double sigma=0.0025;
  // deal with input arguments
  std::map<std::string,int> ArgMap;
  ArgMap.insert(std::make_pair("-peak",1));
  ArgMap.insert(std::make_pair("-p",1));
  ArgMap.insert(std::make_pair("-sigma",2));
  ArgMap.insert(std::make_pair("-w",2));
  ArgMap.insert(std::make_pair("-cut",3));
  ArgMap.insert(std::make_pair("-c",3));
  ArgMap.insert(std::make_pair("-runmode",4));
  std::cout<<"argc is "<<argc<<"\n";
  for (int i=1;i<argc;i++){
    std::cout<<"arg "<<i<<" is "<<argv[i]<<"\n";
    std::string tmpstr=argv[i];
    if(tmpstr.find("-")==0 && argc>i){ // if the first char is '-', then it is a specific argument
      std::cout<<"setting arg:"<<&tmpstr[1]<<"\n";
      args[ArgMap[tmpstr]] = atof(argv[i+1]);
      i++;
    }
  }
  for (int i=1;i<=2;i++){
    std::cout<<"arg "<<i<<": "<<args[i]<<"\n";
  }
 
   double pxa,pya,pza,pxb,pyb,pzb;
   m0 = mpeak;
   sigma=msigma;
   double philow=m0-5*sigma;
   double phiup=m0+5*sigma;
   // try to use roofit
   RooRealVar x("x","energy",m0,philow,phiup,"GeV");
   RooRealVar mean("mean","mean of gaussian",m0,philow,phiup);
   RooRealVar sigma1("sigma1","width of gaussian",0.003,0.0001,0.2);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma1);
   RooRealVar co1("co1","coefficient #1",0,-1000.,1000.);
   //RooRealVar co4("co4","coefficient #4",0);
   RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",3000,1,1000000);//event number
   RooRealVar background("background"," ",0,0,100000);
   RooPlot *xframe;
   RooDataHist *data_k;
   RooAddPdf *sum;
 
   int NP=1;
   // split momentum from 0.13 to 1.13 GeV
   double Ps[NP+1];
   for(int i=0;i<=NP;i++){
     Ps[i]=2.0/NP*i;
   }

   double factor,factorlow,factorup;
   double minimum;
   double minx,miny,minz;
   std::string tmpstr;
   TH1D *hp[NP];
   TH1D *hmass[NP];
   for(int i=0;i<NP;i++){
     char name[100];
     sprintf(name,"hp_%1.2f_to_%1.2f_GeV",Ps[i],Ps[i+1]);
     hp[i]=new TH1D(name,name,100,0,2);
     sprintf(name,"hmass_%1.2f_to_%1.2f_GeV",Ps[i],Ps[i+1]);
     hmass[i]=new TH1D(name,name,100,philow,phiup);
   }
   TH2D *hp2 =new TH2D("hp2","p for K+ K-",100,0,2,100,0,2);
   //TH1D *h1 = new TH1D("h1","momentum of kaon",100,0,3.0);
   //h2->SetLineColor(2);
   // ~~~~~~~~~kaon part~~~~~~~~~~
  
   // likelihood method
   //m0 = 1.01946;
   //double mpeak;
   //m0 = 1.019455;
   //sigma=sigma1.getVal();
   width = 10.*sigma;
   mparticle=0.493677;
   //mparticle=0.000511;
   double massold=0;
   double deltaup=0,deltadown=0;
   double deltaup2=0,deltadown2=0;
   if((int)runmode==0)
   for(int part=0;part<NP;part++){
     ofpar<<Ps[part]<<"\t"<<Ps[part+1]<<std::endl;
     // pre fit
     hmass[part]->Reset();
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string name;
       getline(sample,name);
       //if(ngam>0) continue;
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(massold!=0) hmass[part]->Fill(massold);
       massold=mass;
       // if (Cut(ientry) < 0) continue;
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     
     // likelihood
     px1.clear();
     px2.clear();
     py1.clear();
     py2.clear();
     pz1.clear();
     pz2.clear();
     hp2->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       std::string tmpcha;
       getline(sample,tmpcha);
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
        //if(ngam>0) continue;
       double mass;
       //double totpx,totpy,totpz,tote;
       double kapp,kamp,kape,kame;
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(mass>m0-width/2. && mass<m0+width/2.)
       {
         hp[part]->Fill(kapp);
         hp2->Fill(kapp,kamp);
         //hmass[part]->Fill(mass);
         px1.push_back(pxa);
         px2.push_back(pxb);
         py1.push_back(pya);
         py2.push_back(pyb);
         pz1.push_back(pza);
         pz2.push_back(pzb);
       }
    }
     char name[100];
     TCanvas *c2=new TCanvas("c2","likelihood",800,600);

     hp[part]->Draw();
     sprintf(name,"%s/momentum_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
     /*
     //TF2 *likeli=new TF2("likeli",maxlikelihood1,0.95,1.05,0.01,1.0);
     //TF3 *likeli=new TF3("likeli",maxlikelihood1_3,0.95,1.05,0.01,0.99,-100,100);
     //likeli->Draw();
     //sprintf(name,"%s/momentum_%1.2f_%1.2f_2D.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     //c2->Print(name);
     //likeli->Draw("surf2");
     //sprintf(name,"%s/momentum_%1.2f_%1.2f_2D2.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     //c2->Print(name);
     //likeli->GetMinimumXY(factor,miny);
     //weight = miny;
     weight = 1;
     //slope = minz;
     */
      
     TF1 *likeli_1=new TF1("likeli_1",maxlikelihood,0.95,1.05);
     likeli_1->Draw();
     sprintf(name,"%s/likeli_%1.2f_%1.2f_1D.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     minimum = likeli_1->GetMinimum(0.95,1.05);
     factor  = likeli_1->GetMinimumX(0.95,1.05);
     factorlow=likeli_1->GetX(minimum+1,0.95,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.05);
     /*
     mparticle2 = mparticle;
     TF2 *likeli = new TF2("likeli",maxlikelihood2_2,0.95,1.05,0.95,1.05,1);
     likeli->SetParameter(0,pcut);
     double f1,f2;
     likeli->Draw("surf");
     sprintf(name,"%s/likeli_kk_2D.eps",outputdir.c_str());
     c2->Print(name);
     likeli->GetMinimumXY(f1,f2);
     ofpar<<"\t"<<f1<<"\t"<<f2<<"\t"<<"\t\t"<<weight<<"\t"<<slope<<std::endl;*/
     
     ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
     //ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t\t"<<weight<<"\t"<<slope<<std::endl;
     //detail<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
     //detail<<"signal weight is "<<weight<<" best factor  "<<likeli_1->GetMinimumX(0.99,1.01)<<std::endl;
     delete c2;

     // using the factor to fit
     hmass[part]->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     massold=0;
     //factor = 0.996196;
     //TH2D *massshift= new TH2D("massshift","mass shift",100,philow,phiup,100,-0.01,0.0);
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string tmpcha;
       getline(sample,tmpcha);
      
       //if(ngam>0) continue;
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb,1,&factor);
       if(massold!=0) hmass[part]->Fill(massold);
       massold = mass;
       //ofpar<<"mass: "<<mass<<std::endl;
       // if (Cut(ientry) < 0) continue;
       
       //double masstmp;
       //masstmp = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       //if (mass<m0+sigma && masstmp>m0+sigma) deltaup++;
       //else if (mass<m0-sigma && masstmp>m0-sigma) deltadown++;
       //else if (mass>m0+sigma && masstmp<m0+sigma) deltadown2++;
       //else if (mass>m0-sigma && masstmp<m0-sigma) deltaup2++;
       //massshift->Fill(masstmp,mass-masstmp);
     }
     //massshift->Write();
     //std::cout<<"shift to high value: "<<deltaup<<" , shift to low value: "<<deltadown<<std::endl;
     //std::cout<<"shift to high value: "<<deltaup2<<" , shift to low value: "<<deltadown2<<std::endl;

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     hp[part]->Write();
     hp2->Write();
     //likeli->Write();
     //c2->Write();
   }
   else if ((int)runmode==1) // direct mass mode
   for(int part=0;part<NP;part++){
     ofpar<<Ps[part]<<"\t"<<Ps[part+1]<<std::endl;
     // pre fit
     hmass[part]->Reset();
     massold=0;
     while(!sample.eof()){
       double mass;
       sample>>mass;
       std::string name;
       getline(sample,name);
       if(massold!=0) hmass[part]->Fill(massold);
       massold=mass;
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     //xframe->Write();
     
     // likelihood
     mass1.clear();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       std::string tmpcha;
       getline(sample,tmpcha);
       double mass;
       sample>>mass;
        //if(ngam>0) continue;
       if(mass>m0-width/2. && mass<m0+width/2.)
       {
         mass1.push_back(mass);
       }
     }
  
     char name[100];
     TCanvas *c2=new TCanvas("c2","likelihood",800,600);
      
     TF1 *likeli_1=new TF1("likeli_1",maxlikelihood0,0.95,1.05);
     //likeli_1->SetNpx(1000);
     likeli_1->Draw();
     sprintf(name,"%s/likeli_%1.2f_%1.2f_1D.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     minimum = likeli_1->GetMinimum(0.95,1.05);
     factor  = likeli_1->GetMinimumX(0.95,1.05);
     factorlow=likeli_1->GetX(minimum+1,0.95,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.05);
     
     ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t\t"<<weight<<"\t"<<slope<<std::endl;
     //detail<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
     //detail<<"signal weight is "<<weight<<" best factor  "<<likeli_1->GetMinimumX(0.99,1.01)<<std::endl;
     likeli_1->Write();
     delete c2;

     // using the factor to fit
     hmass[part]->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     massold=0;
     while(!sample.eof()){
       double mass;
       sample>>mass;
       std::string tmpcha;
       getline(sample,tmpcha);
      
       if(massold!=0) hmass[part]->Fill(factor*massold);
       massold = mass;
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     //likeli->Write();
     //ac2->Write();
   } 
   else if ((int)runmode==2) // multimomentum correction factor mode
   for(int part=0;part<NP;part++){
     ofpar<<Ps[part]<<"\t"<<Ps[part+1]<<std::endl;
     // pre fit
     hmass[part]->Reset();
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string name;
       getline(sample,name);
       //if(ngam>0) continue;
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(massold!=0) hmass[part]->Fill(massold);
       massold=mass;
       // if (Cut(ientry) < 0) continue;
     }
     
     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     
     // likelihood
     px1.clear();
     px2.clear();
     py1.clear();
     py2.clear();
     pz1.clear();
     pz2.clear();
     hp2->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       std::string tmpcha;
       getline(sample,tmpcha);
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
        //if(ngam>0) continue;
       double mass;
       //double totpx,totpy,totpz,tote;
       double kapp,kamp,kape,kame;
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(mass>m0-width/2. && mass<m0+width/2.)
       {
         hp[part]->Fill(kapp);
         hp2->Fill(kapp,kamp);
         //hmass[part]->Fill(mass);
         px1.push_back(pxa);
         px2.push_back(pxb);
         py1.push_back(pya);
         py2.push_back(pyb);
         pz1.push_back(pza);
         pz2.push_back(pzb);
       }
     }
  
     char name[100];
     TCanvas *c2=new TCanvas("c2","likelihood",800,600);
     hp[part]->Draw();
     sprintf(name,"%s/momentum_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
     int ndim = 3;
     MultiDimFunction *likeli_multi = new MultiDimFunction(ndim,0.5,1.5);// n part, p start, p end
     likeli_multi->ShowParameters();
     double facs[ndim];
     likeli_multi->GetMinimumX(facs);
     for (int i=0;i<ndim;i++){
       ofpar<<"\t"<<facs[i];
     }
     ofpar<<std::endl;
     delete c2;
     /*
     TF1 *likeli_1=new TF1("likeli_1",maxlikelihood,0.95,1.05);
     //likeli_1->SetNpx(1000);
     likeli_1->Draw();
     sprintf(name,"%s/likeli_%1.2f_%1.2f_1D.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     minimum = likeli_1->GetMinimum(0.95,1.05);
     factor  = likeli_1->GetMinimumX(0.95,1.05);
     factorlow=likeli_1->GetX(minimum+1,0.95,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.05);
     
     ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t\t"<<weight<<"\t"<<slope<<std::endl;
     */

     // using the factor to fit
     hmass[part]->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     massold=0;
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string tmpcha;
       getline(sample,tmpcha);
      
       //if(ngam>0) continue;
       double mass;
       double par[]={0.5,1.5};
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb,ndim,facs,par);
       if(massold!=0) hmass[part]->Fill(massold);
       massold = mass;
       // if (Cut(ientry) < 0) continue;
     }
     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     
     hp[part]->Write();
     hp2->Write();
     //likeli->Write();
     //c2->Write();
   }
   else if ((int)runmode==3)
   for(int part=0;part<NP;part++){
     ofpar<<Ps[part]<<"\t"<<Ps[part+1]<<std::endl;
     // pre fit
     hmass[part]->Reset();
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string name;
       getline(sample,name);
       //if(ngam>0) continue;
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(massold!=0) hmass[part]->Fill(massold);
       massold=mass;
       // if (Cut(ientry) < 0) continue;
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     
     // likelihood
     px1.clear();
     px2.clear();
     py1.clear();
     py2.clear();
     pz1.clear();
     pz2.clear();
     hp2->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       std::string tmpcha;
       getline(sample,tmpcha);
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
        //if(ngam>0) continue;
       double mass;
       //double totpx,totpy,totpz,tote;
       double kapp,kamp,kape,kame;
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(mass>m0-width/2. && mass<m0+width/2.)
       {
         hp[part]->Fill(kapp);
         hp2->Fill(kapp,kamp);
         //hmass[part]->Fill(mass);
         px1.push_back(pxa);
         px2.push_back(pxb);
         py1.push_back(pya);
         py2.push_back(pyb);
         pz1.push_back(pza);
         pz2.push_back(pzb);
       }
    }
     char name[100];
     TCanvas *c2=new TCanvas("c2","likelihood",800,600);

     hp[part]->Draw();
     sprintf(name,"%s/momentum_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
      
     TF1 *likeli_1=new TF1("likeli_1",BiasCoe,0.95,1.05,1);
     likeli_1->SetParameter(0,m0);
     likeli_1->Draw();
     likeli_1->Write();
     c2->Write();
     sprintf(name,"%s/bias_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     minimum = likeli_1->GetMinimum(0.95,1.05);
     factor  = likeli_1->GetMinimumX(0.95,1.05);
     factorlow=likeli_1->GetX(minimum+1,0.95,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.05);
     
     ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
     //ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t\t"<<weight<<"\t"<<slope<<std::endl;
     delete c2;

     // using the factor to fit
     hmass[part]->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     massold=0;
     //factor = 0.996196;
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string tmpcha;
       getline(sample,tmpcha);
      
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb,1,&factor);
       if(massold!=0) hmass[part]->Fill(massold);
       massold = mass;
       //ofpar<<"mass: "<<mass<<std::endl;
       // if (Cut(ientry) < 0) continue;
       
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     hp[part]->Write();
     hp2->Write();
     //likeli->Write();
   }
   else if ((int)runmode==4)// 3 dimension p1 p2 theta
   for(int part=0;part<NP;part++){
     ofpar<<Ps[part]<<"\t"<<Ps[part+1]<<std::endl;
     // pre fit
     hmass[part]->Reset();
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string name;
       getline(sample,name);
       //if(ngam>0) continue;
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(massold!=0) hmass[part]->Fill(massold);
       massold=mass;
       // if (Cut(ientry) < 0) continue;
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     
     
     // 
     TH3D *hmass3D = new TH3D("hmass3D","hmass_2p_theta",100,0,2,100,0,2,100,0,TMath::Pi());
     hp2->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       std::string tmpcha;
       getline(sample,tmpcha);
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
        //if(ngam>0) continue;
       double mass;
       //double totpx,totpy,totpz,tote;
       double kapp,kamp,theta;
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       theta = acos((pxa*pxb+pya*pyb+pza*pzb)/(kapp*kamp));
       hmass3D->Fill(kapp,kamp,theta);
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb);
       if(mass>m0-width/2. && mass<m0+width/2.)
       {
         hp[part]->Fill(kapp);
         hp2->Fill(kapp,kamp);
         //hmass[part]->Fill(mass);
       }
     }
     char name[100];
     TCanvas *c2=new TCanvas("c2","likelihood",800,600);

     hp[part]->Draw();
     sprintf(name,"%s/momentum_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     std::cout<<"m0 is "<<m0<<std::endl;
      
     TF3 *dis = new TF3("dis",distribution,0,2.0,0,2.0,0.99,1.01,2);
     dis->SetParameter(0,2);
     dis->SetParameter(1,1.0);
     hmass3D->Draw();
     hmass3D->Fit(dis);
     hmass3D->Write();
     dis->Write();
     c2->Write();
     sprintf(name,"%s/dis_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     factor = dis->GetParameter(1);
     
     ofpar<<"\t"<<factor<<std::endl;
     delete c2;

     // using the factor to fit
     hmass[part]->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     massold=0;
     //factor = 0.996196;
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string tmpcha;
       getline(sample,tmpcha);
      
       double mass;
       mass = CalInvMass(mparticle,pxa,pya,pza, mparticle,pxb,pyb,pzb,1,&factor);
       if(massold!=0) hmass[part]->Fill(massold);
       massold = mass;
       //ofpar<<"mass: "<<mass<<std::endl;
       // if (Cut(ientry) < 0) continue;
       
     }

     FitAndSave(hmass, part, 
              xframe, data_k, sum, gaus, bkg, 
	    x, mean, sigma1, co1, signal, background,
	    philow, phiup, Ps, ofpar, f);
     hp[part]->Write();
     hp2->Write();
     //likeli->Write();
   }
   
   // ~~~~~~~~~kaon part end~~~~~~~~~~

   f->Close();
   ofpar.close();
   detail.close();

}

//#######define a fit and save function######
void FitAndSave(TH1D **hmass, 
              int &part,
              RooPlot *&xframe,
              RooDataHist* &data_k, 
              RooAddPdf* &sum,
              RooGaussian &gaus,
              RooChebychev &bkg,
              RooRealVar &x,
              RooRealVar &mean,
              RooRealVar &sigma1,
              RooRealVar &co1,
              RooRealVar &signal,
              RooRealVar &background,
              double &philow,
              double &phiup,
              double *Ps,
              ofstream &ofpar,
              TFile *&f
              )
{
   char name[100];
   TCanvas *c2=new TCanvas("c2","likelihood",800,600);
   sprintf(name,"mass_k_%1.2f_to_%1.2f",Ps[part],Ps[part+1]);
   data_k = new RooDataHist(name,"data_k",x,hmass[part]);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
   mean.setVal(m0);
   //sigma.setVal(0.035);
   signal.setVal(3000);
   background.setVal(0);
   co1.setVal(0);
   sum->fitTo(*data_k,Range(philow,phiup));
   //mpeak = mean.getVal();
   //sigma=sigma1.getVal();;
   xframe = x.frame(Title("fit kaon"));
   data_k->plotOn(xframe);
   sum->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
   sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
   xframe->Draw();
   sprintf(name,"%s/mass_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
   c2->Print(name);
   ofpar<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<"\t"<<sigma1.getVal()<<"\t"<<sigma1.getError()<<std::endl;
   //ofpar<<"\t"<<hmass[part]->GetMean()<<std::endl;
   ofpar<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError();
   ofpar<<"\t"<<signal.getVal()/(signal.getVal()+background.getVal())<<std::endl;
   delete data_k;
   //delete xframe;
   delete sum;
   delete c2;

   xframe->Write();
   return;
}
//########### saving function end

