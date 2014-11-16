#include <iostream>
#include <fstream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
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
extern double m0;
extern double mparticle,mparticle2,mparticle3,mparticle4;
extern double sigma;
extern double width;
extern double weight;
extern double slope;

int main(int argc,char** argv)
{
   std::string outputdir=".";

   ifstream sample("sample.txt");
   ofstream ofpar;
   ofpar.open("parkk.txt",std::ios::app);
   ofstream detail;
   detail.open("detailkk.txt",std::ios::app);
   detail<<"kk algrithm: will give factors for kaon"<<std::endl;
   
   double pxa,pya,pza,pxb,pyb,pzb;
   double philow=1.0;
   double phiup=1.05;
   // try to use roofit
   RooRealVar x("x","energy",1.020,philow,phiup,"GeV");
   RooRealVar mean("mean","mean of gaussian",1.020,philow,phiup);
   RooRealVar sigma1("sigma1","width of gaussian",0.003,0.0001,0.005);
   //RooRealVar sigma2("sigma2","width of gaussian",0.02,0.005,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma1);
   //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar co1("co1","coefficient #1",0,-1000.,1000.);
   //RooRealVar co4("co4","coefficient #4",0);
   RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",300,10,1000000);//event number
   //RooRealVar signal2("signal2"," ",1200,10,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
   RooPlot *xframe;
   RooDataHist *data_k;
   RooAddPdf *sum;
 
   int NP=1;
   // split momentum from 0.13 to 1.13 GeV
   double Ps[NP+1];
   for(int i=0;i<=NP;i++){
     Ps[i]=1.0/NP*i+0.13;
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
     hp[i]=new TH1D(name,name,100,Ps[i],Ps[i+1]);
     sprintf(name,"hmass_%1.2f_to_%1.2f_GeV",Ps[i],Ps[i+1]);
     hmass[i]=new TH1D(name,name,100,philow,phiup);
   }
   //TH1D *h1 = new TH1D("h1","momentum of kaon",100,0,3.0);
   //h2->SetLineColor(2);
   // ~~~~~~~~~kaon part~~~~~~~~~~
  
   // likelihood method
   //m0 = 1.01946;
   double mpeak;
   m0 = 1.019455;
   //sigma=0.0024;
   sigma=sigma1.getVal();
   width = 20.*sigma;
   mparticle=0.493677;
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
       double totpx,totpy,totpz,tote;
       double ke[2];
       double kapp,kamp;
       // total invariant mass
       totpx=(pxa+pxb);
       totpy=(pya+pyb);
       totpz=(pza+pzb);
       ke[0]=TMath::Sqrt(mparticle*mparticle + 
             (pxa*pxa+pya*pya+pza*pza));
       ke[1]=TMath::Sqrt(mparticle*mparticle + 
             (pxb*pxb+pyb*pyb+pzb*pzb));
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       if(!(kapp>Ps[part] && kapp<Ps[part+1] && kamp>Ps[part] && kamp<Ps[part+1])) continue;
       tote=ke[0]+ke[1];
       mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
       hmass[part]->Fill(mass);
       // if (Cut(ientry) < 0) continue;
     }

     char tmpchr[100];
     char name[100];
     TCanvas *c2=new TCanvas("c2","likelihood",800,600);
     sprintf(tmpchr,"mass_k_%1.2f_to_%1.2f",Ps[part],Ps[part+1]);
     data_k = new RooDataHist(tmpchr,"data_k",x,hmass[part]);
     sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
     mean.setVal(m0);
     //sigma.setVal(0.035);
     signal.setVal(300);
     background.setVal(200);
     co1.setVal(0);
     sum->fitTo(*data_k,Range(philow,phiup));
     mpeak = mean.getVal();
     sigma=sigma1.getVal();;
     xframe = x.frame(Title("fit kaon"));
     data_k->plotOn(xframe);
     sum->plotOn(xframe);
     sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
     //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
     sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
     xframe->Draw();
     sprintf(name,"%s/mass_%1.2f_%1.2f_pre.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     ofpar<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<"\t"<<sigma1.getVal()<<"\t"<<sigma1.getError()<<std::endl;
     ofpar<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError();
     ofpar<<"\t"<<signal.getVal()/(signal.getVal()+background.getVal())<<std::endl;
     delete data_k;
     //delete xframe;
     delete sum;
     // likelihood
     px1.clear();
     px2.clear();
     py1.clear();
     py2.clear();
     pz1.clear();
     pz2.clear();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string tmpcha;
       getline(sample,tmpcha);
        //if(ngam>0) continue;
       double mass;
       double totpx,totpy,totpz,tote;
       double kapp,kamp,kape,kame;
       // total invariant mass
       totpx=(pxa+pxb);
       totpy=(pya+pyb);
       totpz=(pza+pzb);
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       kape=TMath::Sqrt(mparticle*mparticle+kapp*kapp);
       kame=TMath::Sqrt(mparticle*mparticle+kamp*kamp);
       tote=kape+kame;
       if(!(kapp>Ps[part] && kapp<Ps[part+1] && kamp>Ps[part] && kamp<Ps[part+1])) continue;
       //h1->Fill(kapp);
       mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
       if(mass>mpeak-width/2. && mass<mpeak+width/2.){
         hp[part]->Fill(kapp);
	 //hmass[part]->Fill(mass);
         px1.push_back(pxa);
         px2.push_back(pxb);
         py1.push_back(pya);
         py2.push_back(pyb);
         pz1.push_back(pza);
         pz2.push_back(pzb);
       }
     }
  
     hp[part]->Draw();
     sprintf(name,"%s/momentum_%1.2f_%1.2f.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
     TF2 *likeli=new TF2("likeli",maxlikelihood1,0.95,1.05,0.01,0.5);
     //TF3 *likeli=new TF3("likeli",maxlikelihood1_3,0.95,1.05,0.01,0.99,-100,100);
     //likeli->Draw();
     //sprintf(name,"%s/momentum_%1.2f_%1.2f_2D.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     //c2->Print(name);
     //likeli->Draw("surf2");
     //sprintf(name,"%s/momentum_%1.2f_%1.2f_2D2.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     //c2->Print(name);
     likeli->GetMinimumXY(factor,miny);
     weight = miny;
     //slope = minz;
     TF1 *likeli_1=new TF1("likeli_1",maxlikelihood1_1,0.95,1.05);
     likeli_1->Draw();
     sprintf(name,"%s/momentum_%1.2f_%1.2f_1D.eps",outputdir.c_str(),Ps[part],Ps[part+1]);
     c2->Print(name);
     minimum = likeli_1->GetMinimum(0.98,1.02);
     factorlow=likeli_1->GetX(minimum+1,0.98,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.02);
     ofpar<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t\t"<<weight<<"\t"<<slope<<std::endl;
     detail<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
     detail<<"signal weight is "<<weight<<" best factor  "<<likeli_1->GetMinimumX(0.99,1.01)<<std::endl;

     // using the factor to fit
     hmass[part]->Reset();
     sample.clear();
     sample.seekg(0,std::ios::beg);
     while(!sample.eof()){
       sample>>pxa>>pya>>pza>>pxb>>pyb>>pzb;
       std::string tmpcha;
       getline(sample,tmpcha);
      
       //if(ngam>0) continue;
       double mass;
       double totpx,totpy,totpz,tote;
       double ke[2];
       double kapp,kamp;
       // total invariant mass
       totpx=factor*(pxa+pxb);
       totpy=factor*(pya+pyb);
       totpz=factor*(pza+pzb);
       ke[0]=TMath::Sqrt(mparticle*mparticle + 
             factor*factor*(pxa*pxa+pya*pya+pza*pza));
       ke[1]=TMath::Sqrt(mparticle*mparticle + 
             factor*factor*(pxb*pxb+pyb*pyb+pzb*pzb));
       kapp=TMath::Sqrt(pxa*pxa+pya*pya+pza*pza);
       kamp=TMath::Sqrt(pxb*pxb+pyb*pyb+pzb*pzb);
       if(!(kapp>Ps[part] && kapp<Ps[part+1] && kamp>Ps[part] && kamp<Ps[part+1])) continue;
       tote=ke[0]+ke[1];
       mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
       hmass[part]->Fill(mass);
       // if (Cut(ientry) < 0) continue;
     }

     sprintf(tmpchr,"mass_k_%1.2f_to_%1.2f",Ps[part],Ps[part+1]);
     data_k = new RooDataHist(tmpchr,"data_k",x,hmass[part]);
     sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
     mean.setVal(m0);
     //sigma.setVal(0.035);
     signal.setVal(300);
     background.setVal(200);
     co1.setVal(0);
     sum->fitTo(*data_k,Range(philow,phiup));
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
     ofpar<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError();
     ofpar<<"\t"<<signal.getVal()/(signal.getVal()+background.getVal())<<std::endl;
     delete data_k;
     delete xframe;
     delete sum;
     delete c2;
   }
   // ~~~~~~~~~kaon part end~~~~~~~~~~

   ofpar.close();
   detail.close();

}
