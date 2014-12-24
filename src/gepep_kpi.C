#define gepep_kpi_cxx
#include "gepep_kpi.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TGaxis.h"
#include "TPad.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include "function.h"
#include "TF1.h"
#include "TF2.h"
//#include <fstream>
#include "Pars.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
//#include <iostream>
extern std::string outputdir;
using RooFit::Title;
using RooFit::Name;
using RooFit::Components;
using RooFit::LineStyle;
using RooFit::LineColor;
using RooFit::Range;

extern std::string outputdir;
extern std::vector<double> px1,py1,pz1,px2,py2,pz2;
extern double m0;
extern double mparticle,mparticle2,mparticle3,mparticle4;
extern double sigma;
extern double weight,width;
extern double factor2;

void FitAndSave(TH1D *hmass, 
              int &parti,
	    int &partj,
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
              TFile *&f,
	    std::string suffix=""
              );
void FitThisLoop(TH1D *hmass, 
              int &parti,
	    int &partj,
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
              double *Ps
              );
void SaveThisLoop( 
              int &parti,
	    int &partj,
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
              TFile *&f,
	    std::string suffix
              );
void ResetVars(
              RooRealVar &x,
              RooRealVar &mean,
              RooRealVar &sigma1,
              RooRealVar &co1,
              RooRealVar &signal,
              RooRealVar &background)
{
   x.setVal(1.865);
   mean.setVal(1.865);
   sigma1.setVal(0.0068);
   co1.setVal(0);
   signal.setVal(1200);//event number
   background.setVal(200);
   return;
}

void gepep_kpi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_kpi.C
//      Root > gepep_kpi t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   ofstream ofpar;
   ofpar.open("parkpi.txt",std::ios::app);
   ofstream detail;
   detail.open("detailkpi.txt",std::ios::app);
   detail<<"k- pi+ algrithm: will give factors for pion"<<std::endl;
   ofstream purepar;
   purepar.open("par");
   TFile *f=new TFile("plot_kpi.root","RECREATE");
   // for saving the fit result

   Long64_t nbytes = 0, nb = 0;
   double factor,factorlow,factorup;
   double minimum;
   double minx,miny;
   //std::string tmpstr;
   //ParMap kmap("parkk.txt");
   TH1D *h1 = new TH1D("h1","momentum of kaon",100,0,3.0);
   TH1D *h2 = new TH1D("h2","momentum of kaon",100,0,3.0);
   h2->SetLineColor(2);
   double D0low=1.831;
   double D0up=1.899;
   const int Npart=20;
   double pcut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
   for(int i=0;i<Npart+1;i++){
     double start=0.0;
     double stop =2.0;
     pcut[i] = (stop-start)/Npart*i+start;
   }
   double mk=0.493677;
   double mpi=0.13957018;
   double peakvalue=1.86486;// mD0
   std::vector<std::pair <int,double> > facmap;//it is't a c++ map, but a vector
   std::vector<std::pair <int,int> > partmap;

   // try to use roofit
   RooRealVar x("x","energy",1.865,D0low,D0up,"GeV");
   RooRealVar mean("mean","mean of gaussian",1.865,D0low,D0up);
   RooRealVar sigma1("sigma1","width of gaussian",0.0068,0.006,0.008);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma1);
   RooRealVar co1("co1","coefficient #1",0,-100000.,100000.);
   RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",12000,0,1000000);//event number
   RooRealVar background("background"," ",2000,0,100000);
   RooPlot *xframe;
   RooDataHist *data_kpi;
   RooAddPdf *sum;
 
   TH1D *hmass = new TH1D("hmass","k- pi+ invariant mass",100,D0low,D0up);
   TH1D *hm1   = new TH1D("hm1"  ,"k- pi+ invariant mass",100,D0low,D0up);
   TH1D *hm2   = new TH1D("hm2"  ,"k- pi+ invariant mass",100,D0low,D0up);
   TH2D *h2p   = new TH2D("h2p"  ,"k- pi+ momentum" ,200,0,2,200,0,2);
   h2p->GetXaxis()->SetTitle("pion momentum(GeV)");
   h2p->GetYaxis()->SetTitle("kaon momentum(GeV)");
   TH3D *hthetadis = new TH3D("hthetadis","",100,0,2,100,0,2,100,0,TMath::Pi());
   hthetadis->GetXaxis()->SetTitle("p1");
   hthetadis->GetYaxis()->SetTitle("p2");
   hthetadis->GetZaxis()->SetTitle("#theta");
   TH1D *hthedis = new TH1D("hthedis","",100,0,TMath::Pi());
   hthedis->GetXaxis()->SetTitle("#theta");
   TH1D *hthedis1 = new TH1D("hthedis1","",100,0,TMath::Pi());
   hthedis->GetXaxis()->SetTitle("#theta");
   TH1D *hthedis2 = new TH1D("hthedis2","",100,0,TMath::Pi());
   hthedis->GetXaxis()->SetTitle("#theta");
   //TCanvas *c1= new TCanvas("","",800,600);
   TCanvas *c2= new TCanvas("","",800,600);

   m0 = peakvalue;//mD0
   sigma=0.00687;
   width = 10.*sigma;

   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   TH1D *hmp[Npart][Npart];
   for (int parti=0;parti<Npart;parti++){
   for (int partj=0;partj<Npart;partj++){
     char name[100];
     sprintf(name,"mass_part%d_part%d",parti,partj);
     hmp[parti][partj] = new TH1D(name,name,100,D0low,D0up);
   }
   }

   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
  
     int parti,partj;
     double mass;
     double p1,p2;
     int besidx=0;
     double tmpdeltaold=100;
     double tmpmass;
     for(int i=0; i<npip;i++){
       tmpmass = CalInvMass(mpi,pippx[i],pippy[i],pippz[i],mk,kampx[0],kampy[0],kampz[0]);
       if(fabs(tmpmass-peakvalue)<tmpdeltaold){
         tmpdeltaold = fabs(tmpmass-peakvalue);
         besidx = i;
       }
     }
     //std::cout<<"best index is "<<besidx<<std::endl;
     // total invariant mass, D0 -> k- pi+
     p1=TMath::Sqrt(pippx[besidx]*pippx[besidx]+pippy[besidx]*pippy[besidx]+pippz[besidx]*pippz[besidx]);
     p2=TMath::Sqrt(kampx[0]*kampx[0]+kampy[0]*kampy[0]+kampz[0]*kampz[0]);
     parti = (int)(p1/2.0*Npart);
     partj = (int)(p2/2.0*Npart);
     if (parti>=Npart || partj>=Npart || parti<0 || partj<0) continue;
     mass = CalInvMass(mpi,pippx[besidx],pippy[besidx],pippz[besidx],mk,kampx[0],kampy[0],kampz[0]);
     if(p1>pcut[parti]&&p1<pcut[parti+1]&&p2>pcut[partj]&&p2<pcut[partj+1])
       if (mass>m0-width/2. && mass<m0+width/2.){
         h2p->Fill(p1,p2);
         hmp[parti][partj]->Fill(mass);
       }
     // if (Cut(ientry) < 0) continue;
   }
   h2p->Write();
   for (int parti=0;parti<Npart;parti++)
   for (int partj=0;partj<Npart;partj++){
     hmp[parti][partj]->Write();
     std::cout<<"processed part "<<parti<<", part "<<partj<<std::endl;
     if (hmp[parti][partj]->GetEntries() > 200){
       partmap.push_back(std::make_pair(parti,partj));
     }
   }
   // ~~~~~~~~ draw end

   // ~~~~~~~~~kaon part~~~~~~~~~~
   
   for (int loopi=0;loopi<20;loopi++){
     int parti;
     int partj;
     if (facmap.size()==0 ){
       parti=9;
       partj=9;
     }
     else {
       parti=9;
       partj=loopi;
     }
     if (hmp[parti][partj]->GetEntries()<200) continue;

     double factori=0,factorj=0;
     for (int loopk=0;loopk<facmap.size();loopk++){
       if (facmap.at(loopk).first == parti) factori=facmap.at(loopk).second;
       //if (facmap.at(loopk).first == partj) factorj=facmap.at(loopk).second;
     }
     if(factori==0 && factorj==0 && parti!=partj) continue;
     //if(parti==partj && factori!=0) continue;
     if(factori!=0 && factorj!=0) continue;
     if(factori==0){
       mparticle = mpi;//mpion
       mparticle2= mk;//mkaon
     }
     else {
       mparticle = mk;//
       mparticle2= mpi;//
     }
     FitThisLoop(hmp[parti][partj], parti,partj,data_kpi, sum,gaus,bkg,
              x,mean,sigma1,co1,signal,background,D0low,D0up,pcut);
     double sigNo=signal.getVal();
     double bckNo=width/(D0up-D0low)*background.getVal();
     weight = sigNo/(sigNo+bckNo);
     if (weight<0.2) continue;
     ofpar<<"part "<<parti<<", part "<<partj<<std::endl;
     SaveThisLoop( parti,partj,xframe,data_kpi, sum,gaus,bkg,
              x,mean,sigma1,co1,signal,background,D0low,D0up,pcut,
              ofpar,f,"pre");
     ResetVars(x,mean,sigma1,co1,signal,background);
     
     char tmpchr[100];
     
     //likelihood method
     //factor2 = 0.99751;// factor for kaon
     px1.clear();
     px2.clear();
     py1.clear();
     py2.clear();
     pz1.clear();
     pz2.clear();
     hmass->Reset();
     hthetadis->Reset();
     hthedis->Reset();
     sprintf(tmpchr,"hthetadis_part%d_part%d",parti,partj);
     hthetadis->SetName(tmpchr);
     sprintf(tmpchr,"hthedis_part%d_part%d",parti,partj);
     hthedis->SetName(tmpchr);
     //nb = fChain->GetEntry(1);   nbytes += nb;
     for (Long64_t  jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
    
       double mass;
       double p1,p2;
       int besidx=0;
       double tmpdeltaold=100;
       double tmpmass;
       double theta;
       // search for best pi+
       for(int i=0; i<npip;i++){
         tmpmass = CalInvMass(mpi,pippx[i],pippy[i],pippz[i],mk,kampx[0],kampy[0],kampz[0]);
         if(fabs(tmpmass-m0)<tmpdeltaold){
           tmpdeltaold = fabs(tmpmass-m0);
           besidx = i;
         }
       }
       // total invariant mass, D0 -> k- pi+
       p1=TMath::Sqrt(pippx[besidx]*pippx[besidx]+pippy[besidx]*pippy[besidx]+pippz[besidx]*pippz[besidx]);
       p2=TMath::Sqrt(kampx[0]*kampx[0]+kampy[0]*kampy[0]+kampz[0]*kampz[0]);
       mass = CalInvMass(mpi,pippx[besidx],pippy[besidx],pippz[besidx],mk,kampx[0],kampy[0],kampz[0]);
       if(p1>pcut[parti]&&p1<pcut[parti+1] && p2>pcut[partj]&&p2<pcut[partj+1] ){
         if (mass>m0-width/2. && mass<m0+width/2.){
	 //h2p->Fill(p1,p2);
           theta = acos((pippx[besidx]*kampx[0]
	              +pippx[besidx]*kampx[0]
		    +pippx[besidx]*kampx[0])
		    /(p1*p2));
	 hthetadis->Fill(p1,p2,theta);
	 hthedis->Fill(theta);
	 if (parti==partj){
             px1.push_back(pippx[besidx]);
             px2.push_back(kampx[0]);
             py1.push_back(pippy[besidx]);
             py2.push_back(kampy[0]);
             pz1.push_back(pippz[besidx]);
             pz2.push_back(kampz[0]);
	 }
	 else if (factori==0){
             px1.push_back(pippx[besidx]);
             px2.push_back(kampx[0]*factorj);
             py1.push_back(pippy[besidx]);
             py2.push_back(kampy[0]*factorj);
             pz1.push_back(pippz[besidx]);
             pz2.push_back(kampz[0]*factorj);
	 }
	 else{
             px2.push_back(pippx[besidx]*factori);
             px1.push_back(kampx[0]);
             py2.push_back(pippy[besidx]*factori);
             py1.push_back(kampy[0]);
             pz2.push_back(pippz[besidx]*factori);
             pz1.push_back(kampz[0]);
	 }
         }
       }
     }
     hthetadis->Write();
     hthedis->Write();

     //h2p->Draw();
     //sprintf(tmpchr,"%s/momentum_kpi_part%d.eps",outputdir.c_str(),part);
     //c2->Print(tmpchr);
     std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
     //TF2 *likeli=new TF2("likeli",maxlikelihood2_0,0.95,1.05,0.01,0.99);
     //likeli->SetNpy(100);
     //likeli->Draw("surf1");
     //tmpstr=outputdir+"/likelikpi_2D.eps";
     //c2->Print(tmpstr.c_str());
     //minimum = likeli->GetMinimum(0.98,1.02);
     //factor=likeli->GetMinimumX(0.98,1.02);
     //likeli->GetMinimumXY(factor,miny);
     //weight = miny;
     
     TF1 *likeli_1;
     if (parti==partj)
       likeli_1=new TF1("likeli_1",maxlikelihood2_1,0.95,1.05);
     else{
       factor2=1;
       likeli_1=new TF1("likeli_1",maxlikelihood2,0.95,1.05);
     }
     likeli_1->Draw();
     sprintf(tmpchr,"%s/likelikpi_part%d_part%d.eps",outputdir.c_str(),parti,partj);
     c2->Print(tmpchr);
     minimum = likeli_1->GetMinimum(0.98,1.02);
     factor = likeli_1->GetMinimumX(0.98,1.02);
     factorlow=likeli_1->GetX(minimum+1,0.98,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.02);
     ofpar<<run<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t"<<weight<<std::endl;
     detail<<"weight is "<<miny<<", factor is "<<likeli_1->GetMinimumX(0.98,1.01)<<std::endl;
     //detail<<"minimum 2D  "<<likeli->GetMinimum()<<", minimum 1D "<<likeli_1->GetMinimum()<<std::endl;
     detail<<"best factor  "<<likeli_1->GetMinimumX(0.99,1.01)<<std::endl;
     purepar<<factor<<"\t"<<(factorup-factorlow)/2<<"\t";

     double index;
     if (parti!=partj)
       factori==0? index=parti : index=partj;
     else {index = parti;}
     facmap.push_back(std::make_pair(index,factor));

     // use the factor to refit
     hmass->Reset();
     double factors[2]={ factor, 1};
     if (parti==partj)
       factors[1]=factor;
     for (int vi=0; vi<px1.size();vi++){
       double mass;
       mass = CalInvMass(mparticle,px1.at(vi),py1.at(vi),pz1.at(vi),mparticle2,px2.at(vi),py2.at(vi),pz2.at(vi),-2,factors);
       hmass->Fill(mass);
     }
     FitAndSave(hmass, parti, partj, 
              xframe, data_kpi, sum, gaus, bkg,
	    x, mean, sigma1, co1, signal, background,
	    D0low, D0up, pcut, ofpar, f,"re");
   }
   
   // ~~~~~~~~~kaon part end~~~~~~
   // ~~~~~~~~~pion part~~~~~~~~~~~~
   /*
   for (int loopi=0;loopi<partmap.size();loopi++){
   for (int loopj=0;loopj<partmap.size();loopj++){
     int parti=partmap.at(loopj).first;
     int partj=partmap.at(loopj).second;
     double factori=0,factorj=0;
     for (int loopk=0;loopk<facmap.size();loopk++){
       if (facmap.at(loopk).first == parti) factori=facmap.at(loopk).second;
       if (facmap.at(loopk).first == partj) factorj=facmap.at(loopk).second;
     }
     if(factori==0 && factorj==0 && parti!=partj) continue;
     if(parti==partj && factori!=0) continue;
     if(factori!=0 && factorj!=0) continue;
     if(factori==0){
       mparticle = mpi;//mpion
       mparticle2= mk;//mkaon
     }
     else {
       mparticle = mk;//
       mparticle2= mpi;//
     }
     FitThisLoop(hmp[parti][partj], parti,partj,data_kpi, sum,gaus,bkg,
              x,mean,sigma1,co1,signal,background,D0low,D0up,pcut);
     double sigNo=signal.getVal();
     double bckNo=width/(D0up-D0low)*background.getVal();
     weight = sigNo/(sigNo+bckNo);
     if (weight<0.2) continue;
     ofpar<<"part "<<parti<<", part "<<partj<<std::endl;
     SaveThisLoop( parti,partj,xframe,data_kpi, sum,gaus,bkg,
              x,mean,sigma1,co1,signal,background,D0low,D0up,pcut,
              ofpar,f,"pre");
     ResetVars(x,mean,sigma1,co1,signal,background);
     
     char tmpchr[100];
     
     //likelihood method
     //factor2 = 0.99751;// factor for kaon
     px1.clear();
     px2.clear();
     py1.clear();
     py2.clear();
     pz1.clear();
     pz2.clear();
     hmass->Reset();
     hthetadis->Reset();
     hthedis->Reset();
     sprintf(tmpchr,"hthetadis_part%d_part%d",parti,partj);
     hthetadis->SetName(tmpchr);
     sprintf(tmpchr,"hthedis_part%d_part%d",parti,partj);
     hthedis->SetName(tmpchr);
     //nb = fChain->GetEntry(1);   nbytes += nb;
     for (Long64_t  jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
    
       double mass;
       double p1,p2;
       int besidx=0;
       double tmpdeltaold=100;
       double tmpmass;
       double theta;
       // search for best pi+
       for(int i=0; i<npip;i++){
         tmpmass = CalInvMass(mpi,pippx[i],pippy[i],pippz[i],mk,kampx[0],kampy[0],kampz[0]);
         if(fabs(tmpmass-m0)<tmpdeltaold){
           tmpdeltaold = fabs(tmpmass-m0);
           besidx = i;
         }
       }
       // total invariant mass, D0 -> k- pi+
       p1=TMath::Sqrt(pippx[besidx]*pippx[besidx]+pippy[besidx]*pippy[besidx]+pippz[besidx]*pippz[besidx]);
       p2=TMath::Sqrt(kampx[0]*kampx[0]+kampy[0]*kampy[0]+kampz[0]*kampz[0]);
       mass = CalInvMass(mpi,pippx[besidx],pippy[besidx],pippz[besidx],mk,kampx[0],kampy[0],kampz[0]);
       if(p1>pcut[parti]&&p1<pcut[parti+1] && p2>pcut[partj]&&p2<pcut[partj+1] ){
         if (mass>m0-width/2. && mass<m0+width/2.){
	 //h2p->Fill(p1,p2);
           theta = acos((pippx[besidx]*kampx[0]
	              +pippx[besidx]*kampx[0]
		    +pippx[besidx]*kampx[0])
		    /(p1*p2));
	 hthetadis->Fill(p1,p2,theta);
	 hthedis->Fill(theta);
	 if (parti==partj){
             px1.push_back(pippx[besidx]);
             px2.push_back(kampx[0]);
             py1.push_back(pippy[besidx]);
             py2.push_back(kampy[0]);
             pz1.push_back(pippz[besidx]);
             pz2.push_back(kampz[0]);
	 }
	 else if (factori==0){
             px1.push_back(pippx[besidx]);
             px2.push_back(kampx[0]*factorj);
             py1.push_back(pippy[besidx]);
             py2.push_back(kampy[0]*factorj);
             pz1.push_back(pippz[besidx]);
             pz2.push_back(kampz[0]*factorj);
	 }
	 else{
             px2.push_back(pippx[besidx]*factori);
             px1.push_back(kampx[0]);
             py2.push_back(pippy[besidx]*factori);
             py1.push_back(kampy[0]);
             pz2.push_back(pippz[besidx]*factori);
             pz1.push_back(kampz[0]);
	 }
         }
       }
     }
     hthetadis->Write();
     hthedis->Write();

     //h2p->Draw();
     //sprintf(tmpchr,"%s/momentum_kpi_part%d.eps",outputdir.c_str(),part);
     //c2->Print(tmpchr);
     std::cout<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
     //TF2 *likeli=new TF2("likeli",maxlikelihood2_0,0.95,1.05,0.01,0.99);
     //likeli->SetNpy(100);
     //likeli->Draw("surf1");
     //tmpstr=outputdir+"/likelikpi_2D.eps";
     //c2->Print(tmpstr.c_str());
     //minimum = likeli->GetMinimum(0.98,1.02);
     //factor=likeli->GetMinimumX(0.98,1.02);
     //likeli->GetMinimumXY(factor,miny);
     //weight = miny;
     
     TF1 *likeli_1;
     if (parti==partj)
       likeli_1=new TF1("likeli_1",maxlikelihood2_1,0.95,1.05);
     else{
       factor2=1;
       likeli_1=new TF1("likeli_1",maxlikelihood2,0.95,1.05);
     }
     likeli_1->Draw();
     sprintf(tmpchr,"%s/likelikpi_part%d_part%d.eps",outputdir.c_str(),parti,partj);
     c2->Print(tmpchr);
     minimum = likeli_1->GetMinimum(0.98,1.02);
     factor = likeli_1->GetMinimumX(0.98,1.02);
     factorlow=likeli_1->GetX(minimum+1,0.98,factor);
     factorup =likeli_1->GetX(minimum+1,factor,1.02);
     ofpar<<run<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<"\t"<<weight<<std::endl;
     detail<<"weight is "<<miny<<", factor is "<<likeli_1->GetMinimumX(0.98,1.01)<<std::endl;
     //detail<<"minimum 2D  "<<likeli->GetMinimum()<<", minimum 1D "<<likeli_1->GetMinimum()<<std::endl;
     detail<<"best factor  "<<likeli_1->GetMinimumX(0.99,1.01)<<std::endl;
     purepar<<factor<<"\t"<<(factorup-factorlow)/2<<"\t";

     double index;
     if (parti!=partj)
       factori==0? index=parti : index=partj;
     else {index = parti;}
     facmap.push_back(std::make_pair(index,factor));

     // use the factor to refit
     hmass->Reset();
     double factors[2]={ factor, 1};
     if (parti==partj)
       factors[1]=factor;
     for (int vi=0; vi<px1.size();vi++){
       double mass;
       mass = CalInvMass(mparticle,px1.at(vi),py1.at(vi),pz1.at(vi),mparticle2,px2.at(vi),py2.at(vi),pz2.at(vi),-2,factors);
       hmass->Fill(mass);
     }
     FitAndSave(hmass, parti, partj, 
              xframe, data_kpi, sum, gaus, bkg,
	    x, mean, sigma1, co1, signal, background,
	    D0low, D0up, pcut, ofpar, f,"re");
   }
   }*/
   // ~~~~~~~~~pion part end~~~~~~~~~~
   //TFile f("plot.root","recreate");
   f->Close();
   ofpar.close();
   detail.close();

}

//#######define a fit and save function######
void FitAndSave(TH1D *hmass, 
              int &parti,
	    int &partj,
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
              TFile *&f,
	    std::string suffix
              )
{
  FitThisLoop(hmass, parti,partj,data_k, sum,gaus,bkg,
              x,mean,sigma1,co1,signal,background,philow,phiup,Ps);
  SaveThisLoop( parti,partj,xframe,data_k, sum,gaus,bkg,
              x,mean,sigma1,co1,signal,background,philow,phiup,Ps,
              ofpar,f,suffix);
  return;
}

void FitThisLoop(TH1D *hmass, 
              int &parti,
	    int &partj,
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
              double *Ps
              )
{
   char name[100];
   sprintf(name,"mass_kpi_part%d_part%d",parti,partj);
   data_k = new RooDataHist(name,"data_k",x,hmass);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
   mean.setVal(m0);
   co1.setVal(0);
   sum->fitTo(*data_k,Range(philow,phiup));

   return;
}

void SaveThisLoop( 
              int &parti,
	    int &partj,
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
              TFile *&f,
	    std::string suffix
              )
{
   if (suffix != "") suffix = "_"+ suffix;
   char name[100];
   TCanvas *c2=new TCanvas("c2","likelihood",800,600);

   sprintf(name,"part%d_part%d%s",parti,partj,suffix.c_str());
   xframe = x.frame(Title("fit K pi"));
   xframe->SetName(name);
   data_k->plotOn(xframe);
   sum->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
   xframe->Draw();
   sprintf(name,"%s/masskpi_part%d_part%d%s.eps",outputdir.c_str(),parti,partj,suffix.c_str());
   c2->Print(name);
   ofpar<<suffix<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<"\t"<<sigma1.getVal()<<"\t"<<sigma1.getError()<<std::endl;
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

#ifdef gepep_kpi_cxx
gepep_kpi::gepep_kpi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/RValue_kpi_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/RValue_kpi_3850.root");
      }
      f->GetObject("gepep_kpi",tree);

   }
   Init(tree);
}

gepep_kpi::~gepep_kpi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_kpi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_kpi::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gepep_kpi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("rec", &rec, &b_rec);
   fChain->SetBranchAddress("evttag", &evttag, &b_evttag);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("nkap", &nkap, &b_nkap);
   fChain->SetBranchAddress("kappx", kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", kappz, &b_kappz);
   fChain->SetBranchAddress("kape", kape, &b_kape);
   fChain->SetBranchAddress("nkam", &nkam, &b_nkam);
   fChain->SetBranchAddress("kampx", kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", kampz, &b_kampz);
   fChain->SetBranchAddress("kame", kame, &b_kame);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_kpi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_kpi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_kpi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_kpi_cxx
