#define gepep_fastpipill_cxx
#include "gepep_fastpipill.h"
#include <TH1D.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
//#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include "function.h"
#include "bes3plotstyle.h"
#include "TLegend.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
extern std::string outputdir;
using RooFit::Title;
using RooFit::Components;
using RooFit::LineStyle;
using RooFit::LineColor;
using RooFit::Range;
extern std::string outputdir;
extern std::vector<double> px1,py1,pz1,px2,py2,pz2;
extern std::vector<double> px3,py3,pz3,px4,py4,pz4;
extern std::vector<double> le1,le2;
extern double m0;
extern double mparticle,mparticle2,mparticle3,mparticle4;
extern double factor2,factor3,factor4;
extern double sigma,weight,width;
 
bool gepep_fastpipill::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fastpipill.C
//      Root > gepep_fastpipill t
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
   if (fChain == 0) return false;
   Long64_t nentries = fChain->GetEntriesFast();

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   ofstream ofpare;
   ofpare.open("parpipille.txt",std::ios::app);
   ofstream ofparmu;
   ofparmu.open("parpipillmu.txt",std::ios::app);
   ofstream ofparpi;
   ofparpi.open("parpipillpi.txt",std::ios::app);
   ofstream detail;
   detail.open("detail.txt",std::ios::app);
   detail<<"fastpipi algrithm: will give factors for e,mu,pi"<<std::endl;
   ofstream purepar;
   purepar.open("par");
   
   double jlow=3.0;
   double jup=3.2;
   double psilow=3.675;
   double psiup=3.695;
   // try to use roofit
   RooRealVar x("x","energy",3.097,jlow,jup,"GeV");
   RooRealVar mean("mean","mean of gaussian",3.097,jlow,jup);
   RooRealVar sigma1("sigma1","width of gaussian",0.003,0.0001,0.05);
   //RooRealVar sigma2("sigma2","width of gaussian",0.02,0.005,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma1);
   //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar co1("co1","coefficient #1",0,-1000.,1000.);
   //RooRealVar co4("co4","coefficient #4",0);
   RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",1200,10,1000000);//event number
   //RooRealVar signal2("signal2"," ",1200,10,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
   RooPlot *xframe;
   RooDataHist *data;
   RooAddPdf *sum;
   
   Long64_t nbytes = 0, nb = 0;
   double factor,factorlow,factorup;
   double minimum;
   double minx,miny;
   TH1D *hmass1 = new TH1D("hmassje","mass j",40,jlow,jup);
   TH1D *hmass2 = new TH1D("hmassjmu","mass j",40,jlow,jup);
   TH1D *hmass3 = new TH1D("hmasspsi","mass psi",50,psilow,psiup);
   TH1D *hp     = new TH1D("hp","momentum distrubution",200,0,2);

   std::string tmpstr;

   // ~~~~~~~~~electron part~~~~~~~~~~
   m0 = 3.096916;
   sigma=0.0169;
   mparticle=0.000511;
   weight = 1.0;
   width = 10.*sigma;
   // pre fit
   ofpare<<"prefit"<<"\n";
 
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      if(cos(angle4)>0.90) continue; // cut bhabha 
      if(decay_ee==1){
        double mass;
        double totpx,totpy,totpz,tote;
        double lee[2];
        double le1p,le2p;
        // total invariant mass
        totpx=(lepx4[0]+lepx4[1]);
        totpy=(lepy4[0]+lepy4[1]);
        totpz=(lepz4[0]+lepz4[1]);
        lee[0]=TMath::Sqrt(mparticle*mparticle + 
               (lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
        lee[1]=TMath::Sqrt(mparticle*mparticle +
               (lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
        le1p=TMath::Sqrt(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]);
        le2p=TMath::Sqrt(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]);
        tote=lee[0]+lee[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	if(mass>m0-width/2. && mass<m0+width/2.)
	{
	  hmass1->Fill(mass);
	}
      }
     // total invariant mass
     // if (Cut(ientry) < 0) continue;
    }
    
    TCanvas *c2=new TCanvas("c2","likelihood",800,600);
    char name[100];
    sprintf(name,"mass_jpsi_e");
    data = new RooDataHist(name,"data_e",x,hmass1);
    mean.setVal(m0);
    sigma1.setVal(sigma);
    signal.setVal(120);
    background.setVal(10);
    co1.setVal(0);
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
    sum->fitTo(*data,Range(jlow,jup));
    xframe = x.frame(Title("fit e"));
    data->plotOn(xframe);
    sum->plotOn(xframe);
    sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
    //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
    sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
    xframe->Draw();
    sprintf(name,"%s/mass_jpsi_e_pre.eps",outputdir.c_str());
    c2->Print(name);
    ofpare<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
    ofpare<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError()<<std::endl;
    delete data;
    delete xframe;
    delete sum;
    //delete c2;

  // likelihood
   px1.clear();
   px2.clear();
   py1.clear();
   py2.clear();
   pz1.clear();
   pz2.clear();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(cos(angle4)>0.90) continue; // cut bhabha 
      if(decay_ee==1){
        double mass;
        double totpx,totpy,totpz,tote;
        double p1,p2;
        double lee[2];
        // total invariant mass
        totpx=(lepx4[0]+lepx4[1]);
        totpy=(lepy4[0]+lepy4[1]);
        totpz=(lepz4[0]+lepz4[1]);
        p1=TMath::Sqrt(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]);
        p2=TMath::Sqrt(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]);
        lee[0]=TMath::Sqrt(mparticle*mparticle + p1*p1);
        lee[1]=TMath::Sqrt(mparticle*mparticle + p2*p2);
        tote=lee[0]+lee[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
        if(mass>m0-width/2. && mass<m0+width/2.){
          px1.push_back(lepx4[0]);
          px2.push_back(lepx4[1]);
          py1.push_back(lepy4[0]);
          py2.push_back(lepy4[1]);
          pz1.push_back(lepz4[0]);
          pz2.push_back(lepz4[1]);
        }
      }
      // if (Cut(ientry) < 0) continue;
   }
   
   //factor =0;
   detail<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
   //TF2 *likeli=new TF2("likeli",maxlikelihood1,0.95,1.05,0.1,0.99);
   //likeli->Draw("surf1");
   //sprintf(name,"%s/likelie_2D.eps",outputdir.c_str());
   //c2->Print(name);
   //likeli->GetMinimumXY(factor,weight);
   double sigNo=signal.getVal() ;
   double bckNo=width/(jup-jlow)*background.getVal();
   weight = sigNo/(sigNo+bckNo);
   TF1 *likeli_1=new TF1("likeli_1",maxlikelihood1_1,0.95,1.05);
   likeli_1->Draw();
   sprintf(name,"%s/likelie_1D.eps",outputdir.c_str());
   c2->Print(name);
   factor = likeli_1->GetMinimumX(0.98,1.02);
   minimum = likeli_1->GetMinimum(0.98,1.02);
   factorlow=likeli_1->GetX(minimum+1,0.98,factor);
   factorup =likeli_1->GetX(minimum+1,factor,1.02);
   ofpare<<run<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
   ofpare<<"\t"<<factor<<"\t"<<weight<<"\t"<<likeli_1->GetMinimumX(0.98,1.02)<<std::endl;
   detail<<"weight is "<<weight<<", factor is "<<likeli_1->GetMinimumX(0.98,1.02)<<std::endl;
   detail<<"best factor  "<<factor<<std::endl;
   purepar<<factor<<"\t"<<(factorup-factorlow)/2<<"\t";

   // using the factor to fit
   hmass1->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      if(cos(angle4)>0.90) continue; // cut bhabha 
      if(decay_ee==1){
        double mass;
        double totpx,totpy,totpz,tote;
        double lee[2];
        double le1p,le2p;
        // total invariant mass
        totpx=factor*(lepx4[0]+lepx4[1]);
        totpy=factor*(lepy4[0]+lepy4[1]);
        totpz=factor*(lepz4[0]+lepz4[1]);
        lee[0]=TMath::Sqrt(mparticle*mparticle + 
               factor*factor*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
        lee[1]=TMath::Sqrt(mparticle*mparticle +
               factor*factor*(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
        le1p=TMath::Sqrt(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]);
        le2p=TMath::Sqrt(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]);
        tote=lee[0]+lee[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	if(mass>m0-width/2. && mass<m0+width/2.){
	  hmass1->Fill(mass);
	}
      }
     // total invariant mass
     // if (Cut(ientry) < 0) continue;
    }

    sprintf(name,"mass_jpsi_e");
    data = new RooDataHist(name,"data_e",x,hmass1);
    mean.setVal(m0);
    sigma1.setVal(sigma);
    signal.setVal(120);
    background.setVal(10);
    co1.setVal(0);
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
    sum->fitTo(*data,Range(jlow,jup));
    xframe = x.frame(Title("fit e"));
    data->plotOn(xframe);
    sum->plotOn(xframe);
    sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
    //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
    sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
    xframe->Draw();
    sprintf(name,"%s/mass_jpsi_e.eps",outputdir.c_str());
    c2->Print(name);
    ofpare<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
    ofpare<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError()<<std::endl;
    delete data;
    delete xframe;
    delete sum;
    //delete c2;

   // ~~~~~~~~~electron part end~~~~~~~~~~

   // ~~~~~~~~~muon part~~~~~~~~~~
   m0 = 3.096916;
   sigma=0.0166;
   mparticle=0.105658;
   weight = 0.9;
   width = 10.*sigma;
  
   // pre fit
   ofparmu<<"prefit"<<"\n";
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      if(cos(angle4)>0.90) continue; // cut bhabha 
      if(decay_ee==0){
        double mass;
        double totpx,totpy,totpz,tote;
        double lee[2];
        double le1p,le2p;
        // total invariant mass
        totpx=(lepx4[0]+lepx4[1]);
        totpy=(lepy4[0]+lepy4[1]);
        totpz=(lepz4[0]+lepz4[1]);
        lee[0]=TMath::Sqrt(mparticle*mparticle + 
               (lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
        lee[1]=TMath::Sqrt(mparticle*mparticle +
               (lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
        le1p=TMath::Sqrt(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]);
        le2p=TMath::Sqrt(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]);
        tote=lee[0]+lee[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	if(mass>m0-width/2. && mass<m0+width/2.){
	  hmass2->Fill(mass);
	}
      }
     // total invariant mass
     // if (Cut(ientry) < 0) continue;
    }

    sprintf(name,"mass_jpsi_e");
    data = new RooDataHist(name,"data_e",x,hmass2);
    mean.setVal(m0);
    sigma1.setVal(sigma);
    signal.setVal(120);
    background.setVal(20);
    co1.setVal(0);
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
    sum->fitTo(*data,Range(jlow,jup));
    xframe = x.frame(Title("fit mu"));
    data->plotOn(xframe);
    sum->plotOn(xframe);
    sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
    //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
    sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
    xframe->Draw();
    sprintf(name,"%s/mass_jpsi_mu_pre.eps",outputdir.c_str());
    c2->Print(name);
    ofparmu<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
    ofparmu<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError()<<std::endl;
    delete data;
    delete xframe;
    delete sum;
    //delete c2;

  px1.clear();
   px2.clear();
   py1.clear();
   py2.clear();
   pz1.clear();
   pz2.clear();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(cos(angle4)>0.90) continue; // cut bhabha 
      if(decay_ee==0){
        double mass;
        double totpx,totpy,totpz,tote;
        double lee[2];
        // total invariant mass
        totpx=(lepx4[0]+lepx4[1]);
        totpy=(lepy4[0]+lepy4[1]);
        totpz=(lepz4[0]+lepz4[1]);
        lee[0]=TMath::Sqrt(mparticle*mparticle + 
               (lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
        lee[1]=TMath::Sqrt(mparticle*mparticle +
               (lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
        tote=lee[0]+lee[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	if(mass>3.0 && mass<3.2){
	  px1.push_back(lepx4[0]);
	  px2.push_back(lepx4[1]);
	  py1.push_back(lepy4[0]);
	  py2.push_back(lepy4[1]);
	  pz1.push_back(lepz4[0]);
	  pz2.push_back(lepz4[1]);
	}
      }
      // if (Cut(ientry) < 0) continue;
   }
   
   detail<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
   TF2 *likeli2=new TF2("likeli2",maxlikelihood1,0.95,1.05,0.01,0.99);
   likeli2->Draw("surf1");
   //tmpstr = outputdir+"/likelimu_2D.eps";
   sprintf(name,"%s/likelimu_2D.eps",outputdir.c_str());
   c2->Print(name);
   minimum = likeli2->GetMinimum(0.98,1.02);
   //factor=likeli2->GetMinimumX(0.98,1.02);
   likeli2->GetMinimumXY(factor,weight);
   //weight = miny;
   TF1 *likeli2_1=new TF1("likeli2_1",maxlikelihood1_1,0.95,1.05);
   likeli2_1->Draw();
   //tmpstr = outputdir+"/likelimu_1D.eps";
   sprintf(name,"%s/likelimu_1D.eps",outputdir.c_str());
   c2->Print(name);
   factorlow=likeli2_1->GetX(minimum+1,0.98,factor);
   //factorlow=likeli2->GetX(miny,0.98,factor);
   factorup =likeli2_1->GetX(minimum+1,factor,1.02);
   ofparmu<<run<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
   ofparmu<<"\t"<<factor<<"\t"<<weight<<"\t"<<likeli2_1->GetMinimumX(0.98,1.02)<<std::endl;
   detail<<"weight is "<<weight<<", factor is "<<likeli2_1->GetMinimumX(0.98,1.02)<<std::endl;
   detail<<"best factor  "<<factor<<std::endl;
   purepar<<factor<<"\t"<<(factorup-factorlow)/2<<"\t";
 
   // using the factor to fit
   hmass2->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      if(cos(angle4)>0.90) continue; // cut bhabha 
      if(decay_ee==0){
        double mass;
        double totpx,totpy,totpz,tote;
        double lee[2];
        double le1p,le2p;
        // total invariant mass
        totpx=factor*(lepx4[0]+lepx4[1]);
        totpy=factor*(lepy4[0]+lepy4[1]);
        totpz=factor*(lepz4[0]+lepz4[1]);
        lee[0]=TMath::Sqrt(mparticle*mparticle + 
               factor*factor*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
        lee[1]=TMath::Sqrt(mparticle*mparticle +
               factor*factor*(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
        le1p=TMath::Sqrt(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]);
        le2p=TMath::Sqrt(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]);
        tote=lee[0]+lee[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	if(mass>m0-width/2. && mass<m0+width/2.){
	  hmass2->Fill(mass);
	}
      }
     // total invariant mass
     // if (Cut(ientry) < 0) continue;
    }

    sprintf(name,"mass_jpsi_e");
    data = new RooDataHist(name,"data_e",x,hmass2);
    mean.setVal(m0);
    sigma1.setVal(sigma);
    signal.setVal(120);
    background.setVal(20);
    co1.setVal(0);
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
    sum->fitTo(*data,Range(jlow,jup));
    xframe = x.frame(Title("fit mu"));
    data->plotOn(xframe);
    sum->plotOn(xframe);
    sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
    //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
    sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
    xframe->Draw();
    sprintf(name,"%s/mass_jpsi_mu.eps",outputdir.c_str());
    c2->Print(name);
    ofparmu<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
    ofparmu<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError()<<std::endl;
    delete data;
    delete xframe;
    delete sum;
    //delete c2;

   // ~~~~~~~~~muon part end~~~~~~~~~~

   // ~~~~~~~~~pion part~~~~~~~~~~
   m0 = 3.686109;
   sigma =0.00258;
   width = 10*sigma;
   mparticle =0.13957018;//pi
   mparticle2=0.000511;//electron
   mparticle3=0.105658;//muon
  
   // pre fit
   ofparpi<<"pre fit"<<"\n";
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     if(cos(angle4)>0.90) continue; // cut bhabha 
     double mass,massjpsi;
     double totpx,totpy,totpz,tote;
     double jpsipx,jpsipy,jpsipz,jpsie;
     double lee[2],pie[2];
     // total invariant mass
     if(cos(angle4)>0.90) continue; // cut bhabha 
     totpx=(lepx4[0]+lepx4[1])+(pipx4[0]+pipx4[1]);
     totpy=(lepy4[0]+lepy4[1])+(pipy4[0]+pipy4[1]);
     totpz=(lepz4[0]+lepz4[1])+(pipz4[0]+pipz4[1]);
     jpsipx=(lepx4[0]+lepx4[1]);
     jpsipy=(lepy4[0]+lepy4[1]);
     jpsipz=(lepz4[0]+lepz4[1]);
     pie[0]=TMath::Sqrt(mparticle*mparticle + 
            (pipx4[0]*pipx4[0]+pipy4[0]*pipy4[0]+pipz4[0]*pipz4[0]) );
     pie[1]=TMath::Sqrt(mparticle*mparticle +
            (pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]+pipz4[1]*pipz4[1]) );
     tote=lee4[0]+lee4[1]+pie[0]+pie[1];
     jpsie=lee4[0]+lee4[1];
     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
     massjpsi=TMath::Sqrt(jpsie*jpsie-jpsipx*jpsipx-jpsipy*jpsipy-jpsipz*jpsipz);
     mass = mass-massjpsi+3.096916;
     //if(mass>m0-width/2. && mass<m0+width/2.)
     {
       hmass3->Fill(mass);
     }
   }
   // total invariant mass
   // if (Cut(ientry) < 0) continue;

   sprintf(name,"mass_psi_pi");
   x.setRange(psilow,psiup);
   x.setVal(m0);
   data = new RooDataHist(name,"data_e",x,hmass3);
   mean.setRange(psilow,psiup);
   mean.setVal(m0);
   sigma1.setVal(sigma);
   sigma1.setRange(0.9*sigma,1.1*sigma);
   signal.setVal(120);
   background.setVal(0);
   co1.setVal(0);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
   sum->fitTo(*data,Range(psilow,psiup));
   xframe = x.frame(Title("fit pi"));
   data->plotOn(xframe);
   sum->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
   xframe->Draw();
   sprintf(name,"%s/mass_psi_pi_pre.eps",outputdir.c_str());
   c2->Print(name);
   ofparpi<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofparpi<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError()<<std::endl;
   delete data;
   delete xframe;
   delete sum;

   px1.clear();
   px2.clear();
   py1.clear();
   py2.clear();
   pz1.clear();
   pz2.clear();
   px3.clear();
   px4.clear();
   py3.clear();
   py4.clear();
   pz3.clear();
   pz4.clear();
   le1.clear();
   le2.clear();
   hp->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(cos(angle4)>0.90) continue; // cut bhabha 
      double mass,massjpsi;
      double totpx,totpy,totpz,tote;
      double p1,p2;
      double jpsipx,jpsipy,jpsipz,jpsie;
      double lee[2],pie[2];
      // total invariant mass
      // if(cos(angle4)>0.90) continue; // cut bhabha 
      totpx=(lepx4[0]+lepx4[1])+(pipx4[0]+pipx4[1]);
      totpy=(lepy4[0]+lepy4[1])+(pipy4[0]+pipy4[1]);
      totpz=(lepz4[0]+lepz4[1])+(pipz4[0]+pipz4[1]);
      p1=TMath::Sqrt(pipx4[0]*pipx4[0]+pipy4[0]*pipy4[0]+pipz4[0]*pipz4[0]);
      p2=TMath::Sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]+pipz4[1]*pipz4[1]);
      tote=lee4[0]+lee4[1]+pie4[0]+pie4[1];
      jpsipx=(lepx4[0]+lepx4[1]);
      jpsipy=(lepy4[0]+lepy4[1]);
      jpsipz=(lepz4[0]+lepz4[1]);
      //tote=lee4[0]+lee4[1]+pie4[0]+pie4[1];
      //tote=lee[0]+lee[1]+pie[0]+pie[1];
      jpsie=lee4[0]+lee4[1];
      hp->Fill(p1);
      hp->Fill(p2);
      mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
      massjpsi=TMath::Sqrt(jpsie*jpsie-jpsipx*jpsipx-jpsipy*jpsipy-jpsipz*jpsipz);
      mass = mass-massjpsi+3.096916;
      if (mass>m0-width/2. && mass<m0+width/2.){
        px1.push_back(pipx4[0]);
        px2.push_back(pipx4[1]);
        py1.push_back(pipy4[0]);
        py2.push_back(pipy4[1]);
        pz1.push_back(pipz4[0]);
        pz2.push_back(pipz4[1]);
        px3.push_back(lepx4[0]);
        px4.push_back(lepx4[1]);
        py3.push_back(lepy4[0]);
        py4.push_back(lepy4[1]);
        pz3.push_back(lepz4[0]);
        pz4.push_back(lepz4[1]);
        le1.push_back(lee4[0]);
        le2.push_back(lee4[1]);
      }
      // if (Cut(ientry) < 0) continue;
   }
   hp->Draw();
   sprintf(name,"%s/momentumpi_pipill.eps",outputdir.c_str());
   c2->Print(name);
   
   detail<<"m0 is "<<m0<<", data size is "<<px1.size()<<std::endl;
   //TF2 *likeli3=new TF2("likeli3",maxlikelihood4_0,0.95,1.05,0.1,0.99);
   //TCanvas *c2=new TCanvas("c2","likelihood",800,600);
   //likeli3->Draw("surf1");
   //tmpstr = outputdir +"/likelipi_2D.eps";
   //c2->Print(tmpstr.c_str());
   //minimum = likeli3->GetMinimum(0.98,1.02);
   //likeli3->GetMinimumXY(factor,miny);
   //weight = miny;
   sigNo=signal.getVal() ;
   bckNo=width/(psiup-psilow)*background.getVal();
   weight = sigNo/(sigNo+bckNo);
   TF1 *likeli3_1=new TF1("likeli3_1",maxlikelihood4_1,0.95,1.05);
   likeli3_1->Draw();
   tmpstr = outputdir +"/likelipi_1D.eps";
   c2->Print(tmpstr.c_str());
   factor = likeli3_1->GetMinimumX(0.98,1.02);
   minimum = likeli3_1->GetMinimum(0.98,1.02);
   factorlow=likeli3_1->GetX(minimum+1,0.98,factor);
   factorup =likeli3_1->GetX(minimum+1,factor,1.02);
   ofparpi<<run<<"\t"<<factor<<"\t"<<factorlow<<"\t"<<factorup<<std::endl;
   ofparpi<<"\t"<<factor<<"\t"<<weight<<"\t"<<likeli3_1->GetMinimumX(0.98,1.02)<<std::endl;
   detail<<"weight is "<<miny<<", factor is "<<likeli3_1->GetMinimumX(0.98,1.02)<<std::endl;
   detail<<"best factor  "<<factor<<std::endl;
   purepar<<factor<<"\t"<<(factorup-factorlow)/2;
   
   // using the factor to fit
   hmass3->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     if(cos(angle4)>0.90) continue; // cut bhabha 
     double mass,massjpsi;
     double totpx,totpy,totpz,tote;
     double jpsipx,jpsipy,jpsipz,jpsie;
     double lee[2],pie[2];
     // total invariant mass
     if(cos(angle4)>0.90) continue; // cut bhabha 
     totpx=(lepx4[0]+lepx4[1])+factor*(pipx4[0]+pipx4[1]);
     totpy=(lepy4[0]+lepy4[1])+factor*(pipy4[0]+pipy4[1]);
     totpz=(lepz4[0]+lepz4[1])+factor*(pipz4[0]+pipz4[1]);
     jpsipx=(lepx4[0]+lepx4[1]);
     jpsipy=(lepy4[0]+lepy4[1]);
     jpsipz=(lepz4[0]+lepz4[1]);
     pie[0]=TMath::Sqrt(mparticle*mparticle + 
            factor*factor*(pipx4[0]*pipx4[0]+pipy4[0]*pipy4[0]+pipz4[0]*pipz4[0]) );
     pie[1]=TMath::Sqrt(mparticle*mparticle +
            factor*factor*(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]+pipz4[1]*pipz4[1]) );
     tote=lee4[0]+lee4[1]+pie[0]+pie[1];
     jpsie=lee4[0]+lee4[1];
     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
     massjpsi=TMath::Sqrt(jpsie*jpsie-jpsipx*jpsipx-jpsipy*jpsipy-jpsipz*jpsipz);
     mass = mass-massjpsi+3.096916;
     if(mass>m0-width/2. && mass<m0+width/2.){
       hmass3->Fill(mass);
     }
   }
   // total invariant mass
   // if (Cut(ientry) < 0) continue;

   sprintf(name,"mass_psi_pi");
   x.setRange(psilow,psiup);
   x.setVal(m0);
   data = new RooDataHist(name,"data_e",x,hmass3);
   mean.setRange(psilow,psiup);
   mean.setVal(m0);
   sigma1.setVal(sigma);
   sigma1.setRange(0.9*sigma,1.1*sigma);
   signal.setVal(120);
   background.setVal(0);
   co1.setVal(0);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
   sum->fitTo(*data,Range(psilow,psiup));
   xframe = x.frame(Title("fit pi"));
   data->plotOn(xframe);
   sum->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
   xframe->Draw();
   sprintf(name,"%s/mass_psi_pi.eps",outputdir.c_str());
   c2->Print(name);
   ofparpi<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofparpi<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError()<<std::endl;
   delete data;
   delete xframe;
   delete sum;
   delete c2;


     // ~~~~~~~~~pion part end~~~~~~~~~~

   ofpare.close();
   ofparmu.close();
   ofparpi.close();
   detail.close();
   return true;
}


gepep_fastpipill::gepep_fastpipill(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_Rvalue_pipill_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_Rvalue_pipill_1.root");
      }
      f->GetObject("gepep_fastpipill",tree);

   }
   Init(tree);
}

gepep_fastpipill::~gepep_fastpipill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fastpipill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fastpipill::LoadTree(Long64_t entry)
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

void gepep_fastpipill::Init(TTree *tree)
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
   fChain->SetBranchAddress("mcgpx", mcgpx, &b_mcgpx);
   fChain->SetBranchAddress("mcgpy", mcgpy, &b_mcgpy);
   fChain->SetBranchAddress("mcgpz", mcgpz, &b_mcgpz);
   fChain->SetBranchAddress("mcge", mcge, &b_mcge);
   fChain->SetBranchAddress("mcgid", mcgid, &b_mcgid);
   fChain->SetBranchAddress("nGammatch", &nGammatch, &b_nGammatch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("ntot", &ntot, &b_ntot);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ngam", &ngam, &b_ngam);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("netap2", &netap2, &b_netap2);
   fChain->SetBranchAddress("delang", delang, &b_delang);
   fChain->SetBranchAddress("delphi", delphi, &b_delphi);
   fChain->SetBranchAddress("delthe", delthe, &b_delthe);
   fChain->SetBranchAddress("npart", npart, &b_npart);
   fChain->SetBranchAddress("nemchits", nemchits, &b_nemchits);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("dx", dx, &b_dx);
   fChain->SetBranchAddress("dy", dy, &b_dy);
   fChain->SetBranchAddress("dz", dz, &b_dz);
   fChain->SetBranchAddress("dtheta", dtheta, &b_dtheta);
   fChain->SetBranchAddress("dphi", dphi, &b_dphi);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("dE", dE, &b_dE);
   fChain->SetBranchAddress("eSeed", eSeed, &b_eSeed);
   fChain->SetBranchAddress("nSeed", nSeed, &b_nSeed);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("e5x5", e5x5, &b_e5x5);
   fChain->SetBranchAddress("secondMoment", secondMoment, &b_secondMoment);
   fChain->SetBranchAddress("latMoment", latMoment, &b_latMoment);
   fChain->SetBranchAddress("a20Moment", a20Moment, &b_a20Moment);
   fChain->SetBranchAddress("a42Moment", a42Moment, &b_a42Moment);
   fChain->SetBranchAddress("getTime", getTime, &b_getTime);
   fChain->SetBranchAddress("getEAll", getEAll, &b_getEAll);
   fChain->SetBranchAddress("mpi0", mpi0, &b_mpi0);
   fChain->SetBranchAddress("chisq1cpi0", chisq1cpi0, &b_chisq1cpi0);
   fChain->SetBranchAddress("ig1pi0", ig1pi0, &b_ig1pi0);
   fChain->SetBranchAddress("ig2pi0", ig2pi0, &b_ig2pi0);
   fChain->SetBranchAddress("chisq6c", &chisq6c, &b_chisq6c);
   fChain->SetBranchAddress("chisq4c", &chisq4c, &b_chisq4c);
   fChain->SetBranchAddress("chisq3c", &chisq3c, &b_chisq3c);
   fChain->SetBranchAddress("gpx4", gpx4, &b_gpx4);
   fChain->SetBranchAddress("gpy4", gpy4, &b_gpy4);
   fChain->SetBranchAddress("gpz4", gpz4, &b_gpz4);
   fChain->SetBranchAddress("ge4", ge4, &b_ge4);
   fChain->SetBranchAddress("lepx4", lepx4, &b_lepx4);
   fChain->SetBranchAddress("lepy4", lepy4, &b_lepy4);
   fChain->SetBranchAddress("lepz4", lepz4, &b_lepz4);
   fChain->SetBranchAddress("lee4", lee4, &b_lee4);
   fChain->SetBranchAddress("pipx4", pipx4, &b_pipx4);
   fChain->SetBranchAddress("pipy4", pipy4, &b_pipy4);
   fChain->SetBranchAddress("pipz4", pipz4, &b_pipz4);
   fChain->SetBranchAddress("pie4", pie4, &b_pie4);
   fChain->SetBranchAddress("ggm4", &ggm4, &b_ggm4);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   fChain->SetBranchAddress("kpi01m4", &kpi01m4, &b_kpi01m4);
   fChain->SetBranchAddress("kpi02m4", &kpi02m4, &b_kpi02m4);
   fChain->SetBranchAddress("kkpi0m4", &kkpi0m4, &b_kkpi0m4);
   fChain->SetBranchAddress("kkpipim4", &kkpipim4, &b_kkpipim4);
   fChain->SetBranchAddress("kkpipipi0m4", &kkpipipi0m4, &b_kkpipipi0m4);
   fChain->SetBranchAddress("Recoilmass", &Recoilmass, &b_Recoilmass);
   fChain->SetBranchAddress("llm4", &llm4, &b_llm4);
   fChain->SetBranchAddress("pipim4", &pipim4, &b_pipim4);
   fChain->SetBranchAddress("llpipi1m4", &llpipi1m4, &b_llpipi1m4);
   fChain->SetBranchAddress("llpipi2m4", &llpipi2m4, &b_llpipi2m4);
   fChain->SetBranchAddress("llpipi3m4", &llpipi3m4, &b_llpipi3m4);
   fChain->SetBranchAddress("llpipi4m4", &llpipi4m4, &b_llpipi4m4);
   fChain->SetBranchAddress("gpipim4", &gpipim4, &b_gpipim4);
   fChain->SetBranchAddress("decay_ee", &decay_ee, &b_decay_ee);
   fChain->SetBranchAddress("hepp", &hepp, &b_hepp);
   fChain->SetBranchAddress("hepm", &hepm, &b_hepm);
   fChain->SetBranchAddress("emcp", &emcp, &b_emcp);
   fChain->SetBranchAddress("emcm", &emcm, &b_emcm);
   fChain->SetBranchAddress("mucp", &mucp, &b_mucp);
   fChain->SetBranchAddress("mucm", &mucm, &b_mucm);
   fChain->SetBranchAddress("zcpm4", &zcpm4, &b_zcpm4);
   fChain->SetBranchAddress("zcmm4", &zcmm4, &b_zcmm4);
   fChain->SetBranchAddress("ecms", &ecms, &b_ecms);
   fChain->SetBranchAddress("angle1", &angle1, &b_angle1);
   fChain->SetBranchAddress("angle2", &angle2, &b_angle2);
   fChain->SetBranchAddress("angle3", &angle3, &b_angle3);
   fChain->SetBranchAddress("angle4", &angle4, &b_angle4);
   fChain->SetBranchAddress("psipm4", &psipm4, &b_psipm4);
   Notify();
}

Bool_t gepep_fastpipill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fastpipill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fastpipill::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
