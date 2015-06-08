#define gepep_kk_cxx
#include "gepep_kk.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
//#include <iostream>
extern std::string outputdir;
using namespace RooFit;
namespace D0COR{

  void FitSpectrum(TTree *&dataraw, char* namesfx, bool out=false);
}

void gepep_kk::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_kk.C
//      Root > gepep_kk t
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
   int nBins=200;
   double factorstart=0.99;
   double mk=0.493677;
   // D0 -> K K
   double philow=1.82;
   double phiup=1.90;
   double peakvalue=1.86484;// mphi
   // phi -> K K
   //double philow=0.995;
   //double phiup=1.045;
   //double peakvalue=1.019455;// mphi
   
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,philow,phiup,"GeV");
   //RooRealVar x("x","energy",peakvalue,1.015,1.025,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,philow,phiup);
   RooRealVar mean2("mean2","mean of gaussian2",peakvalue,philow,phiup);
   //RooRealVar sigma("sigma","width of gaussian",0.0023,0.0015,0.0040);//phi version
   //RooRealVar sigma2("sigma2","width of gaussian",0.005,0.003,0.008);
   RooRealVar sigma("sigma","width of gaussian",0.007,0.003,0.0075);//D0 version
   RooRealVar sigma2("sigma2","width of gaussian",0.02,0.008,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   
   RooRealVar co1("co1","coefficient #1",0,-0.5,0.5);
   RooRealVar co2("co2","coefficient #2",0,-0.01,0.5);
   RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
   
   RooRealVar alpha("alpha","#alpha",1.16,-5,5);
   RooRealVar nnn("nnn","n",50,1,200);
   RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);

   //RooBreitWigner brewig("brewig","brewig",x,mean,sigma);
   //RooLandau landau("landau","landau",x,mean,sigma);

   RooRealVar a0("a0","coefficient #0",100,100,100000);
   RooRealVar a1("a1","coefficient #1",-50,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooRealVar signal("signal"," ",1200,10,1000000);//event number
   RooRealVar signal2("signal2"," ",1200,0,1000000);//event number
   RooRealVar background("background"," ",200,0,1000000);
   
   RooPlot *xframe;
   //RooDataHist *data_k;
   RooDataSet *dataset;
   RooAddPdf *sum;
   
   //RooDataSet *dataset = new RooDataSet("dataset","data",dataraw,x);
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   int Npar;
   char fname[1000];
   ofstream ofpar;
   sprintf(fname,"%s/pars.txt",outputdir.c_str());
   ofpar.open(fname,std::ios::app);
   ofstream purepar;
   sprintf(fname,"%s/parspure.txt",outputdir.c_str());
   purepar.open(fname,std::ios::app);
   sprintf(fname,"%s/plot_kk.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");
   TTree *dataraw=new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");
   TTree *vars = new TTree("vars","vars");
   double phi1,phi2;
   double costheta,costheta1,costheta2;
   double p1,p2;
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta" ,&costheta ,"costheta/D" );
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"p1/D");
   vars->Branch("p2",&p2,"p2/D");
   vars->Branch("mass",&mass,"mass/D");
   double chi2;
   vars->Branch("chi2",&chi2,"chi2/D");

   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1    = new TH1D("h1","2 kaon invariant mass",nBins,philow,phiup);
   TCanvas *c1 = new TCanvas("","",800,600);

   const int Npart=1;
   double m0=peakvalue;
   double sigma_m=0.0024;//0.0024 for phi,
   double width = 10.*sigma_m;
   double mparticle=0.493677;
   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;

   int realsize=0;
   double partid[Npart];
   double parter[Npart];
   double corfac[Npart];
   double corerr[Npart];
   
   double pcut[Npart+1];
   double facv[Npart];
   double facev[Npart];


     pcut[0]=0; facv[0] = 1.000610; facev[0]=8.6e-5;
     pcut[1] = 2.0;
// pcut[0] =0.0 ;    facv[0] =1.0;       facev[0] =1.0;  
// pcut[1] =0.10;    facv[1] =1.0     ;  facev[1] =1.0        ;  
// pcut[2] =0.20;    facv[2] =1.0     ;  facev[2] =1.0        ;
// pcut[3] =0.30;    facv[3] =1.0     ;  facev[3] =1.0        ;
// pcut[4] =0.40;    facv[4] =1.0     ;  facev[4] =1.0        ;
// pcut[5] =0.50;    facv[5] =1.00338 ;  facev[5] =1.0        ;
// pcut[6] =0.60;    facv[6] =1.00484 ;  facev[6] =1.0        ;
// pcut[7] =0.70;    facv[7] =1.0024  ;  facev[7] =1.0        ;
// pcut[8] =0.80;    facv[8] =0.999958;  facev[8] =1.0        ;
// pcut[9] =0.90;    facv[9] =0.997556;  facev[9] =1.0        ;
// pcut[10]=1.00;    facv[10]=0.997925;  facev[10]=1.0        ;
// pcut[11]=1.10;    facv[11]=0.999218;  facev[11]=1.0        ;
// pcut[12]=1.20;    facv[12]=1.00027 ;  facev[12]=1.0        ;
// pcut[13]=1.30;    facv[13]=1.0     ;  facev[13]=1.0        ;
// pcut[14]=1.40;    facv[14]=1.0     ;  facev[14]=1.0        ;
// pcut[15]=1.50;    facv[15]=1.0     ;  facev[15]=1.0        ;
// pcut[16]=1.60;    facv[16]=1.0     ;  facev[16]=1.0        ;
// pcut[17]=1.70;    facv[17]=1.0     ;  facev[17]=1.0        ;
// pcut[18]=1.80;    facv[18]=1.0;       facev[18]=1.0;  
// pcut[19]=1.90;    facv[19]=1.0;       facev[19]=1.0;  
// pcut[20]=2.00;           

   char name[100];
   TH1D *hmD0[Npart];
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%0d",partj);
     hmD0[partj] = new TH1D(name,name,100,philow,phiup);
   }
 
   // loop data
   TH1D *hp = new TH1D("hp","hp",200,0,2);
   TH2D *hp2= new TH2D("hp2","hp2",100,0,2, 100,0,2);
   Event evt(mk);
   std::vector<Event> evts;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
	    evt.SetVal(kappx,kappy,kappz,kampx,kampy,kampz);
        mass = evt.InvMass();
        p1 = evt.GetP1();
        p2 = evt.GetP2();
        costheta  = evt.GetCostheta();
        costheta1 = evt.GetCostheta1();
		costheta2 = evt.GetCostheta2();
        phi1 = evt.GetPhi1();
        phi2 = evt.GetPhi2();
		chi2 = kkm4;
        vars->Fill();
        //if (costheta1 > 0.8 && costheta2 <-0.8) continue;
		//if (p1<0.5 || p1>1.3) continue;

        //if ( partj>=Npart || partj<0 ) continue;
        for (int partj=0;partj<Npart;partj++){
          if (p1<pcut[partj] || p1>pcut[partj+1]) continue;
          if (p2<pcut[partj] || p2>pcut[partj+1]) continue;
          if (mass>philow-0.02 && mass<phiup+0.02){
            hmD0[partj]->Fill(mass);
			evts.push_back(evt);
	      }
          break;
		}
   
   }
   std::cout<<"fffffffffffffffff"<<std::endl;

  
     xframe = x.frame(Title("fit kaon"));

      dataraw->Reset();
      std::cout<<"factor is "<<1<<std::endl;
      for (Long64_t jentry=0; jentry<evts.size();jentry++) {
          mass = evts.at(jentry).InvMass();
          if (mass>philow && mass<phiup){
            dataraw->Fill();
          }
         // if (Cut(ientry) < 0) continue;
      }
      //dataraw->Write();

      dataset = new RooDataSet("dataset","data",dataraw,x);
      char tmpchr[100];
      sprintf(tmpchr,"data_k");
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      mean.setVal(peakvalue);
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
      sum->plotOn(xframe);
      xframe->Draw();
      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(4000);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.035);
      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      pt->Draw();
      c1->Write("raw_spec");
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;


 
     xframe = x.frame(Title("fit kaon"));

     double factor = 1.000610;
      //h1->Reset();
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evts.size();jentry++) {
          mass = evts.at(jentry).InvMass(factor,factor);
          if (mass>philow && mass<phiup){
            dataraw->Fill();
          }
         // if (Cut(ientry) < 0) continue;
      }
      //dataraw->Write();

      dataset = new RooDataSet("dataset","data",dataraw,x);
      //char tmpchr[100];
      sprintf(tmpchr,"data_k");
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
      sum->plotOn(xframe);
      xframe->Draw();
      //TPaveText *
	  pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(4000);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.035);
      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      pt->Draw();
      c1->Write("cor_spec");
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;


   f->Close();
   return;
}

#ifdef gepep_kk_cxx
gepep_kk::gepep_kk(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/RValue_kk_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/RValue_kk_3850.root");
      }
      f->GetObject("gepep_kk",tree);

   }
   Init(tree);
}

gepep_kk::~gepep_kk()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_kk::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_kk::LoadTree(Long64_t entry)
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

void gepep_kk::Init(TTree *tree)
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
   fChain->SetBranchAddress("kappx", &kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", &kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", &kappz, &b_kappz);
   fChain->SetBranchAddress("kape", &kape, &b_kape);
   fChain->SetBranchAddress("kampx", &kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", &kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", &kampz, &b_kampz);
   fChain->SetBranchAddress("kame", &kame, &b_kame);
   //fChain->SetBranchAddress("mphi", &mphi, &b_mphi);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_kk::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_kk::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_kk::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_kk_cxx
