#define gepep_kpi_cxx
#include "gepep_kpi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
//#include <iostream>
extern std::string outputdir;
using RooFit::Title;
using RooFit::Components;
using RooFit::LineStyle;
using RooFit::LineColor;
using RooFit::Range;


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
   int nBins=40;
   double factorstart=0.99;
   double D0low=1.82;
   double D0up=1.90;
   double mk=0.493677;
   double mpi=0.13957018;
   double peakvalue=1.86486;// mD0
   int pointNo=20;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   double factor=factorstart;
   double factorstep=(1.-factor)*2/pointNo;
   double factork=1;//0.998881;
   
   // try to use roofit
   RooRealVar x("x","energy",1.865,D0low,D0up,"GeV");
   RooRealVar mean("mean","mean of gaussian",1.865,D0low,D0up);
   RooRealVar sigma("sigma","width of gaussian",0.0068,0.005,0.007);
   RooRealVar sigma2("sigma2","width of gaussian",0.0075,0.007,0.010);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   
   RooRealVar alpha("alpha","#alpha",1.32,-5,5);
   RooRealVar nnn("nnn","n",5,1,200);
   RooCBShape cbshape("cbshape","cbshape",x,mean,sigma,alpha,nnn);

   RooRealVar co1("co1","coefficient #1",0,-1000.,1000.);
   RooRealVar co2("co2","coefficient #2",0,-1000.,0.1);
   RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
   
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooRealVar signal("signal"," ",12000,10,1000000);//event number
   RooRealVar signal2("signal2"," ",12000,10,1000000);//event number
   RooRealVar background("background"," ",2000,0,100000);
   RooPlot *xframe;
   RooDataHist *data_kpi;
   RooDataSet *dataset;
   RooAddPdf *sum;
   int Npar;
 
   ofstream ofpar;
   ofpar.open("parkpi.txt",std::ios::app);
   ofpar<<"k- pi+ algrithm: will give factors for pion"<<std::endl;
   ofstream ofpardetail;
   ofpardetail.open("detail.txt",std::ios::app);
   ofstream purepar;
   purepar.open("par");
   char fname[100];
   sprintf(fname,"%s/plot_kpi.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");

   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1    = new TH1D("h1","k- pi+ invariant mass",nBins,D0low,D0up);
   TCanvas *c1 = new TCanvas("","",800,600);
   TTree *dataraw = new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");

 
   Long64_t nbytes = 0, nb = 0;

   TH2D *h2p = new TH2D("h2p","#pi K momentum",200,0,2,200,0,2);
   TH2D *thedis = new TH2D("thedis","#theta",200,0,2,100,0,TMath::Pi());
   TH2D *thedis2= new TH2D("thedis2","#theta",100,0,TMath::Pi(),100,0,TMath::Pi());
   int Npart=1;
   double factorpi=1.00132;
   double thetacut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
   double start=0.0;
   double stop =TMath::Pi();
   for(int i=0;i<Npart+1;i++){
     thetacut[i] = (stop-start)/Npart*i+start;
   }
   double m0=peakvalue;
   double sigma_m=0.0068;//0.0024 for phi,
   double width = 10.*sigma_m;
   //double mparticle=mk;
   std::vector<std::pair<int,int> > partmap;
   std::vector<std::pair<int,double> > facmap;

   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   TH1D *hmtheta[Npart][Npart];
   for (int parti=0;parti<Npart;parti++){
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%d_part%d",parti,partj);
     hmtheta[parti][partj] = new TH1D(name,name,100,D0low,D0up);
   }
   }

   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
  
     int parti,partj;
     //double mass;
     double p1,p2;
	 double theta1,theta2;
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
     mass = CalInvMass(mpi,pippx[besidx],pippy[besidx],pippz[besidx],mk,kampx[0],kampy[0],kampz[0]);
     
     // total invariant mass, D0 -> k- pi+
     p1 = CalMom(pippx[besidx],pippy[besidx],pippz[besidx]);
	 p2 = CalMom(kampx[0],kampy[0],kampz[0]);
     theta1= acos(pippz[besidx]/p1);
     theta2= acos(kampz[0]/p2);
	 parti = (int)(theta1/stop*Npart);
     partj = (int)(theta2/stop*Npart);
     if (parti>=Npart || partj>=Npart || parti<0 || partj<0) continue;
     if (theta1>thetacut[parti]&&theta1<thetacut[parti+1]
	   &&theta2>thetacut[partj]&&theta2<thetacut[partj+1])
       if (mass>m0-width/2. && mass<m0+width/2.)
	   {
         h2p->Fill(p1,p2);
		 thedis->Fill(p1,theta1);
		 thedis->Fill(p2,theta2);
		 thedis2->Fill(theta1,theta2);
         hmtheta[parti][partj]->Fill(mass);
       }
     // if (Cut(ientry) < 0) continue;
   }
   h2p->Write();
   thedis->Write();
   thedis2->Write();
   for (int parti=0;parti<Npart;parti++)
   for (int partj=0;partj<Npart;partj++){
     hmtheta[parti][partj]->Write();
     std::cout<<"processed part "<<parti<<", part "<<partj<<std::endl;
     if (hmtheta[parti][partj]->GetEntries() > 100
	   &&hmtheta[parti][partj]->GetMaximumBin()>30
	   &&hmtheta[parti][partj]->GetMaximumBin()<70){
       partmap.push_back(std::make_pair(parti,partj));
     }
   }
   // ~~~~~~~~ draw end

   for (int loopi=0;loopi<2;loopi++){
   for (int loopj=0;loopj<partmap.size();loopj++){
     int parti=partmap.at(loopj).first;
     int partj=partmap.at(loopj).second;
	 if (loopi==0 && parti!=partj) continue;
	 else if(loopi!=0 && parti==partj) continue;
     double factori=0,factorj=0;
     for (int loopk=0;loopk<facmap.size();loopk++){
       if (facmap.at(loopk).first == parti) factori=facmap.at(loopk).second;
       if (facmap.at(loopk).first == partj) factorj=facmap.at(loopk).second;
     }
     if(factori==0 && factorj==0 && parti!=partj) continue;
     if(parti==partj && factori!=0) continue;
     if(factori!=0 && factorj!=0) continue;
   ofpar<<"parti: "<<parti<<"--"<<factori<<", partj: "<<partj<<"--"<<factorj<<std::endl;
   std::cout<<"parti: "<<parti<<"--"<<factori<<", partj: "<<partj<<"--"<<factorj<<std::endl;

   // for saving the fit result
   char fitepsname[100];
   char fiteps_start[100];
   char fiteps_stop[100];
   sprintf(fitepsname,"%s/fitkpi_part%d_part%d.eps",outputdir.c_str(),parti,partj);
   sprintf(fiteps_start,"%s[",fitepsname);
   sprintf(fiteps_stop,"%s]",fitepsname);
   c1->Print(fiteps_start);
 
   factor = factorstart;
   int fittimes = 0;
   for (int i=0;i<pointNo;i++){
	  xframe = x.frame(Title("fit k pi"));

      h1->Reset();
	  dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
	     //double mass;
		 double p1,p2;
		 double theta1,theta2;
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
 		 double factors[2]={1,1};
		 factori==0? factors[0]=factor : factors[0]=factori;
		 //factorj==0? factors[1]=factor : factors[1]=factorj;
 	     factors[1]=factorpi;
		 p1=CalMom(pippx[besidx],pippy[besidx],pippz[besidx]);
         p2=CalMom(kampx[0],kampy[0],kampz[0]);
         theta1=acos(pippz[besidx]/p1);
		 theta2=acos(kampz[0]/p2);
         if (theta1>thetacut[parti]&&theta1<thetacut[parti+1]
		   &&theta2>thetacut[partj]&&theta2<thetacut[partj+1] ){
           mass = CalInvMass(mpi,pippx[besidx],pippy[besidx],pippz[besidx],mk,kampx[0],kampy[0],kampz[0],-2,factors);
	       h1->Fill(mass);
		 }
		 if (mass>D0low && mass<D0up){
		   dataraw->Fill();
		 }
         // if (Cut(ientry) < 0) continue;
      }

	  char tmpchr[100];
	  sprintf(tmpchr,"data_kpi_%02d",fittimes);
      //data_kpi = new RooDataHist(tmpchr,"data_kpi",x,h1);
      dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
	  sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
      mean.setVal(peakvalue+1.4*(factor-1.0));
	  Npar=8;
	  //sigma.setVal(0.035);
	  //signal.setVal(1200);
	  //background.setVal(500);
	  //co1.setVal(0);
	  //co2.setVal(-0.2);
      sum->fitTo(*dataset,Range(D0low,D0up));
	  dataset->plotOn(xframe);
	  sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	  sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
	  sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
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
      c1->Update();
      c1->Print(fitepsname);
	  delete dataset;
	  delete xframe;
	  delete sum;

	  // save pars
	  factors[i]=factor;
	  factorserr[i]=0;
	  deltapeaks[i] = mean.getVal() - peakvalue;
	  deltapeakserr[i] = mean.getError();

	  fittimes++;
	  factor += factorstep;
   }
   std::cout<<"entry is "<<nentries<<std::endl;
   c1->Print(fiteps_stop);
   c1->Clear();
   
   TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   factor =facfit->GetParameter(0);
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   double factorerr=deltapeakserr[9]/facfit->GetParameter(1);
   ofpar<<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<std::endl;
   ofpar<<signal.getVal()<<"\t"<<signal.getError()<<std::endl;
   purepar<<"\t"<<factor<<"\t"<<factorerr<<std::endl;
   //std::cout<<"fit factor: "<<factor1<<", error is "<<factor1err<<std::endl;
   sprintf(name,"%s/factorkpi_part%d_part%d.eps",outputdir.c_str(),parti,partj);
   c1->Print(name);
   sprintf(name,"factors_kpi_part%d_part%d",parti,partj);
   graph1->SetName(name);
   graph1->Write();
   int index;
   factori==0? index=parti : index=partj;
   facmap.push_back(std::make_pair(index,factor));

   // draw the best fit
   xframe = x.frame(Title("fit k pi"));
   h1->Reset();
   dataraw->Reset();
   std::cout<<"factor is "<<factor<<std::endl;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
	 //double mass;
	 double p1,p2;
	 double theta1,theta2;
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
 	 double factors[2]={1,1};
	 factori==0? factors[0]=factor : factors[0]=factori;
	 factorj==0? factors[1]=factor : factors[1]=factorj;
 	 p1=CalMom(pippx[besidx],pippy[besidx],pippz[besidx]);
     p2=CalMom(kampx[0],kampy[0],kampz[0]);
	 theta1=acos(pippz[besidx]/p1);
	 theta2=acos(kampz[0]/p2);
     if (theta1>thetacut[parti]&&theta1<thetacut[parti+1]
	   &&theta2>thetacut[partj]&&theta2<thetacut[partj+1] ){
       mass = CalInvMass(mpi,pippx[besidx],pippy[besidx],pippz[besidx],mk,kampx[0],kampy[0],kampz[0],-2,factors);
	   h1->Fill(mass);
	 }
	 if (mass>D0low && mass<D0up){
	   dataraw->Fill();
	 }
     // if (Cut(ientry) < 0) continue;
   }
	char tmpchr[100];
	sprintf(tmpchr,"data_kpi");
    //data_kpi = new RooDataHist(tmpchr,"data_kpi",x,h1);
    dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
    mean.setVal(peakvalue+0.05*(factor-1.0));
	Npar=8;
	//sigma.setVal(0.035);
	//signal.setVal(1200);
	//background.setVal(200);
	//co1.setVal(0);
	//co2.setVal(-0.2);
    sum->fitTo(*dataset,Range(D0low,D0up));
	dataset->plotOn(xframe);
	sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
	sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
	sum->plotOn(xframe);
    xframe->Draw();
    TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
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
    c1->Update();
    sprintf(name,"fitkpi_best_part%d_part%d.eps",parti,partj);
	c1->Print(name);
	xframe->SetName(name);
	xframe->Write();
	delete dataset;
	delete xframe;
	delete sum;
   
   }
   }

   f->Close();
   ofpar.close();
   ofpardetail.close();

}

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
