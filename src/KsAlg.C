#define KsAlg_cxx
#include "KsAlg.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "function.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include <vector>
extern std::string outputdir;
using namespace RooFit;

void KsAlg::Loop(double cutd)
{
//   In a ROOT session, you can do:
//      Root > .L KsAlg.C
//      Root > KsAlg t
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

   std::vector<double> ppx,ppy,ppz;
   std::vector<double> mpx,mpy,mpz;

   TH1D *hmass = new TH1D("hmass","Ks mass",500,0.4,0.6);
   //TH2D *h2pos = new TH2D("h2pos","pos",10, 0,5, 10,0,5);
   TH2D *h2p = new TH2D("h2p","P",200, 0,2, 200,0,2);
   int nBins=100;
   double mpi=0.13957;
   double mparticle=0.497614;
   double kslow=0.43;
   double ksup =0.57;
   TF1 *facfit = new TF1("facfit",line2,kslow,ksup,2);
   char fname[100];
   sprintf(fname,"%s/plot_ks2pi.root",outputdir.c_str());
   TFile *f=new TFile(fname,"RECREATE");
   TCanvas *c1 = new TCanvas();
   int Npar;

   // roofit variables and functions
   RooRealVar x("x","energy",mparticle,kslow,ksup,"GeV");
   RooRealVar mean("mean","mean of gaussian",  kslow,ksup);
   RooRealVar mean2("mean2","mean of gaussian",kslow,ksup);
   RooRealVar sigma("sigma","width of gaussian",0.0017,0.0010,0.0050);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.01);
   RooRealVar brewid("brewid","width of breit wigner",0.0023,0.0010,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
 
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
 
   RooRealVar signal("signal"," ",500,10,1000000);//event number
   RooRealVar signal2("signal2"," ",100,1,1000000);//event number
   RooRealVar background("background"," ",10,0,1000000);
 
   RooAddPdf *sum;
   RooDataHist *datahist;
   RooDataSet *dataset;
   RooPlot *xframe;  
  
   TTree *dataraw = new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");
   TTree *vars = new TTree("vars","vars");
   double phi,phi1,phi2;
   double costheta,costheta1,costheta2;
   double p1,p2;
   vars->Branch("phi",&phi,"phi/D");
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta",&costheta,"costheta/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"p1/D");
   vars->Branch("p2",&p2,"p2/D");

   int Npart=10;
   double pcut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
   double start=0;
   double stop =2.0;
   for(int i=0;i<Npart+1;i++){
     pcut[i] = (stop-start)/Npart*i+start;
   }

   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;
  
   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different range
   TH1D *hmKs[Npart];
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%d",partj);
     hmKs[partj] = new TH1D(name,name,100,kslow,ksup);
   }
 
   // loop data
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	  
	  int ipip=0,ipim=0;
	  double chi2=10000000;
	  if (npip+npim>2){
	  for (int i=0; i<npip;i++)
	  for (int j=0; j<npim;j++){
	    if (Ks_decayL[i*npim+j]<cutd) continue;
	    //hmass->Fill(Ks_mass[i]);
		//double mass = CalInvMass(mpi,pippx[i],pippy[i],pippz[i],
		//						 mpi,pimpx[j],pimpy[j],pimpz[j],0,0,0);
		if (Ks_chi2[i*npim+j]<chi2){
		  chi2 = Ks_chi2[i*npim+j];
		  ipip=i;
		  ipim=j;
		}
	  }
      mass = CalInvMass(mpi,pippx[ipip],pippy[ipip],pippz[ipip],mpi,pimpx[ipim],pimpy[ipim],pimpz[ipim]);
	  if (mass>kslow && mass<ksup){
	    p1=CalMom(pippx[ipip],pippy[ipip],pippz[ipip]);
	    p2=CalMom(pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        costheta1 = pippz[ipip]/p1;
		costheta2 = pimpz[ipim]/p2;
		if (pippy[ipip]>0)  phi1 = acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
		else  phi1 = 2*TMath::Pi()-acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
		if (pimpy[ipim]>0)  phi2 = acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
		else  phi2 = 2*TMath::Pi()-acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
		vars->Fill();
	  }

	  if (mass>kslow-0.02 && mass<ksup+0.02){
	    ppx.push_back(pippx[ipip]);
	    ppy.push_back(pippy[ipip]);
	    ppz.push_back(pippz[ipip]);
	    mpx.push_back(pimpx[ipim]);
	    mpy.push_back(pimpy[ipim]);
	    mpz.push_back(pimpz[ipim]);
	  }
      //if (fabs(p1-p2)<0.0001) continue;
	  //hmass->Fill(Ks_mass[ipip*npim+ipim]);
	  //h2pos->Fill(ipip,ipim);
	  int parti,partj;
      partj = (int)((p2-start)/(stop-start)*Npart);
      if ( partj>=Npart || partj<0 ) continue;
      if ( p2>pcut[partj] && p2<pcut[partj+1] )
      {
        hmKs[partj]->Fill(mass);
      }

	  //if (Ks_mass[ipip*npim+ipim]>0.48&&Ks_mass[ipip*npim+ipim]<0.52)
	  //  h2p->Fill(p1,p2);
	  //hmass->Fill(mass);
	  }
   }
   std::cout<<"fffffffffffffffff"<<std::endl;
   vars->Write();
   for (int partj=0;partj<Npart;partj++){
     hmKs[partj]->Write();
     if (hmKs[partj]->GetEntries() > 100
       &&hmKs[partj]->GetMaximumBin()> 30 //check the the peak position,
       &&hmKs[partj]->GetMaximumBin()< 70 
       ){
       partmap.push_back(partj);
     }
   }

   //hmass->Draw();
   //return hmass;
   //TCanvas *c2 = new TCanvas();
   //h2p->Draw();
   //h2pos->Draw();
  
   const int pointNo=20;
   double factor=0.995;
   double factorstep=(1.-factor)*2/pointNo;
   int fittimes=0;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   double peakvalue = mparticle;

   std::cout<<"part map size is "<<partmap.size()<<std::endl;
   for (int loopj=0;loopj<partmap.size();loopj++){
     int partj = partmap.at(loopj);
   std::cout<<"part is "<<partj<<std::endl;
	 double factori=1.0;
	 factor = 0.995;
     fittimes=0;

	 for (fittimes=0; fittimes<pointNo;fittimes++){
	   xframe = x.frame(Title("fit Ks"));
	   dataraw->Reset();
       std::cout<<"factor is "<<factor<<std::endl;

       //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   std::cout<<"ppx size is "<<ppx.size()<<std::endl;
       for (Long64_t jin=0; jin<ppx.size();jin++) {
          
          mass = CalInvMass(mpi,factori*ppx[jin],factori*ppy[jin],factori*ppz[jin],mpi,factor*mpx[jin],factor*mpy[jin],factor*mpz[jin]);
          p2=CalMom(mpx[jin],mpy[jin],mpz[jin]);
		  if (partj<0 || partj>Npart) continue;
		  if (p2 < pcut[partj] || p2> pcut[partj+1]) continue;
		  if (mass>kslow && mass<ksup){
		    dataraw->Fill();
		  }
       }// loop data end
       // fit data

   std::cout<<"dataraw size is "<<dataraw->GetEntries()<<std::endl;
	   sprintf(name,"data_Ks_%02d",fittimes);
       dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*dataraw));
       sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
       Npar=8;
	   sum->fitTo(*dataset,Range(kslow,ksup));
       dataset->plotOn(xframe,Binning(nBins));
       sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
       sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
       sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
       sum->plotOn(xframe);
	   xframe->Draw();
       TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
       pt->SetBorderSize(0);
       pt->SetFillStyle(4000);
       pt->SetTextAlign(12);
       pt->SetTextFont(42);
       pt->SetTextSize(0.035);
       sprintf(name,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
       pt->AddText(name);
       sprintf(name,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
       pt->AddText(name);
       sprintf(name,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
       pt->AddText(name);
       sprintf(name,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
       pt->AddText(name);
       sprintf(name,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
       pt->AddText(name);
       sprintf(name,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
       pt->AddText(name);
       sprintf(name,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
       pt->AddText(name);
       pt->Draw();
       //c1->Update();
	   sprintf(name,"part%d_fitFor_%dth_time",partj,fittimes);
	   c1->SetName(name);
	   c1->Write();

       factors[fittimes]=factor;
       factorserr[fittimes]=0;
       deltapeaks[fittimes] = mean.getValV() - peakvalue;
       deltapeakserr[fittimes] = mean.getError();
      
	   factor += factorstep;
       
	   delete sum;
       delete dataset;
       delete xframe;
 
	 }// n point end
	 //will fit linear, get factor needed
    c1->Clear();
    TGraphErrors *graph = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
    graph->SetTitle("delta peak");
    graph->Draw("AP");
    graph->SetMarkerStyle(5);
    gStyle->SetOptFit(1111);
    facfit->SetParameters(1,0.4);
    facfit->SetParNames("factor","slope");
    graph->Fit(facfit,"","",factors[0],factors[pointNo-1]);
    factor = facfit->GetParameter(0);
    TPaveText *pt1 = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(4000);
    pt1->SetTextAlign(12);
    pt1->SetTextFont(42);
    pt1->SetTextSize(0.035);
    sprintf(name,"factor = %1.6f #pm %1.6f",facfit->GetParameter(0),facfit->GetParError(0));
    pt1->AddText(name);
    sprintf(name,"slope = %1.6f #pm %1.6f",facfit->GetParameter(1),facfit->GetParError(1));
    pt1->AddText(name);
	pt1->Draw();
 
    sprintf(name,"part%d_factors",partj);
	c1->SetName(name);
	c1->Write();

    // try to use the best factor to fit
	xframe = x.frame(Title("fit Ks"));
	dataraw->Reset();
    std::cout<<"factor is "<<factor<<std::endl;
       
	for (Long64_t jin=0; jin<ppx.size();jin++) {
       mass = CalInvMass(mpi,factori*ppx[jin],factori*ppy[jin],factori*ppz[jin],mpi,factor*mpx[jin],factor*mpy[jin],factor*mpz[jin]);
       p2 = CalMom(mpx[jin],mpy[jin],mpz[jin]);
	   if (partj<0 || partj>Npart) continue;
	   if (p2 < pcut[partj] || p2> pcut[partj+1]) continue;
	   if (mass>kslow && mass<ksup){
	     dataraw->Fill();
	   }
    }// loop data end

    // fit data
	sprintf(name,"data_Ks_part%d",partj);
    dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*dataraw));
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
    Npar=8;
	sum->fitTo(*dataset,Range(kslow,ksup));
    dataset->plotOn(xframe,Binning(nBins));
    sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
    sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
    sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
    sum->plotOn(xframe);
	xframe->Draw();
    TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(4000);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    sprintf(name,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
    pt->AddText(name);
    sprintf(name,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
    pt->AddText(name);
    sprintf(name,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
    pt->AddText(name);
    sprintf(name,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
    pt->AddText(name);
    sprintf(name,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
    pt->AddText(name);
    sprintf(name,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
    pt->AddText(name);
    sprintf(name,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
    pt->AddText(name);
    double factor4err=TMath::Sqrt(TMath::Power(mean.getError()/facfit->GetParameter(1),2) + TMath::Power(facfit->GetParError(0),2));
    sprintf(name,"factor = %.6f #pm %.6f",factor,factor4err);
    pt->AddText(name);

    pt->Draw();
    c1->Update();
	sprintf(name,"part%d_final_fit",partj);
	c1->SetName(name);
	c1->Write();
   
	delete sum;
    delete dataset;
    delete xframe;
 
   }// n part end

}

#ifdef KsAlg_cxx
KsAlg::KsAlg(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("multipi_001_131210_0034011.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("multipi_001_131210_0034011.root");
      }
      f->GetObject("gepep_2pi",tree);

   }
   Init(tree);
}

KsAlg::KsAlg(std::string filename) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TFile *f = new TFile(filename.c_str());
   TTree *tree;
   f->GetObject("gepep_2pi",tree);

   Init(tree);
}

KsAlg::~KsAlg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t KsAlg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t KsAlg::LoadTree(Long64_t entry)
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

void KsAlg::Init(TTree *tree)
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
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("nKsCan", &nKsCan, &b_nKsCan);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   fChain->SetBranchAddress("Ks_mass", Ks_mass, &b_Ks_mass);
   fChain->SetBranchAddress("Ks_ratio", Ks_ratio, &b_Ks_ratio);
   fChain->SetBranchAddress("Ks_decayL", Ks_decayL, &b_Ks_decayL);
   fChain->SetBranchAddress("Ks_chi2", Ks_chi2, &b_Ks_chi2);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t KsAlg::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void KsAlg::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t KsAlg::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef KsAlg_cxx
