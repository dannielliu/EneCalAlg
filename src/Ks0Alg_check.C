#define Ks0Alg_cxx
#include "Ks0Alg.h"
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
#include "RooMsgService.h"
#include <vector>
#include <fstream>
extern std::string outputdir;
using namespace RooFit;


void Ks0Alg::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Ks0Alg.C
//      Root > Ks0Alg t
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

   //Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //   Long64_t ientry = LoadTree(jentry);
   //   if (ientry < 0) break;
   //   nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   //}


   std::cout<<"loop"<<std::endl;

   TH1D *hmass  = new TH1D("hmass" ,"Ks mass",200,0.45,0.55);
   TH1D *hmassU = new TH1D("hmassU","Ks mass",200,0.45,0.55);
   TH1D *hmassc = new TH1D("hmassc","Ks mass",200,0.45,0.55);
   //TH2D *h2pos = new TH2D("h2pos","pos",10, 0,5, 10,0,5);
   TH2D *h2p = new TH2D("h2p","P",200, 0,2, 200,0,2);
   int nBins=100;
   double mpi=0.13957;
   double mparticle=0.497614;
   double kslow=0.45;
   double ksup =0.55;
   TF1 *facfit = new TF1("facfit",line2,kslow,ksup,2);
   char fname[100];
   sprintf(fname,"%s/plot_ks.root",outputdir.c_str());
   TFile *f=new TFile(fname,"RECREATE");
   
   //sprintf(fname,"%s/pars.txt",outputdir.c_str());
   //ofstream ofpar(fname,std::ios::app);
   //sprintf(fname,"%s/parspur.txt",outputdir.c_str());
   //ofstream ofparpur(fname,std::ios::app);

   TCanvas *c1 = new TCanvas();
   int Npar;

   // roofit variables and functions
   RooRealVar x("x","energy",mparticle,kslow,ksup,"GeV");
   RooRealVar mean("mean","mean of gaussian",mparticle,  kslow,ksup);
   //RooRealVar mean2("mean2","mean of gaussian",kslow,ksup);
   RooRealVar sigma("sigma","width of gaussian",0.0017,0.0010,0.0050);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooRealVar brewid("brewid","width of breit wigner",0.0023,0.0010,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
 
   RooRealVar a0("a0","coefficient #0",100,-1000000,1000000);
   RooRealVar a1("a1","coefficient #1",-1,-1000000,1000000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
 
   RooRealVar signal("signal"," ",500,10,10000000);//event number
   RooRealVar signal2("signal2"," ",100,1,10000000);//event number
   RooRealVar background("background"," ",10,0,10000000);
 
   RooAddPdf *sum;
   RooDataHist *datahist;
   RooDataSet *dataset;
   RooPlot *xframe;  
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   TTree *dataraw = new TTree("dataraw" ,"dataraw");
   TTree *datarawl= new TTree("datarawl" ,"dataraw");
   TTree *datarawu= new TTree("datarawu" ,"dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");
   datarawl->Branch("x",&mass,"x/D");
   datarawu->Branch("x",&mass,"x/D");
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
   vars->Branch("mass",&mass,"mass/D");

   int Npart=20;
   double start=0.1;
   double stop =1.0;
   double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   //pcut[0]=0.1;
   //pcut[1]=0.4;
   //pcut[2]=0.9;
   double facmap[Npart];
   double facemap[Npart];

//  set normal factor in (0.1, 0.4) to 1.00082, get factor in different range
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00263 ;  facemap[2] =0.000106524;
   pcut[3] =0.15;  facmap[3] =1.00126 ;  facemap[3] =5.65503e-05;
   pcut[4] =0.20;  facmap[4] =1.00087 ;  facemap[4] =5.36283e-05;
   pcut[5] =0.25;  facmap[5] =1.00060 ;  facemap[5] =4.96854e-05;
   pcut[6] =0.30;  facmap[6] =1.00063 ;  facemap[6] =5.4565e-05 ;
   pcut[7] =0.35;  facmap[7] =1.00034 ;  facemap[7] =6.70578e-05;
   pcut[8] =0.40;  facmap[8] =0.999913;  facemap[8] =8.50691e-05;
   pcut[9] =0.45;  facmap[9] =0.999925;  facemap[9] =0.000103853;
   pcut[10]=0.50;  facmap[10]=0.999614;  facemap[10]=0.000101297;
   pcut[11]=0.60;  facmap[11]=0.999918;  facemap[11]=9.75854e-05;
   pcut[12]=0.70;  facmap[12]=1.000260;  facemap[12]=0.000228099;
   pcut[13]=0.80;  facmap[13]=1.000340;  facemap[13]=0.000356495;
   pcut[14]=0.90;  facmap[14]=1.00056 ;  facemap[14]=0.000554875;
   pcut[15]=1.00;  facmap[15]=1.00185 ;  facemap[15]=0.000650392;
   pcut[16]=1.20;  facmap[16]=0.994443;  facemap[16]=0.00288715 ;
   pcut[17]=1.40;  facmap[17]=0.994232;  facemap[17]=0.00301238 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00263 ;  facemap[2] =0.000104711;
   pcut[3] =0.15;  facmap[3] =1.00127 ;  facemap[3] =5.66523e-05;
   pcut[4] =0.20;  facmap[4] =1.00089 ;  facemap[4] =5.40327e-05;
   pcut[5] =0.25;  facmap[5] =1.00061 ;  facemap[5] =4.96881e-05;
   pcut[6] =0.30;  facmap[6] =1.00063 ;  facemap[6] =5.45699e-05;
   pcut[7] =0.35;  facmap[7] =1.00034 ;  facemap[7] =6.70568e-05;
   pcut[8] =0.40;  facmap[8] =0.999916;  facemap[8] =8.51944e-05;
   pcut[9] =0.45;  facmap[9] =0.999916;  facemap[9] =0.00010567 ;
   pcut[10]=0.50;  facmap[10]=0.99955 ;  facemap[10]=8.6275e-05 ;
   pcut[11]=0.60;  facmap[11]=0.999922;  facemap[11]=0.000112065;
   pcut[12]=0.70;  facmap[12]=1.00024 ;  facemap[12]=0.000229724;
   pcut[13]=0.80;  facmap[13]=1.00037 ;  facemap[13]=0.000353896;
   pcut[14]=0.90;  facmap[14]=1.00057 ;  facemap[14]=0.000551636;
   pcut[15]=1.00;  facmap[15]=1.00208 ;  facemap[15]=0.000718673;
   pcut[16]=1.20;  facmap[16]=0.994435;  facemap[16]=0.00286605 ;
   pcut[17]=1.40;  facmap[17]=0.98542 ;  facemap[17]=0.00410906 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;

   //for(int i=0;i<Npart+1;i++){
   //  pcut[i] = (stop-start)/Npart*i+start;
   //}
 
   char name[100];
   double factori=1;
   double factorj=1;
   // ~~~~~~~~ draw nxn histogram, m distribution in different range
 
   // loop data
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
      int ipip=0,ipim=0;
      for (int ipip=0; ipip<npip;ipip++)
      for (int ipim=0; ipim<npim;ipim++){
        mass = CalInvMass(mpi,pippx[ipip],pippy[ipip],pippz[ipip],mpi,pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        double chi2 = m_chisq0[ipip*npim+ipim]+m_chisqvtxnd[ipip*npim+ipim];
        if (chi2 > 100) continue;
        double ratio = m_decayL[ipip*npim+ipim]/m_decayLerr[ipip*npim+ipim];
        if (ratio<2) continue;
        if (Ks_id[ipip*npim+ipim]==0){
          hmassU->Fill(mass);
          continue;
        }
        hmass->Fill(Mpippim[ipip*npim+ipim]);
        hmassc->Fill(mass);
        p1=CalMom(pippx[ipip],pippy[ipip],pippz[ipip]);
        p2=CalMom(pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        if (p1<0.1 || p1>1.0) continue;
        if (p2<0.1 || p2>1.0) continue;
        costheta1 = pippz[ipip]/p1;
        costheta2 = pimpz[ipim]/p2;
        if (pippy[ipip]>0)  phi1 = acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
        else  phi1 = 2*TMath::Pi()-acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
        if (pimpy[ipim]>0)  phi2 = acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
        else  phi2 = 2*TMath::Pi()-acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
        if (mass>kslow && mass<ksup){
          vars->Fill();
        }
        // allocate factor according to p, average factor
        for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            factori = facmap[i];
            break;
          }
        }

        for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            factorj = facmap[i];
            break;
          }
        }
        mass = CalInvMass(mpi,factori*pippx[ipip],factori*pippy[ipip],factori*pippz[ipip],
			  mpi,factorj*pimpx[ipim],factorj*pimpy[ipim],factorj*pimpz[ipim]);
        if (mass>kslow && mass<ksup){
          dataraw->Fill();
        }
 
        // allocate factor according to p, factor at low edge
        for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            factori = facmap[i]-facemap[i];
            break;
          }
        }

       for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            factorj = facmap[i]-facemap[i];
            break;
          }
        }
        mass = CalInvMass(mpi,factori*pippx[ipip],factori*pippy[ipip],factori*pippz[ipip],
			  mpi,factorj*pimpx[ipim],factorj*pimpy[ipim],factorj*pimpz[ipim]);
        if (mass>kslow && mass<ksup){
          datarawl->Fill();
        }

        // allocate factor according to p, factor at up edge
        for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            factori = facmap[i]+facemap[i];
            break;
          }
        }

       for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            factorj = facmap[i]+facemap[i];
            break;
          }
        }
        mass = CalInvMass(mpi,factori*pippx[ipip],factori*pippy[ipip],factori*pippz[ipip],
			  mpi,factorj*pimpx[ipim],factorj*pimpy[ipim],factorj*pimpz[ipim]);
        if (mass>kslow && mass<ksup){
          datarawu->Fill();
        }
       //int parti,partj;
        //partj = (int)((p2-start)/(stop-start)*Npart);

      }
      
   }// loop data end
   std::cout<<"fffffffffffffffff"<<std::endl;
   vars->Write();
   
   hmass->Write();
   hmassU->Write();
   hmassc->Write();
   //return hmass;
   
   // average factor part
   xframe = x.frame(Title("fit Ks"));
   // fit data
   sprintf(name,"data_Ks");
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
   TPaveText *pt = new TPaveText(0.65,0.50,0.95,0.95,"BRNDC");
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
   c1->Update();
   sprintf(name,"final_fit");
   c1->SetName(name);
   c1->Write();

   delete sum;
   delete dataset;
   delete xframe;
   //average part end
 
   // low edge part
   xframe = x.frame(Title("fit Ks"));
   // fit data
   sprintf(name,"data_Ks_low");
   dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*datarawl));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
   Npar=8;
   sum->fitTo(*dataset,Range(kslow,ksup));
   dataset->plotOn(xframe,Binning(nBins));
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   pt = new TPaveText(0.65,0.50,0.95,0.95,"BRNDC");
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
   c1->Update();
   sprintf(name,"final_fit_low");
   c1->SetName(name);
   c1->Write();

   delete sum;
   delete dataset;
   delete xframe;
   // low part end

//  up part 
 
   xframe = x.frame(Title("fit Ks"));
   // fit data
   sprintf(name,"data_Ks_up");
   dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*datarawu));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
   Npar=8;
   sum->fitTo(*dataset,Range(kslow,ksup));
   dataset->plotOn(xframe,Binning(nBins));
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   pt = new TPaveText(0.65,0.50,0.95,0.95,"BRNDC");
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
   c1->Update();
   sprintf(name,"final_fit_up");
   c1->SetName(name);
   c1->Write();

   delete sum;
   delete dataset;
   delete xframe;
// up part end



}

#ifdef Ks0Alg_cxx
Ks0Alg::Ks0Alg(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ksto2pi_583_140226_0035862.root");
      if (!f) {
         f = new TFile("Ksto2pi_583_140226_0035862.root");
      }
      tree = (TTree*)gDirectory->Get("Ks_info");

   }
   Init(tree);
}

Ks0Alg::~Ks0Alg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Ks0Alg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ks0Alg::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ks0Alg::Init(TTree *tree)
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

   fChain->SetBranchAddress("rec_truth_Mks", &rec_truth_Mks, &b_rec_truth_Mks);
   fChain->SetBranchAddress("rec_truth_Pks", &rec_truth_Pks, &b_rec_truth_Pks);
   fChain->SetBranchAddress("rec_truth_Eks", &rec_truth_Eks, &b_rec_truth_Eks);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("nKsCan", &nKsCan, &b_nKsCan);
   fChain->SetBranchAddress("m_chisq0", m_chisq0, &b_m_chisq0);
   fChain->SetBranchAddress("m_chisqvtxnd", m_chisqvtxnd, &b_m_chisqvtxnd);
   fChain->SetBranchAddress("m_decayL", m_decayL, &b_m_decayL);
   fChain->SetBranchAddress("m_decayLerr", m_decayLerr, &b_m_decayLerr);
   fChain->SetBranchAddress("m_ctau", m_ctau, &b_m_ctau);
   fChain->SetBranchAddress("Ks_id", Ks_id, &b_Ks_id);
   fChain->SetBranchAddress("Mpippim", Mpippim, &b_Mpippim);
   fChain->SetBranchAddress("Ppippim", Ppippim, &b_Ppippim);
   fChain->SetBranchAddress("Epippim", Epippim, &b_Epippim);
   fChain->SetBranchAddress("Thepippim", Thepippim, &b_Thepippim);
   fChain->SetBranchAddress("Phipippim", Phipippim, &b_Phipippim);
   fChain->SetBranchAddress("Ppip", Ppip, &b_Ppip);
   fChain->SetBranchAddress("Ppim", Ppim, &b_Ppim);
   fChain->SetBranchAddress("Thepip", Thepip, &b_Thepip);
   fChain->SetBranchAddress("Thepim", Thepim, &b_Thepim);
   fChain->SetBranchAddress("Phipip", Phipip, &b_Phipip);
   fChain->SetBranchAddress("Phipim", Phipim, &b_Phipim);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   Notify();
}

Bool_t Ks0Alg::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ks0Alg::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ks0Alg::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Ks0Alg_cxx
