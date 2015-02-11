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

   std::vector<double> ppx,ppy,ppz;
   std::vector<double> mpx,mpy,mpz;

   TH1D *hmass  = new TH1D("hmass" ,"Ks mass",200,0.45,0.55);
   TH1D *hmassU = new TH1D("hmassU","Ks mass",200,0.45,0.55);
   TH1D *hmassc = new TH1D("hmassc","Ks mass",200,0.45,0.55);
   //TH2D *h2pos = new TH2D("h2pos","pos",10, 0,5, 10,0,5);
   TH2D *h2p = new TH2D("h2p","P",200, 0,2, 200,0,2);
   int nBins=100;
   double mpi=0.13957;
   double mparticle=0.497614;
   double kslow=0.43;
   double ksup =0.57;
   TF1 *facfit = new TF1("facfit",line2,kslow,ksup,2);
   char fname[100];
   sprintf(fname,"%s/plot_ks.root",outputdir.c_str());
   TFile *f=new TFile(fname,"RECREATE");
   
   sprintf(fname,"%s/pars.txt",outputdir.c_str());
   ofstream ofpar(fname,std::ios::app);
   sprintf(fname,"%s/parspur.txt",outputdir.c_str());
   ofstream ofparpur(fname,std::ios::app);

   TCanvas *c1 = new TCanvas();
   int Npar;

   // roofit variables and functions
   RooRealVar x("x","energy",mparticle,kslow,ksup,"GeV");
   RooRealVar mean("mean","mean of gaussian",  kslow,ksup);
   RooRealVar mean2("mean2","mean of gaussian",kslow,ksup);
   RooRealVar sigma("sigma","width of gaussian",0.0017,0.0010,0.0050);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
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
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

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

   int Npart=1; // 20
   int realsize=0;
   double partid[Npart];
   double parter[Npart];
   double corfac[Npart];
   double corerr[Npart];
   double start=0;
   double stop =2.0;
   double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   pcut[0]=0.1;
   pcut[1]=0.4;
   //pcut[2]=0.9;
 //pcut[0] =0.0;
 //pcut[1] =0.05;
 //pcut[2] =0.10;
 //pcut[3] =0.15;
 //pcut[4] =0.20;
 //pcut[5] =0.25;
 //pcut[6] =0.30;
 //pcut[7] =0.35;
 //pcut[8] =0.40;
 //pcut[9] =0.45;
 //pcut[10]=0.50;
 //pcut[11]=0.60;
 //pcut[12]=0.70;
 //pcut[13]=0.80;
 //pcut[14]=0.90;
 //pcut[15]=1.00;
 //pcut[16]=1.20;
 //pcut[17]=1.40;
 //pcut[18]=1.60;
 //pcut[19]=1.80;
 //pcut[20]=2.00;
   //for(int i=0;i<Npart+1;i++){
   //  pcut[i] = (stop-start)/Npart*i+start;
   //}

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
      //double chi2=10000000;
      //if (npip+npim<3) continue;
      for (int ipip=0; ipip<npip;ipip++)
      for (int ipim=0; ipim<npim;ipim++){
        mass = CalInvMass(mpi,pippx[ipip],pippy[ipip],pippz[ipip],mpi,pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        double chi2 = m_chisq0[ipip*npim+ipim]+m_chisqvtxnd[ipip*npim+ipim];
        if (chi2 > 100) continue;
        double ratio = m_decayL[ipip*npim+ipim]/m_decayLerr[ipip*npim+ipim];
        if (ratio<3) continue;
        if (Ks_id[ipip*npim+ipim]==0){
          hmassU->Fill(mass);
          continue;
        }
        hmass->Fill(Mpippim[ipip*npim+ipim]);
        hmassc->Fill(mass);
        p1=CalMom(pippx[ipip],pippy[ipip],pippz[ipip]);
        p2=CalMom(pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        if (p1+p2<0.4 || p1+p2 > 0.6) continue;
        costheta1 = pippz[ipip]/p1;
        costheta2 = pimpz[ipim]/p2;
        if (pippy[ipip]>0)  phi1 = acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
        else  phi1 = 2*TMath::Pi()-acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
        if (pimpy[ipim]>0)  phi2 = acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
        else  phi2 = 2*TMath::Pi()-acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
        if (mass>kslow && mass<ksup){
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
        //int parti,partj;
        //partj = (int)((p2-start)/(stop-start)*Npart);
        //if ( partj>=Npart || partj<0 ) continue;
        for ( int partj=0;partj<Npart;partj++){
          if (p1>pcut[partj] && p1<pcut[partj+1] && p2>pcut[partj] && p2<pcut[partj+1] )
          {
            if (mass>kslow-0.02 && mass<ksup+0.02)
              hmKs[partj]->Fill(mass);
            break;
          }
        }

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

   hmass->Write();
   hmassU->Write();
   hmassc->Write();
   //return hmass;
   //TCanvas *c2 = new TCanvas();
   //h2p->Draw();
   //h2pos->Draw();
  
   const int pointNo=10;
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
     factori=factor;
     fittimes=0;

     for (fittimes=0; fittimes<pointNo;fittimes++){
       xframe = x.frame(Title("fit Ks"));
       dataraw->Reset();
       std::cout<<"factor is "<<factor<<std::endl;

       //for (Long64_t jentry=0; jentry<nentries;jentry++) {
       std::cout<<"ppx size is "<<ppx.size()<<std::endl;
       for (Long64_t jin=0; jin<ppx.size();jin++) {
          
          mass = CalInvMass(mpi,factori*ppx[jin],factori*ppy[jin],factori*ppz[jin],mpi,factor*mpx[jin],factor*mpy[jin],factor*mpz[jin]);
          p1=CalMom(ppx[jin],ppy[jin],ppz[jin]);
          p2=CalMom(mpx[jin],mpy[jin],mpz[jin]);
          if (partj<0 || partj>Npart) continue;
          if (p1 < pcut[partj] || p1> pcut[partj+1] || p2 < pcut[partj] || p2> pcut[partj+1]) continue;
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
       mean.setVal(peakvalue+0.165*(factor-1));
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
       deltapeaks[fittimes] = mean.getVal() - peakvalue;
       deltapeakserr[fittimes] = mean.getError();
      
       factor += factorstep;
       factori = factor;

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
    factori=factor;
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
      p1 = CalMom(ppx[jin],ppy[jin],ppz[jin]);
      p2 = CalMom(mpx[jin],mpy[jin],mpz[jin]);
      if (partj<0 || partj>Npart) continue;
      if (p1 < pcut[partj] || p1> pcut[partj+1] || p2 < pcut[partj] || p2> pcut[partj+1]) continue;
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
    double factor4err=TMath::Sqrt(TMath::Power(mean.getError()/facfit->GetParameter(1),2) + TMath::Power(facfit->GetParError(0),2));
    sprintf(name,"factor = %.6f #pm %.6f",factor,factor4err);
    pt->AddText(name);

    pt->Draw();
    c1->Update();
    sprintf(name,"part%d_final_fit",partj);
    c1->SetName(name);
    c1->Write();

    partid[loopj] = pcut[partj]+(pcut[partj+1]-pcut[partj])/2;
    parter[loopj] = 0;
    corfac[loopj] = factor;
    corerr[loopj] = factor4err;
   
    delete sum;
    delete dataset;
    delete xframe;
 
  }// n part end

  realsize = partmap.size();
  TGraphErrors *graph = new TGraphErrors(realsize,partid,corfac,parter,corerr);
  graph->SetTitle("compare factors");
  graph->Draw("AP");
  graph->SetMarkerStyle(5);
  gStyle->SetOptFit(1111);
  c1->SetName("compare_factors");
  c1->Write();

  for (int i=0;i<realsize;i++){
    ofpar<<"p="<<partid[i]<<"\tfactor: "<<corfac[i]<<"\t +/- \t"<< corerr[i]<<std::endl;
    ofparpur<<partid[i]<<"\t"<<corfac[i]<<"\t"<< corerr[i]<<std::endl;
  }





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
