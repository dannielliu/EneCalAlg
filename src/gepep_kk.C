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
//#include "RooLandau.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
//#include <iostream>
extern std::string outputdir;
//using RooFit::Title;
//using RooFit::Components;
//using RooFit::LineStyle;
//using RooFit::LineColor;
//using RooFit::Range;
using namespace RooFit;

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
   double peakvalue=1.86486;// mphi
   // phi -> K K
   //double philow=0.995;
   //double phiup=1.045;
   //double peakvalue=1.019455;// mphi
   int pointNo=10;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   double factor=factorstart;
   double factorstep=(1.-factor)*2/pointNo;
   
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

   RooBreitWigner brewig("brewig","brewig",x,mean,sigma);
   //RooLandau landau("landau","landau",x,mean,sigma);

   RooRealVar tao("tao","#tao",-200.,-1000.,-100.);
   RooExponential expo("expo","expo",x,tao);

   RooRealVar a0("a0","coefficient #0",100,100,100000);
   RooRealVar a1("a1","coefficient #1",-50,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooRealVar signal("signal"," ",1200,10,1000000);//event number
   RooRealVar signal2("signal2"," ",1200,0,1000000);//event number
   RooRealVar background("background"," ",200,0,1000000);
   RooRealVar background2("background2"," ",200,0,1000000);
   
   RooPlot *xframe;
   RooDataHist *data_k;
   RooAddPdf *sum;
   int Npar;
   //sum->selfNormalized();
   //RooFitResult *result;
   
   ofstream ofpar;
   ofpar.open("parkk.txt",std::ios::app);
   ofpar<<"kk algrithm: will give factors for kaon"<<std::endl;
   ofstream ofpardetail;
   ofpardetail.open("detail.txt",std::ios::app);
   ofstream purepar;
   purepar.open("par");
   char fname[100];
   sprintf(fname,"%s/plot_kk.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");
   TTree *dataraw=new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");
   RooDataSet *dataset = new RooDataSet("dataset","data",dataraw,x);

   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1    = new TH1D("h1","2 kaon invariant mass",nBins,philow,phiup);
   TCanvas *c1 = new TCanvas("","",800,600);

   TH2D *h2p = new TH2D("h2p","KK momentum",200,0,2,200,0,2);
   TH3D *the = new TH3D("the","K momentum",200,0,2,200,0,2,100,0,TMath::Pi());
   TH2D *thedis = new TH2D("thedis","K momentum",200,0,2,100,0,TMath::Pi());
   TH2D *thedis2 = new TH2D("thedis2","K momentum",100,0,TMath::Pi(),100,0,TMath::Pi());
   int Npart=1;
   double thetacut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
   double start=0.0;
   double stop =TMath::Pi();;
   for(int i=0;i<Npart+1;i++){
     thetacut[i] = (stop-start)/Npart*i+start;
   }
   double m0=peakvalue;
   double sigma_m=0.0024;//0.0024 for phi,
   double width = 10.*sigma_m;
   double mparticle=0.493677;
   std::vector<std::pair<int,int> > partmap;
   std::vector<std::pair<int,double> > facmap;

   Long64_t nbytes = 0, nb = 0;
   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   TH1D *hmtheta[Npart][Npart];
   for (int parti=0;parti<Npart;parti++){
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%d_part%d",parti,partj);
     hmtheta[parti][partj] = new TH1D(name,name,100,philow,phiup);
   }
   }

   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
  
     int parti,partj;
     double p1,p2;
     double theta,theta1,theta2;
     mass = CalInvMass(mk,kappx,kappy,kappz,mk,kampx,kampy,kampz);
     
     // total invariant mass, D0 -> k- pi+
     p1 = CalMom(kappx,kappy,kappz);
     p2 = CalMom(kampx,kampy,kampz);
     theta1 = acos(kappz/p1);
     theta2 = acos(kampz/p2);
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
         theta = acos((kappx*kampx+kappy*kampy+kappz*kampz)/(p1*p2));
         the->Fill(p1,p2,theta);
         hmtheta[parti][partj]->Fill(mass);
       }
     // if (Cut(ientry) < 0) continue;
   }
   h2p->Write();
   thedis->Write();
   thedis2->Write();
   the->Write();
   for (int parti=0;parti<Npart;parti++)
   for (int partj=0;partj<Npart;partj++){
     hmtheta[parti][partj]->Write();
     std::cout<<"processed part "<<parti<<", part "<<partj<<std::endl;
     if (hmtheta[parti][partj]->GetEntries() > 100
       &&hmtheta[parti][partj]->GetMaximumBin()> 30 //check the the peak position,
       &&hmtheta[parti][partj]->GetMaximumBin()< 70 //make sure there is a peak in the region
       ){
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

   // for saving the fit result
   char fitepsname[100];
   char fiteps_start[100];
   char fiteps_stop[100];
   sprintf(fitepsname,"%s/fitkk_part%d_part%d.eps",outputdir.c_str(),parti,partj);
   sprintf(fiteps_start,"%s[",fitepsname);
   sprintf(fiteps_stop,"%s]",fitepsname);
   c1->Print(fiteps_start);
   
   int fittimes = 0;
   factor = factorstart;
   for (int i=0;i<pointNo;i++){
      //xframe->Clear();
      xframe = x.frame(Title("fit kaon"));

      //h1->Reset();
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
         //if(ngam>0) continue;
         double p1,p2;
         double theta1,theta2;
         double factors[2]={1,1};
         factori==0? factors[0]=factor : factors[0]=factori;
         factorj==0? factors[1]=factor : factors[1]=factorj;
         p1 = CalMom(kappx,kappy,kappz);
         p2 = CalMom(kampx,kampy,kampz);
         theta1 = acos(kappz/p1);
         theta2 = acos(kampz/p2);
         if (theta1>thetacut[parti]&&theta1<thetacut[parti+1]
           &&theta2>thetacut[partj]&&theta2<thetacut[partj+1]){
           mass = CalInvMass(mk,kappx,kappy,kappz,mk,kampx,kampy,kampz,-2,factors);
           //h1->Fill(mass);
           if (mass>philow && mass<phiup){
             dataraw->Fill();
           }
         }
         // if (Cut(ientry) < 0) continue;
      }
      //dataraw->Write();

      dataset = new RooDataSet("dataset","data",dataraw,x);
      char tmpchr[100];
      sprintf(tmpchr,"data_k_%02d",fittimes);
      //data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
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
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;

      //sprintf(tmpchr,"data_k_%d.eps",fittimes);
      //h1->Draw();
      //c1->Print(tmpchr);
      
      // save pars
      factors[i]=factor;
      factorserr[i]=0;
      deltapeaks[i] = mean.getValV() - peakvalue;
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
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   factor = facfit->GetParameter(0);
   double factorerr=deltapeakserr[9]/facfit->GetParameter(1);
   ofpar<<facfit->GetParameter(0)<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<std::endl;
   ofpar<<signal.getValV()<<"\t"<<signal.getError()<<std::endl;
   purepar<<"\t"<<factor<<"\t"<<factorerr<<std::endl;
   //std::cout<<"fit factor: "<<factor1<<", error is "<<factor1err<<std::endl;
   sprintf(name,"factors_kk_part%d_part%d.eps",parti,partj);
   c1->Print(name);
   sprintf(name,"factors_kk_part%d_part%d",parti,partj);
   graph1->SetName(name);
   graph1->Write();
   int index;
   factori==0? index=parti : index=partj;
   facmap.push_back(std::make_pair(index,factor));

   // draw the best fitting
      xframe = x.frame(Title("fit kaon"));
      //h1->Reset();
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
         //if(ngam>0) continue;
         double p1,p2;
         double theta1,theta2;
         double factors[2]={1,1};
         factori==0? factors[0]=factor : factors[0]=factori;
         factorj==0? factors[1]=factor : factors[1]=factorj;
         p1 = CalMom(kappx,kappy,kappz);
         p2 = CalMom(kampx,kampy,kampz);
         theta1 = acos(kappz/p1);
         theta2 = acos(kampz/p2);
         if(theta1>thetacut[parti]&&theta1<thetacut[parti+1]
          &&theta2>thetacut[partj]&&theta2<thetacut[partj+1]){
           mass = CalInvMass(mk,kappx,kappy,kappz,mk,kampx,kampy,kampz,-2,factors);
           //h1->Fill(mass);
           if (mass>philow && mass<phiup){
             dataraw->Fill();
           }
         }
         // if (Cut(ientry) < 0) continue;
      }
      dataraw->Write();

      char tmpchr[100];
      sprintf(tmpchr,"data_k");
      //data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      dataset = new RooDataSet("dataset","data",dataraw,x);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      //sum->plotOn(xframe,Components(expo),LineStyle(3),LineColor(3));
      sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
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
      sprintf(name,"fit_kk_best_part%d_part%d.eps",parti,partj);
      c1->Print(name);
      sprintf(name,"mass_kk_part%d_part%d",parti,partj);
      xframe->SetName(name);
      xframe->Write();
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;
   }
   }

   f->Close();
   ofpar.close();
   ofpardetail.close();
   

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
   fChain->SetBranchAddress("mphi", &mphi, &b_mphi);
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
