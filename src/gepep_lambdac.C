#define gepep_lambdac_cxx
#include "gepep_lambdac.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "function.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
  #include "RooChebychev.h"
//#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include <vector>
#include <iomanip>
#include <sstream>
#include "EventClass.h"

extern std::string outputdir;
using namespace RooFit;
using namespace std;

namespace EEtoLambdac{
  void FitSpe(std::vector<Lambdac> &evts, double beame, const char* namesfx=0);
  void FitSpectrum(TTree *&dataraw, double beame, const char* namesfx=0, double *dm=0, double *dmerr=0);
  void FitSpectrum2(TTree *&dataraw, double beame, const char* namesfx=0); // fit cme with gaussian function
  void FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx=0); // fit cme with crystal ball function
  //double GetEnergy(int runNo);
}

void gepep_lambdac::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_lambdac.C
//      Root > gepep_lambdac t
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
   
   fChain->GetEntry(1);
   double beamene = GetEnergy(run);
   std::cout<<"current beam energy is "<<beamene<<", run id "<<run<<std::endl;
   if (beamene < 0.1){
     std::cout<<"can not get a suitable beam energy!"<<std::endl;
	 return;
   }
   
   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=100;
   //double peakvalue=beamene;// mbeam
   //double beamlow=beamene-0.25;
   //double beamup=beamene+0.25;
   double mlambdac = 2.28646;
   double mpi=0.13957;
   double mp = 0.938272;
   double mka=0.493677;
   double totpcut = sqrt((beamene/2)*(beamene/2)-2.28646*2.28646)*2;
   std::cout<<"INFO: Momentum cut is "<< totpcut << std::endl;
    
  char fname[1000];
  sprintf(fname,"%s/plot_lambdac.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");

  char name[100];
  TCanvas *c1=new TCanvas("c1","",800,600);
  TTree *vars = new TTree("vars","vars");
  double mass, pp, pk, ppi,pall,cme;
  int nlam=0;
  int nlambar=0;
  vars->Branch("pp",&pp,"pp/D");
  vars->Branch("pk",&pk,"pk/D");
  vars->Branch("ppi",&ppi,"ppi/D");
  // near threshold, total momentum should be very small, 
  // this can be used to cut most background! pall < 0.2
  vars->Branch("pall",&pall,"pall/D");
  vars->Branch("mass",&mass,"mass/D");
  vars->Branch("CME",&cme,"cme/D");
  vars->Branch("nlam",&nlam,"nlam/I");
  vars->Branch("nlambar",&nlambar,"nlambar/I");

  Lambdac evt;
  Lambdac evt1, evt2;
  vector<Lambdac> evts;
  vector<HepLorentzVector> doubletags;
 
   Long64_t nbytes = 0, nb = 0;
   // tag number of all particles
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      nlam = 0;
      if (npp+nkam+npip >=3){
        if (npp<1) continue;
	if (nkam<1) continue;
	if (npip<1) continue;
        for (int i=0;i<npp;i++) {
          for (int j=0;j<nkam;j++) {
            for (int k=0;k<npip;k++) {
               HepLorentzVector pro, kam, pip;
               pro.setVectM(Hep3Vector(protonpx[i],protonpy[i],protonpz[i]),mp);
               kam.setVectM(Hep3Vector(kampx[j],kampy[j],kampz[j]),mka);
               pip.setVectM(Hep3Vector(pippx[k],pippy[k],pippz[k]),mpi);
 
               evt.set(pro,kam,pip);
               evt1.set(pro,kam,pip);
               mass = evt.m();
               pp = evt.proton.rho();
               pk = evt.kam.rho();
               ppi= evt.pip.rho();
               pall=(pro+kam+pip).rho();
               
	       //if (pall > 0.2) continue;
	       if (pall > totpcut) continue;
	       if (mass > 2.4) continue;
	       if (mass < 2.2) continue;
               nlam ++;
               evts.push_back(evt);
               vars->Fill();
            } // end loop k
	  } // end loop j
	} // end loop i
      }
      nlambar = 0;
      if (npm+nkap+npim >=3){
        if (npm<1) continue;
	if (nkap<1) continue;
	if (npim<1) continue;
        for (int i=0;i<npm;i++) {
          for (int j=0;j<nkap;j++) {
            for (int k=0;k<npim;k++) {
              HepLorentzVector pbar, kap, pim;
              pbar.setVectM(Hep3Vector(pbarpx[i],pbarpy[i],pbarpz[i]),mp);
              kap.setVectM(Hep3Vector(kappx[j],kappy[j],kappz[j]),mka);
              pim.setVectM(Hep3Vector(pimpx[k],pimpy[k],pimpz[k]),mpi);
   
              evt.set(pbar,kap,pim);
              evt2.set(pbar,kap,pim);
              mass = evt.m();
              pp = pbar.rho();
              pk = kap.rho();
              ppi= pim.rho();
              pall=(pbar+kap+pim).rho();
              
	      //if (pall > 0.2) continue;
	      if (pall > totpcut) continue;
	      if (mass > 2.4) continue;
	      if (mass < 2.2) continue;
              nlambar ++;
              evts.push_back(evt);
              vars->Fill();
	    } // end loop k
	  } // end loop j
	} // end loop i
      }

      if (nlam ==1 && nlambar==1){
        std::cout<<"evt id "<< jentry<<std::endl;
	doubletags.push_back(evt1.get4P()+evt2.get4P());
      }
      
   }
   
   /* version 1, with out tagging number of proton
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (nkam<1) continue;
      if (npip<1) continue;
      for (int i=0;i<nkam;i++)
      for (int j=0;j<npip;j++){
        HepLorentzVector pro, kam, pip;
	pro.setVectM(Hep3Vector(protonpx,protonpy,protonpz),mp);
	kam.setVectM(Hep3Vector(kampx[i],kampy[i],kampz[i]),mka);
	pip.setVectM(Hep3Vector(pippx[j],pippy[j],pippz[j]),mpi);

	evt.set(pro,kam,pip);
	mass = evt.m();
	pp = evt.proton.rho();
	pk = evt.kam.rho();
	ppi= evt.pip.rho();
	pall=(pro+kam+pip).rho();
	
	cme = (pro+kam+pip).e()*2;

	evts.push_back(evt);
	vars->Fill();
      }
   }
   */
   vars->Write();


   EEtoLambdac::FitSpe(evts,beamene,0);


/*
  TTree *dataraw = new TTree("dataraw","dataraw");
  //double mass;
  dataraw->Branch("x",&mass,"x/D");
  
  for (Long64_t ievt=0; ievt<evts.size();ievt++) {
    Lambdac evt = evts.at(ievt);
    mass = evt.m();
    dataraw->Fill();
  }

  double dm=0;
  EEtoLambdac::FitSpectrum(dataraw, mlambdac, "lambdac", &dm);
   

   ////////////////////////////////////
   TTree *dataraw1 = new TTree("dataraw1","dataraw1");
   //double cme;
   dataraw1->Branch("x",&cme,"x/D");
   std::cout<<"######### peak shift by "<<dm<<std::endl;
   HepLorentzVector lamc1,lamc2,lamctmp;
   std::vector<HepLorentzVector> singletags;
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     lamc1.setVectM( lamctmp.vect(),lamctmp.m()-dm);
     lamc2.setVectM(-lamctmp.vect(),lamctmp.m()-dm);
     //singletags.push_back(lamc);
     
     cme = (lamc1+lamc2).m();
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beamene,"doublee");
   
   // check p correction
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     lamc1.setVectM( lamctmp.vect()*1.01,lamctmp.m()-dm);
     lamc2.setVectM(-lamctmp.vect()*1.01,lamctmp.m()-dm);
     //singletags.push_back(lamc);
     
     //cme = lamc.e()*2;
     cme = (lamc1+lamc2).m();
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beamene,"fp1.01");
   
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     lamc1.setVectM( lamctmp.vect()*0.99,lamctmp.m()-dm);
     lamc2.setVectM(-lamctmp.vect()*0.99,lamctmp.m()-dm);
     //singletags.push_back(lamc);
     
     //cme = lamc.e()*2;
     cme = (lamc1+lamc2).m();
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beamene,"fp0.99");
   
   // check px shift
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     Hep3Vector totp;
     totp.setX(4.575*sin(0.011)-lamctmp.vect().getX());
     totp.setY(-lamctmp.vect().getY());
     totp.setZ(-lamctmp.vect().getZ());
     lamc1.setVectM(totp,lamctmp.m()-dm);
     lamc2.setVectM(lamctmp.vect(),lamctmp.m()-dm);
     //singletags.push_back(lamc);
     
     cme = (lamc1+lamc2).m();
     //cme = lamc.e()*2;
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beamene,"pxcon");

   // tag 2 particles
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<doubletags.size();ievt++) {
     cme = doubletags.at(ievt).m();
     //cme = lamc.e()*2;
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beamene,"2tag");
*/
}

void EEtoLambdac::FitSpe(std::vector<Lambdac> &evts, double beame, const char* namesfx)
{
  TTree *dataraw = new TTree("dataraw","dataraw");
  double mass;
  dataraw->Branch("x",&mass,"x/D");
  
  // get correction factor for momentum
  char tmpchr[100];
  const int np = 20;
  double factors[np];
  double factorserr[np];
  for (int i=0; i<np; i++){
    factors[i] = (1.01-0.99)/np*i+0.99;
    factorserr[i] = 0;
  }
  double deltapeaks[np];
  double peakerrors[np];

  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (Long64_t ievt=0; ievt<evts.size();ievt++) {
      Lambdac evt = evts.at(ievt);
      evt.setCorrectionFactors(factors[i]);
      mass = evt.m();
      dataraw->Fill();
    }
    sprintf(tmpchr,"%.3f",factors[i]);
    EEtoLambdac::FitSpectrum(dataraw, beame, tmpchr, &deltapeaks[i], &peakerrors[i]);
  }

   TCanvas *c1 = new TCanvas("c1_1","c1");
   TGraphErrors *graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   TF1 *facfit = new TF1("facfit","[0]*(x-[1])");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","factor");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   double factor = facfit->GetParameter(1);
   double fiterr = facfit->GetParError(1);
   double peslop = peakerrors[np/2]/facfit->GetParameter(0);
   double factorerr = sqrt(TMath::Power(fiterr,2)+TMath::Power(peslop,2));
   graph1->SetName("PsVF");
   graph1->Draw();
   c1->Write();
 

   ////////////////////////////////////
   TTree *dataraw1 = new TTree("dataraw1","dataraw1");
   double cme;
   dataraw1->Branch("x",&cme,"x/D");
   std::cout<<"######### average correction factor is "<< factor <<std::endl;
   std::cout<<"######### average correction factor error is "<< factorerr <<std::endl;
   
   HepLorentzVector lamc1,lamc2,lamctmp;
   std::vector<HepLorentzVector> singletags;
   
   // correct p
   dataraw->Reset();
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     evt.setCorrectionFactors(factor);
     lamctmp = evt.get4P();
     lamc1.setVectM( lamctmp.vect(),lamctmp.m());
     lamc2.setVectM(-lamctmp.vect(),lamctmp.m());
     //singletags.push_back(lamc);
     
     std::cout<<"id:"<<ievt<< "\tE:"<<lamc1.e()<<"\tp:"<<lamc1.rho()<<"\t"<<lamc1.m()<<std::endl;
     mass = evt.m();
     dataraw->Fill();
     cme = (lamc1+lamc2).m();
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum(dataraw, beame, "corp_inv_m");
   EEtoLambdac::FitSpectrum2(dataraw1,beame,"corp_doubleE");
   EEtoLambdac::FitSpectrum3(dataraw1,beame,"corp_doubleE");
   
   // correct m
   //double dm = 1.413e-3;
   double dm = - deltapeaks[np/2];
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     lamc1.setVectM( lamctmp.vect(),lamctmp.m()+dm);
     lamc2.setVectM(-lamctmp.vect(),lamctmp.m()+dm);
     //singletags.push_back(lamc);
     
     std::cout<<"id:"<<ievt<< "\tE:"<<lamc1.e()<<"\tp:"<<lamc1.rho()<<"\t"<<lamc1.m()<<std::endl;
     cme = (lamc1+lamc2).m();
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beame,"shiftm_doubleE");
   EEtoLambdac::FitSpectrum3(dataraw1,beame,"shiftm_doubleE");
  
   // just fit E
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     lamc1.setVectM( lamctmp.vect(),lamctmp.m()+dm);
     lamc2.setVectM(-lamctmp.vect(),lamctmp.m()+dm);
     //singletags.push_back(lamc);
     
     //std::cout<<"id:"<<ievt<< "\tE:"<<lamc1.e()<<"\tp:"<<lamc1.rho()<<std::endl;
     mass = lamc1.e();
     dataraw->Fill();
   }
   EEtoLambdac::FitSpectrum(dataraw,beame,"shiftm_E");
     
 //// check p correction
 //dataraw1->Reset();
 //for (Long64_t ievt=0; ievt<evts.size();ievt++) {
 //  Lambdac evt = evts.at(ievt);
 //  lamctmp = evt.get4P();
 //  lamc1.setVectM( lamctmp.vect()*1.01,lamctmp.m()-dm);
 //  lamc2.setVectM(-lamctmp.vect()*1.01,lamctmp.m()-dm);
 //  //singletags.push_back(lamc);
 //  
 //  //cme = lamc.e()*2;
 //  cme = (lamc1+lamc2).m();
 //  dataraw1->Fill();
 //}
 //EEtoLambdac::FitSpectrum2(dataraw1,beame,"fp1.01");
   
   // check px shift
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     evt.setCorrectionFactors(factor);
     lamctmp = evt.get4P();
     Hep3Vector totp;
     totp.setX(beame*sin(0.011)-lamctmp.vect().getX());
     totp.setY(-lamctmp.vect().getY());
     totp.setZ(-lamctmp.vect().getZ());
     lamc1.setVectM(totp,lamctmp.m());
     lamc2.setVectM(lamctmp.vect(),lamctmp.m());
     //singletags.push_back(lamc);
     
     cme = (lamc1+lamc2).m();
     //cme = lamc.e()*2;
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beame,"pxboost");
   EEtoLambdac::FitSpectrum3(dataraw1,beame,"pxboost");
   
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     evt.setCorrectionFactors(factor-factorerr);
     lamctmp = evt.get4P();
     Hep3Vector totp;
     totp.setX(beame*sin(0.011)-lamctmp.vect().getX());
     totp.setY(-lamctmp.vect().getY());
     totp.setZ(-lamctmp.vect().getZ());
     lamc1.setVectM(totp,lamctmp.m());
     lamc2.setVectM(lamctmp.vect(),lamctmp.m());
     //singletags.push_back(lamc);
     
     cme = (lamc1+lamc2).m();
     //cme = lamc.e()*2;
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beame,"pxboost_flow");
   EEtoLambdac::FitSpectrum3(dataraw1,beame,"pxboost_flow");
    
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     evt.setCorrectionFactors(factor+factorerr);
     lamctmp = evt.get4P();
     Hep3Vector totp;
     totp.setX(beame*sin(0.011)-lamctmp.vect().getX());
     totp.setY(-lamctmp.vect().getY());
     totp.setZ(-lamctmp.vect().getZ());
     lamc1.setVectM(totp,lamctmp.m());
     lamc2.setVectM(lamctmp.vect(),lamctmp.m());
     //singletags.push_back(lamc);
     
     cme = (lamc1+lamc2).m();
     //cme = lamc.e()*2;
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beame,"pxboost_fup");
   EEtoLambdac::FitSpectrum3(dataraw1,beame,"pxboost_fup");
 
   // mass shift
   dataraw1->Reset();
   for (Long64_t ievt=0; ievt<evts.size();ievt++) {
     Lambdac evt = evts.at(ievt);
     lamctmp = evt.get4P();
     Hep3Vector totp;
     totp.setX(beame*sin(0.011)-lamctmp.vect().getX());
     totp.setY(-lamctmp.vect().getY());
     totp.setZ(-lamctmp.vect().getZ());
     lamc1.setVectM(totp,lamctmp.m()+dm);
     lamc2.setVectM(lamctmp.vect(),lamctmp.m()+dm);
     //singletags.push_back(lamc);
     
     cme = (lamc1+lamc2).m();
     //cme = lamc.e()*2;
     dataraw1->Fill();
   }
   EEtoLambdac::FitSpectrum2(dataraw1,beame,"pxboost_shiftm");
   EEtoLambdac::FitSpectrum3(dataraw1,beame,"pxboost_shiftm");


}

void EEtoLambdac::FitSpectrum(TTree *&dataraw, double beame, const char* namesfx, double *dm, double *dmerr)
{
   int nBins=50;
   int Npar;
   double peakvalue = 2.28646;//beame;
   double beamlow=peakvalue-0.04;
   double beamup=peakvalue+0.04;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue -0.002,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.004,0.002,0.006);
   //if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.005);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   
   RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
   RooRealVar co2("co2","coefficient #2",-0.2,-100.,100.);
   RooRealVar co3("co3","coefficient #3",0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);
   
   double factor = atof(namesfx);
   if (factor !=0) mean.setVal(peakvalue + 1*(factor - 1.001) );

   char tmpchr[100];
   sprintf(tmpchr,"data_lambdac_%s",namesfx);
   xframe = x.frame(Title("fit lambdac"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe,Binning(nBins));
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"mass_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();

   if (dm != 0)  *dm = mean.getVal() - peakvalue;
   if (dmerr!=0) *dmerr = mean.getError();
  std::cout<<"############ffff"<<std::endl;

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
} 

void EEtoLambdac::FitSpectrum2(TTree *&dataraw, double beame, const char* namesfx)
{ 
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.1;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue /*-0.002*/,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.003,0.001,0.018);
   //if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.005);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   
   RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
   RooRealVar co2("co2","coefficient #2",-0.2,-100.,100.);
   RooRealVar co3("co3","coefficient #3",0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("c1_2","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_lambdac_bar_%s",namesfx);
   xframe = x.frame(Title("fit lambdac and bar"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe,Binning(nBins));
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"cme_gaus_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();

   delete xframe;
   delete dataset;
   delete sum;
   return;
}

void EEtoLambdac::FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx)
{ 
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.1;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue /* -0.002*/,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.003,0.001,0.018);
   //if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.005);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   
   RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
   RooRealVar co2("co2","coefficient #2",-0.2,-100.,100.);
   RooRealVar co3("co3","coefficient #3",0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",140,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("c1_2","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_lambdac_bar_%s",namesfx);
   xframe = x.frame(Title("fit lambdac and bar"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe,Binning(nBins));
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"cme_crystal_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();

   delete xframe;
   delete dataset;
   delete sum;
   return;
}

