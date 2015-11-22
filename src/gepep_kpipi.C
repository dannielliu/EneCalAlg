// this is from a copy of fkkpipi algorithm
#define gepep_kpipi_cxx
#include "gepep_kpipi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "EventClass.h"
#include <iomanip>
extern std::string outputdir;
using namespace RooFit;
using namespace std;
namespace KPIPI{
  void FitSpe(std::vector<D2KPiPi> &evts, double beame,  char* namesfx);
  void FitSpectrum(TTree *&dataraw,double beame, char* namesfx, double &peak, double &peakerror);
  //double GetEnergy(int runNo);
}

void gepep_kpipi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_kpipi.C
//      Root > gepep_kpipi t
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
//
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   fChain->GetEntry(1);
   double beamene = GetEnergy(run);
   //beamene = 4.26;
   std::cout<<"current beam energy is "<<beamene<<", run id "<<run<<std::endl;
   if (beamene < 0.1){
     std::cout<<"can not get a suitable beam energy!"<<std::endl;
	 return;
   }

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=100;
// double peakvalue=beamene;// mbeam
// double beamlow=beamene-0.1;
// double beamup=beamene+0.1;
   double peakvalue=1.86962;// mbeam
   double beamlow=peakvalue-0.03;
   double beamup=peakvalue+0.03;
   double mpi=0.13957;
   double mk = 0.493677;
    
  char fname[1000];
  sprintf(fname,"%s/plot_k2pi_AR.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");

  char name[100];
  TCanvas *c1=new TCanvas("c1","",800,600);
  TTree *vars = new TTree("vars","vars");
  double mass;
  //double phi,phi1,phi2;
  //double costheta,costheta1,costheta2;
  double ppi1,ppi2,pk1,pk2;
  vars->Branch("ppi1",&ppi1,"ppi1/D");
  vars->Branch("ppi2",&ppi2,"ppi2/D");
  vars->Branch("pk1",&pk1,"pk1/D");
  //vars->Branch("pk2",&pk2,"pk2/D");
  vars->Branch("mass",&mass,"mass/D");
  double chisq;
  vars->Branch("chi2",&chisq,"chi2/D");
   
   D2KPiPi evt;
   std::vector<D2KPiPi> evts;
   std::vector<D2KPiPi> evts2;
   
   
   TH1D *hppi = new TH1D("hppi","pt_{#pi}",200,0,2);
   TH1D *hpka = new TH1D("hpka","pt_{K}",200,0,2);
   TH1D *hmD  = new TH1D("hmD" ,"M_{K#pi#pi}",200,1,2);

   // select useful events
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;	  
	  
	  HepLorentzVector kam,pip1,pip2;
	  if (nkam==1 && npip==2) {
		  kam.setVectM(Hep3Vector(kampx[0],kampy[0],kampz[0]),mk);
		  pip1.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
		  pip2.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mpi);
	  }
	  else if (nkap==1 && npim==2) {
		  kam.setVectM(Hep3Vector(kappx[0],kappy[0],kappz[0]),mk);
		  pip1.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
		  pip2.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mpi);
	  }
	  else continue;
	  evt.set(kam,pip1,pip2);
	  double costheta1 = cos(kam.vect().theta());
	  double costheta2 = cos(pip1.vect().theta());
	  double costheta3 = cos(pip2.vect().theta());
////  if (fabs(costheta1)>0.5) continue;
////  if (fabs(costheta2)>0.5) continue;
////  if (fabs(costheta3)>0.5) continue;
          mass = evt.m();
	  ppi1 = pip1.rho();
	  ppi2 = pip2.rho();
	  pk1 = kam.rho();
	  double ptpi1 = pip1.perp();
	  double ptpi2 = pip2.perp();
	  double ptk1 = kam.perp();
	 //pk2 = kap.vect().mag();
	  chisq = kkm4;

	  
	  if (mass>beamlow-0.01 && mass<beamup+0.01){
            //evts2.push_back(evt);
	    //if (pk1<0.6 || pk1>0.8) continue;
            evts.push_back(evt);
	 
	    //vars->Fill();
	  }

	  //if (ptk1<0.4 || ptk1>0.8) continue;
	  //if (pk1<0.6 || pk1>0.8) continue;
	  hmD->Fill(mass);
	  if (mass>peakvalue-0.01 && mass<peakvalue+0.01){
	    hppi->Fill(ptpi1);
	    hppi->Fill(ptpi2);
	    hpka->Fill(ptk1);
	  }
   }//select end
   //vars->Write();
   //sprintf(name,"%f",peakvalue);
 
// TFile *ftmp = new TFile("P_cmp.root","update");
// ftmp->WriteTObject(hpka,"hptka_DKpipi");
// ftmp->WriteTObject(hppi,"hptpi_DKpipi");
// ftmp->WriteTObject(hmD ,"hmD_DKpipi");
// ftmp->Close();
// delete ftmp;
// f->cd();
// return ;

   double pka_tot = hpka->GetMean();
   double pkasig_tot = hpka->GetRMS();
   double ppi_tot = hppi->GetMean();
   double ppisig_tot = hppi->GetRMS();
   ofstream ofpar("parkpipi.txt",std::ios::app);
   ofpar<<"p_K = " << pka_tot << " +/- "<< pkasig_tot <<"\t";
   ofpar<<"p_pi = " << ppi_tot << " +/- "<< ppisig_tot <<endl;
   ofpar.close();
   
   sprintf(name,"cut");
   KPIPI::FitSpe(evts,beamene,name);
// sprintf(name,"nocut");
// KPIPI::FitSpe(evts2,beamene,name);
   return;

}

void KPIPI::FitSpe(std::vector<D2KPiPi> &evts, double beame, char *namesfx)
{
  double peak = 1.86962;
  double beamlow=peak-0.1;
  double beamup=peak+0.1;
  // for factor fit
 
  TTree *datarawo = new TTree("datarawo","dataraw");
  TTree *dataraw = new TTree("dataraw","dataraw");
//TTree *datarawl = new TTree("datarawl","dataraw");
//TTree *datarawu = new TTree("datarawu","dataraw");
  double mass;
  datarawo->Branch("x",&mass,"x/D");
  dataraw->Branch("x",&mass,"x/D");
//datarawl->Branch("x",&mass,"x/D");
//datarawu->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
 
  char tmpchr[100];

  //~~~~~~~~~~part start~~~~~~~~

//for (Long64_t jentry=0; jentry<evts.size();jentry++) {
//   // total invariant mass
//   D2KPiPi evt = evts.at(jentry);
//   // without correction
//   mass = evt.m();
//   if (mass>beamlow-0.001 && mass<beamup+0.001) datarawo->Fill();
// //evt.setCorrectionFactors(1.000714,1.000902);
// //mass = evt.m();
// //if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
//}
////dataraw->Write();
//// no correction
//sprintf(tmpchr,"raw_%s",namesfx);
//KPIPI::FitSpectrum(datarawo,beame,tmpchr);
  
  int np = 20;
  double factors[np];
  double fpi=1.00;
  double fk =1.00;
  double factorserr[np];
  for (int i=0; i<np; i++){
    factors[i] = (1.005-0.995)/np*i+0.995;
    factorserr[i] = 0;
  }
  double peaks[np];
  double deltapeaks[np];
  double peakerrors[np];


  // dm Vs fk
  fpi = 1.0;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D2KPiPi evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fk_%f",factors[i]);
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
    deltapeaks[i] = peaks[i] - peak;
    if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   TCanvas *c1 = new TCanvas("c1_fk","c1_fpi1.0");
   TGraphErrors *graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   TF1 *facfit = new TF1("facfit","[0]*(x-1)+[1]");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","interupt");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   double a = facfit->GetParameter(0);
   double ae= facfit->GetParError(0);
   double b = facfit->GetParameter(1);
   double be= facfit->GetParError(1);
   double factor = 1 - b/a;
   double factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);

   graph1->SetName("PsVF");
   graph1->Draw();

  TPaveText *pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fk");

  ofstream ofpar("parkpipi.txt",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b
   << ", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np+1)/2] 
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D2KPiPi evt = evts.at(evtid);
      evt.setCorrectionFactors(factor,fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fk_%f",factor);
    double peakt,errt;
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

 // dm Vs fk
  fpi = 1.001;
  /*
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D2KPiPi evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
	deltapeaks[i] = peaks[i] - peak;
	if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   c1 = new TCanvas("c1_fk","c1_fpi1.001");
   delete graph1;
   //TGraphErrors *;
   graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   //TF1 *facfit = new TF1("facfit","[0]*(x-1)+[1]");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","interupt");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   a = facfit->GetParameter(0);
   ae= facfit->GetParError(0);
   b = facfit->GetParameter(1);
   be= facfit->GetParError(1);
   factor = 1 - b/a;
   factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);

   graph1->SetName("PsVF");
   graph1->Draw();

  pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fk");

  //ofstream ofpar("parkpipi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b
   << ", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np+1)/2] 
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D2KPiPi evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


  // dm Vs fk
  fpi = 0.999;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D2KPiPi evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
	deltapeaks[i] = peaks[i] - peak;
	if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   c1 = new TCanvas("c1_fk","c1_fpi0.999");
   //TGraphErrors *
   delete graph1;
   graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   //TF1 *facfit = new TF1("facfit","[0]*(x-1)+[1]");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","interupt");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   a = facfit->GetParameter(0);
   ae= facfit->GetParError(0);
   b = facfit->GetParameter(1);
   be= facfit->GetParError(1);
   factor = 1 - b/a;
   factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);

   graph1->SetName("PsVF");
   graph1->Draw();

  pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fk");

  //ofstream ofpar("parkpipi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b
   << ", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np+1)/2] 
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D2KPiPi evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
*/


  // dm Vs fpi
  //fk = 0.999830;
  fk = 1.0;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
       D2KPiPi evt = evts.at(evtid);
       evt.setCorrectionFactors(fk,factors[i]);
       mass = evt.m();
       if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
    deltapeaks[i] = peaks[i] - peak;
    if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   //TCanvas *c1 = new TCanvas("c1_1","c1");
   delete c1;
   c1 = new TCanvas("c1_fpi","c1_fk1.00");
   c1->Clear();
   //TGraphErrors *
   delete graph1;
   graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   //TF1 *
   delete facfit;
   facfit = new TF1("facfit","[0]*(x-1)+[1]");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","interupt");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   a = facfit->GetParameter(0);
   ae= facfit->GetParError(0);
   b = facfit->GetParameter(1);
   be= facfit->GetParError(1);
   factor = 1 - b/a;
   factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);
   graph1->SetName("PsVF");
   graph1->Draw();

  //TPaveText *
  delete pt;
  pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fpi");

  //ofstream ofpar("parkpipi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b 
   <<", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np-1)/2]
   <<resetiosflags(ios::fixed)<<std::endl;
    
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D2KPiPi evt = evts.at(evtid);
      evt.setCorrectionFactors(1,factor);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    //double factor =1;
    sprintf(tmpchr,"%s_fpi_%f",namesfx,factor);
    //double peakt,errt;
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

//  dataraw->Reset();
//  TF1 ff("ff","1.00065+0.000630*x",0,2);
//  for (int evtid=0; evtid<evts.size();evtid++){
//    D2KPiPi evt = evts.at(evtid);
//    double f1 = ff.Eval(evt.pip.rho());
//    double f2 = ff.Eval(evt.pim.rho());
//    evt.setCorrectionFactors(fk,f1,f2);
//    mass = evt.m();
//    if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
//  }
//  //double 
//  factor =0;
//  sprintf(tmpchr,"factor_%f",factor);
//  //double peakt,errt;
//  KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

/*
  // dm Vs fpi
  fk = 1.001;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D2KPiPi evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factors[i]);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
	deltapeaks[i] = peaks[i] - peak;
	if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   //TCanvas *c1 = new TCanvas("c1_1","c1");
   c1 = new TCanvas("c1_fpi","c1_fk1.001");
   //TGraphErrors *
   c1->Clear();
   graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   //TF1 *
   facfit = new TF1("facfit","[0]*(x-1)+[1]");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","interupt");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   a = facfit->GetParameter(0);
   ae= facfit->GetParError(0);
   b = facfit->GetParameter(1);
   be= facfit->GetParError(1);
   factor = 1 - b/a;
   factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);
   graph1->SetName("PsVF");
   graph1->Draw();

  //TPaveText *
  pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fpi");

  //ofstream ofpar("parkpipi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b 
   <<", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np-1)/2]
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D2KPiPi evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


  // dm Vs fpi
  fk = 0.999;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D2KPiPi evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factors[i]);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
	deltapeaks[i] = peaks[i] - peak;
	if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   //TCanvas *c1 = new TCanvas("c1_1","c1");
   c1 = new TCanvas("c1_fpi","c1_fk1.001");
   //TGraphErrors *
   c1->Clear();
   graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   //TF1 *
   facfit = new TF1("facfit","[0]*(x-1)+[1]");
   facfit->SetParameters(0.3,1);
   facfit->SetParNames("slope","interupt");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   a = facfit->GetParameter(0);
   ae= facfit->GetParError(0);
   b = facfit->GetParameter(1);
   be= facfit->GetParError(1);
   factor = 1 - b/a;
   factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);
   graph1->SetName("PsVF");
   graph1->Draw();

  //TPaveText *
  pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fpi");

  //ofstream ofpar("parkpipi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b 
   <<", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np-1)/2]
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D2KPiPi evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
*/


  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void KPIPI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx, double &peak, double &peakerror)
{
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = 1.86962;
   double beamlow = peakvalue-0.03;
   double beamup  = peakvalue+0.03;
   // try to use roofit
   RooRealVar x("x","M(K#pi#pi)",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.003,0.001,0.005);
   RooRealVar sigma2("sigma2","width of gaussian",0.006,0.005,0.01);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   //RooRealVar co1("co1","coefficient #1",0,-100.,100.);
   //RooRealVar co4("co4","coefficient #4",0);
   //RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",200,0,10000000);//event number
   RooRealVar sigfra("sigfra"," ",0.5,0.3,1.0);//event number
   RooRealVar signal2("signal2"," ",200,0,10000000);//event number
   RooRealVar background("background"," ",20,0,1000000);
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooAddPdf sig("sig","signal",RooArgList(gaus,gaus2),sigfra);
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
   TH1D *hmass = new TH1D("hmass","M(K #pi #pi)",100, beamlow, beamup);
   RooDataHist *datahist;
   if (largesample) {
     dataraw->Draw("x>>hmass");
   }
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_kpipi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("M(K#pi#pi)"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   if (largesample)  datahist = new RooDataHist(tmpchr,"data",RooArgSet(x),hmass);
   sum = new RooAddPdf("sum","sum",RooArgList(sig,ground),RooArgList(signal,background));
   Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   if (largesample) {sum->fitTo(*datahist);datahist->plotOn(xframe);}
   else {
     sum->fitTo(*dataset,Range(beamlow,beamup));
     dataset->plotOn(xframe);
   }
   //sum->fitTo(*dataset,Range(beamlow,beamup));
   //dataset->plotOn(xframe);
   //sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   //sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   //dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //if (dataraw->GetEntries()>2000) sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(4));
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
//if (largesample){
//  sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
//  pt->AddText(tmpchr);
//}
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
//if (largesample){
//  sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
//  pt->AddText(tmpchr);
//}
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("kpipi",std::ios::app);
   outf<<beame<<'\t'<< /* namesfx<<"\t"<< */ mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   peak = mean.getVal();
   peakerror = mean.getError();
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   if (largesample) delete hmass;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}






#ifdef gepep_kpipi_cxx
gepep_kpipi::gepep_kpipi(TTree *tree) : fChain(0) 
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

gepep_kpipi::~gepep_kpipi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_kpipi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_kpipi::LoadTree(Long64_t entry)
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

void gepep_kpipi::Init(TTree *tree)
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

Bool_t gepep_kpipi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_kpipi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_kpipi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_kpipipi_cxx
