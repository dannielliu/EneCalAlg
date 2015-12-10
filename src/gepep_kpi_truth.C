#define gepep_kpi_cxx
#include "gepep_kpi.h"
#include "mctruth.h"
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
#include "EventClass.h"
#include <iomanip>
extern std::string outputdir;
extern mctruth *truthalg;
using namespace RooFit;
using namespace std;
namespace KPI{
  void FitSpe(std::vector<D0KPI> &evts, double beame,  char* namesfx);
  void FitSpectrum(TTree *&dataraw,double beame, char* namesfx, double &peak, double &peakerror);
  //double GetEnergy(int runNo);
}


void gepep_kpi::Loop()
{
//   In a ROOT session, you can do:
//	  Root > .L gepep_kpi.C
//	  Root > gepep_kpi t
//	  Root > t.GetEntry(12); // Fill t data members with entry number 12
//	  Root > t.Show();	   // Show values of entry 12
//	  Root > t.Show(16);	 // Read and show values of entry 16
//	  Root > t.Loop();	   // Loop on all entries
//

//	 This is the loop skeleton where:
//	jentry is the global entry number in the chain
//	ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//	jentry for TChain::GetEntry
//	ientry for TTree::GetEntry and TBranch::GetEntry
//
//	   To read only selected branches, Insert statements like:
// METHOD1:
//	fChain->SetBranchStatus("*",0);  // disable all branches
//	fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//	fChain->GetEntry(jentry);	   //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   fChain->GetEntry(1);
   double beamene = GetEnergy(run);

   std::cout<<"Toral entry is "<<nentries<< "@ "<< beamene << " GeV." <<std::endl;
// int nBins=40;
// double factorstart=0.995;
   double D0low=1.82;
   double D0up=1.90;
   double mk=0.493677;
   double mpi=0.13957018;
   double peakvalue=1.86484;// mD0
   int pointNo=20;

   ofstream ofpar;
   ofpar.open("parkpi",std::ios::app);
   //ofpar<<"k- pi+ algrithm: will give factors for pion"<<std::endl;
   ofstream ofpardetail;
   ofpardetail.open("detail.txt",std::ios::app);
   ofstream purepar;
   purepar.open("par");
   char fname[100];
   sprintf(fname,"%s/plot_kpi_AR.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");
  
   double mass;
 
   TTree *vars = new TTree("vars","vars");
   double phi1,phi2;
   double costheta1,costheta2;
   double p1,p2;
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"pk/D");
   vars->Branch("p2",&p2,"ppi/D");
   vars->Branch("mass",&mass,"mass/D");


   Long64_t nbytes = 0, nb = 0;

   TH1D *hppires = new TH1D("hppires","p_{#pi} (rec-truth)",200,-0.1,0.1);
   TH1D *hpkares = new TH1D("hpkares","p_{#ka} (rec-truth)",200,-0.1,0.1);
   TH1D *hppi = new TH1D("hppi","p_{#pi}",200,0,2);
   TH1D *hpka = new TH1D("hpka","p_{K}",200,0,2);
   TH1D *hmD0 = new TH1D("hmD0","M_{K #pi}",200,1,2);
   
   TH2D *h2p = new TH2D("h2p","#pi K momentum",200,0,2,200,0,2);
   TH2D *thedis = new TH2D("thedis","#theta",200,0,2,100,0,TMath::Pi());
   TH2D *thedis2= new TH2D("thedis2","#theta",100,0,TMath::Pi(),100,0,TMath::Pi());
   const  int Npart=17;
 //double thetacut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
 //double start=0.0;
 //double stop =TMath::Pi();
 //for(int i=0;i<Npart+1;i++){
 //  thetacut[i] = (stop-start)/Npart*i+start;
 //}
 //double coscut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
 //double start=-1;
 //double stop =1;
 //for(int i=0;i<Npart+1;i++){
 //      coscut[i] = (stop-start)/Npart*i+start;
 //}
   double pcut[Npart+1]={0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                         1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};

   double m0=peakvalue;
   double sigma_m = 0.0068;//0.0024 for phi,
   double width = 20.*sigma_m;
   //double mparticle=mk;
   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;
   std::vector<D0KPI> evts_set[Npart];
   std::vector<D0KPI> evts_set2;
   std::vector<D0KPI> evts_setmc;

   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   TH1D *hmtheta[Npart];
   TH1D *hpkaires[Npart];
   TH1D *hppiires[Npart];
   for (int partj=0;partj<Npart;partj++){
	 sprintf(name,"mass_part%d",partj);
	 hmtheta[partj] = new TH1D(name,name,100,D0low,D0up);
	 sprintf(name,"pkares_part%d",partj);
	 hpkaires[partj] = new TH1D(name,name,200,-0.1,0.1);
	 sprintf(name,"ppires_part%d",partj);
	 hppiires[partj] = new TH1D(name,name,200,-0.1,0.1);
   }

   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	 Long64_t ientry = LoadTree(jentry);
	 if (ientry < 0) break;
	 nb = fChain->GetEntry(jentry);   nbytes += nb;
	 //if (run >34028) break;

	 truthalg->fChain->GetEntry(jentry);

	 //int parti,partj;
	 //double mass;
	 double theta1,theta2;
	 int besidx=0;
	 double tmpdeltaold=100;
	 double tmpmass;
	 D0KPI evt, evtda,evtmc;
         std::vector<D0KPI> tmpv,tmpvmc;
	 HepLorentzVector kam, pip;
	 //if (npip==1){
	 //for(int i=0; i<1;i++){
	 //int i=0;
	 if (npim>=1 && nkap==1){
	   kam.setVectM(Hep3Vector(kappx[0],kappy[0],kappz[0]),mk);
	   pip.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
	   evtda.set(kam,pip);
	   for (int i=1;i<npim;i++){
	     pip.setVectM(Hep3Vector(pimpx[i],pimpy[i],pimpz[i]),mpi);
	     D0KPI tmpevt;
	     tmpevt.set(kam,pip);
	     if (fabs(tmpevt.m()-m0)<fabs(evtda.m()-m0)) evtda = tmpevt;
	   }
	   tmpv.push_back(evtda);
	   kam.setVectM(Hep3Vector(truthalg->truthpkap[0],
	                           truthalg->truthpkap[1],
				   truthalg->truthpkap[2]),
		        mk);
	   pip.setVectM(Hep3Vector(truthalg->mcpimpx[0],
	                           truthalg->mcpimpy[0],
	                           truthalg->mcpimpz[0]),
                        mpi);
	   evtmc.set(kam,pip);
	   tmpvmc.push_back(evtmc);
	   //hppires->Fill(evtda.pip.rho()-evtmc.pip.rho());
	 }
	 if (npip>=1 && nkam==1){
	   kam.setVectM(Hep3Vector(kampx[0],kampy[0],kampz[0]),mk);
	   pip.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
	   evtda.set(kam,pip);
	   for (int i=1;i<npip;i++){
	     pip.setVectM(Hep3Vector(pippx[i],pippy[i],pippz[i]),mpi);
	     D0KPI tmpevt;
	     tmpevt.set(kam,pip);
	     if (fabs(tmpevt.m()-m0)<fabs(evtda.m()-m0)) evtda = tmpevt;
	   }
	   tmpv.push_back(evtda);
	   kam.setVectM(Hep3Vector(truthalg->truthpkam[0],
	                           truthalg->truthpkam[1],
				   truthalg->truthpkam[2]),
		        mk);
	   pip.setVectM(Hep3Vector(truthalg->mcpippx[0],
	                           truthalg->mcpippy[0],
	                           truthalg->mcpippz[0]),
                        mpi);
	   evtmc.set(kam,pip);
	   tmpvmc.push_back(evtmc);
	   //hppires->Fill(evtda.pip.rho()-evtmc.pip.rho());
	 }
	   
	 for (int ie=0;ie<tmpv.size();ie++){
	   evt = tmpv.at(ie);
	   evtmc = tmpvmc.at(ie);
	 
/////////costheta1 = cos(evt.kam.theta());
/////////costheta2 = cos(evt.pip.theta());
/////////phi1 = evt.kam.phi();
/////////phi2 = evt.pip.phi();
	 //if (fabs(costheta1)>0.5) continue;
	 //if (fabs(costheta2)>0.5) continue;
	 
	 // total invariant mass, D0 -> k- pi+
	 p1 = evt.kam.rho();
	 p2 = evt.pip.rho();
	 double pt1 = evt.kam.perp();
	 double pt2 = evt.pip.perp();
//////// theta1= evt.kam.theta();
//////// theta2= evt.pip.theta();
	 mass = evt.m();
	 //if (parti>=Npart || partj>=Npart || parti<0 || partj<0) continue;
	 for (int parti=0; parti<Npart; parti++){
	    if (evtmc.pip.rho()>pcut[parti]&&evtmc.pip.rho()<pcut[parti+1]) {
	      hmtheta[parti]->Fill(mass);
	    //if (costheta1>coscut[parti]&&costheta1<coscut[parti+1])
	      //&&theta2>thetacut[partj]&&theta2<thetacut[partj+1])
	      if (mass>m0-width/2.-0.02 && mass<m0+width/2.+0.02)
	      {
	            evts_set2.push_back(evt);
	            evts_setmc.push_back(evtmc);
	    ////////h2p->Fill(p1,p2);
	    ////////thedis->Fill(p1,theta1);
	    ////////thedis->Fill(p2,theta2);
	    ////////thedis2->Fill(theta1,theta2);
	            evts_set[parti].push_back(evt);
	            //vars->Fill();
	      }
	      if (mass>m0-sigma_m && mass<m0+sigma_m){
	        //    hpkares->Fill(evt.kam.rho()-evtmc.kam.rho());
	            hppires->Fill(evt.pip.rho()-evtmc.pip.rho());
	          //  hpkaires[parti]->Fill(evt.kam.rho()-evtmc.kam.rho());
	            hppiires[parti]->Fill((evt.pip.rho()-evtmc.pip.rho())/evtmc.pip.rho());
	      }
	    }
	    if (evtmc.kam.rho()>pcut[parti]&&evtmc.kam.rho()<pcut[parti+1]) {
	      if (mass>m0-sigma_m && mass<m0+sigma_m){
	            hpkares->Fill(evt.kam.rho()-evtmc.kam.rho());
	            //hppires->Fill(evt.pip.rho()-evtmc.pip.rho());
	            hpkaires[parti]->Fill((evt.kam.rho()-evtmc.kam.rho())/evtmc.kam.rho());
	            //hppiires[parti]->Fill(evt.pip.rho()-evtmc.pip.rho());
	      }
	    }
	 }

	// fill momentum spectrum
	 //if (pt1<0.4 || pt1>0.8) continue;
	 //if (p1<0.6 || p1>0.8) continue;
	 hmD0->Fill(mass);
	 if (mass>m0-sigma_m && mass<m0+sigma_m){
	   hpka->Fill(p1);
	   hppi->Fill(p2);
	 }
	 // if (Cut(ientry) < 0) continue;
       }
   }

   TFile *ftmp = new TFile("P_cmp.root","recreate");
   double xxka[Npart];// = {0.5,0.7,0.9,1.1,1.3};
   double xeka[Npart];// = {0.1,0.1,0.1,0.1,0.1};
   double yyka[Npart],yeka[Npart];
   double yykaave, yekaave;
   double yypi[Npart],yepi[Npart];
   double yypiave, yepiave;
   double yymm[Npart],yemm[Npart];
   double yymmave, yemmave;
   ftmp->WriteTObject(hpka,"hpka_DKpi");
   ftmp->WriteTObject(hppi,"hppi_DKpi");
   hpkares->Fit("gaus","R","",-0.01,0.01);
   yykaave = hpkares->GetFunction("gaus")->GetParameter(1);
   yekaave = hpkares->GetFunction("gaus")->GetParError(1);
   ftmp->WriteTObject(hpkares,"hpkares_DKpi");
   hppires->Fit("gaus","R","",-0.01,0.01);
   yypiave = hppires->GetFunction("gaus")->GetParameter(1);
   yepiave = hppires->GetFunction("gaus")->GetParError(1);
   ftmp->WriteTObject(hppires,"hppires_DKpi");
   ftmp->WriteTObject(hmD0,"hmD0_DKpi");
   for (int parti=0; parti<Npart;parti++){
     xxka[parti] = (pcut[parti]+pcut[parti+1])/2.0;
     xeka[parti] = (pcut[parti+1]-pcut[parti])/2.0;
     hpkaires[parti]->Fit("gaus","R","",-0.01,0.01);
     yyka[parti] = hpkaires[parti]->GetFunction("gaus")->GetParameter(1);
     yeka[parti] = hpkaires[parti]->GetFunction("gaus")->GetParError(1);
     ftmp->WriteTObject(hpkaires[parti]);
     hppiires[parti]->Fit("gaus","R","",-0.01,0.01);
     yypi[parti] = hppiires[parti]->GetFunction("gaus")->GetParameter(1);
     yepi[parti] = hppiires[parti]->GetFunction("gaus")->GetParError(1);
     ftmp->WriteTObject(hppiires[parti]);
     hmtheta[parti]->Fit("gaus","R","",m0-sigma_m,m0+sigma_m);
     yymm[parti] = hmtheta[parti]->GetFunction("gaus")->GetParameter(1);
     yemm[parti] = hmtheta[parti]->GetFunction("gaus")->GetParError(1);
     ftmp->WriteTObject(hmtheta[parti]);
   }
   TGraphErrors* greshka = new TGraphErrors(Npart,xxka,yyka,xeka,yeka);
   greshka->SetName("pkresolutionshift_pk");
   greshka->GetXaxis()->SetTitle("p_{K} (GeV/c)");
   greshka->GetYaxis()->SetTitle("average (p_{K}(rec)-p_{K}(truth))/p_{K}(truth) (GeV/c)");
   greshka->Write();
   delete greshka;
   TGraphErrors* greshpi = new TGraphErrors(Npart,xxka,yypi,xeka,yepi);
   greshpi->SetName("ppiresolutionshift_pka");
   greshpi->GetXaxis()->SetTitle("p_{K} (GeV/c)");
   greshpi->GetYaxis()->SetTitle("average (p_{#pi}(rec)-p_{#pi}(truth))/p_{#pi}(truth) (GeV/c)");
   greshpi->Write();
   delete greshpi;
   TGraphErrors* greshmm = new TGraphErrors(Npart,xxka,yymm,xeka,yemm);
   greshmm->SetName("massshift_pka");
   greshmm->GetXaxis()->SetTitle("p_{K} (GeV/c)");
   greshmm->GetYaxis()->SetTitle("mass shift(MC) (GeV/c)");
   greshmm->Write();
   delete greshmm;

   ftmp->Close();
   delete ftmp;
   f->cd();
   return ;

   //vars->Write();
   //h2p->Write();
   //thedis->Write();
   //thedis2->Write();
   for (int parti=0;parti<Npart;parti++){
	 hmtheta[parti]->Write();
	 std::cout<<"processed part "<<parti<<std::endl;
	 std::cout<<"entry is "<<hmtheta[parti]->GetEntries()<<'\t'<<hmtheta[parti]->GetMaximumBin()<<std::endl;
	 if (hmtheta[parti]->GetEntries() > 50
	   &&hmtheta[parti]->GetMaximumBin()>20
	   &&hmtheta[parti]->GetMaximumBin()<70){
	   partmap.push_back(parti);
	 }
   }
   // ~~~~~~~~ draw end
   
     char namesfx[100];
     sprintf(namesfx,"nocut");
     KPI::FitSpe(evts_set2,beamene,namesfx);
     sprintf(namesfx,"nocut_mc");
     KPI::FitSpe(evts_setmc,beamene,namesfx);
   for (int ip=0;ip<partmap.size();ip++){
     int ipart = partmap.at(ip);
     sprintf(namesfx,"%d",ipart);
     KPI::FitSpe(evts_set[ipart],beamene,namesfx);
   }

}


void KPI::FitSpe(std::vector<D0KPI> &evts, double beame, char *namesfx)
{
  double peak = 1.86484;
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

  TCanvas *c1;
/*
  // dm Vs fk
  fpi = 1.00;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0KPI evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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

  ofstream ofpar("parkpipi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b
   << ", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np+1)/2] 
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

 // dm Vs fk
  fpi = 1.001;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0KPI evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
    deltapeaks[i] = peaks[i] - peak;
    if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   c1 = new TCanvas("c1_fk","c1_fpi1.001");
   //TGraphErrors *;
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
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


  // dm Vs fk
  fpi = 0.999;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0KPI evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
  */

  // dm Vs fpi
  //fk = 0.999830;
//for (int i=0; i<np; i++){
//  dataraw->Reset();
//  for (int evtid=0; evtid<evts.size();evtid++){
//        D0KPI evt = evts.at(evtid);
//        evt.setCorrectionFactors(fk,factors[i]);
//        mass = evt.m();
//    if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
//      }
//  sprintf(tmpchr,"factor_%f",factors[i]);
//  KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
//      deltapeaks[i] = peaks[i] - peak;
//      if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
//}

// //TCanvas *c1 = new TCanvas("c1_1","c1");
// c1 = new TCanvas("c1_fpi","c1_fk1.00");
// c1->Clear();
// TGraphErrors *
// graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
// graph1->SetTitle("delta peak");
// graph1->SetMarkerStyle(5);
// graph1->Draw("AP");
// gStyle->SetOptFit(1111);
// TF1 *
// facfit = new TF1("facfit","[0]*(x-1)+[1]");
// facfit->SetParameters(0.3,1);
// facfit->SetParNames("slope","interupt");
// graph1->Fit(facfit,"","",factors[0],factors[np-1]);
// double a = facfit->GetParameter(0);
// double ae= facfit->GetParError(0);
// double b = facfit->GetParameter(1);
// double be= facfit->GetParError(1);
// double factor = 1 - b/a;
// double factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);
// graph1->SetName("PsVF");
// graph1->Draw();

//TPaveText *
//pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
//sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
//pt->AddText(tmpchr);
//sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
//pt->AddText(tmpchr);
//pt->Draw();
//c1->Write("de_fpi");

//ofstream ofpar("parkpipi",std::ios::app);
//ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
// << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
// <<"dE = "<< a <<"*(x-1) + "<< b 
// <<", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np-1)/2]
// <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0KPI evt = evts.at(evtid);
      evt.setCorrectionFactors(1,1);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    double factor =1;
    sprintf(tmpchr,"%s_factor_%f",namesfx,factor);
    double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
   
 // dataraw->Reset();
 // TF1 ff("ff","1.00065+0.000630*x",0,2);
 // for (int evtid=0; evtid<evts.size();evtid++){
 //   D0KPI evt = evts.at(evtid);
 //   double fpi = ff.Eval(evt.pip.rho());
 //   evt.setCorrectionFactors(fk,fpi);
 //   mass = evt.m();
 //   if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
 // }
 // //double 
 // factor =0;
 // sprintf(tmpchr,"factor_%f",factor);
 // //double peakt,errt;
 // KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


/*
  // dm Vs fpi
  fk = 1.001;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factors[i]);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


  // dm Vs fpi
  fk = 0.999;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factors[i]);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factors[i]);
    KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
*/


  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void KPI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx, double &peak, double &peakerror)
{
   int nBins=100;
   bool largesample = false;
   //if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = 1.86484;
   double beamlow = 1.82;
   double beamup  = 1.90;
   // try to use roofit
   RooRealVar x("x","M(K#pi)",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.0050,0.0025,0.0075);
   RooRealVar sigma2("sigma2","width of gaussian",0.0075,0.0075,0.010);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   //RooRealVar co1("co1","coefficient #1",0,-100.,100.);
   //RooRealVar co4("co4","coefficient #4",0);
   //RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",5000,0,10000000);//event number
   //RooRealVar signal2("signal2"," ",200,0,10000000);//event number
   RooRealVar sigfra("sigfra"," ",0.5,0.2,1.0);//event number
   RooRealVar background("background"," ",1000,0,1000000);
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooAddPdf sig("sig","signal",RooArgList(gaus,gaus2),sigfra);
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_kpi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("M(K #pi)"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //if (!largesample) {
     sum = new RooAddPdf("sum","sum",RooArgList(sig,ground),RooArgList(signal,background));
     Npar = 8;
   //}
   //else {
   //  sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
//	 Npar=8;
  // }
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe);
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
   delete xframe;
   delete dataset;
   delete sum;
   return;
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
   fChain = 0;
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
