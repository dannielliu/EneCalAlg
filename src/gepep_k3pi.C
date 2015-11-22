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
namespace K3PI{
  void FitSpe(std::vector<D0K3PI> &evts, double beame,  char* namesfx);
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

   std::cout<<"Toral entry is "<<nentries<< " @ "<< beamene << " GeV." <<std::endl;
// int nBins=40;
// double factorstart=0.995;
   double D0low=1.82;
   double D0up=1.90;
   double mk=0.493677;
   double mpi=0.13957018;
   double peakvalue=1.86484;// mD0
   //int pointNo=20;
 
   //ofpar<<"k- pi+ algrithm: will give factors for pion"<<std::endl;
 //ofstream ofpardetail;
 //ofpardetail.open("detail.txt",std::ios::app);
 //ofstream purepar;
 //purepar.open("par");
   char fname[100];
   sprintf(fname,"%s/plot_k3pi_AR.root",outputdir.c_str());
   TFile *f = new TFile(fname,"update");
  
   //TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
// TF1 *facfit = new TF1("facfit","[0]*(x-1)+[1]",0.9,1.1);
// facfit->SetParNames("slope","interupt");
// TH1D *h1	= new TH1D("h1","k- pi+ invariant mass",nBins,D0low,D0up);
// TCanvas *c1 = new TCanvas("c1","c1",800,600);
// TTree *dataraw = new TTree("dataraw","dataraw");
   double mass;
// dataraw->Branch("x",&mass,"x/D");
 
   TTree *vars = new TTree("vars","vars");
   double phi1,phi2;
   double costheta1,costheta2;
   double p1,p2,p3,p4;
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"pk/D");
   vars->Branch("p2",&p2,"ppi/D");
   vars->Branch("mass",&mass,"mass/D");


   Long64_t nbytes = 0, nb = 0;

   TH1D *hppires = new TH1D("hppires","p_{#pi} (rec-truth)",200,-0.1,0.1);
   TH1D *hpkares = new TH1D("hpkares","p_{K} (rec-truth)",200,-0.1,0.1);
   TH1D *hppi = new TH1D("hppi","p_{#pi}",200,0,2);
   TH1D *hpka = new TH1D("hpka","p_{K}",200,0,2);
   TH1D *hmD  = new TH1D("hmD" ,"M(K#pi#pi#pi)",200,1,2);
   TH1D *hmDc = new TH1D("hmDc","M(K#pi#pi#pi)",200,1,2);
   TH1D *htotp = new TH1D("htotp","total p",200,0,2);

   TH2D *h2p = new TH2D("h2p","#pi K momentum",200,0,2,200,0,2);
   TH2D *thedis = new TH2D("thedis","#theta",200,0,2,100,0,TMath::Pi());
   TH2D *thedis2= new TH2D("thedis2","#theta",100,0,TMath::Pi(),100,0,TMath::Pi());
   const  int Npart=1;
   double pcut[Npart+1]={0., 4.6};
   //double pcut[Npart+1]={0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
   //double coscut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
// double start=-1;
// double stop =1;
// for(int i=0;i<Npart+1;i++){
//       coscut[i] = (stop-start)/Npart*i+start;
// }
   double m0=peakvalue;
   double sigma_m = 0.0068;//0.0024 for phi,
   double width = 20.*sigma_m;
   //double mparticle=mk;
   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;
   std::vector<D0K3PI> evts_set[Npart];
   std::vector<D0K3PI> evts_set2;

   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   TH1D *hmtheta[Npart];
   TH1D *hpkaires[Npart];
   TH1D *hppiires[Npart];
   for (int partj=0;partj<Npart;partj++){
	 sprintf(name,"mass_part%d",partj);
	 hmtheta[partj] = new TH1D(name,name,100,D0low,D0up);
	 sprintf(name,"pkares_%d",partj);
	 hpkaires[partj] = new TH1D(name,name,200,-0.1,0.1);
	 sprintf(name,"ppires_%d",partj);
	 hppiires[partj] = new TH1D(name,name,200,-0.1,0.1);
   }

   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	 Long64_t ientry = LoadTree(jentry);
	 if (ientry < 0) break;
	 nb = fChain->GetEntry(jentry);   nbytes += nb;
	 //if (run >34028) break;

	 //truthalg->fChain->GetEntry(jentry);
         //cout<<nkap<<nkam<<npip<<npim<<"\t";
         //cout<<truthalg->mcnpip<<truthalg->mcnpim<<"\t";
	 //int parti,partj;
	 //double mass
	 double theta1,theta2,theta3,theta4;
	 int besidx=0;
	 double tmpdeltaold=100;
	 double tmpmass;
	 D0K3PI evtda,evt,evtmc;
	 D0K3PI tmpevt;
	 std::vector<D0K3PI> tmpv,tmpvmc;
	 tmpv.clear();
	 tmpvmc.clear();
	 HepLorentzVector kam, pip1, pip2, pim;
	 //if (npip==1){
	 //for(int i=0; i<1;i++){
	 //int i=0;
	 if (npim>=1 && nkam==1 && npip ==2){
	   kam.setVectM(Hep3Vector(kampx[0],kampy[0],kampz[0]),mk);
	   pip1.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
	   pip2.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mpi);
	   pim.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
	   evtda.set(kam,pip1,pip2,pim);
	   for (int i=1;i<npim;i++){
	     pim.setVectM(Hep3Vector(pimpx[i],pimpy[i],pimpz[i]),mpi);
	     tmpevt.set(kam,pip1,pip2,pim);
	     if (fabs(tmpevt.m()-m0)<fabs(evtda.m()-m0)) evtda = tmpevt;
	   }
	   tmpv.push_back(evtda);
         
	 /*
	 //kam.setVectM(Hep3Vector(truthalg->truthpkam[0],
	 //                        truthalg->truthpkam[1],
	 //     		   truthalg->truthpkam[2]),
	 //             mk);
	 //pip1.setVectM(Hep3Vector(truthalg->mcpippx[0],
	 //                         truthalg->mcpippy[0],
	 //     		    truthalg->mcpippz[0]),
	 //             mpi);
	 //pip2.setVectM(Hep3Vector(truthalg->mcpippx[1],
	 //                         truthalg->mcpippy[1],
	 //     		    truthalg->mcpippz[1]),
	 //             mpi);
	 //pim.setVectM(Hep3Vector(truthalg->mcpimpx[0],
	 //                        truthalg->mcpimpy[0],
	 //     		   truthalg->mcpimpz[0]),
	 //             mpi);
         //evtmc.set(kam,pip1,pip2,pim);
         //
         //for (int i=0;i<3;i++){
	 //  pim.setVectM(Hep3Vector(truthalg->mcpimpx[i],
	 //                          truthalg->mcpimpy[i],
	 //                          truthalg->mcpimpz[i]),
	 //               mpi);
         //  tmpevt.set(kam,pip1,pip2,pim);
	 //  if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
	 //}

	 //pip2.setVectM(Hep3Vector(truthalg->mcpippx[2],
	 //                         truthalg->mcpippy[2],
	 //     		    truthalg->mcpippz[2]),
	 //             mpi);
         //tmpevt.set(kam,pip1,pip2,pim);
	 //if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
         //for (int i=0;i<3;i++){
	 //  pim.setVectM(Hep3Vector(truthalg->mcpimpx[i],
	 //                          truthalg->mcpimpy[i],
	 //                          truthalg->mcpimpz[i]),
	 //               mpi);
         //  tmpevt.set(kam,pip1,pip2,pim);
	 //  if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
	 //}
	 //
	 //pip1.setVectM(Hep3Vector(truthalg->mcpippx[1],
	 //                         truthalg->mcpippy[1],
	 //     		    truthalg->mcpippz[1]),
	 //             mpi);
         //tmpevt.set(kam,pip1,pip2,pim);
	 //if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
         //for (int i=0;i<3;i++){
	 //  pim.setVectM(Hep3Vector(truthalg->mcpimpx[i],
	 //                          truthalg->mcpimpy[i],
	 //                          truthalg->mcpimpz[i]),
	 //               mpi);
         //  tmpevt.set(kam,pip1,pip2,pim);
	 //  if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
	 //}
	 //
	 */
	 //tmpvmc.push_back(evtmc);
	   //evt = evtda;
	 }
	 if (npip>=1 && nkap==1 && npim ==2){
	   kam.setVectM(Hep3Vector(kappx[0],kappy[0],kappz[0]),mk);
	   pip1.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
	   pip2.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mpi);
	   pim.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
	   evtda.set(kam,pip1,pip2,pim);
	   for (int i=1;i<npip;i++){
	     pim.setVectM(Hep3Vector(pippx[i],pippy[i],pippz[i]),mpi);
	     D0K3PI tmpevt;
	     tmpevt.set(kam,pip1,pip2,pim);
	     if (fabs(tmpevt.m()-m0)<fabs(evtda.m()-m0)) evtda = tmpevt;
	   }
	   tmpv.push_back(evtda);
	   
	   /*
	   kam.setVectM(Hep3Vector(truthalg->truthpkap[0],
	                           truthalg->truthpkap[1],
				   truthalg->truthpkap[2]),
	                mk);
	   pip1.setVectM(Hep3Vector(truthalg->mcpimpx[0],
	                            truthalg->mcpimpy[0],
				    truthalg->mcpimpz[0]),
	                mpi);
	   pip2.setVectM(Hep3Vector(truthalg->mcpimpx[1],
	                            truthalg->mcpimpy[1],
				    truthalg->mcpimpz[1]),
	                mpi);
	   pim.setVectM(Hep3Vector(truthalg->mcpippx[0],
	                           truthalg->mcpippy[0],
				   truthalg->mcpippz[0]),
	                mpi);
	   evtmc.set(kam,pip1,pip2,pim);
	   
           for (int i=0;i<3;i++){
	     pim.setVectM(Hep3Vector(truthalg->mcpippx[i],
	                             truthalg->mcpippy[i],
		                     truthalg->mcpippz[i]),
	                  mpi);
             tmpevt.set(kam,pip1,pip2,pim);
	     if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
	   }

	   pip2.setVectM(Hep3Vector(truthalg->mcpimpx[2],
	                            truthalg->mcpimpy[2],
				    truthalg->mcpimpz[2]),
	                mpi);
           tmpevt.set(kam,pip1,pip2,pim);
	   if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
           for (int i=0;i<3;i++){
	     pim.setVectM(Hep3Vector(truthalg->mcpippx[i],
	                             truthalg->mcpippy[i],
		                     truthalg->mcpippz[i]),
	                  mpi);
             tmpevt.set(kam,pip1,pip2,pim);
	     if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
	   }
	   
	   pip1.setVectM(Hep3Vector(truthalg->mcpimpx[1],
	                            truthalg->mcpimpy[1],
				    truthalg->mcpimpz[1]),
	                mpi);
           tmpevt.set(kam,pip1,pip2,pim);
	   if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
           for (int i=0;i<3;i++){
	     pim.setVectM(Hep3Vector(truthalg->mcpippx[i],
	                             truthalg->mcpippy[i],
		                     truthalg->mcpippz[i]),
	                  mpi);
             tmpevt.set(kam,pip1,pip2,pim);
	     if (fabs(tmpevt.m()-m0)<fabs(evtmc.m()-m0)) evtmc = tmpevt;
	   }
	   tmpvmc.push_back(evtmc);
	   */
	   //evt = evtda;
	 }
	
//////// // total invariant mass, D0 -> k- pi+
//////// p1 = evt.kam.rho();
//////// p2 = evt.pip1.rho();
//////// p3 = evt.pip2.rho();
//////// p4 = evt.pim.rho();
//////// double pt1 = evt.kam.perp();
//////// double pt2 = evt.pip1.perp();
//////// double pt3 = evt.pip2.perp();
//////// double pt4 = evt.pim.perp();
//////// theta1 = evt.kam.theta();
//////// theta2 = evt.pip1.theta();
//////// mass = evt.m();
	 //if (parti>=Npart || partj>=Npart || parti<0 || partj<0) continue;
       for (int ie=0;ie<tmpv.size();ie++){
	 evt = tmpv.at(ie);
	 //evtmc = tmpvmc.at(ie);
	 mass = evt.m();
         p1 = evt.kam.rho();
         p2 = evt.pip1.rho();
         p3 = evt.pip2.rho();
         p4 = evt.pim.rho();
	 for (int parti=0; parti<Npart; parti++){
	   if (p1>pcut[parti]&&p1<pcut[parti+1]){
	   //&&theta2>thetacut[partj]&&theta2<thetacut[partj+1])
	     hmtheta[parti]->Fill(mass);
	     if (mass>m0-width/2.-0.02 && mass<m0+width/2.+0.02)
	     {
	   //  evts_set2.push_back(evt);
	       evts_set[parti].push_back(evt);
	   //  hpkaires[parti]->Fill(evt.kam.rho()-evtmc.kam.rho());
	   //  hppiires[parti]->Fill(evt.pim.rho()-evtmc.pim.rho());
           //  hpkares->Fill(evt.kam.rho()-evtmc.kam.rho());
           //  hppires->Fill(evt.pim.rho()-evtmc.pim.rho());
	   //  htotp->Fill(evtmc.Get4P().rho());
	       //vars->Fill();
	     }
	   }
	 }

	 //if (pt1<0.6 || pt1>0.8) continue;
	 //if (p1<0.6 || p1>0.8) continue;
	 hmD->Fill(mass);
	 //fill momentum spectrum
         if (mass>m0-2*sigma_m && mass<m0+2*sigma_m){
           hpka->Fill(evt.kam.rho());
           hppi->Fill(evt.pim.rho());

//         hpkares->Fill(evtda.kam.rho()-evtmc.kam.rho());
//         hppires->Fill(evtda.pim.rho()-evtmc.pim.rho());
         }
	
       }
   }
   //vars->Write();
   //h2p->Write();
   //thedis->Write();
   //thedis2->Write();

/*
   TFile *ftmp = new TFile("P_cmp.root","update");
   double xxka[5] = {0.3,0.5,0.7,0.9,1.1};
   double xeka[5] = {0.1,0.1,0.1,0.1,0.1};
   double yyka[5],yeka[5];
   double yykaave, yekaave;
   double yypi[5],yepi[5];
   double yypiave, yepiave;
   double yymm[5],yemm[5];
   double yymmave, yemmave;
   ftmp->WriteTObject(hpka,"hpka_DKpipipi");
   ftmp->WriteTObject(hppi,"hppi_DKpipipi");
   ftmp->WriteTObject(hmD,"hmD_DKpipipi");
   ftmp->WriteTObject(htotp,"htotp_DKpipipi");
   hpkares->Fit("gaus","R","",-0.01,0.01);
   yykaave = hpkares->GetFunction("gaus")->GetParameter(1);
   yekaave = hpkares->GetFunction("gaus")->GetParError(1);
   ftmp->WriteTObject(hpkares,"pkares_K3pi");
   hppires->Fit("gaus","R","",-0.01,0.01);
   yypiave = hppires->GetFunction("gaus")->GetParameter(1);
   yepiave = hppires->GetFunction("gaus")->GetParError(1);
   ftmp->WriteTObject(hppires,"ppires_K3pi");
   for (int ip=0;ip<Npart;ip++){
     hpkaires[ip]->Fit("gaus","R","",-0.01,0.01);
     yyka[ip] = hpkaires[ip]->GetFunction("gaus")->GetParameter(1);
     yeka[ip] = hpkaires[ip]->GetFunction("gaus")->GetParError(1);
     ftmp->WriteTObject(hpkaires[ip]);
     hppiires[ip]->Fit("gaus","R","",-0.01,0.01);
     yypi[ip] = hppiires[ip]->GetFunction("gaus")->GetParameter(1);
     yepi[ip] = hppiires[ip]->GetFunction("gaus")->GetParError(1);
     ftmp->WriteTObject(hppiires[ip]);
     hmtheta[ip]->Fit("gaus","R","",m0-sigma_m,m0+sigma_m);
     yymm[ip] = hmtheta[ip]->GetFunction("gaus")->GetParameter(1);
     yemm[ip] = hmtheta[ip]->GetFunction("gaus")->GetParError(1);
     ftmp->WriteTObject(hmtheta[ip]);
   }
   TGraphErrors* greshka = new TGraphErrors(5,xxka,yyka,xeka,yeka);
   greshka->SetName("pkresolutionshift_pk");
   greshka->GetXaxis()->SetTitle("p_{K} (GeV/c)");
   greshka->GetYaxis()->SetTitle("average p_{K}(data)-p_{K}(MC) (GeV/c)");
   greshka->Write();
   delete greshka;
   TGraphErrors* greshpi = new TGraphErrors(5,xxka,yypi,xeka,yepi);
   greshpi->SetName("ppiresolutionshift_pka");
   greshpi->GetXaxis()->SetTitle("p_{K} (GeV/c)");
   greshpi->GetYaxis()->SetTitle("average p_{#pi}(data)-p_{#pi}(MC) (GeV/c)");
   greshpi->Write();
   delete greshpi;
   TGraphErrors* greshmm = new TGraphErrors(5,xxka,yymm,xeka,yemm);
   greshmm->SetName("massshift_pka");
   greshmm->GetXaxis()->SetTitle("p_{K} (GeV/c)");
   greshmm->GetYaxis()->SetTitle("M(D^{0}) (GeV/c)");
   greshmm->Write();
   delete greshmm;
   ftmp->Close();
   delete ftmp;
   f->cd();
   return ;
*/

   // for several parts, check the size first.
   for (int parti=0;parti<Npart;parti++){
        // hmtheta[parti]->Write();
         std::cout<<"processed part "<<parti<<std::endl;
         std::cout<<"entry is "<<hmtheta[parti]->GetEntries()<<'\t'<<hmtheta[parti]->GetMaximumBin()<<std::endl;
         if (hmtheta[parti]->GetEntries() > 50
           &&hmtheta[parti]->GetMaximumBin()>20
           &&hmtheta[parti]->GetMaximumBin()<70){
           partmap.push_back(parti);
         }
   }
   // ~~~~~~~~ draw end
   
   hpka->Write();
   hppi->Write();
   //return ;
   double pka_tot = hpka->GetMean();
   double pkasig_tot = hpka->GetRMS();
   double ppi_tot = hppi->GetMean();
   double ppisig_tot = hppi->GetRMS();
   ofstream outf("park3pi.txt",std::ios::app);
   outf<<"p_K = " << pka_tot << " +/- "<< pkasig_tot <<"\t";
   outf<<"p_pi = " << ppi_tot << " +/- "<< ppisig_tot <<endl;
   outf.close();


   for (int ip=0;ip<partmap.size();ip++){
     int ipart = partmap.at(ip);
     char namesfx[100];
     sprintf(namesfx,"%d",ipart);
     K3PI::FitSpe(evts_set[ipart],beamene,namesfx);
//   sprintf(namesfx,"%d_nocut",ipart);
//   K3PI::FitSpe(evts_set2,beamene,namesfx);
   }
}


void K3PI::FitSpe(std::vector<D0K3PI> &evts, double beame, char *namesfx)
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

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0K3PI evt = evts.at(evtid);
      evt.setCorrectionFactors(1.00,1.00);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"%s_factor_0",namesfx);
    double peakt,errt;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


 // dataraw->Reset();
 // for (int evtid=0; evtid<evts.size();evtid++){
 //       D0K3PI evt = evts.at(evtid);
 //       // with f from R scan
 //       evt.setCorrectionFactors(1.000170,1.000785);
 //       // with f from 4230
 //       //evt.setCorrectionFactors(1.000686,1.000488);
 //       mass = evt.m();
 //   if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
 // }
 // sprintf(tmpchr,"factor_170785");
 // //double peakt,errt;
 // K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);



  // dm Vs fk
  fpi = 1.0;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0K3PI evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fk_%f",factors[i]);
    peaks[i] = peak + 0.3*(factors[i]-1) - 0.001;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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

  ofstream ofpar("park3pi.txt",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
   <<"dE = "<< a <<"*(x-1) + "<< b
   << ", ae = "<< ae << " be = "<< be  << " perr  = " << peakerrors[(np+1)/2] 
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fk_%f",factor);
    //double peakt,errt;
    peakt = peak + 0.3*(factor-1) -0.001;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


/*
 // dm Vs fk
  fpi = 1.001;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0K3PI evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
*/
/*
  // dm Vs fk
  fpi = 0.999;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0K3PI evt = evts.at(evtid);
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor,fpi);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
*/ 
  
 
  // dm Vs fpi
  //fk = 0.999830;
  fk = 1.0;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      D0K3PI evt = evts.at(evtid);
      evt.setCorrectionFactors(fk,factors[i]);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fpi_%f",factors[i]);
    peaks[i] = peak + 1.02*(factors[i]-1) -0.001;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fpi_%f",factor);
    //double peakt,errt;
    peakt = peak + 1.02*(factor-1) -0.001;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


  //dataraw->Reset();
  //TF1 ff("ff","1.00065+0.000630*x",0,2);
  //double f1 = fk;
  //for (int evtid=0; evtid<evts.size();evtid++){
  //      D0K3PI evt = evts.at(evtid);
  //      double f2 = ff.Eval(evt.pip1.rho());
  //      double f3 = ff.Eval(evt.pip2.rho());
  //      double f4 = ff.Eval(evt.pim.rho());
  //      //evt.setCorrectionFactors(fk,factor);
  //      evt.setCorrectionFactors(f1,f2,f3,f4);
  //      mass = evt.m();
  //  if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
  //}
  //double factor = 0;
  //sprintf(tmpchr,"factor_%f",factor);
  ////double peakt,errt;
  //K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


/*
  // dm Vs fpi
  fk = 1.001;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factors[i]);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factors[i]);
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);


  // dm Vs fpi
  fk = 0.999;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factors[i]);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    sprintf(tmpchr,"factor_%f",factors[i]);
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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
	  D0K3PI evt = evts.at(evtid);
	  evt.setCorrectionFactors(fk,factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    //double peakt,errt;
    K3PI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
*/


  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void K3PI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx, double &peak, double &peakerror)
{
   int entries = dataraw->GetEntries();
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = 1.86484;
   double beamlow = 1.82;
   double beamup  = 1.90;

   TH1D *hmass = new TH1D("hmass","M(K 3#pi)",100, beamlow, beamup);
   RooDataHist *datahist;
   if (largesample) {
     dataraw->Draw("x>>hmass");
   }
   // try to use roofit
   if (fabs(peak-0)>1e-3) peakvalue = peak;
   RooRealVar x("x","M(K 3#pi)",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.004,0.0025,0.005);
   RooRealVar sigma2("sigma2","width of gaussian",0.01,0.005,0.020);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   //RooRealVar co1("co1","coefficient #1",0,-100.,100.);
   //RooRealVar co4("co4","coefficient #4",0);
   //RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",15000,0,10000000);//event number
   RooRealVar sigfra("sigfra"," ",0.5,0.3,1.0);//event number
   //RooRealVar signal2("signal2"," ",2000,0,10000000);//event number
   RooRealVar background("background"," ",30000,0,1000000);
  
   RooRealVar a0("a0","coefficient #0",1,-100000,100000);
   RooRealVar a1("a1","coefficient #1",1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
  
   RooAddPdf sig("sig","signal",RooArgList(gaus,gaus2),sigfra);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_kpipi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("M(K^{-} #pi^{+} #pi^{+} #pi^{-})"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   if (largesample)  datahist = new RooDataHist(tmpchr,"data",RooArgSet(x),hmass);
   //if (!largesample) {
     sum = new RooAddPdf("sum","sum",RooArgList(sig,ground),RooArgList(signal,background));
     Npar = 8;
   //}
   //else {
   //  sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
//	 Npar=8;
  // }
   if (largesample) {sum->fitTo(*datahist);datahist->plotOn(xframe);}
   else {
     sum->fitTo(*dataset,Range(beamlow,beamup));
     dataset->plotOn(xframe);
   }
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
  sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
  pt->AddText(tmpchr);
//if (largesample){
//  sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
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

   ofstream outf("park3pi.txt",std::ios::app);
   outf<<beame<<'\t'<< /* namesfx<<"\t"<< */ mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   peak = mean.getVal();
   peakerror = mean.getError();
   outf.close();
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   if (largesample) delete hmass;
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
