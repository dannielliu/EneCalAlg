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
#include "EventClass.h"
#include <iomanip>
extern std::string outputdir;
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
// ofstream ofpardetail;
// ofpardetail.open("detail.txt",std::ios::app);
// ofstream purepar;
// purepar.open("par");
   char fname[100];
   sprintf(fname,"%s/plot_kpi_AR.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");
  
   double mass;
 
   TTree *vars = new TTree("vars","vars");
   double phi1,phi2;
   double costheta1,costheta2;
   double p1,p2;
   vars->Branch("indexmc", &indexmc, "indexmc/I");
   vars->Branch("pdgid", pdgid, "pdgid[indexmc]/I");
   vars->Branch("motheridx", motheridx, "motheridx[indexmc]/I");
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"pk/D");
   vars->Branch("p2",&p2,"ppi/D");
   vars->Branch("mass",&mass,"mass/D");


   Long64_t nbytes = 0, nb = 0;

   TH1D *hppi = new TH1D("hppi","pt_{#pi}",200,0,2);
   TH1D *hpka = new TH1D("hpka","pt_{K}",200,0,2);
   TH1D *hmD0 = new TH1D("hmD0","M_{K #pi}",200,1,2);
   
   TH2D *h2p = new TH2D("h2p","#pi K momentum",200,0,2,200,0,2);
   TH2D *thedis = new TH2D("thedis","#theta",200,0,2,100,0,TMath::Pi());
   TH2D *thedis2= new TH2D("thedis2","#theta",100,0,TMath::Pi(),100,0,TMath::Pi());
   
   const  int Npart=1;
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
 //double pcut[Npart+1]={0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   double pcut[Npart+1]={0., 4.6};

   double m0=peakvalue;
   double sigma_m = 0.0068;//0.0024 for phi,
   double width = 20.*sigma_m;
   //double mparticle=mk;
   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;
   std::vector<D0KPI> evts_set[Npart];
   std::vector<D0KPI> evts_set2;

   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   TH1D *hmtheta[Npart];
   for (int partj=0;partj<Npart;partj++){
	 sprintf(name,"mass_part%d",partj);
	 hmtheta[partj] = new TH1D(name,name,100,D0low,D0up);
   }

   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	 Long64_t ientry = LoadTree(jentry);
	 if (ientry < 0) break;
	 nb = fChain->GetEntry(jentry);   nbytes += nb;
	 //if (run >34028) break;

	 //int parti,partj;
	 //double mass;
	 double theta1,theta2;
	 int besidx=0;
	 double tmpdeltaold=100;
	 double tmpmass;
	 D0KPI evt;
	 HepLorentzVector kam, pip;
	 //if (npip==1){
	 //for(int i=0; i<1;i++){
	 int i=0;
	 if (npim==1 && nkap==1){
	   kam.setVectM(Hep3Vector(kappx[0],kappy[0],kappz[0]),mk);
	   pip.setVectM(Hep3Vector(pimpx[i],pimpy[i],pimpz[i]),mpi);
	   evt.set(kam,pip);
	 }
	 if (npip==1 && nkam==1){
	   kam.setVectM(Hep3Vector(kampx[0],kampy[0],kampz[0]),mk);
	   pip.setVectM(Hep3Vector(pippx[i],pippy[i],pippz[i]),mpi);
	   evt.set(kam,pip);
	 }
	   
////////   tmpmass = evt.m();
////////   if(fabs(tmpmass-peakvalue)<tmpdeltaold){
////////         tmpdeltaold = fabs(tmpmass-peakvalue);
////////         besidx = i;
////////   }
//////// }
//////// 
//////// pip.setVectM(Hep3Vector(pippx[besidx],pippy[besidx],pippz[besidx]),mpi);
//////// evt.set(kam,pip);
	 
	 costheta1 = cos(evt.kam.theta());
	 costheta2 = cos(evt.pip.theta());
	 phi1 = evt.kam.phi();
	 phi2 = evt.pip.phi();
	 //if (fabs(costheta1)>0.5) continue;
	 //if (fabs(costheta2)>0.5) continue;
	 
	 // total invariant mass, D0 -> k- pi+
	 p1 = evt.kam.rho();
	 p2 = evt.pip.rho();
	 double pt1 = evt.kam.perp();
	 double pt2 = evt.pip.perp();
	 theta1= evt.kam.theta();
	 theta2= evt.pip.theta();
	 mass = evt.m();
	 //if (parti>=Npart || partj>=Npart || parti<0 || partj<0) continue;
	 for (int parti=0; parti<Npart; parti++){
	 if (p1>pcut[parti]&&p1<pcut[parti+1])
	 //if (costheta1>coscut[parti]&&costheta1<coscut[parti+1])
	   //&&theta2>thetacut[partj]&&theta2<thetacut[partj+1])
	   if (mass>m0-width/2.-0.02 && mass<m0+width/2.+0.02)
	   {
		 evts_set2.push_back(evt);
	      	 //if (p1<0.6 || p1>0.8) continue;
	/////////h2p->Fill(p1,p2);
	/////////thedis->Fill(p1,theta1);
	/////////thedis->Fill(p2,theta2);
	/////////thedis2->Fill(theta1,theta2);
	         hmtheta[parti]->Fill(mass);
		 evts_set[parti].push_back(evt);
	   }
	   if (mass>m0-2*sigma_m && mass<m0+2*sigma_m){
	     hpka->Fill(p1);
	     hppi->Fill(p2);
	     //vars->Fill();
	   }
	 }

	// fill momentum spectrum
	 //if (pt1<0.4 || pt1>0.8) continue;
	 //if (p1<0.6 || p1>0.8) continue;
	 hmD0->Fill(mass);
/////////if (mass>m0-2*sigma_m && mass<m0+2*sigma_m){
/////////  hpka->Fill(p1);
/////////  hppi->Fill(p2);
/////////}
	 // if (Cut(ientry) < 0) continue;
   }

   TFile *ftmp = new TFile("P_cmp.root","update");
   ftmp->WriteTObject(hpka,"hptka_DKpi");
   ftmp->WriteTObject(hppi,"hptpi_DKpi");
   ftmp->WriteTObject(hmD0,"hmD0_DKpi");
   ftmp->Close();
   delete ftmp;
   f->cd();
     return ;

   //vars->Write();
   //h2p->Write();
   //thedis->Write();
   //thedis2->Write();
   for (int parti=0;parti<Npart;parti++){
	 //hmtheta[parti]->Write();
	 std::cout<<"processed part "<<parti<<std::endl;
	 std::cout<<"entry is "<<hmtheta[parti]->GetEntries()<<'\t'<<hmtheta[parti]->GetMaximumBin()<<std::endl;
	 if (hmtheta[parti]->GetEntries() > 50
	   &&hmtheta[parti]->GetMaximumBin()>20
	   &&hmtheta[parti]->GetMaximumBin()<70){
	   partmap.push_back(parti);
	 }
   }
   // ~~~~~~~~ draw end
   
   double pka_tot = hpka->GetMean();
   double pkasig_tot = hpka->GetRMS();
   double ppi_tot = hppi->GetMean();
   double ppisig_tot = hppi->GetRMS();
   ofstream outf("parkpi.txt",std::ios::app);
   outf<<"p_K = " << pka_tot << " +/- "<< pkasig_tot <<"\t";
   outf<<"p_pi = " << ppi_tot << " +/- "<< ppisig_tot <<endl;
   outf.close();

   char namesfx[100];
 //sprintf(namesfx,"nocut");
 //KPI::FitSpe(evts_set2,beamene,namesfx);
   for (int ip=0;ip<partmap.size();ip++){
     int ipart = partmap.at(ip);
     sprintf(namesfx,"%d",ipart);
     KPI::FitSpe(evts_set[ipart],beamene,namesfx);
   }
   
   /*
   // find delta E relation with fk
// for (int loopi=0;loopi<partmap.size();loopi++){
//       int parti=partmap.at(loopi);

//       factor = factorstart;
//       int fittimes = 0;
//       double fpi = 1.00;
//       for (int i=0;i<pointNo;i++){
//        xframe = x.frame(Title("fit k pi"));

//        h1->Reset();
//        dataraw->Reset();
//        std::cout<<"factor is "<<factor<<std::endl;
//        int evtNo = evts_set[parti].size();
//        for (Long64_t evtid=0; evtid<evtNo;evtid++) {
//      	 D0KPI evt = evts_set[parti].at(evtid);
//      	 evt.setCorrectionFactors(factor,fpi);
//      	 mass = evt.m();
//      	 if (mass>D0low && mass<D0up){
//      	   dataraw->Fill();
//      	 }
//      	 // if (Cut(ientry) < 0) continue;
//        }

//        char tmpchr[100];
//        sprintf(tmpchr,"data_kpi_%02d",fittimes);
//        //data_kpi = new RooDataHist(tmpchr,"data_kpi",x,h1);
//        dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
//        sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
//        mean.setVal(peakvalue+0.8*(factor-1.002));
//        Npar=8;
//        //sigma.setVal(0.035);
//        //signal.setVal(1200);
//        //background.setVal(500);
//        //co1.setVal(0);
//        //co2.setVal(-0.2);
//        sum->fitTo(*dataset,Range(D0low,D0up));
//        dataset->plotOn(xframe);
//        sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
//        sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
//        sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
//        sum->plotOn(xframe);
//        xframe->Draw();
//        TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
//      //pt->SetBorderSize(0);
//      //pt->SetFillStyle(4000);
//      //pt->SetTextAlign(12);
//      //pt->SetTextFont(42);
//      //pt->SetTextSize(0.035);
//        sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"factor = %.6f",factor);
//        pt->AddText(tmpchr);
//        pt->Draw();
//        sprintf(tmpchr,"spe_at_%.6f",factor);
//        c1->Write(tmpchr);
//        //c1->Update();
//        //c1->Print(fitepsname);
//        delete dataset;
//        delete xframe;
//        delete sum;

//        // save pars
//        factors[i]=factor;
//        factorserr[i]=0;
//        deltapeaks[i] = mean.getVal() - peakvalue;
//        deltapeakserr[i] = mean.getError();
//        //if (deltapeakserr[i]<5e-5) deltapeakserr[i] = 1e-3;

//        fittimes++;
//        factor += factorstep;
// }// loop n point
// std::cout<<"entry is "<<nentries<<std::endl;
// //c1->Print(fiteps_stop);
// c1->Clear();
// 
// TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
// graph1->SetTitle("delta peak");
// graph1->SetMarkerStyle(5);
// graph1->Draw("AP");
// gStyle->SetOptFit(1111);
// facfit->SetParameters(1,1.);
// graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);

// //double factorerr=deltapeakserr[10]/facfit->GetParameter(1);
// double a = facfit->GetParameter(0);
// double ae= facfit->GetParError(0);
// double b = facfit->GetParameter(1);
// double be= facfit->GetParError(1);
// factor = 1 - b/a;
// double factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);
// ofpar<<beamene<<'\t'<<setiosflags(ios::fixed)<<setprecision(6)<<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
// <<"dE = "<< a <<"*(x-1) + "<< b <<", ae = "<< ae << " be = "<< be << " perr  = " << deltapeakserr[(pointNo-1)/2]<<resetiosflags(ios::fixed)<<std::endl;
// 
// std::cout<<"parameter@"<<beamene<<'\t'<<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<std::endl;
// purepar<<"\t"<<factor<<"\t"<<factorerr<<std::endl;
// sprintf(name,"factors_kpi_part%d",parti);
// graph1->SetName(name);
// graph1->Write("de_fk");

// // draw the best fit
// xframe = x.frame(Title("fit k pi"));
// h1->Reset();
// dataraw->Reset();
// std::cout<<"factor is "<<factor<<std::endl;
// int evtNo = evts_set[parti].size();
// for (Long64_t jentry=0; jentry<evtNo;jentry++) {
//       D0KPI evt = evts_set[parti].at(jentry);

//       evt.setCorrectionFactors(factor,fpi);
//       mass = evt.m();
//       if (mass>D0low && mass<D0up){
//         dataraw->Fill();
//       }
//       // if (Cut(ientry) < 0) continue;
// }
//      char tmpchr[100];
//      sprintf(tmpchr,"data_kpi");
//      //data_kpi = new RooDataHist(tmpchr,"data_kpi",x,h1);
//      dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
//      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
//      mean.setVal(peakvalue+0.05*(factor-1.0));
//      Npar=8;
//      //sigma.setVal(0.035);
//      //signal.setVal(1200);
//      //background.setVal(200);
//      //co1.setVal(0);
//      //co2.setVal(-0.2);
//      sum->fitTo(*dataset,Range(D0low,D0up));
//      dataset->plotOn(xframe);
//      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
//      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
//      sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
//      sum->plotOn(xframe);
//      xframe->Draw();
//      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
//      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"factor = %.6f #pm %.6f",factor,factorerr);
//      pt->AddText(tmpchr);
//      pt->Draw();
//      c1->Write("correctk");
//      xframe->SetName(name);
//      xframe->Write();
//      delete dataset;
//      delete xframe;
//      delete sum;
// 
// }// loop i end


// // find dE Vs fpi
// for (int loopi=0;loopi<partmap.size();loopi++){
//       int parti=partmap.at(loopi);

//       factor = factorstart;
//       int fittimes = 0;
//       double fk = 1.00;
//       for (int i=0;i<pointNo;i++){
//        xframe = x.frame(Title("fit k pi"));

//        h1->Reset();
//        dataraw->Reset();
//        std::cout<<"factor is "<<factor<<std::endl;
//        int evtNo = evts_set[parti].size();
//        for (Long64_t evtid=0; evtid<evtNo;evtid++) {
//      	 D0KPI evt = evts_set[parti].at(evtid);
//      	 evt.setCorrectionFactors(fk, factor);
//      	 mass = evt.m();
//      	 if (mass>D0low && mass<D0up){
//      	   dataraw->Fill();
//      	 }
//      	 // if (Cut(ientry) < 0) continue;
//        }

//        char tmpchr[100];
//        sprintf(tmpchr,"data_kpi_%02d",fittimes);
//        //data_kpi = new RooDataHist(tmpchr,"data_kpi",x,h1);
//        dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
//        sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
//        mean.setVal(peakvalue+0.8*(factor-1.002));
//        Npar=8;
//        //sigma.setVal(0.035);
//        //signal.setVal(1200);
//        //background.setVal(500);
//        //co1.setVal(0);
//        //co2.setVal(-0.2);
//        sum->fitTo(*dataset,Range(D0low,D0up));
//        dataset->plotOn(xframe);
//        sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
//        sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
//        sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
//        sum->plotOn(xframe);
//        xframe->Draw();
//        TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
//      //pt->SetBorderSize(0);
//      //pt->SetFillStyle(4000);
//      //pt->SetTextAlign(12);
//      //pt->SetTextFont(42);
//      //pt->SetTextSize(0.035);
//        sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
//        pt->AddText(tmpchr);
//        sprintf(tmpchr,"factor = %.6f",factor);
//        pt->AddText(tmpchr);
//        pt->Draw();
//        sprintf(tmpchr,"spe_at_%.6f",factor);
//        c1->Write(tmpchr);
//        //c1->Update();
//        //c1->Print(fitepsname);
//        delete dataset;
//        delete xframe;
//        delete sum;

//        // save pars
//        factors[i]=factor;
//        factorserr[i]=0;
//        deltapeaks[i] = mean.getVal() - peakvalue;
//        deltapeakserr[i] = mean.getError();
//        //if (deltapeakserr[i]<5e-5) deltapeakserr[i] = 1e-3;

//        fittimes++;
//        factor += factorstep;
// }// loop n point
// std::cout<<"entry is "<<nentries<<std::endl;
// c1->Clear();
// 
// TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
// graph1->SetTitle("delta peak");
// graph1->SetMarkerStyle(5);
// graph1->Draw("AP");
// gStyle->SetOptFit(1111);
// facfit->SetParameters(1,1.);
// graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);
// 
// double a = facfit->GetParameter(0);
// double ae= facfit->GetParError(0);
// double b = facfit->GetParameter(1);
// double be= facfit->GetParError(1);
// factor = 1 - b/a;
// double factorerr = sqrt(TMath::Power(ae/a,2)+TMath::Power(be/b,2))*fabs(b/a);
// //factor =facfit->GetParameter(0);
// //double factorerr=deltapeakserr[9]/facfit->GetParameter(1);
// ofpar<<beamene<<'\t'<<setiosflags(ios::fixed)<<setprecision(6)<<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<'\t'
// <<"dE = "<< a <<"*(x-1) + "<< b<<", ae = "<< ae << " be = "<< be  << " perr  = " << deltapeakserr[(pointNo-1)/2]<<resetiosflags(ios::fixed)<<std::endl;
// std::cout<<"parameter@"<<beamene<<'\t'<<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(0)<<std::endl;
// std::cout<<" deltapeakerr "<<deltapeakserr[9]<<" slope "<<facfit->GetParameter(1)<<std::endl;
// purepar<<"\t"<<factor<<"\t"<<factorerr<<std::endl;
// sprintf(name,"factors_kpi_part%d",parti);
// graph1->SetName(name);
// graph1->Write("de_fpi");

// // draw the best fit
// xframe = x.frame(Title("fit k pi"));
// h1->Reset();
// dataraw->Reset();
// std::cout<<"factor is "<<factor<<std::endl;
// int evtNo = evts_set[parti].size();
// for (Long64_t jentry=0; jentry<evtNo;jentry++) {
//       D0KPI evt = evts_set[parti].at(jentry);

//       evt.setCorrectionFactors(fk, factor);
//       mass = evt.m();
//       if (mass>D0low && mass<D0up){
//         dataraw->Fill();
//       }
//       // if (Cut(ientry) < 0) continue;
// }
//      char tmpchr[100];
//      sprintf(tmpchr,"data_kpi");
//      //data_kpi = new RooDataHist(tmpchr,"data_kpi",x,h1);
//      dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
//      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
//      mean.setVal(peakvalue+0.05*(factor-1.0));
//      Npar=8;
//      //sigma.setVal(0.035);
//      //signal.setVal(1200);
//      //background.setVal(200);
//      //co1.setVal(0);
//      //co2.setVal(-0.2);
//      sum->fitTo(*dataset,Range(D0low,D0up));
//      dataset->plotOn(xframe);
//      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
//      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
//      sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
//      sum->plotOn(xframe);
//      xframe->Draw();
//      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
//      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
//      pt->AddText(tmpchr);
//      sprintf(tmpchr,"factor = %.6f #pm %.6f",factor,factorerr);
//      pt->AddText(tmpchr);
//      pt->Draw();
//      c1->Write("correctpi");
//      xframe->SetName(name);
//      xframe->Write();
//      delete dataset;
//      delete xframe;
//      delete sum;
// 
// }// loop i end
// f->Close();
// ofpar.close();
// ofpardetail.close();
*/

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

  //TCanvas *c1;

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
    sprintf(tmpchr,"fk_%f",factors[i]);
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

   graph1->SetName("PsVFk");
   graph1->Draw();

  TPaveText *pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fk");

  ofstream ofpar("parkpi.txt",std::ios::app);
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
    sprintf(tmpchr,"fk_%f",factor);
    double peakt,errt;
    KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

 // dm Vs fk
  fpi = 1.001;
  /*
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
*/

  // dm Vs fk
  fpi = 0.999;
  /* 
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
  fk = 1.0;
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
       D0KPI evt = evts.at(evtid);
       evt.setCorrectionFactors(fk,factors[i]);
       mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"fpi_%f",factors[i]);
    KPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
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

  //dataraw->Reset();
  //for (int evtid=0; evtid<evts.size();evtid++){
  //  D0KPI evt = evts.at(evtid);
  //  evt.setCorrectionFactors(1,1);
  //  mass = evt.m();
  //  if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
  //}
  //double factor =1;
  //sprintf(tmpchr,"%s_factor_%f",namesfx,factor);
  //double peakt,errt;
  //KPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);
   
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
  ofpar.close();
  return;
}

void KPI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx, double &peak, double &peakerror)
{
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = 1.86484;
   double beamlow = 1.82;
   double beamup  = 1.90;
   // try to use roofit
   RooRealVar x("x","M(K#pi)",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue+0.0005,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.0050,0.0045,0.0075);
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
   TH1D *hmass = new TH1D("hmass","M(K #pi)",100, beamlow, beamup);
   RooDataHist *datahist;
   if (largesample) {
     dataraw->Draw("x>>hmass");
   }
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_kpi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("M(K #pi)"));
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
   //sum->fitTo(*dataset,Range(beamlow,beamup));
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

 //ofstream outf("parkpi.txt",std::ios::app);
 //outf<<beame<<'\t'<< /* namesfx<<"\t"<< */ mean.getVal()<<"\t"<<mean.getError()<<std::endl;
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
