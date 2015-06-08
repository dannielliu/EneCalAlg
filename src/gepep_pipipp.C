#define gepep_pipipp_cxx
#include "gepep_pipipp.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TPaveText.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "function.h"
#include "EventClass.h"
using namespace RooFit;
extern std::string outputdir;

void gepep_pipipp::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_pipipp.C
//      Root > gepep_pipipp t
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
   fChain->GetEntry(1);
   
   double mpi = 0.13957;
   double mp  = 0.938272;
   double beamene = 0;
   double beamlow, beamup;
   beamene = GetEnergy(run);
   //beamene = 3.097;
   beamlow = beamene - 0.05;
   beamup  = beamene + 0.05;
   double peakvalue = beamene;
   int nBins = 100;
   char name[1000];
   sprintf(name,"%s/plot_pipipp.root",outputdir.c_str());
   TFile *f = new TFile(name,"recreate");
   ofstream ofpeak("result_pipipp.txt",std::ios::app);

  double mass;
  TTree *dataraw = new TTree("dataraw","dataraw");
  dataraw->Branch("x",&mass,"x/D");
  TTree *datarawc = new TTree("datarawc","datarawc");
  datarawc->Branch("x",&mass,"x/D");
  TTree *vars = new TTree("vars","vars");
  double phi,phi1,phi2;
  double costheta,costheta1,costheta2;
  double p1,p2,p3,p4;
  vars->Branch("phi",&phi,"phi/D");
  vars->Branch("phi1",&phi1,"phi1/D");
  vars->Branch("phi2",&phi2,"phi2/D");
  vars->Branch("costheta",&costheta,"costheta/D");
  vars->Branch("costheta1",&costheta1,"costheta1/D");
  vars->Branch("costheta2",&costheta2,"costheta2/D");
  vars->Branch("p1",&p1,"p1/D");
  vars->Branch("p2",&p2,"p2/D");
  vars->Branch("p3",&p3,"p3/D");
  vars->Branch("p4",&p4,"p4/D");
  vars->Branch("mass",&mass,"mass/D");


  PiPiPP evt;
  std::vector<PiPiPP> evts_set;
   
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  HepLorentzVector pip,pim,pp,pm;
	  pip.setVectM(Hep3Vector(pippx,pippy,pippz),mpi);
	  pim.setVectM(Hep3Vector(pimpx,pimpy,pimpz),mpi);
	  pp.setVectM(Hep3Vector(protonpx,protonpy,protonpz),mp);
	  pm.setVectM(Hep3Vector(pbarpx,pbarpy,pbarpz),mp);
	  evt.set(pip,pim,pp,pm);
	  p1 = evt.pip.rho();
	  p2 = evt.pim.rho();
	  p3 = evt.pp.rho();
	  p4 = evt.pm.rho();

	  mass = evt.m();
	  if (mass>beamlow-0.002 && mass< beamup+0.002) evts_set.push_back(evt);
	  vars->Fill();
   }
   vars->Write();


   double factor = 1.000769;
   double fpi1,fpi2, fp1,fp2;
   fp1 = 1.000;
   fp2 = 1.000;
   int count[2] = {0,0};
   for (int ievt=0; ievt<evts_set.size(); ievt++){
         evt = evts_set.at(ievt);
	 mass = evt.m();
	 if (mass>beamlow && mass<beamup) {
	 	dataraw->Fill();
		count[0]++;
	}
	 p1 = evt.pip.rho();
	 p2 = evt.pim.rho();
	 p1<0.4 ? fpi1 = 1.000902 : fpi1 = factor;
	 p2<0.4 ? fpi2 = 1.000902 : fpi2 = factor;
	 evt.setCorrectionFactors(fpi1,fpi2,fp1,fp2);
	 mass = evt.m();
	 if (mass>beamlow && mass<beamup) {
	 	datarawc->Fill();
		count[1]++;
	 }
   }
   std::cout<<"raw size "<<count[0]<<std::endl;
   std::cout<<"cor size "<<count[1]<<std::endl;
   
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.009,0.001,0.02);
   RooRealVar sigma2("sigma2","width of gaussian",0.022,0.02,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar signal("signal"," ",100,0,10000000);//event number
   RooRealVar signal2("signal2"," ",200,0,10000000);//event number
   RooRealVar background("background"," ",200,0,1000000);
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
 

   TCanvas *c1 = new TCanvas();

// raw data part ===============
   sprintf(name,"data_pipipp_raw");
   xframe = x.frame(Title("fit pipipp"));
   dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();

  TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
  sprintf(name,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
  pt->AddText(name);
  sprintf(name,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
  pt->AddText(name);
  sprintf(name,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(name);
  sprintf(name,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(name);
  sprintf(name,"#chi^{2}/(%d-%d) = %5.6f",nBins,6,xframe->chiSquare(6));
  pt->AddText(name);
  pt->Draw();
  sprintf(name,"mass_spectrum_raw");
  c1->SetName(name);
  c1->Write();

  delete xframe;
  delete sum;
  delete dataset;
//==========raw part end
// save result
  ofpeak<<beamene<<'\t'<<mean.getVal()<<'\t'<<mean.getError()<<std::endl;

// ======= correction======
   sprintf(name,"data_pipipp_cor");
   xframe = x.frame(Title("fit pipipp"));
   dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*datarawc));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();

  pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
  sprintf(name,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
  pt->AddText(name);
  sprintf(name,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
  pt->AddText(name);
  sprintf(name,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(name);
  sprintf(name,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(name);
  sprintf(name,"#chi^{2}/(%d-%d) = %5.6f",nBins,6,xframe->chiSquare(6));
  pt->AddText(name);
  pt->Draw();
  sprintf(name,"mass_spectrum_cor");
  c1->SetName(name);
  c1->Write();

  delete xframe;
  delete sum;
  delete dataset;
// ======= correction end===

  ofpeak<<beamene<<'\t'<<mean.getVal()<<'\t'<<mean.getError()<<std::endl;

}
