#define mumu_cxx
#include "mumu.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "function.h"
#include "bes3plotstyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <fstream>
#include <TFile.h>
#include <TFolder.h>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "EventClass.h"
#include <iomanip>
extern std::string outputdir;
using namespace RooFit;
using namespace std;

namespace MUMU{
  void FitSpe(std::vector<APair> &evts, double beame,  char* namesfx);
  void FitSpectrum(TTree *&dataraw,double beame, char* namesfx, double &peak, double &peakerror);
}


void mumu::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L mumu.C
//      Root > mumu t
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
  
  std::cout<<"Total entry is "<<nentries<<std::endl;
  //int nBins=12;
  int nBins=100;
 
  double mmu=0.105658;

  double jpsilow = 2.9;
  double jpsiup  = 3.3;
  
  char fname[1000];
  sprintf(fname,"%s/plot_mumu.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");
  
  TCanvas *c1=new TCanvas("c1","",800,600);
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
  vars->Branch("mass",&mass,"mass/D");

  APair evt;
  std::vector<APair> evts;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      HepLorentzVector mup,mum;
      mup.setVectM(Hep3Vector(lepx4[0],lepy4[0],lepz4[0]),mmu);
      mum.setVectM(Hep3Vector(lepx4[1],lepy4[1],lepz4[1]),mmu);
      evt.set(mup,mum);
      p1 = evt.partilep.rho();
      p2 = evt.partilem.rho();
      mass = evt.m();
      if (mass > jpsilow && mass<jpsiup){
        vars->Fill();
	evts.push_back(evt);
      }
      
   }
   vars->Write();

   MUMU::FitSpe(evts, beamene, 0);
}

void MUMU::FitSpe(std::vector<APair> &evts, double beame,  char* namesfx)
{
  double peakv = 3.097;
  double beamlow= peakv-0.3;
  double beamup = peakv+0.3;
  // for factor fit
 
  TTree *dataraw = new TTree("dataraw","dataraw");
  double mass;
  dataraw->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
  
  int np = 20;
  double factors[np];
  //double fpi=1.00;
  double factorserr[np];
  for (int i=0; i<np; i++){
    factors[i] = (1.005-0.995)/np*i+0.995;
    factorserr[i] = 0;
  }
  double peaks[np];
  double deltapeaks[np];
  double peakerrors[np];

  char tmpchr[100];

  //~~~~~~~~~~part start~~~~~~~~
  for (int id=0; id < np; id++){
    dataraw->Reset();
    for (Long64_t i=0; i<evts.size(); i++) {
       // total invariant mass
       APair evt = evts.at(i);
       evt.setCorrectionFactors(factors[id]);
       mass = evt.m();
       if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"f_%.4f_%s",factors[id],namesfx);
    MUMU::FitSpectrum(dataraw,peakv-0.002+3.0*(factors[id]-1),tmpchr,deltapeaks[id],peakerrors[id]);
    deltapeaks[id] -= 3.096916;
    if (peakerrors[id]<2e-5) peakerrors[id] = 2e-3;
  }
   
   TCanvas *c1 = new TCanvas("c1_1","c1");
   TGraphErrors *graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   TF1 *facfit = new TF1("facfit","[0]*(x-[1])");
   facfit->SetParLimits(0,0.,10);
   facfit->SetParLimits(1,0.9,1.1);
   facfit->SetParameters(1,1);
   facfit->SetParNames("slope","factor");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   double factor = facfit->GetParameter(1);
   double factorerr = peakerrors[np/2]/facfit->GetParameter(0);
   graph1->SetName("PsVF");
   graph1->Draw();

  TPaveText *pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"slope = %1.6f #pm %1.6f",facfit->GetParameter(0),facfit->GetParError(0));
  pt->AddText(tmpchr);
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fk");

  ofstream ofpar("parmumu",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed) << setprecision(6)
   << '\t' << factor << "\t" << factorerr << "\t" << facfit->GetParError(0) << '\t'
   << resetiosflags(ios::fixed) << std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  APair evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor);
	  mass = evt.m();
        if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%.4f",factor);
    double peakt,errt;
    MUMU::FitSpectrum(dataraw,peakv,tmpchr,peakt,errt);


  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void MUMU::FitSpectrum(TTree *&dataraw,double beame, char* namesfx, double &peak, double &peakerror)
{
   int nBins=100;
   bool largesample = false;
   int Npar;
   double peakvalue = beame;//3.096916;
   double beamlow = peakvalue-0.15;
   double beamup  = peakvalue+0.15;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.014,0.01,0.02);
   RooRealVar sigma2("sigma2","width of gaussian",0.006,0.005,0.01);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
     RooRealVar co1("co1","coefficient #1",0.1,-100.,100.);
     RooRealVar co2("co2","coefficient #2",0,-100.,100.);
     RooRealVar co4("co4","coefficient #4",0);
     RooChebychev ground("ground","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",70,0,10000000);//event number
   RooRealVar signal2("signal2"," ",70,0,10000000);//event number
   RooRealVar background("background"," ",240,0,1000000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
     
   RooRealVar alpha1("alpha1","#alpha",3.5,3.2,4.4);
   RooRealVar nnn1("n1","n",100,50,150);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("c2","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_mumu_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit mumu"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
     sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
     Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   dataset->plotOn(xframe,Binning(nBins));
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
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

  // ofstream outf("parmumu",std::ios::app);
  // outf<<beame<<'\t'<< /* namesfx<<"\t"<< */ mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   peak = mean.getVal();
   peakerror = mean.getError();
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}
