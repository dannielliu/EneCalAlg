#define gepep_ppi_cxx
#include "gepep_ppi.h"
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
class Lambda
{
public:
   HepLorentzVector pp;
   HepLorentzVector pim;
public:
   Lambda() {};
   ~Lambda() {};

   void set(HepLorentzVector fpp, HepLorentzVector fpim)
   {
     pp = fpp;
     pim = fpim;
   }
   double m()
   {
     return (pp+pim).m();
   }
   void setCorrectionFactors(double fp, double fpi)
   {
     double mp = pp.m();
     double mpi = pim.m();
     pp.setVectM(pp.vect()*fp, mp);
     pim.setVectM(pim.vect()*fpi, mpi);
   }

};
namespace PPI{
  void FitSpe(std::vector<Lambda> &evts, double beame,  char* namesfx);
  void FitSpectrum(TTree *&dataraw,double beame, char* namesfx, double &peak, double &peakerror);
}

void gepep_ppi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_ppi.C
//      Root > gepep_ppi t
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

   fChain->GetEntry(0);
   double beamene = GetEnergy(run);
   std::cout<<"current beam energy is "<<beamene<<", run id "<<run<<std::endl;
   if (beamene < 0.1){
     std::cout<<"can not get a suitable beam energy!"<<std::endl;
      return;
   }

   double mpi = 0.13957;
   double mp  = 0.938272;
  
   char fname[1000];
   sprintf(fname,"%s/plot_ppi.root",outputdir.c_str());
   
   TFile *f=new TFile(fname,"RECREATE");
   double mass;
   double pproton, ppi;
   double thetap, thetapi;
   double phip, phipi;

   TTree *vars = new TTree("vars","vars");
   vars->Branch("pp", &pproton, "pp/D" );
   vars->Branch("ppi", &ppi,    "ppi/D");
   vars->Branch("thetap", &thetap, "thetap/D");
   vars->Branch("thetapi",&thetapi,"thetapi/D");
   vars->Branch("phip" ,&phip,"phip/D");
   vars->Branch("phipi",&phipi,"phipi/D");
   vars->Branch("mass", &mass, "mass/D");

   Lambda evt;
   std::vector<Lambda> evts;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (npp != 1 || npim != 1) continue;
      HepLorentzVector pp, pim;
      pp.setVectM(Hep3Vector(pppx[0],pppy[0],pppz[0]), mp);
      pim.setVectM(Hep3Vector(pimpx,pimpy,pimpz), mpi);
      evt.set(pp, pim);

      pproton = evt.pp.rho();
      ppi = evt.pim.rho();
      thetap = evt.pp.theta();
      thetapi = evt.pim.theta();
      phip  = evt.pp.phi();
      phipi = evt.pim.phi();
      mass = evt.m();

      vars->Fill();
      evts.push_back(evt);
   }
   vars->Write();

   PPI::FitSpe(evts,beamene, 0);
}


void PPI::FitSpe(std::vector<Lambda> &evts, double beame, char* namesfx)
{
  double peak = 1.115683;
  double beamlow=peak-0.1;
  double beamup=peak+0.1;
  // for factor fit
 
  TTree *dataraw = new TTree("dataraw","dataraw");
  double mass;
  dataraw->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
  char tmpchr[100];
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
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      Lambda evt = evts.at(evtid);
      if (evt.pim.rho()<0.4) fpi = 1.000902;
      else fpi = 1.000769;
      evt.setCorrectionFactors(factors[i],fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factors[i]);
    PPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
    deltapeaks[i] = peaks[i] - peak;
    if (peakerrors[i]<2e-5) peakerrors[i] = 1e-2;
  }

   TCanvas *c1 = new TCanvas("c1_1","c1");
   TGraphErrors *graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   TF1 *facfit = new TF1("facfit","[0]*(x-[1])");
   facfit->SetParameters(-0.01,1);
   facfit->SetParLimits(0,-2,2);
   facfit->SetParNames("slope","factor");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   double factor = facfit->GetParameter(1);
   double factorerr = deltapeaks[np/2]/facfit->GetParameter(0);
   graph1->SetName("PsVF");
   graph1->Draw();

  TPaveText *pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"slope = %1.6f #pm %1.6f", facfit->GetParameter(0), facfit->GetParError(0));
  pt->AddText(tmpchr);
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  pt->Draw();
  c1->Write("de_fk");

  ofstream ofpar("parppi",std::ios::app);
  ofpar<< beame << setiosflags(ios::fixed)<<setprecision(6)
   << '\t' <<factor<<"\t"<<factorerr<<"\t"<<facfit->GetParError(1)<<'\t'
   <<resetiosflags(ios::fixed)<<std::endl;

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
      Lambda evt = evts.at(evtid);
      if (evt.pim.rho()<0.4) fpi = 1.000902;
      else fpi = 1.000769;
      evt.setCorrectionFactors(factor,fpi);
      mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%f",factor);
    double peakt,errt;
    PPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void PPI::FitSpectrum(TTree *&dataraw,double beame, char* namesfx, double &peak, double &peakerror)
{
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = 1.115683;
   double beamlow = peakvalue-0.012;
   double beamup  = peakvalue+0.012;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.002,0.001,0.005);
   RooRealVar sigma2("sigma2","width of gaussian",0.006,0.005,0.01);
   RooGaussian fsig("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
     RooRealVar co1("co1","coefficient #1",0,-100.,100.);
     RooRealVar co2("co2","coefficient #2",0,-100.,100.);
     RooRealVar co3("co3","coefficient #3",0,-100.,100.);
     RooRealVar co4("co4","coefficient #4",0);
     RooChebychev ground("ground","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",100,0,10000000);//event number
   RooRealVar signal2("signal2"," ",200,0,10000000);//event number
   RooRealVar background("background"," ",400,0,1000000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));

 //RooRealVar alpha1("alpha1","#alpha",1.32,-5,5);
 //RooRealVar nnn1("n1","n",5,1,200);
 //RooCBShape fsig("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_ppi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit ppi"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(fsig,ground),RooArgList(signal,background));
   Npar = 6;
 //if (!largesample) {
 //  sum = new RooAddPdf("sum","sum",RooArgList(fsig,ground),RooArgList(signal,background));
 //  Npar = 6;
 //}
 //else {
 //  sum = new RooAddPdf("sum","sum",RooArgList(fsig,gaus2,ground),RooArgList(signal,signal2,background));
 //      Npar=8;
 //}
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(fsig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
  // if (dataraw->GetEntries()>2000) sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(4));
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

   ofstream outf("ppi",std::ios::app);
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



