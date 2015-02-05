#define gepep_fastpipill_cxx
#include "gepep_fastpipill.h"
#include <TH1D.h>
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
extern std::string outputdir;
using namespace RooFit;
namespace PIPILL{
  void FitAndSave(TH1D *&hmass, RooPlot *&xframe, RooDataHist *&data,
          RooRealVar &x,RooRealVar &mean, RooRealVar &sigma, RooRealVar &co1,
    RooRealVar &co2, RooRealVar &signal, RooRealVar &background,
    RooAddPdf *&sum, RooGaussian &gaus, RooChebychev &bkg,
    int idx,double factor,double psiplow,double psipup,std::string suffix="");
  void ResetVars(
          RooRealVar &mean, RooRealVar &sigma, RooRealVar &co1,
    RooRealVar &co2, RooRealVar &signal, RooRealVar &background,
    double factor=1, double meanv=3.686, double sigmav=0.0023,
    double signalv=60, double backgroundv=10, double co1v=0,
    double co2v=0
          );
}

bool gepep_fastpipill::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fastpipill.C
//      Root > gepep_fastpipill t
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
  ofstream logfile("log.txt");
  std::cout.rdbuf(logfile.rdbuf());
  std::clog.rdbuf(logfile.rdbuf());
  // std::cerr.rdbuf(logfile.rdbuf());
  
  if (fChain == 0) return false;
  Long64_t nentries = fChain->GetEntriesFast();

  std::cout<<"Total entry is "<<nentries<<std::endl;
  //int nBins=12;
  int nBins=100;
  double NBratio=((double)nentries)/nBins;
  double me=0.000511;
  double mmu=0.105658;
  double mpi=0.13957;
  double psilow=3.0;
  double psiup=3.2;
  //double psiplow=3.67;
  //double psipup=3.70;
  double psiplow=3.665;
  double psipup=3.705;
  double factorstart=0.995;
  // for factor fit
  TF1 *facfit = new TF1("facfit",line2,3.0,3.2,2);
  char fname[100];
  sprintf(fname,"%s/plot_pipill.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");

  // try to use roofit
  RooRealVar x("x","energy",3.097,3.0,3.2,"GeV");
  RooRealVar mean("mean","mean of gaussian",3.0,3.8);
  RooRealVar mean2("mean2","mean of gaussian",3.0,3.8);
  RooRealVar sigma("sigma","width of gaussian",0.0017,0.001,0.002);
  RooRealVar sigma2("sigma2","width of gaussian",0.0025,0.002,0.004);
  RooRealVar brewid("brewid","width of breit wigner",0.0023,0.0010,0.05);
  RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
  RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
  
  RooRealVar co1("co1","coefficient #1",0,-300.,300.);
  RooRealVar co2("co2","coefficient #2",-1,-500.,500.);
  RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
  
  RooRealVar alpha("alpha","#alpha",1.32,-5,5);
  RooRealVar nnn("nnn","n",5,1,200);
  RooCBShape cbshape("cbshape","crystal ball",x,mean,sigma,alpha,nnn);
  
  RooBreitWigner brewig("brewig","brewig",x,mean,brewid);
  
  RooRealVar signal("signal"," ",500,10,1000000);//event number
  RooRealVar signal2("signal2"," ",100,1,1000000);//event number
  RooRealVar background("background"," ",10,0,1000000);
  //RooAddPdf sum("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
  
  RooRealVar a0("a0","coefficient #0",100,-100000,100000);
  RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
  RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
 
  // Roofit part
  RooAddPdf *sum;
  RooDataHist *data_e;
  RooDataHist *data_mu;
  RooDataHist *data_pi;
  RooDataSet *dataset;
  RooPlot *xframe;  
  
  double par1[6];// for e
  double parerr1[6];
  double par2[6];// for mu
  double parerr2[6];
  double par3[6];// for pi
  double parerr3[6];
  
  int Npar;
  ofstream ofpar;
  ofpar.open("parpipill.txt",std::ios::app);
  ofpar<<"fastpipi algrithm: will give factors for e,mu,pi"<<std::endl;
  ofstream ofpardetail;
  ofpardetail.open("detail.txt",std::ios::app);
  ofstream purepar;
  purepar.open("par");
 
  TCanvas *c1=new TCanvas("","",800,600);
  //TH1D *h1   = new TH1D("h1","2 electron invariant mass",nBins,psilow,psiup);
  //TH1D *h2   = new TH1D("h2","2 muon invariant mass",nBins,psilow,psiup);
  TH1D *h3   = new TH1D("h3","total invariant mass(e)",nBins,psiplow,psipup);
  TH1D *h4   = new TH1D("h4","total invariant mass(mu)",nBins,psiplow,psipup);
  TH1D *h5   = new TH1D("h5","total invariant mass",nBins,psiplow,psipup);
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

  // for initial spectrum
  Long64_t nbytes = 0, nb = 0;

  int Npart=10;
  // cut phi
  //double phicut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
  //double start=0.0;
  //double stop =2*TMath::Pi();
  //for(int i=0;i<Npart+1;i++){
  //  phicut[i] = (stop-start)/Npart*i+start;
  //}
  //
  // cut cos theta
  //double costhecut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
  //double start=-1;
  //double stop =1;
  //for(int i=0;i<Npart+1;i++){
  //  costhecut[i] = (stop-start)/Npart*i+start;
  //}
  // cut p
  double pcut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
  double start=0;
  double stop =0.5;
  for(int i=0;i<Npart+1;i++){
    pcut[i] = (stop-start)/Npart*i+start;
  }

  //double m0=peakvalue;
  //double sigma_m=0.0017;//0.0024 for phi,
  //double width = 10.*sigma_m;
  double mparticle=mpi;
  //std::vector<std::pair<int,int> > partmap;
  //std::vector<std::pair<int,double> > facmap;
  std::vector<int> partmap;
  std::vector<std::pair<int,double> > facmap;
  
  char name[100];
  // ~~~~~~~~ draw nxn histogram, m distribution in different range
  TH1D *hmphi[Npart];
  for (int partj=0;partj<Npart;partj++){
    sprintf(name,"mass_part%d",partj);
    hmphi[partj] = new TH1D(name,name,100,psiplow,psipup);
  }
  
  // loop the data
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

     double mlepton;
     double massjpsi;
     double totpx,totpy,totpz,tote;
     double lee[2],pie[2];
     // total invariant mass
     if(cos(angle4)>0.90) continue; // cut bhabha 
     if (decay_ee == 0)    mlepton = mmu;
     else if (decay_ee==1) mlepton = me;
     else {
       std::cout<<"can not identify lepton. "<<std::endl;
       continue;;
     } 
     totpx=(lepx4[0]+lepx4[1])+(pipx4[0]+pipx4[1]);
     totpy=(lepy4[0]+lepy4[1])+(pipy4[0]+pipy4[1]);
     totpz=(lepz4[0]+lepz4[1])+(pipz4[0]+pipz4[1]);
     lee[0]=CalEne(mlepton,lepx4[0],lepy4[0],lepz4[0]);
     lee[1]=CalEne(mlepton,lepx4[1],lepy4[1],lepz4[1]);
     massjpsi = CalInvMass(mlepton,lepx4[0],lepy4[0],lepz4[0],mlepton,lepx4[1],lepy4[1],lepz4[1]);
     pie[0]=CalEne(mpi,pipx4[0],pipy4[0],pipz4[0]);
     pie[1]=CalEne(mpi,pipx4[1],pipy4[1],pipz4[1]);
     tote=lee[0]+lee[1]+pie[0]+pie[1];
     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
     mass = mass-massjpsi+3.096916;
     // if (Cut(ientry) < 0) continue;
     int parti, partj;
     if ( mass>psiplow && mass<psipup ) {
       p1 = CalMom(pipx4[0],pipy4[0],pipz4[0]);
       p2 = CalMom(pipx4[1],pipy4[1],pipz4[1]);
       costheta1 = pipz4[0]/p1;
       costheta2 = pipz4[1]/p2;
       if (pipy4[0]>0) phi1 = acos(pipx4[0]/sqrt(pipx4[0]*pipx4[0]+pipy4[0]*pipy4[0]));
       else phi1=2*TMath::Pi()-acos(pipx4[0]/sqrt(pipx4[0]*pipx4[0]+pipy4[0]*pipy4[0]));
       if (pipy4[1]>0) phi2 = acos(pipx4[1]/sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]));
       else phi2=2*TMath::Pi()-acos(pipx4[1]/sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]));
       vars->Fill();
     }
     //parti = (int)(phi1/stop*Npart);
     //partj = (int)((costheta2-start)/(stop-start)*Npart);
     partj = (int)((p2-start)/(stop-start)*Npart);
     if ( partj>=Npart || partj<0 ) continue;
     //if ( costheta2>costhecut[partj] && costheta2<costhecut[partj+1] )
     if ( p2>pcut[partj] && p2<pcut[partj+1] )
     {
       hmphi[partj]->Fill(mass);
     }
  }
  vars->Write();
  for (int partj=0;partj<Npart;partj++){
    hmphi[partj]->Write();
    std::cout<<"processed part "<<partj<<std::endl;
    if (hmphi[partj]->GetEntries() > 100
      &&hmphi[partj]->GetMaximumBin()> 30 //check the the peak position,
      &&hmphi[partj]->GetMaximumBin()< 70 //make sure there is a peak in the region
      ){
      partmap.push_back(partj);
    }
  }
  // ~~~~~~~~ draw end

  // try to correct the spectrum
  // iniialize the fit function
  // double factor1,factor1err;// for e+,e-
  // double factor2,factor2err;// for mu+,mu-
  double factor3,factor3err;// for pi
  double factor4,factor4err;// for pi

  // psi 3.097
  const int pointNo=20;
  double factor=factorstart;
  double factorstep=(1.-factor)*2/pointNo;
  double peakvalue=3.096916;
  int fittimes=0;
  double factors[pointNo];
  double factorserr[pointNo];
  double deltapeaks[pointNo];
  double deltapeakserr[pointNo];
  double errsum=0;
  std::string fitepsname;
  std::string fiteps_start;
  std::string fiteps_stop;
  std::string tmpstr;
  TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
     
  char tmpchr[100];


//~~~~~~~~~~pion part start~~~~~~~~

  std::cout<<"pion part start"<<std::endl;
  peakvalue=3.686109;
  // iniialize the fit function
  // psi(2S) 3.686
  
  //TFolder *folder;//(name,name);
  for (int loopj=0;loopj<partmap.size();loopj++){
    int partj = partmap.at(loopj);
    //double factori=1.00090;
    double factori=1.0000;
    factor=factorstart;
    //folder = new TFolder(name,name);
    
	sprintf(name,"part%d",partj);
	//f->Cd(name);

    fittimes=0;
 
    // for saving the fit result
    std::string fitepsname2=outputdir+"/fitpi.eps";
    std::string fiteps2_start=fitepsname2+"[";
    std::string fiteps2_stop =fitepsname2+"]";
    //c1->Print(fiteps2_start.c_str());
    //
    //x.setVal(3.686);
    x.setRange(psiplow,psipup);
    mean.setRange(psiplow,psipup);
    mean2.setRange(3.684,3.688);
    mean.setVal(3.686);
    mean2.setVal(3.686);
    //sigma.setRange(0.002,0.004);
    sigma.setRange(0.001,0.0025);
    sigma.setVal(0.0015);
    //sigma2.setRange(0.003,0.08);
    //sigma2.setVal(0.004);
    
    double factors2[pointNo];
    double factorserr2[pointNo];
    double deltapeaks2[pointNo];
    double deltapeakserr2[pointNo];
    for (int i=0;i<pointNo;i++){
      //delete xframe;
      xframe = x.frame(Title("fit pi"));
      //h5->Reset();
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
 
        //if(ngam>0) continue;
        double mlepton;
        double massjpsi;
        double totpx,totpy,totpz,tote;
        double lee[2],pie[2];
        // total invariant mass
        if(cos(angle4)>0.90) continue; // cut bhabha 
        if (decay_ee == 0)       mlepton = mmu;
        else if (decay_ee==1) mlepton = me;
        else {
          std::cout<<"can not identify lepton. "<<std::endl;
          return false;
        } 
        //if (pipy4[1]>0) phi2 = acos(pipx4[1]/sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]));
        //else phi2=2*TMath::Pi()-acos(pipx4[1]/sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]));
        p2 = CalMom(pipx4[1],pipy4[1],pipz4[1]);
        //costheta2 = pipz4[1]/p2;
        //if (!(costheta2>costhecut[partj] && costheta2<costhecut[partj+1])) continue;
        if (!(p2>pcut[partj] && p2<pcut[partj+1])) continue;
        totpx=(lepx4[0]+lepx4[1])+factori*pipx4[0]+factor*pipx4[1];
        totpy=(lepy4[0]+lepy4[1])+factori*pipy4[0]+factor*pipy4[1];
        totpz=(lepz4[0]+lepz4[1])+factori*pipz4[0]+factor*pipz4[1];
        lee[0]=CalEne(mlepton,lepx4[0],lepy4[0],lepz4[0]);
        lee[1]=CalEne(mlepton,lepx4[1],lepy4[1],lepz4[1]);
        massjpsi = CalInvMass(mlepton,lepx4[0],lepy4[0],lepz4[0],mlepton,lepx4[1],lepy4[1],lepz4[1]);
        pie[0]=CalEne(mpi,pipx4[0],pipy4[0],pipz4[0],factori);
        pie[1]=CalEne(mpi,pipx4[1],pipy4[1],pipz4[1],factor);
        tote=lee[0]+lee[1]+pie[0]+pie[1];
        mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
        mass = mass-massjpsi+3.096916;
        //h5->Fill(mass);

        if (mass>psiplow && mass<psipup){
          dataraw->Fill();
        }
        // if (Cut(ientry) < 0) continue;
      }
      //dataraw->Write();
      
      sprintf(tmpchr,"data_pi_%2d",fittimes);
      //data_pi = new RooDataHist(tmpchr,"data_pi",x,h4);
      //dataset = new RooDataSet(tmpchr,"data",dataraw,x);
      dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
      mean.setVal(3.686+0.44*(factor-1.0));
      signal.setError(10);
      signal2.setError(10);
      background.setError(10);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
      Npar=8;
      sum->fitTo(*dataset,Range(psiplow,psipup));
      //data_pi->plotOn(xframe);
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
      sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      pt->Draw();
      c1->Update();
      sprintf(tmpchr,"mass_pi_%2d",fittimes);
      xframe->SetName(tmpchr);
      //xframe->Write();
	  //folder->Add(xframe);
      // save pars
      factors2[i]=factor;
      factorserr2[i]=0;
      deltapeaks2[i] = mean.getVal() - peakvalue;
      deltapeakserr2[i] = mean.getError();
      if (deltapeakserr2[i]<1e-5) deltapeakserr2[i]=5e-4;
      //c1->Print(fitepsname2.c_str());
      sprintf(name,"part%d_fitFor%dTimes",partj,fittimes);
      c1->SetName(name);
      c1->Write();
	  //folder->Add(c1);
	  //folder->Write();
	  //delete folder;
      delete sum;
      delete dataset;
      delete xframe;
 
      fittimes++;
      factor += factorstep;
    }
    //c1->Print(fiteps2_stop.c_str());
 
    c1->Clear();
    TGraphErrors *graph3 = new TGraphErrors(pointNo,factors2,deltapeaks2,factorserr2,deltapeakserr2);
    graph3->SetTitle("delta peak");
    graph3->Draw("AP");
    graph3->SetMarkerStyle(5);
    gStyle->SetOptFit(1111);
    facfit->SetParameters(1,0.4);
    facfit->SetParNames("factor","slope");
    graph3->Fit(facfit,"","",factors2[0],factors2[pointNo-1]);
    factor4=facfit->GetParameter(0);
 
    tmpstr=outputdir+"/factorpi.eps";
    //c1->Print(tmpstr.c_str());
	sprintf(name,"factors_pi_part%d",fittimes);
    graph3->SetName(name);
    graph3->Write();
 
    //delete xframe;
    xframe = x.frame(Title("fit pi"));
    dataraw->Reset();
    factor = factor4;
    //factor = 1;
    std::cout<<"factor is "<<factor<<std::endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Long64_t ientry = LoadTree(jentry);
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       //if(ngam>0) continue;
       double massjpsi;
       double mlepton;
       //double p1,p2;
       //double theta1,theta2;
       double totpx,totpy,totpz,tote;
       double lee[2],pie[2];
       // total invariant mass
       if (cos(angle4)>0.90) continue; // cut bhabha 
       if (decay_ee == 0)    mlepton = mmu;
       else if(decay_ee==1)  mlepton = me;
       else {
         std::cout<<"can not identify lepton. "<<std::endl;
         return false;
       }
       //if (pipy4[1]>0) phi2 = acos(pipx4[1]/sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]));
       //else phi2=2*TMath::Pi()-acos(pipx4[1]/sqrt(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]));
       //if (!(phi2>phicut[partj] && phi2<phicut[partj+1])) continue;
       p2 = CalMom(pipx4[1],pipy4[1],pipz4[1]);
       //costheta2 = pipz4[1]/p2;
       //if (!(costheta2>costhecut[partj] && costheta2<costhecut[partj+1])) continue;
       if (!(p2>pcut[partj] && p2<pcut[partj+1])) continue;
       totpx=(lepx4[0]+lepx4[1])+factori*pipx4[0]+factor*pipx4[1];
       totpy=(lepy4[0]+lepy4[1])+factori*pipy4[0]+factor*pipy4[1];
       totpz=(lepz4[0]+lepz4[1])+factori*pipz4[0]+factor*pipz4[1];
       lee[0]= CalEne(mlepton,lepx4[0],lepy4[0],lepz4[0]);
       lee[1]= CalEne(mlepton,lepx4[1],lepy4[1],lepz4[1]);
       massjpsi = CalInvMass(mlepton,lepx4[0],lepy4[0],lepz4[0],mlepton,lepx4[1],lepy4[1],lepz4[1]);
 
       pie[0]= CalEne(mpi,pipx4[0],pipy4[0],pipz4[0],factori);
       pie[1]= CalEne(mpi,pipx4[1],pipy4[1],pipz4[1],factor);
 
       tote=lee[0]+lee[1]+pie[0]+pie[1];
       mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
       //h5->Fill(mass-massjpsi+3.096916);
       mass = mass-massjpsi+3.096916;
       if (mass>psiplow && mass<psipup)
         dataraw->Fill();
    }
    //h2p->Write();
    //thedis2->Write();
    dataraw->Write();
    //vars->Write();
 
    sprintf(tmpchr,"data_pi_%2d",fittimes);
    //data_pi = new RooDataHist(tmpchr,"data_pi",x,h4);
    dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
    mean.setVal(3.686+0.45*(factor-1.0));
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
    Npar=8;
    sum->fitTo(*dataset,Range(psiplow,psipup));
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
    sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
    pt->AddText(tmpchr);
    factor4err=TMath::Sqrt(TMath::Power(mean.getError()/facfit->GetParameter(1),2)
               +TMath::Power(facfit->GetParError(0),2));
    sprintf(tmpchr,"factor = %.6f #pm %.6f",factor,factor4err);
    pt->AddText(tmpchr);
	pt->Draw();
    c1->Update();
    tmpstr=outputdir+"/fitpi_best.eps";
    //c1->Print(tmpstr.c_str());
    xframe->SetName("mass_pi_best");
    //xframe->Write();
    sprintf(name,"part%d_correct",partj);
    c1->SetName(name);
    c1->Write();
    delete sum;
 
    // get factor and its error
    factor4=facfit->GetParameter(0);
    ofpar<<factor4<<"\t"<<factor4err<<std::endl;
    ofpar<<"\tChisqure "<<xframe->chiSquare()<<" "<<xframe->chiSquare(Npar)<<std::endl;
    ofpar<<signal.getVal()<<"\t"<<signal.getError()<<std::endl;
    purepar<<"\t"<<factor4<<"\t"<<factor4err<<std::endl;
    delete dataset;
    delete xframe;

	//folder.Write();
  
  }

//~~~~~~~~~~pion part end~~~~~~~~
  std::cout<<"dataset pointer:"<<dataset<<std::endl;
  std::cout<<"dataraw pointer:"<<dataraw<<std::endl;
  std::cout<<"xframe  pointer:"<<xframe <<std::endl;
  std::cout<<"legend  pointer:"<<legend <<std::endl;
  std::cout<<"TFile f pointer:"<<f      <<std::endl;
  std::cout<<"dataset pointer:"<<dataset<<std::endl;
  ofpar.close();
  ofpardetail.close();
  //delete h1;
  //delete h2;
  delete legend;
  delete dataraw;
  delete c1;
  return true;
}

void PIPILL::FitAndSave( TH1D *&hmass, RooPlot *&xframe, RooDataHist *&data,
       RooRealVar &x, RooRealVar &mean, RooRealVar &sigma, RooRealVar &co1,
       RooRealVar &co2, RooRealVar &signal, RooRealVar &background,
       RooAddPdf *&sum, RooGaussian &gaus, RooChebychev &bkg,
       int idx,double factor,double psiplow,double psipup,std::string suffix)
{
  delete xframe;
  xframe = x.frame(Title("fit pi"));
  char tmpchr[100];
  sprintf(tmpchr,"data_pi_%2d",idx);
  data = new RooDataHist(tmpchr,"data_pi",x,hmass);
  mean.setVal(3.686-0.001+0.45*(factor-1.0));
  sigma.setRange(0.002,0.004);
  sigma.setVal(0.0023);
  signal.setVal(60);
  background.setVal(10);
  co1.setVal(0);
  sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
  sum->fitTo(*data,Range(psiplow,psipup));
  data->plotOn(xframe);
  sum->plotOn(xframe);
  sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
  //sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
  sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
  xframe->Draw();
  sprintf(tmpchr,"mass_pi_e_%2d",idx);
  xframe->SetName(tmpchr);
  xframe->Write();

  return;
}
void PIPILL::ResetVars(  RooRealVar &mean, RooRealVar &sigma, RooRealVar &co1,
       RooRealVar &co2, RooRealVar &signal, RooRealVar &background,
       double factor, double meanv, double sigmav,
       double signalv, double backgroundv, double co1v,
       double co2v)
{
  mean.setVal(meanv-0.001+0.45*(factor-1.0));
  sigma.setRange(0.002,0.004);
  sigma.setVal(sigmav);
  //sigma2.setRange(0.005,0.03);
  //sigma2.setVal(0.014);
  signal.setVal(signalv);
  background.setVal(backgroundv);
  co1.setVal(co1v);

  return;
}

gepep_fastpipill::gepep_fastpipill(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_Rvalue_pipill_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_Rvalue_pipill_1.root");
      }
      f->GetObject("gepep_fastpipill",tree);

   }
   Init(tree);
}

gepep_fastpipill::~gepep_fastpipill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fastpipill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fastpipill::LoadTree(Long64_t entry)
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

void gepep_fastpipill::Init(TTree *tree)
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
   fChain->SetBranchAddress("mcgpx", mcgpx, &b_mcgpx);
   fChain->SetBranchAddress("mcgpy", mcgpy, &b_mcgpy);
   fChain->SetBranchAddress("mcgpz", mcgpz, &b_mcgpz);
   fChain->SetBranchAddress("mcge", mcge, &b_mcge);
   fChain->SetBranchAddress("mcgid", mcgid, &b_mcgid);
   fChain->SetBranchAddress("nGammatch", &nGammatch, &b_nGammatch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("ntot", &ntot, &b_ntot);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ngam", &ngam, &b_ngam);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("netap2", &netap2, &b_netap2);
   fChain->SetBranchAddress("delang", delang, &b_delang);
   fChain->SetBranchAddress("delphi", delphi, &b_delphi);
   fChain->SetBranchAddress("delthe", delthe, &b_delthe);
   fChain->SetBranchAddress("npart", npart, &b_npart);
   fChain->SetBranchAddress("nemchits", nemchits, &b_nemchits);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("dx", dx, &b_dx);
   fChain->SetBranchAddress("dy", dy, &b_dy);
   fChain->SetBranchAddress("dz", dz, &b_dz);
   fChain->SetBranchAddress("dtheta", dtheta, &b_dtheta);
   fChain->SetBranchAddress("dphi", dphi, &b_dphi);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("dE", dE, &b_dE);
   fChain->SetBranchAddress("eSeed", eSeed, &b_eSeed);
   fChain->SetBranchAddress("nSeed", nSeed, &b_nSeed);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("e5x5", e5x5, &b_e5x5);
   fChain->SetBranchAddress("secondMoment", secondMoment, &b_secondMoment);
   fChain->SetBranchAddress("latMoment", latMoment, &b_latMoment);
   fChain->SetBranchAddress("a20Moment", a20Moment, &b_a20Moment);
   fChain->SetBranchAddress("a42Moment", a42Moment, &b_a42Moment);
   fChain->SetBranchAddress("getTime", getTime, &b_getTime);
   fChain->SetBranchAddress("getEAll", getEAll, &b_getEAll);
   fChain->SetBranchAddress("mpi0", mpi0, &b_mpi0);
   fChain->SetBranchAddress("chisq1cpi0", chisq1cpi0, &b_chisq1cpi0);
   fChain->SetBranchAddress("ig1pi0", ig1pi0, &b_ig1pi0);
   fChain->SetBranchAddress("ig2pi0", ig2pi0, &b_ig2pi0);
   fChain->SetBranchAddress("chisq6c", &chisq6c, &b_chisq6c);
   fChain->SetBranchAddress("chisq4c", &chisq4c, &b_chisq4c);
   fChain->SetBranchAddress("chisq3c", &chisq3c, &b_chisq3c);
   fChain->SetBranchAddress("gpx4", gpx4, &b_gpx4);
   fChain->SetBranchAddress("gpy4", gpy4, &b_gpy4);
   fChain->SetBranchAddress("gpz4", gpz4, &b_gpz4);
   fChain->SetBranchAddress("ge4", ge4, &b_ge4);
   fChain->SetBranchAddress("lepx4", lepx4, &b_lepx4);
   fChain->SetBranchAddress("lepy4", lepy4, &b_lepy4);
   fChain->SetBranchAddress("lepz4", lepz4, &b_lepz4);
   fChain->SetBranchAddress("lee4", lee4, &b_lee4);
   fChain->SetBranchAddress("pipx4", pipx4, &b_pipx4);
   fChain->SetBranchAddress("pipy4", pipy4, &b_pipy4);
   fChain->SetBranchAddress("pipz4", pipz4, &b_pipz4);
   fChain->SetBranchAddress("pie4", pie4, &b_pie4);
   fChain->SetBranchAddress("ggm4", &ggm4, &b_ggm4);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   fChain->SetBranchAddress("kpi01m4", &kpi01m4, &b_kpi01m4);
   fChain->SetBranchAddress("kpi02m4", &kpi02m4, &b_kpi02m4);
   fChain->SetBranchAddress("kkpi0m4", &kkpi0m4, &b_kkpi0m4);
   fChain->SetBranchAddress("kkpipim4", &kkpipim4, &b_kkpipim4);
   fChain->SetBranchAddress("kkpipipi0m4", &kkpipipi0m4, &b_kkpipipi0m4);
   fChain->SetBranchAddress("Recoilmass", &Recoilmass, &b_Recoilmass);
   fChain->SetBranchAddress("llm4", &llm4, &b_llm4);
   fChain->SetBranchAddress("pipim4", &pipim4, &b_pipim4);
   fChain->SetBranchAddress("llpipi1m4", &llpipi1m4, &b_llpipi1m4);
   fChain->SetBranchAddress("llpipi2m4", &llpipi2m4, &b_llpipi2m4);
   fChain->SetBranchAddress("llpipi3m4", &llpipi3m4, &b_llpipi3m4);
   fChain->SetBranchAddress("llpipi4m4", &llpipi4m4, &b_llpipi4m4);
   fChain->SetBranchAddress("gpipim4", &gpipim4, &b_gpipim4);
   fChain->SetBranchAddress("decay_ee", &decay_ee, &b_decay_ee);
   fChain->SetBranchAddress("hepp", &hepp, &b_hepp);
   fChain->SetBranchAddress("hepm", &hepm, &b_hepm);
   fChain->SetBranchAddress("emcp", &emcp, &b_emcp);
   fChain->SetBranchAddress("emcm", &emcm, &b_emcm);
   fChain->SetBranchAddress("mucp", &mucp, &b_mucp);
   fChain->SetBranchAddress("mucm", &mucm, &b_mucm);
   fChain->SetBranchAddress("zcpm4", &zcpm4, &b_zcpm4);
   fChain->SetBranchAddress("zcmm4", &zcmm4, &b_zcmm4);
   fChain->SetBranchAddress("ecms", &ecms, &b_ecms);
   fChain->SetBranchAddress("angle1", &angle1, &b_angle1);
   fChain->SetBranchAddress("angle2", &angle2, &b_angle2);
   fChain->SetBranchAddress("angle3", &angle3, &b_angle3);
   fChain->SetBranchAddress("angle4", &angle4, &b_angle4);
   fChain->SetBranchAddress("psipm4", &psipm4, &b_psipm4);
   Notify();
}

Bool_t gepep_fastpipill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fastpipill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fastpipill::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
