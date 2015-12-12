#define gepep_kk_cxx
#include "gepep_kk.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include <fstream>

#include <stdio.h>
#include "unistd.h"
#include <sys/types.h>

extern std::string outputdir;
using namespace std;

namespace GEPEP_KK{
  void FitSpecPhi(std::vector<Event> &evts, double beame, char* namesfx);
  void FitSpectrumPhi(TTree*& dataraw, double beame, char* namesfx, double &peak, double &peakerror);
  void FitSpecD(std::vector<Event> &evts, double beame, char* namesfx);
  void FitSpectrumD(TTree*& dataraw, double beame, char* namesfx, double &peak, double &peakerror);
  void FitSpectrumD1(TH1D*& dataraw, double beame, char* namesfx, double &peak, double &peakerror, int id);
}

class AlgElements{
private:
  vector<Event> *evts;
  TFile* file;
  ofstream* file2;
  TTree* result;
  char* flag;
  TF1*  function;

  double *phipeak;
  double *phipeake;
  double *Dpeak;
  double *Dpeake;
public:
  int id; 

public:
  AlgElements():evts(0),result(0),flag(0){}
  ~AlgElements(){}
  void setData(vector<Event> *vevt){ evts   = vevt;  }
  void setOutputFile(TFile* fil)   { file   = fil;   }
  void setOutputFile2(ofstream* fil){ file2 = fil;   }
  void setOutputTree(TTree* res)   { result = res;   }
  void setFlag(char* fl)           { flag   = fl;    }
  void setFunction(TF1* func)      { function = func;}
  void connectResult(double *mphi,double *mphie, double *mD, double *mDe)  
  { 
    phipeak = mphi;   
    phipeake = mphie;   
    Dpeak = mD;   
    Dpeake = mDe;   
  }

  vector<Event>* getData(){ return evts;}
  TFile* getOutputFile()  { return file;}
  ofstream* getOutputFile2()  { return file2;}
  TTree* getOutputTree()  { return result;}
  char* getFlag()         { return flag;}
  TF1* getFunction()      { return function;}
  double *getmphi() {return phipeak;}
  double *getmphie() {return phipeake;}
  double *getmD() {return Dpeak;}
  double *getmDe() {return Dpeake;}
};

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
   fChain->GetEntry(1);
   double beamene = GetEnergy(run);

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=100;
   double factorstart=0.99;
   double mk=0.493677;
   // D0 -> K K
   double D0low=1.82;
   double D0up=1.90;
   double peakvalue=1.86484;// mD0
   // phi -> K K
     double philow=1.005;
     double phiup=1.05;
     double phipeak=1.019461;// mphi
  

   char fname[1000];
// ofstream ofpar;
// sprintf(fname,"%s/parskk.txt",outputdir.c_str());
// ofpar.open(fname,std::ios::app);
// ofstream purepar;
// sprintf(fname,"parkk");
// purepar.open(fname,std::ios::app);
// ofstream purepar2;
// sprintf(fname,"parkk_phi");
// purepar2.open(fname,std::ios::app);
   sprintf(fname,"%s/plot_kk.root",outputdir.c_str());
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
   vars->Branch("p1",&p1,"p1/D");
   vars->Branch("p2",&p2,"p2/D");
   vars->Branch("mass",&mass,"mass/D");


   const int Npart =1 ;
   double m0=peakvalue;
   //double peakest[Npart];
   double sigma_m=0.0024;//0.0024 for phi,
   double width = 10.*sigma_m;
   double mparticle=0.493677;
// std::vector<int> partmap;
// std::vector<std::pair<int,double> > facmap;

   char name[100];
   TH1D* hpKD0 = new TH1D("hpKD0","p_{K} (D0)",200,0,2);
   TH1D* hpKphi = new TH1D("hpKphi","p_{K} (phi)",200,0,2);
 
   // loop data
   Event evt(mk);
   std::vector<Event> evts[Npart];
   std::vector<Event> evtsphi[Npart];
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
	evt.SetVal(kappx,kappy,kappz,kampx,kampy,kampz);
        mass = evt.InvMass();
        p1 = evt.GetP1();
        p2 = evt.GetP2();
        costheta1 = evt.GetCostheta1();
	costheta2 = evt.GetCostheta2();
        phi1 = evt.GetPhi1();
        phi2 = evt.GetPhi2();
		
	//if (fabs(costheta1)>0.7) continue;
	//if (fabs(costheta2)>0.7) continue;
      //if (mass>D0low && mass<D0up){
      //  vars->Fill();
      //}
      //if (mass>philow && mass<phiup){
      //  vars->Fill();
      //}
        //if (costheta1 > 0.8 && costheta2 <-0.8) continue;
	//if (p1<0.5 || p1>1.3) continue;
	//if (p1<0.4 || p1>1.3) continue;
	//if (p2<0.4 || p2>1.3) continue;

        //if ( partj>=Npart || partj<0 ) continue;
        for (int partj=0;partj<Npart;partj++){
          //if (p1<pcut[partj] || p1>pcut[partj+1]) continue;
          //if (p2<pcut[partj] || p2>pcut[partj+1]) continue;
          
          //if ( costheta1>costhecut[partj] && costheta1<costhecut[partj+1] )
	  if (mass>D0low-0.02 && mass<D0up+0.02){
            //hmD0[partj]->Fill(mass);
  	    //evts[partj].push_back(evt);
  	    evtsphi[partj].push_back(evt);
            if (mass > peakvalue-0.007 && mass < peakvalue+0.007) {
	      hpKD0->Fill(p1);
	      hpKD0->Fill(p2);
	    }
	    break;
	  }
          if (mass>philow-0.02 && mass<phiup+0.02){
            //hmphi[partj]->Fill(mass);
	    evtsphi[partj].push_back(evt);
            if (mass > phipeak-0.0024 && mass < phipeak+0.0024) {
	      hpKphi->Fill(p1);
	      hpKphi->Fill(p2);
	    }
            break;
	  }
	}
   
   }
   cout<<" D0: ave pK = "<<hpKD0->GetMean()<<" +/- "<<hpKD0->GetRMS()<<endl;
   cout<<"phi: ave pK = "<<hpKphi->GetMean()<<" +/- "<<hpKphi->GetRMS()<<endl;
   std::cout<<"part No: "<< Npart<<std::endl;
   //vars->Write();

   cout<<"Prepare going into parallel process."<<endl;
   GEPEP_KK::FitSpecPhi(evtsphi[0],0,0);

   // ~~~~~~~~ draw end
   f->Close();
   return;
}

#include "TPaveText.h"
#include "TGraphErrors.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"

#include "RooRelativisticBW.h"
#include "RooRBW_evtgen.h"
#include "RooFFTConvPdf.h"

using namespace RooFit;
pthread_rwlock_t rwl = PTHREAD_RWLOCK_INITIALIZER;
pthread_mutex_t lock_a = PTHREAD_MUTEX_INITIALIZER;

double ffun(double *x, double *par)
{
  return 1./(par[0]*pow(x[0],par[1]))+par[2];
}

void* mythread(void* par)
{
  AlgElements* curpar = (AlgElements*)par;
  
  vector<Event> &evts = *(curpar->getData());
  TF1 &facfun = *(curpar->getFunction());
  TTree *result = curpar->getOutputTree();


  double philow=1.005;
  double phiup=1.05;
  double mphi=1.019461;// mphi
  double D0low=1.82;
  double D0up=1.90;
  double mD=1.86484;// mD0
  
  double para = facfun.GetParameter(0);
  double parb = facfun.GetParameter(1);
  double parc = facfun.GetParameter(2);
//if (pthread_rwlock_rdlock(&rwl) != 0) {
//  perror("pthread_rwlock_rdlock");
//  exit(-1);
//}
  if (pthread_mutex_lock(&lock_a) != 0) {
    perror("lock_reader:pthread_mutex_lock");
    exit(-1);
  }
  cout<<"Current parameters: "<<para << "\t" << parb<<"\t"<<parc<<endl;
  if (pthread_mutex_unlock(&lock_a) != 0) {
    perror("lock_reader:pthread_mutex_unlock");
    exit(-1);
  }
//if (pthread_rwlock_unlock(&rwl) != 0) {
//  perror("pthread_rwlock_unlock");
//  exit(-1);
//}
  //facfun.SetParameters(para,parb,parc);
  double mass;
  char tname[100]; char lname[100];
//TTree *dataraw = new TTree("dataraw","dataraw");
//sprintf(tname,"x_%2d",curpar->id);
//sprintf(lname,"x_%2d/D",curpar->id);
//dataraw->Branch(tname,&mass,lname);
//TTree *datarawf = new TTree("datarawf","dataraw");
//sprintf(tname,"xx_%2d",curpar->id);
//sprintf(lname,"xx_%2d/D",curpar->id);
//datarawf->Branch(tname,&mass,lname);

//dataraw->Reset();
//datarawf->Reset();
  
  sprintf(tname,"hmD_%2d",curpar->id);
  TH1D *hmD    = new TH1D(tname,"2 kaon invariant mass",100,D0low,D0up);
  sprintf(tname,"hmphi_%2d",curpar->id);
  TH1D *hmphi  = new TH1D(tname,"2 kaon invariant mass",100,philow,phiup);
  hmD->Reset();
  hmphi->Reset();
  // D to K K
  for (Long64_t jentry=0; jentry<evts.size();jentry++) {
    double p1 = evts.at(jentry).GetP1();
    double p2 = evts.at(jentry).GetP2();
    double f1 = facfun.Eval(p1);
    double f2 = facfun.Eval(p2);
    mass = evts.at(jentry).InvMass(f1,f2);
    if (mass>D0low && mass<D0up){
      //dataraw->Fill();
      hmD->Fill(mass);
    }
    if (mass>philow && mass<phiup){
      //datarawf->Fill();
      hmphi->Fill(mass);
    }
  }
  //fitmphi=0; fitmphie=0; fitmD=0; fitmDe=0;
  char parset[100];
  sprintf(parset,"a%4f_b%.3f_c%.5f",para,parb,parc);
  //char* parset = curpar.getFlag();
  GEPEP_KK::FitSpectrumD1(hmD,mD,parset,*curpar->getmD(), *curpar->getmDe(), curpar->id);
  //GEPEP_KK::FitSpectrumPhi(hmphi,mphi,parset,*curpar->getmphi(), *curpar->getmphie());
//curpar->getOutputFile()->cd();
//result->Fill();
  delete hmD;
  delete hmphi;
  if (pthread_rwlock_rdlock(&rwl) != 0) {
    perror("pthread_rwlock_rdlock");
    exit(-1);
  }
  (*curpar->getOutputFile2())<< para<<"\t"<<parb<<"\t"<<parc;
  (*curpar->getOutputFile2())<< "\t"<< *curpar->getmD() <<"\t"<<*curpar->getmDe()<<"\t"<<*curpar->getmphi()<<"\t"<<*curpar->getmphie()<<endl;
  if (pthread_rwlock_unlock(&rwl) != 0) {
    perror("pthread_rwlock_unlock");
    exit(-1);
  }
  return NULL;
}

void GEPEP_KK::FitSpecPhi(std::vector<Event> &evts, double beame, char* namesfx)
{
   double para[3], parb[3], parc[3], fitmphi[3], fitmphie[3], fitmD[3], fitmDe[3];
   TFile* files[3];
   ofstream ofile[3];
   TTree* result[3];
   AlgElements parset[3];
   for (int i=0;i<3;i++){
     char name[100];
     sprintf(name,"output_%02d.root",i);
     files[i] = new TFile(name,"recreate");
     result[i] = new TTree("result","result");
     result[i]->Branch("para",&para[i],"para/D");
     result[i]->Branch("parb",&parb[i],"parb/D");
     result[i]->Branch("parc",&parc[i],"parc/D");
     result[i]->Branch("mphi",&fitmphi[i],"mphi/D");
     result[i]->Branch("mphie",&fitmphie[i],"mphie/D");
     result[i]->Branch("mD",&fitmD[i],"mD/D");
     result[i]->Branch("mDe",&fitmDe[i],"mDe/D");

     sprintf(name,"output_%02d.txt",i);
     ofile[i] = ofstream(name);
     
     parset[i].id = i;
     parset[i].setData(&evts);
     //parset[i].setOutputFile(files[i]);
     parset[i].setOutputFile2(&ofile[i]);
     parset[i].setOutputTree(result[i]);
     parset[i].setFlag(namesfx);
   }
   

 //TF1 facfun("facfun",ffun,0,3,3);
 //double pars[3] = {1.9e3, 2.69, -5.7e-4};
 //facfun.SetParameters(pars);

   for (double par1=1000; par1<=1000; par1+=100){
     for (double par2=1.0; par2<=5; par2+=0.1){
       for (double par3=-2e-3; par3<2e-3;){
         TF1* f1[3];
         pthread_t tid[3];
         void *vp[3];
	 char name[100];
	 for (int i=0;i<3;i++){
           para[i] = par1;
           parb[i] = par2;
           parc[i] = par3;
	   sprintf(name,"fun%02d",i);
           f1[i] = new TF1(name,ffun,0,3,3);
           f1[i]->SetParameters(para[i],parb[i],parc[i]);
           parset[i].setFunction(f1[i]);
           fitmphi[i]=0; fitmphie[i]=0; fitmD[i]=0; fitmDe[i]=0;
	   parset[i].connectResult(&fitmphi[i],&fitmphie[i],&fitmD[i],&fitmDe[i]);
           if (pthread_create(&tid[i], NULL, mythread, &parset[i]) != 0) {
             perror("pthread_create");
             exit(-1);
           }
           par3 += 1e-4;
	 }
        
         for (int i=0;i<3;i++)
	 {
	   if (pthread_join(tid[i], &vp[i]) != 0) { // wait 
             perror("pthread_join");
             exit(-1);
           }
	   delete f1[i];
         }

       }
     }
   }

   //result->Write();
   return;
}

void GEPEP_KK::FitSpectrumPhi(TTree*& dataraw, double beame, char* namesfx, double &peak, double &peakerror)
{
   double philow=1.005;
   double phiup=1.05;
   double phipeak=1.019461;// mphi
   TCanvas *c1 = new TCanvas("","",800,600);
   
   if (fabs(peak)>1e-1) phipeak = peak; 

   RooRealVar xx("xx","energy",phipeak,philow,phiup,"GeV");
   RooRealVar meanf("meanf","mean of breitwigner",phipeak,philow,phiup);
   RooRealVar sigmaf("sigmaf","width of breitwiger",0.00426);//phi version
   RooRealVar meang("meang","mean of breitwigner",0);
   RooRealVar sigmag("sigmag","width of gaussian",0.002,0.001,0.03);//phi version
   //RooRealVar sigma2f("sigma2f","width of gaussian",0.004,0.003,0.005);
   RooGaussian gausg("gausg","gauss(x,m,s)",xx,meang,sigmag);
  // RooGaussian gaus2f("gaus2f","gauss(x,m,s)",xx,meanf,sigma2f);
   
   //RooBreitWigner brewig("brewig","breitwigner",xx,meanf,sigmaf);
   //RooRelativisticBW brewig("brewig","breitwigner",xx,meanf,sigmaf);
   RooRBW_evtgen brewig("brewig","breitwigner",xx,meanf,sigmaf);
   //RooRBW_evtgen sig("brewig","breitwigner",xx,meanf,sigmaf);

   RooFFTConvPdf sig("sig","signal function",xx, brewig, gausg);
   
   RooRealVar co1f("co1f","coefficient #1",0,-1000,1000);
   RooRealVar co2f("co2f","coefficient #2",-1,-1000,1000);
   RooChebychev bkgf("bkgf","background",xx,RooArgList(co1f,co2f));
   
   RooRealVar signalf("signalf"," ",1200,10,1000000);//event number
   RooRealVar signal2f("signal2f"," ",1200,0,1000000);//event number
   RooRealVar backgroundf("backgroundf"," ",200,0,1000000);
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

    RooDataSet* dataset = new RooDataSet("datasetf","data",dataraw,xx);
    RooAddPdf *sum = new RooAddPdf("sum","sum",RooArgList(sig,bkgf),RooArgList(signalf,backgroundf));
    
    RooPlot* xframe = xx.frame();
    char tmpchr[100];
    sprintf(tmpchr,"data_kf_%s",namesfx);
    //sum = new RooAddPdf("sumf","sum",RooArgList(gausf,gaus2f,bkgf),RooArgList(signalf,signal2f,backgroundf));
    int Npar=6;
    //sum = new RooAddPdf("sum","sum",RooArgList(sig,bkgf),RooArgList(signalf,backgroundf));
    //meanf.setVal(phipeak+0.06*(factor-1.0));
    sum->fitTo(*dataset,Range(philow,phiup));
    //sum->fitTo(*dataset,Range(phipeak-0.002,phipeak+0.002));
    dataset->plotOn(xframe);
    sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
    //sum->plotOn(xframe,Components(gaus2f),LineStyle(2),LineColor(2));
    sum->plotOn(xframe,Components(bkgf),LineStyle(3),LineColor(3));
    sum->plotOn(xframe);
    xframe->Draw();
    TPaveText *
    pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(4000);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",meanf.getVal(),meanf.getError());
    pt->AddText(tmpchr);
    sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigmag.getVal(),sigmag.getError());
    pt->AddText(tmpchr);
    //sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2f.getVal(),sigma2f.getError());
    //pt->AddText(tmpchr);
    sprintf(tmpchr,"signal = %.2f #pm %.2f",signalf.getVal(),signalf.getError());
    pt->AddText(tmpchr);
    //sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2f.getVal(),signal2f.getError());
    //pt->AddText(tmpchr);
    sprintf(tmpchr,"backNo = %.2f #pm %.2f",backgroundf.getVal(),backgroundf.getError());
    pt->AddText(tmpchr);
    sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
    pt->AddText(tmpchr);
  //sprintf(name,"factor = %.6f",factor);
  //pt->AddText(name);
    pt->Draw();
   char name[100];
    sprintf(name,"phipart%s",namesfx);
    c1->SetName(name);
    c1->Write();
    //delete data_k;
    
    peak = meanf.getVal();
    peakerror = meanf.getError();

    delete dataset;
    delete xframe;
    delete sum;
    delete pt;
    delete c1;
    return;

}


void GEPEP_KK::FitSpecD(std::vector<Event> &evts, double beame, char* namesfx)
{

}

void GEPEP_KK::FitSpectrumD(TTree*& dataraw, double beame, char* namesfx, double &peak, double &peakerror)
{
   double D0low=1.82;
   double D0up=1.90;
   double peakvalue=1.86484;// mD0
   TCanvas *c1 = new TCanvas("","",800,600);
 
   if (fabs(peak)>1e-1) peakvalue = peak; 
   
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,D0low,D0up,"GeV");
   //RooRealVar x("x","energy",peakvalue,1.015,1.025,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,D0low,D0up);
   RooRealVar mean2("mean2","mean of gaussian2",peakvalue,D0low,D0up);
   //RooRealVar sigma("sigma","width of gaussian",0.0023,0.0015,0.0040);//phi version
   //RooRealVar sigma2("sigma2","width of gaussian",0.005,0.003,0.008);
   RooRealVar sigma("sigma","width of gaussian",0.007,0.003,0.0075);//D0 version
   RooRealVar sigma2("sigma2","width of gaussian",0.02,0.008,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   
   RooRealVar co1("co1","coefficient #1",0,-0.5,0.5);
   RooRealVar co2("co2","coefficient #2",0,-0.01,0.5);
   RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
   
 //RooRealVar a0("a0","coefficient #0",100,100,100000);
 //RooRealVar a1("a1","coefficient #1",-50,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooRealVar signal("signal"," ",1200,10,1000000);//event number
   RooRealVar signal2("signal2"," ",1200,0,1000000);//event number
   RooRealVar background("background"," ",200,0,1000000);
   
   RooDataSet *dataset;
   RooAddPdf *sum;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

    RooPlot* xframe = x.frame();
      dataset = new RooDataSet("dataset","data",dataraw,x);
      char tmpchr[100];
      sprintf(tmpchr,"data_k_%s",namesfx);
      //data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      int Npar=8;
      //mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(D0low,D0up));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
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
      char name[100];
      pt->Draw();
      sprintf(name,"partD%s",namesfx);
      c1->SetName(name);
      c1->Write();
      
      
      peak = mean.getVal();
      peakerror = mean.getError();
      delete dataset;
      delete sum;
      delete c1;
      delete pt;

      return;

}
void GEPEP_KK::FitSpectrumD1(TH1D*& dataraw, double beame, char* namesfx, double &peak, double &peakerror, int id)
{  
   double D0low=1.82;
   double D0up=1.90;
   double peakvalue=1.86484;// mD0
   char name[100];
   sprintf(name,"c_%d",id);
   TCanvas *c1 = new TCanvas(name,name,800,600);
   if (dataraw->GetEntries()<1) {
     cout<< "TOO FEW EVENTS!!!!!!!!!!!"<<endl;
     delete c1; return ;
   }
   if (fabs(peak)>1e-1) peakvalue = peak; 
   
   // try to use roofit
   sprintf(name,"x_%d",id);
   RooRealVar x(name,"energy",peakvalue,D0low,D0up,"GeV");
   sprintf(name,"mean_%d",id);
   RooRealVar mean(name,"mean of gaussian",peakvalue,D0low,D0up);
   //RooRealVar mean2("mean2","mean of gaussian2",peakvalue,D0low,D0up);
   //RooRealVar sigma("sigma","width of gaussian",0.0023,0.0015,0.0040);//phi version
   //RooRealVar sigma2("sigma2","width of gaussian",0.005,0.003,0.008);
   sprintf(name,"sigma_%d",id);
   RooRealVar sigma(name,"width of gaussian",0.007,0.003,0.0075);//D0 version
   sprintf(name,"sigma2_%d",id);
   RooRealVar sigma2(name,"width of gaussian",0.02,0.008,0.05);
   sprintf(name,"gaus_%d",id);
   RooGaussian gaus(name,"gauss(x,m,s)",x,mean,sigma);
   sprintf(name,"gaus2_%d",id);
   RooGaussian gaus2(name,"gauss(x,m,s)",x,mean,sigma2);
   
   sprintf(name,"co1_%d",id);
   RooRealVar co1(name,"co_1efficient #1",0,-0.5,0.5);
   sprintf(name,"co2_%d",id);
   RooRealVar co2(name,"co_1efficient #2",0,-0.01,0.5);
   sprintf(name,"bkg_%d",id);
   RooChebychev bkg(name,"background",x,RooArgList(co1,co2));
   
 //RooRealVar a0("a0","coefficient #0",100,100,100000);
 //RooRealVar a1("a1","coefficient #1",-50,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   sprintf(name,"signal_%d",id);
   RooRealVar signal(name," ",1200,10,1000000);//event number
   sprintf(name,"signal2_%d",id);
   RooRealVar signal2(name," ",1200,0,1000000);//event number
   sprintf(name,"background_%d",id);
   RooRealVar background(name," ",200,0,1000000);
   
   //RooDataSet *dataset;
   RooAddPdf *sum;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   sprintf(name,"frame_%d",id);
    RooPlot* xframe = x.frame(Title(name));
   sprintf(name,"dataset_%d",id);
   RooDataHist* dataset = new RooDataHist(name,name,x,dataraw);
      char tmpchr[100];
      sprintf(tmpchr,"data_k_%s",namesfx);
      //data_k = new RooDataHist(tmpchr,"data_k",x,h1);
   sprintf(name,"sum_%d",id);
      sum = new RooAddPdf(name,"sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      int Npar=8;
      //mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(D0low,D0up));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
      sum->plotOn(xframe);
      xframe->Draw();
      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(4000);
      pt->SetTextAlign(12);
      pt->SetTextFont(42);
      pt->SetTextSize(0.035);
      sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
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
      //char name[100];
      pt->Draw();
      sprintf(name,"partD%s",namesfx);
      c1->SetName(name);
      c1->Write();
      
      
      peak = mean.getVal();
      peakerror = mean.getError();
      delete dataset;
      delete sum;
      delete c1;
      delete pt;

      return;
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
