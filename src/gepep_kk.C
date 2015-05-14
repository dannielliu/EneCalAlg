#define gepep_kk_cxx
#include "gepep_kk.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <fstream>
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
extern std::string outputdir;
using namespace RooFit;

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
     double phiup=1.035;
     double phipeak=1.019455;// mphi
   int pointNo=10;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   double factor=factorstart;
   double factorstep=(1.-factor)*2/pointNo;
   
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
   
   RooPlot *xframe;
   //RooDataHist *data_k;
   RooDataSet *dataset;
   RooAddPdf *sum;
   
   //RooDataSet *dataset = new RooDataSet("dataset","data",dataraw,x);
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   int Npar;
   char fname[1000];
   ofstream ofpar;
   sprintf(fname,"%s/pars.txt",outputdir.c_str());
   ofpar.open(fname,std::ios::app);
   ofstream purepar;
   sprintf(fname,"parkk");
   purepar.open(fname,std::ios::app);
   ofstream purepar2;
   sprintf(fname,"parkk_phi");
   purepar2.open(fname,std::ios::app);
   sprintf(fname,"%s/plot_kk.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");
   double mass;
   TTree *dataraw = new TTree("dataraw","dataraw");
   dataraw->Branch("x",&mass,"x/D");
   TTree *datarawf = new TTree("datarawf","dataraw");
   datarawf->Branch("xx",&mass,"xx/D");
   
   
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


   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1    = new TH1D("h1","2 kaon invariant mass",nBins,D0low,D0up);
   TCanvas *c1 = new TCanvas("","",800,600);

   const int Npart = 1;
    double costhecut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
    double start=-1;
    double stop =1;
    for(int i=0;i<Npart+1;i++){
      costhecut[i] = (stop-start)/Npart*i+start;
    }
   double m0=peakvalue;
   double peakest[Npart];
   double sigma_m=0.0024;//0.0024 for phi,
   double width = 10.*sigma_m;
   double mparticle=0.493677;
   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;

   int realsize=0;
   double partid[Npart];
   double parter[Npart];
   double corfac[Npart];
   double corerr[Npart];
   
   double pcut[Npart+1];
   double facv[Npart];
   double facev[Npart];


   //pcut[0]=0;
   //pcut[1] = 2.0;
 //pcut[0] =0.0 ;    facv[0] =1.0;       facev[0] =1.0;  
 //pcut[1] =0.10;    facv[1] =1.0     ;  facev[1] =1.0        ;  
 //pcut[2] =0.20;    facv[2] =1.0     ;  facev[2] =1.0        ;
 //pcut[3] =0.30;    facv[3] =1.0     ;  facev[3] =1.0        ;
 //pcut[4] =0.40;    facv[4] =1.0     ;  facev[4] =1.0        ;
 //pcut[5] =0.50;    facv[5] =1.0     ;  facev[5] =1.0        ;
 //pcut[6] =0.60;    facv[6] =1.0     ;  facev[6] =1.0        ;
 //pcut[7] =0.70;    facv[7] =1.0     ;  facev[7] =1.0        ;
 //pcut[8] =0.80;    facv[8] =1.0     ;  facev[8] =1.0        ;
 //pcut[9] =0.90;    facv[9] =1.0     ;  facev[9] =1.0        ;
 //pcut[10]=1.00;    facv[10]=1.0     ;  facev[10]=1.0        ;
 //pcut[11]=1.10;    facv[11]=1.0     ;  facev[11]=1.0        ;
 //pcut[12]=1.20;    facv[12]=1.0     ;  facev[12]=1.0        ;
 //pcut[13]=1.30;    facv[13]=1.0     ;  facev[13]=1.0        ;
 //pcut[14]=1.40;    facv[14]=1.0     ;  facev[14]=1.0        ;
 //pcut[15]=1.50;    facv[15]=1.0     ;  facev[15]=1.0        ;
 //pcut[16]=1.60;    facv[16]=1.0     ;  facev[16]=1.0        ;
 //pcut[17]=1.70;    facv[17]=1.0     ;  facev[17]=1.0        ;
 //pcut[18]=1.80;    facv[18]=1.0;       facev[18]=1.0;  
 //pcut[19]=1.90;    facv[19]=1.0;       facev[19]=1.0;  
 //pcut[20]=2.00;    facv[20]=1.0;       facev[19]=1.0;  
 //pcut[21]=2.20;           

   char name[100];
   TH1D *hmD0[Npart];
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%0d",partj);
     hmD0[partj] = new TH1D(name,name,100,D0low,D0up);
   }
   TH1D *hmphi[Npart];
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"massphi_part%0d",partj);
     hmphi[partj] = new TH1D(name,name,100,philow,phiup);
   }
 
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
        if (mass>D0low && mass<D0up){
          vars->Fill();
        }
        if (mass>philow && mass<phiup){
          vars->Fill();
        }
        //if (costheta1 > 0.8 && costheta2 <-0.8) continue;
		//if (p1<0.5 || p1>1.3) continue;
		//if (p1<0.4 || p1>1.3) continue;
		//if (p2<0.4 || p2>1.3) continue;

        //if ( partj>=Npart || partj<0 ) continue;
        for (int partj=0;partj<Npart;partj++){
          //if (p1<pcut[partj] || p1>pcut[partj+1]) continue;
          //if (p2<pcut[partj] || p2>pcut[partj+1]) continue;
          
          if ( costheta1>costhecut[partj] && costheta1<costhecut[partj+1] )
		  if (mass>D0low-0.02 && mass<D0up+0.02){
            hmD0[partj]->Fill(mass);
			evts[partj].push_back(evt);
            break;
		  }
          if (mass>philow-0.02 && mass<phiup+0.02){
            hmphi[partj]->Fill(mass);
			evtsphi[partj].push_back(evt);
            break;
		  }
		}
   
   }
   std::cout<<"part No: "<< Npart<<std::endl;
   vars->Write();
   for (int partj=0;partj<Npart;partj++){
     hmD0[partj]->Write();
	 std::cout<<"histogram part "<<partj<<" has entries: "<<hmD0[partj]->GetEntries()<<std::endl;
     if (hmD0[partj]->GetEntries() > 100
       &&hmD0[partj]->GetMaximumBin()> 30 //check the the peak position,
       &&hmD0[partj]->GetMaximumBin()< 70 
       ){
       partmap.push_back(partj);
     	      	
   	////PeakEstimate::SetMResonace(peakvalue);
	////PeakEstimate::SetMFinal(mk);
	////TH1D *hp = new TH1D("hp","hp",200,0,2);
   	//////TH2D *hp2= new TH2D("hp2","hp2",100,0,2, 100,0,2);
	////double plow = pcut[partj];
	////double pup  = pcut[partj+1];
    ////std::cout<<"p range is "<<plow<<", "<<pup<<std::endl;
	////char cuts[100];
	////sprintf(cuts,"p1>%f & p1<%f",plow,pup);
    ////hp->Reset();
	////vars->Draw("p1>>hp",cuts);
	////hp->Fit("gaus","","",plow,pup);
	////double mean = hp->GetFunction("gaus")->GetParameter(1);
	////double sigma= hp->GetFunction("gaus")->GetParameter(2);
    ////hp->Reset();
	////vars->Draw("p2>>hp",cuts);
	////hp->Fit("gaus","","",pcut[0],pcut[Npart]);
	////double mean2 = hp->GetFunction("gaus")->GetParameter(1);
	////double sigma2= hp->GetFunction("gaus")->GetParameter(2);
	////PeakEstimate::SetPdisPar(mean, sigma, mean2, sigma2);
	//////std::cout<<"before est"<<std::endl;
	////double ps = PeakEstimate::peakshift(plow,pup,0.,2.);
	////peakest[partj] = ps;
	////delete hp;
		//std::cout<<"after  est"<<std::endl;
       	//hmD0[partj]->Fit("gaus");
	   	//double mpeak = hmD0[partj]->GetFunction("gaus")->GetParameter(1);
	   	//std::cout<<"fit peak is "<< mpeak << ", est peak is "<<ps<<std::endl;

	   	peakest[partj] = peakvalue;
	 }
   
   }

   // ~~~~~~~~ draw end

   // ~~~~~~~~~~ D0 part ~~~~~~~~~~~~~
   for (int loopj=0;loopj<partmap.size();loopj++){
     int partj=partmap.at(loopj);

   // for saving the fit result
   //double factori = 1.000714;
   int fittimes = 0;
   factor = factorstart;
   for (int i=0;i<pointNo;i++){
      //xframe->Clear();
      xframe = x.frame(Title("fit kaon"));

      //h1->Reset();
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evts[partj].size();jentry++) {
		//p1 = evts[partj].at(jentry).GetP1();
		//p2 = evts[partj].at(jentry).GetP2();
        //for (int i=0;i<Npart;i++){
        //  if (p1>=pcut[i]&&p1<pcut[i+1]){
        //    factori = facv[i];
        //    break;
        //  }
        //}
          mass = evts[partj].at(jentry).InvMass(factor,factor);
          if (mass>D0low && mass<D0up){
            dataraw->Fill();
          }
         // if (Cut(ientry) < 0) continue;
      }
      //dataraw->Write();

      dataset = new RooDataSet("dataset","data",dataraw,x);
      char tmpchr[100];
      sprintf(tmpchr,"data_k_%02d",fittimes);
      //data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      mean.setVal(peakvalue+0.06*(factor-1.0));
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
      sprintf(name,"factor = %.6f",factor);
      pt->AddText(name);
      pt->Draw();
      sprintf(name,"part%d_fitFor_%dth_time",partj,fittimes);
      c1->SetName(name);
      c1->Write();
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;

      //sprintf(tmpchr,"data_k_%d.eps",fittimes);
      //h1->Draw();
      //c1->Print(tmpchr);
      
      // save pars
      factors[i]=factor;
      factorserr[i]=0;
	  deltapeaks[i] = mean.getVal() - peakest[partj];
	  //deltapeaks[i] = mean.getVal() - peakvalue;
      deltapeakserr[i] = mean.getError();

      fittimes++;
      factor += factorstep;
   }
   
   TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   factor = facfit->GetParameter(0);
   sprintf(name,"factors_D0kk_part%d",partj);
   graph1->SetName(name);
   graph1->Write();

   // draw the best fitting
      xframe = x.frame(Title("fit kaon"));
      dataraw->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evts[partj].size();jentry++) {
	/////p1 = evts[partj].at(jentry).GetP1();
	/////p2 = evts[partj].at(jentry).GetP2();
       //for (int i=0;i<Npart;i++){
       //   if (p1>=pcut[i]&&p1<pcut[i+1]){
       //     factori = facv[i];
       //     break;
       //   }
       //}
        // factori = factor;

         mass = evts[partj].at(jentry).InvMass(factor,factor);
         if (mass>D0low && mass<D0up){
           dataraw->Fill();
         }
      }
      //dataraw->Write();

      char tmpchr[100];
      sprintf(tmpchr,"data_k");
      dataset = new RooDataSet("dataset","data",dataraw,x);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      Npar=8;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      mean.setVal(peakvalue+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(D0low,D0up));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
      sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
      sum->plotOn(xframe);
      xframe->Draw();
      TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
   ///pt->SetBorderSize(0);
   ///pt->SetFillStyle(4000);
   ///pt->SetTextAlign(12);
   ///pt->SetTextFont(42);
   ///pt->SetTextSize(0.035);
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
      double factor4err=TMath::Sqrt(TMath::Power(mean.getError()/facfit->GetParameter(1),2) + TMath::Power(facfit->GetParError(0),2));
      sprintf(name,"factor = %.6f #pm %.6f",factor,factor4err);
      pt->AddText(name);

      pt->Draw();
      sprintf(name,"part%d_final_fit",partj);
      c1->SetName(name);
      c1->Write();      
      
      partid[loopj] = pcut[partj]+(pcut[partj+1]-pcut[partj])/2;
      parter[loopj] = 0;
      corfac[loopj] = factor;
      corerr[loopj] = factor4err;
   

	  //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;
   }

  realsize = partmap.size();
   
  for (int i=0;i<realsize;i++){
    ofpar<<"cos="<<costhecut[i]<<"\tfactor: "<<corfac[i]<<"\t +/- \t"<< corerr[i]<<std::endl;
    purepar<<beamene<<"\t"<<corfac[i]<<"\t"<< corerr[i]<<std::endl;
  }
  // ~~~~~~~~~~~~D0 part end ........
//return;
   // ~~~~~~~~~~ phi part ~~~~~~~~~~~~~
   RooRealVar xx("xx","energy",phipeak,philow,phiup,"GeV");
   RooRealVar meanf("meanf","mean of gaussian",phipeak,philow,phiup);
   RooRealVar sigmaf("sigmaf","width of gaussian",0.004,0.001,0.008);//phi version
  // RooRealVar sigmaf("sigmaf","width of gaussian",0.0024,0.001,0.003);//phi version
  // RooRealVar sigma2f("sigma2f","width of gaussian",0.004,0.003,0.005);
  // RooGaussian gausf("gausf","gauss(x,m,s)",xx,meanf,sigmaf);
  // RooGaussian gaus2f("gaus2f","gauss(x,m,s)",xx,meanf,sigma2f);
   
   RooBreitWigner brewig("brewig","breitwigner",xx,meanf,sigmaf);
   
   RooRealVar co1f("co1f","coefficient #1",0,-0.5,0.5);
   RooRealVar co2f("co2f","coefficient #2",0,-0.01,0.5);
   RooChebychev bkgf("bkgf","background",xx,RooArgList(co1f,co2f));
   
   RooRealVar signalf("signalf"," ",1200,10,1000000);//event number
   RooRealVar signal2f("signal2f"," ",1200,0,1000000);//event number
   RooRealVar backgroundf("backgroundf"," ",200,0,1000000);
   
   // for saving the fit result
   factor = factorstart;
   int partj = 0;
   int  fittimes = 0 ;
   double simplex[100];
   double simplexe[100];
   double simpley[100];
   double simpleye[100];
   for (int i=0;i<pointNo;i++){
      //xframe->Clear();
      xframe = xx.frame(Title("fit kaon"));

      datarawf->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evtsphi[partj].size();jentry++) {
		  p1 = evtsphi[partj].at(jentry).GetP1();
		  p2 = evtsphi[partj].at(jentry).GetP2();
        //for (int i=0;i<Npart;i++){
        //  if (p1>=pcut[i]&&p1<pcut[i+1]){
        //    factori = facv[i];
        //    break;
        //  }
        //}
          mass = evtsphi[partj].at(jentry).InvMass(factor,factor);
          if (mass>philow && mass<phiup){
            datarawf->Fill();
          }
         // if (Cut(ientry) < 0) continue;
      }
      //datarawf->Write();
      std::cout<<"Loop end aaaa"<<std::endl;

      dataset = new RooDataSet("datasetf","data",datarawf,xx);
      char tmpchr[100];
      sprintf(tmpchr,"data_kf_%02d",fittimes);
      //sum = new RooAddPdf("sumf","sum",RooArgList(gausf,gaus2f,bkgf),RooArgList(signalf,signal2f,backgroundf));
      Npar=6;
      sum = new RooAddPdf("sum","sum",RooArgList(brewig,bkgf),RooArgList(signalf,backgroundf));
      meanf.setVal(phipeak+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(brewig),LineStyle(2),LineColor(2));
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
      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",meanf.getVal(),meanf.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigmaf.getVal(),sigmaf.getError());
      pt->AddText(tmpchr);
    //sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2f.getVal(),sigma2f.getError());
    //pt->AddText(tmpchr);
      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signalf.getVal(),signalf.getError());
      pt->AddText(tmpchr);
    //sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2f.getVal(),signal2f.getError());
    //pt->AddText(tmpchr);
      sprintf(tmpchr,"backNo = %.2f #pm %.2f",backgroundf.getVal(),backgroundf.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      sprintf(name,"factor = %.6f",factor);
      pt->AddText(name);
      pt->Draw();
      sprintf(name,"part%d_fitFor_%dth_time",partj,fittimes);
      c1->SetName(name);
      c1->Write();
      //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;
   std::cout<<"Fit and save end"<<std::endl;

      //sprintf(tmpchr,"data_k_%d.eps",fittimes);
      //h1->Draw();
      //c1->Print(tmpchr);
      
      // save pars
      factors[i]=factor;
      factorserr[i]=0;
	  //deltapeaks[i] = meanf.getVal() - peakest[partj];
	  deltapeaks[i] = meanf.getVal() - phipeak;
      deltapeakserr[i] = meanf.getError();


	  TH1D hmasstmp("hmasstmp","hmasstmp",100,philow,phiup);
	  std::cout<<"entries "<<datarawf->GetEntries()<<std::endl;
	  datarawf->Draw("xx>>hmasstmp");
	  hmasstmp.Write("rawhist");
	  hmasstmp.Fit("gaus","","",meanf.getVal()-0.003, meanf.getVal()+0.003);
	  simplex[i] = factor;//hmasstmp.GetFunction("gaus")->GetParameter(1);
	  simplexe[i] = 0;// hmasstmp.GetFunction("gaus")->GetParError(1);
	  simpley[i] = hmasstmp.GetFunction("gaus")->GetParameter(1) - phipeak;
	  simpleye[i] = hmasstmp.GetFunction("gaus")->GetParError(1);
	  sprintf(name,"simple_fit_factor_%d",fittimes);
	  hmasstmp.Write(name);

      fittimes++;
      factor += factorstep;
   }
   
   TGraphErrors *graph1f;
   graph1f = new TGraphErrors(pointNo,simplex,simpley,simplexe,simpleye);
   graph1f->SetTitle("delta peak");
   graph1f->SetMarkerStyle(5);
   graph1f->Draw("AP");
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph1f->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   //factor = facfit->GetParameter(0);
   //sprintf(name,"factors_kk_part%d",partj);
   //graph1f->SetName(name);
   graph1f->Write("simpleFItVersion");


   graph1f = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph1f->SetTitle("delta peak");
   graph1f->SetMarkerStyle(5);
   graph1f->Draw("AP");
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.05);
   facfit->SetParNames("factor","slope");
   graph1f->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   factor = facfit->GetParameter(0);
   sprintf(name,"factors_kk_part%d",partj);
   graph1f->SetName(name);
   graph1f->Write("RooFitVersion_graph");
   
   
   // draw the best fitting
      xframe = xx.frame(Title("fit kaon"));
      datarawf->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<evtsphi[partj].size();jentry++) {
		 p1 = evtsphi[partj].at(jentry).GetP1();
		 p2 = evtsphi[partj].at(jentry).GetP2();
       //for (int i=0;i<Npart;i++){
       //   if (p1>=pcut[i]&&p1<pcut[i+1]){
       //     factori = facv[i];
       //     break;
       //   }
       //}
         mass = evtsphi[partj].at(jentry).InvMass(factor,factor);
         if (mass>philow && mass<phiup){
           datarawf->Fill();
         }
      }
      //datarawf->Write();

      char tmpchr[100];
      sprintf(tmpchr,"data_k");
      dataset = new RooDataSet("dataset","data",datarawf,xx);
      //sum = new RooAddPdf("sum","sum",RooArgList(gausf,gaus2f,bkgf),RooArgList(signalf,signal2f,backgroundf));
      sum = new RooAddPdf("sum","sum",RooArgList(brewig,bkgf),RooArgList(signalf,backgroundf));
      //Npar=6;
      //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background2));
      meanf.setVal(phipeak+0.06*(factor-1.0));
      sum->fitTo(*dataset,Range(philow,phiup));
      dataset->plotOn(xframe);
      sum->plotOn(xframe,Components(brewig),LineStyle(2),LineColor(2));
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
      sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",meanf.getVal(),meanf.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigmaf.getVal(),sigmaf.getError());
      pt->AddText(tmpchr);
    //sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2f.getVal(),sigma2f.getError());
    //pt->AddText(tmpchr);
      sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signalf.getVal(),signalf.getError());
      pt->AddText(tmpchr);
    //sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2f.getVal(),signal2f.getError());
    //pt->AddText(tmpchr);
      sprintf(tmpchr,"backNo = %.2f #pm %.2f",backgroundf.getVal(),backgroundf.getError());
      pt->AddText(tmpchr);
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      double factor4err=TMath::Sqrt(TMath::Power(meanf.getError()/facfit->GetParameter(1),2) + TMath::Power(facfit->GetParError(0),2));
      sprintf(name,"factor = %.6f #pm %.6f",factor,factor4err);
      pt->AddText(name);

      pt->Draw();
      sprintf(name,"part%d_final_fit",partj);
      c1->SetName(name);
      c1->Write();      
      
	  //delete data_k;
      delete dataset;
      delete xframe;
      delete sum;

   
    ofpar<<"\tfactor: "<<factor<<"\t +/- \t"<< factor4err<<std::endl;
    purepar2<<beamene<<'\t' << "\t"<<factor<<"\t"<< factor4err<<std::endl;
  
  // ~~~~~~~~~~~~phi part end ........
   
   f->Close();
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
