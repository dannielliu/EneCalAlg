#define Ks0Alg_cxx
#include "Ks0Alg.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "function.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include <vector>
#include <fstream>
extern std::string outputdir;
using namespace RooFit;


void Ks0Alg::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Ks0Alg.C
//      Root > Ks0Alg t
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

   //Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //   Long64_t ientry = LoadTree(jentry);
   //   if (ientry < 0) break;
   //   nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   //}


   std::cout<<"loop"<<std::endl;

   TH1D *hmass  = new TH1D("hmass" ,"Ks mass",200,0.45,0.55);
   TH1D *hmassU = new TH1D("hmassU","Ks mass",200,0.45,0.55);
   TH1D *hmassc = new TH1D("hmassc","Ks mass",200,0.45,0.55);
   //TH2D *h2pos = new TH2D("h2pos","pos",10, 0,5, 10,0,5);
   TH2D *h2p = new TH2D("h2p","P",200, 0,2, 200,0,2);
   int nBins=100;
   double mpi=0.13957;
   double mparticle=0.497614;
   double peakvalue = mparticle;
   //double kslow=0.45;
   //double ksup =0.55;
   double kslow=0.475;
   double ksup =0.525;
   TF1 *facfit = new TF1("facfit",line2,kslow,ksup,2);
   char fname[1000];
   sprintf(fname,"%s/plot_ks.root",outputdir.c_str());
   TFile *f=new TFile(fname,"RECREATE");
   
   sprintf(fname,"%s/pars.txt",outputdir.c_str());
   ofstream ofpar(fname,std::ios::app);
   sprintf(fname,"%s/parspur.txt",outputdir.c_str());
   ofstream ofparpur(fname,std::ios::app);

   TCanvas *c1 = new TCanvas();
   int Npar;

   // roofit variables and functions
   RooRealVar x("x","energy",mparticle,kslow,ksup,"GeV");
   RooRealVar mean("mean","mean of gaussian",  kslow,ksup);
   //RooRealVar mean2("mean2","mean of gaussian",kslow,ksup);
   RooRealVar sigma("sigma","width of gaussian",0.0020,0.0010,0.0022);
   RooRealVar sigma2("sigma2","width of gaussian",0.004,0.0025,0.02);
   //RooRealVar brewid("brewid","width of breit wigner",0.0023,0.0010,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
 
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));

   RooRealVar signal("signal"," ",1000,0,10000000);//event number
   RooRealVar signal2("signal2"," ",10,0,10000000);//event number
   RooRealVar background("background"," ",10,0,1000000);
 
   RooAddPdf *sum;
   //RooDataHist *datahist;
   RooDataSet *dataset;
   RooPlot *xframe;  
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   TTree *dataraw = new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");
   TTree *vars = new TTree("vars","vars");
   double phi,phi1,phi2;
   double costheta,costheta1,costheta2;
   double p1,p2;
   double pt1,pt2;
   vars->Branch("phi",&phi,"phi/D");
   vars->Branch("phi1",&phi1,"phi1/D");
   vars->Branch("phi2",&phi2,"phi2/D");
   vars->Branch("costheta",&costheta,"costheta/D");
   vars->Branch("costheta1",&costheta1,"costheta1/D");
   vars->Branch("costheta2",&costheta2,"costheta2/D");
   vars->Branch("p1",&p1,"p1/D");
   vars->Branch("p2",&p2,"p2/D");
   vars->Branch("pt1",&pt1,"pt1/D");
   vars->Branch("pt2",&pt2,"pt2/D");
   vars->Branch("mass",&mass,"mass/D");

   const int Npart= 10;
   double coscut[Npart+1];
   for (int i=0; i<Npart+1; i++){
     coscut[i] = 2./Npart*i - 1;
   }
   //int Npart1 = 20;
   //int Npartcos = 10;
   int realsize=0;
   double peakest[Npart];
   double partid[Npart];
   //double partjd[Npart1];
   double parter[Npart];
   //double partje[Npart1];
   double comsgm[Npart][2];
   double corfac[Npart];
   double corerr[Npart];
   double start=0;
   double stop =2.0;
   double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   double facv[Npart];
   double facev[Npart];
   //pcut[0]=0.1;
   //pcut[1]=0.4;
   //pcut[1]=2.0;
   //pcut[2]=2.0;
   //pcut[2]=0.9;
// rvalue combine pi+ and pi-, factor in (0.2, 0.3) to 1.00061
   // factors in range (0.15, 0.6)
   // factor got from previous correction
// pcut[0] =0.0 ;    facv[0] =1.0     ;  facev[0] =1.0;  
// pcut[1] =0.05;    facv[1] =1.00201 ;  facev[1] =0.00390561;  
// pcut[2] =0.10;    facv[2] =1.00203 ;  facev[2] =0.0010036 ;
// pcut[3] =0.15;    facv[3] =1.00184 ;  facev[3] =0.00059564;
// pcut[4] =0.20;    facv[4] =1.00044 ;  facev[4] =0.00047191;
// pcut[5] =0.25;    facv[5] =1.0007  ;  facev[5] =0.00041868;
// pcut[6] =0.30;    facv[6] =1.00089 ;  facev[6] =0.00050795;
// pcut[7] =0.35;    facv[7] =1.00082 ;  facev[7] =0.00067682;
// pcut[8] =0.40;    facv[8] =1.00257 ;  facev[8] =0.00119303;
// pcut[9] =0.45;    facv[9] =1.00062 ;  facev[9] =0.00125913;
// pcut[10]=0.50;    facv[10]=0.997233;  facev[10]=0.00532376;
// pcut[11]=0.60;    facv[11]=1.00207 ;  facev[11]=0.00209495;
// pcut[12]=0.70;    facv[12]=1.00469 ;  facev[12]=0.00122247;
// pcut[13]=0.80;    facv[13]=1.0175  ;  facev[13]=0.00354285;
// pcut[14]=0.90;    facv[14]=1.00342 ;  facev[14]=0.00451459;
// pcut[15]=1.00;    facv[15]=1.00781 ;  facev[15]=0.00679501;
// pcut[16]=1.20;    facv[16]=0.993716;  facev[16]=0.00178269;
// pcut[17]=1.40;    facv[17]=1.01757 ;  facev[17]=0.0216559 ;
// pcut[18]=1.60;    facv[18]=1.0;       facev[18]=1.0;  
// pcut[19]=1.80;    facv[19]=1.0;       facev[19]=1.0;  
// pcut[20]=2.20;           
// 
   //for(int i=0;i<Npart+1;i++){
   //  pcut[i] = (stop-start)/Npart*i+start;
   //}

   std::vector<int> partmap;
   std::vector<std::pair<int,double> > facmap;
  
   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different range
   TH1D *hmKs[Npart];
   for (int partj=0;partj<Npart;partj++){
     sprintf(name,"mass_part%0d",partj);
     hmKs[partj] = new TH1D(name,name,100,kslow,ksup);
   }
 
   // loop data
   std::cout<<"aaaaaaaaaaaaaaaaf"<<std::endl;
   Event evt;
   std::vector<Event> evts[Npart];
   Long64_t nbytes = 0, nb = 0;
   int counts[10]={0,0,0,0,0,0,0,0,0,0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      counts[0]++;
	  Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   //std::cout<<"bbbbbbbbbbbbbbbbf"<<std::endl;
     
	  // for new Ks data
	  /*
      for (int idx=0; idx<nKsCan;idx++){
	    evt.SetVal(pippx[idx],pippy[idx],pippz[idx],pimpx[idx],pimpy[idx],pimpz[idx]);
        mass = evt.InvMass();
        double chi2 = m_chisq0[idx]+m_chisqvtxnd[idx];
        if (chi2 > 100) continue;
        double ratio = m_decayL[idx]/m_decayLerr[idx];
        if (ratio<2) continue;
        if (Ks_id[idx]==0){
          hmassU->Fill(mass);
          continue;
        }
        hmass->Fill(Mpippim[idx]);
        hmassc->Fill(mass);
        p1 = evt.GetP1();
        p2 = evt.GetP2();
        //if (p1<0.15||p1>0.6) continue;
        if (p1<0.05||p1>1.4) continue;
        //if (p2<0.1) continue;
        //if (p1+p2<0.4 || p1+p2 > 0.6) continue;
        costheta1 = evt.GetCostheta1();
        costheta2 = evt.GetCostheta2();
        phi1 = evt.GetPhi1();
        phi2 = evt.GetPhi2();
        if (mass>kslow && mass<ksup){
          vars->Fill();
        }

        //if ( partj>=Npart || partj<0 ) continue;
        for (int partj=0;partj<Npart1;partj++){
          if (p2<pcut1[partj] || p2>pcut1[partj+1]) continue;
          //if (p2<pcut[partj] || p2>pcut[partj+1]) continue;
          if (mass>kslow-0.02 && mass<ksup+0.02){
            hmKs[partj]->Fill(mass);
			evts[partj].push_back(evt);
	      }
          break;
		}
       }
*/
      
	  // for old Ks data
	   
	  for (int ipip=0; ipip<npip;ipip++)
      for (int ipim=0; ipim<npim;ipim++){
	    evt.SetVal(pimpx[ipim],pimpy[ipim],pimpz[ipim],pippx[ipip],pippy[ipip],pippz[ipip]);// change + -
        mass = evt.InvMass();
        double chi2 = m_chisq0[ipip*npim+ipim]+m_chisqvtxnd[ipip*npim+ipim];
        if (chi2 > 100) continue;
		counts[1]++;
        double ratio = m_decayL[ipip*npim+ipim]/m_decayLerr[ipip*npim+ipim];
        //if (ratio>2 || ratio<1) continue;//select Ks from ip
        if (ratio<2 ) continue;//select Ks from ip
		counts[2]++;
        if (Ks_id[ipip*npim+ipim]==0){
          hmassU->Fill(mass);
          continue;
        } 
		counts[3]++;
        hmass->Fill(Mpippim[ipip*npim+ipim]);
        hmassc->Fill(mass);
        p1 = evt.GetP1();
        p2 = evt.GetP2();
        pt1 = evt.GetPt1();
        pt2 = evt.GetPt2();
        //if (p1<0.15|| p1>0.6) continue;
        //if (p1<0.1 || p1>0.4) continue;
        //if (p1<0.1) continue;
        //if (p2<0.1) continue;
        //if (p1+p2<0.4 || p1+p2 > 0.6) continue;
        costheta1 = evt.GetCostheta1();
        costheta2 = evt.GetCostheta2();
        phi1 = evt.GetPhi1();
        phi2 = evt.GetPhi2();
        if (mass>kslow && mass<ksup){
		counts[4]++;
          vars->Fill();
        }

        //if ( partj>=Npart || partj<0 ) continue;
        for (int partj=0;partj<Npart;partj++){
          //if (p1<pcut[partj] || p1>pcut[partj+1]) continue;
          //if (pt2<pcut[partj] || pt2>pcut[partj+1]) continue;
          if (costheta1>coscut[partj] && costheta1<coscut[partj+1])
          if (costheta2>coscut[partj] && costheta2<coscut[partj+1])
		  if (mass>kslow-0.02 && mass<ksup+0.02){
            hmKs[partj]->Fill(mass);
			evts[partj].push_back(evt);
            break;
	      }
		}
       }
   }
   std::cout<<"ffffffffff cut flow: "<<counts[0]<<' '<<counts[1]<<' '<<counts[2]<<' '<<counts[3]<<' '<<counts[4]<<std::endl;
   vars->Write();
////TH1D *hp = new TH1D("hp","hp",200,0,2);
////TF1 f2("f2","[0]*(x-[1])*exp(-(x-[1])*(x-[1])/[2])");
////f2.SetParameters(10000, 0, 0.2);
////hp->Reset();
////vars->Draw("p2>>hp");
//////hp->Fit("gaus","","",plow-0.1,pup+0.1);
////hp->Fit("f2","","",pcut[0],pcut[Npart]);
////hp->Fit("f2","","",pcut[0],pcut[Npart]);
////hp->SetTitle("p2_dis");
////hp->Write();
////double fmean2 = hp->GetFunction("f2")->GetParameter(1);
////double fsigma2= hp->GetFunction("f2")->GetParameter(2);
////std::cout<<"mean2, sigma2:  "<<fmean2<<' '<<fsigma2<<std::endl;
   for (int partj=0;partj<Npart;partj++){
     hmKs[partj]->Write();
     if (hmKs[partj]->GetEntries() > 100
       &&hmKs[partj]->GetMaximumBin()> 30 //check the the peak position,
       &&hmKs[partj]->GetMaximumBin()< 70 
       ){ 
       partmap.push_back(partj);
	   peakest[partj]=peakvalue;
     	      	
   	////PeakEstimate::SetMResonace(peakvalue);
	////PeakEstimate::SetMFinal(mpi);
   	//////TH2D *hp2= new TH2D("hp2","hp2",100,0,2, 100,0,2);
	////double plow = pcut[partj];
	////double pup  = pcut[partj+1];
    ////std::cout<<"p range is "<<plow<<", "<<pup<<std::endl;
	////char cuts[100];
	////sprintf(cuts,"p2>%f & p2<%f",plow,pup);
    ////hp->Reset();
	////vars->Draw("p1>>hp",cuts);
	////TF1 f1("f1","[0]*sqrt(x-[1])*exp(-(x-[1])/[2])");
	////f1.SetParameters(10000, 0., 0.2);
	////hp->Fit("f1","","",pcut[0],pcut[Npart]);
	////hp->Fit("f1","","",pcut[0],pcut[Npart]);
	////hp->SetTitle(cuts);
	////hp->Write();
	////double fmean = hp->GetFunction("f1")->GetParameter(1);
	////double fsigma= hp->GetFunction("f1")->GetParameter(2);

	////PeakEstimate::SetPdisPar(fmean, fsigma, fmean2, fsigma2);
	////std::cout<<fmean<<' '<<fsigma<<' '<<fmean2<<' '<<fsigma2<<std::endl;
	////std::cout<<"before est"<<std::endl;
	////double ps = PeakEstimate::peakshift(0.0, 2.0, plow, pup);
	////std::cout<<"after est"<<std::endl;
	////peakest[partj] = ps;
	//////peakest[partj] = peakvalue;
	////if (fabs(ps-peakvalue)>0.01) {
	////	peakest[partj]=peakvalue;
	////	std::cout<<"Warning: not a suitable peak shift at "<< partj<<", ps is "<<ps<<std::endl;
	////}
		//delete hp;

     }
   }

//return;

   hmass->Write();
   hmassU->Write();
   hmassc->Write();
  
   const int pointNo=5;
   double factor=0.995;
   double factorstep=(1.-factor)*2/pointNo;
   int fittimes=0;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];

   std::cout<<"part map size is "<<partmap.size()<<std::endl;
   for (int loopj=0;loopj<partmap.size();loopj++){
     int partj = partmap.at(loopj);
	 //int cosj  = partmap.at(loopj).second;
     std::cout<<"part is "<<partj<<std::endl;
     double factori=1.000  ;
     double factorj=1.000  ;
     factor = 0.995;
     //factori=1.000815;
     //factori=factor;
     fittimes=0;
    
     signal.setVal(10000);
     signal.setError(100);
     signal2.setVal(0);
     background.setVal(0);
     background.setError(10);
     a0.setVal(0);
     a1.setVal(0);
     sigma.setVal(0.002);
     sigma2.setVal(0.004);
     double slope = partj/100.;
     mean.setVal(peakest[partj]+slope*(factor-1));
     mean.setError(0.00001);


     for (fittimes=0; fittimes<pointNo;fittimes++){
       xframe = x.frame(Title("fit Ks"));
       dataraw->Reset();
       std::cout<<"factor is "<<factor<<std::endl;
       
	   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
       std::cout<<"evts size is "<<evts[partj].size()<<std::endl;
       for (Long64_t jin=0; jin<evts[partj].size();jin++) {
		  //factori = factor;
		  //factori = 1.0;
		  factorj = factor;
		//p1 = evts[partj].at(jin).GetP1();
		//p2 = evts[partj].at(jin).GetP2();
	  //  pt1 = evts[partj].at(jin).GetPt1();
	  //  pt2 = evts[partj].at(jin).GetPt2();
		  //if ( p1<0.4 && p1>0.1) factori = 1.0009 ;
		  //else continue;
		  //if ( p2<0.4) factorj = 1.0009 ;
		  //costheta2 = evts.at(jin).GetCostheta2();
          //if (p1<0.15||p1>0.6) continue;
       // for (int i=0;i<Npart;i++){
       //   if (pt1>=pcut[i]&&pt1<pcut[i+1]){
       //     factori = facv[i];
       //     break;
       //   }
       // }
          //if (partj<0 || partj>Npart) continue;
          //if (p1 < pcut[partj] || p1> pcut[partj+1]) continue;
          //if (p2 < pcut1[partj] || p2> pcut1[partj+1]) continue;
          mass = evts[partj].at(jin).InvMass(factor,factor);
		  if (mass>kslow && mass<ksup){
            dataraw->Fill();
          }
       }// loop data end
       // fit data

       std::cout<<"dataraw size is "<<dataraw->GetEntries()<<std::endl;
       sprintf(name,"data_Ks_%02d",fittimes);
       dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*dataraw));
       sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
       Npar=8;
       
  // signal.setVal(10000);
  // signal.setError(100);
  // signal2.setVal(0);
  // background.setVal(0);
  // //background.setError(10);
  // background.setConstant(kTRUE);
  // //a0.setVal(0);
  // //a1.setVal(0);
  // sigma.setVal(0.002);
  // sigma.setError(0.0007);
  // sigma2.setVal(0.004);
  // sigma2.setError(0.0002);
  // double slope = partj/100.;
  // mean.setVal(peakest[partj]+slope*(factor-1));
  // mean.setError(0.00003);

       sum->fitTo(*dataset,Range(kslow,ksup));
       dataset->plotOn(xframe,Binning(nBins));
       sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
       sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
       sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
       sum->plotOn(xframe);
	   double chi2 = xframe->chiSquare(Npar);
       xframe->Draw();
       TPaveText *pt = new TPaveText(0.70,0.60,0.95,0.95,"BRNDC");
       pt->SetBorderSize(0);
       pt->SetFillStyle(4000);
       pt->SetTextAlign(12);
       pt->SetTextFont(42);
       pt->SetTextSize(0.035);
       sprintf(name,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
       pt->AddText(name);
       sprintf(name,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
       pt->AddText(name);
       sprintf(name,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
       pt->AddText(name);
       sprintf(name,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
       pt->AddText(name);
       sprintf(name,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
       pt->AddText(name);
       sprintf(name,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
       pt->AddText(name);
       sprintf(name,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
       pt->AddText(name);
       pt->Draw();
       //c1->Update();
       sprintf(name,"part%d_fitFor_%dth_time",partj,fittimes);
       c1->SetName(name);
       c1->Write();

       factors[fittimes]=factor;
       factorserr[fittimes]=0;
       deltapeaks[fittimes] = mean.getVal() - peakest[partj];
       deltapeakserr[fittimes] = mean.getError();
	   if (chi2>10) deltapeakserr[fittimes] = 1.0;
      
       factor += factorstep;
       //factori = factor;

       delete sum;
       delete dataset;
       delete xframe;

	   if (factor == 1.0){
	     double s1 = sigma.getVal();
		 double s2 = sigma2.getVal();
		 double n1 = signal.getVal();
		 double n2 = signal2.getVal();
	     comsgm[loopj][0] = sqrt((s1*s1*n1+s2*s2*n2)/(n1+n2));
	   } 
    }// n point end
    //will fit linear, get factor needed
    c1->Clear();
    TGraphErrors *graph = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
    graph->SetTitle("delta peak");
    graph->Draw("AP");
    graph->SetMarkerStyle(5);
    gStyle->SetOptFit(1111);
    facfit->SetParameters(1.0,0.4);
    facfit->SetParNames("factor","slope");
    graph->Fit(facfit,"","",factors[0],factors[pointNo-1]);
    factor = facfit->GetParameter(0);
    //factori=factor;
    TPaveText *pat1 = new TPaveText(0.12,0.70,0.5,0.90,"BRNDC");
    pat1->SetBorderSize(0);
    pat1->SetFillStyle(4000);
    pat1->SetTextAlign(12);
    pat1->SetTextFont(42);
    pat1->SetTextSize(0.035);
    sprintf(name,"factor = %1.6f #pm %1.6f",facfit->GetParameter(0),facfit->GetParError(0));
    pat1->AddText(name);
    sprintf(name,"slope = %1.6f #pm %1.6f",facfit->GetParameter(1),facfit->GetParError(1));
    pat1->AddText(name);
    pat1->Draw();
 
    sprintf(name,"part%d_factors",partj);
    c1->SetName(name);
    c1->Write();

    // try to use the best factor to fit
    xframe = x.frame(Title("fit Ks"));
    dataraw->Reset();
    std::cout<<"factor is "<<factor<<std::endl;
       
    for (Long64_t jin=0; jin<evts[partj].size();jin++) {
	  //factori = factor;
	  //factori = 1.0;
	  factorj = factor;
    //p1 = evts[partj].at(jin).GetP1();
    //p2 = evts[partj].at(jin).GetP2();
  //  pt1 = evts[partj].at(jin).GetPt1();
  //  pt2 = evts[partj].at(jin).GetPt2();
	  //if ( p1<0.4 && p1>0.1) factori = 1.0009 ;
	  //else continue;
	  //if ( p2<0.4) factorj = 1.0009 ;
		  //costheta2 = evts.at(jin).GetCostheta2();
          //if (p1<0.15||p1>0.6) continue;
      //  for (int i=0;i<Npart;i++){
      //    if (pt1>=pcut[i]&&pt1<pcut[i+1]){
      //      factori = facv[i];
      //      break;
      //    }
      //  }
          //if (partj<0 || partj>Npart) continue;
          //if (p1 < pcut[partj] || p1> pcut[partj+1]) continue;
          //if (p2 < pcut1[partj] || p2> pcut1[partj+1]) continue;
          mass = evts[partj].at(jin).InvMass(factor,factor);
		  if (mass>kslow && mass<ksup){
            dataraw->Fill();
          }
       }// loop data end
       // fit data


    // fit data
    sprintf(name,"data_Ks_part%d",partj);
    dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*dataraw));
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
    //Npar=8;
	       
  // signal.setVal(10000);
  // signal.setError(100);
  // signal2.setVal(0);
  // //background.setVal(0);
  // //background.setError(10);
  // a0.setVal(0);
  // a1.setVal(0);
  // sigma.setVal(0.002);
  // sigma.setError(0.0007);
  // sigma2.setVal(0.004);
  // sigma2.setError(0.0002);
  // mean.setVal(peakest[partj]);
  // mean.setError(0.00003);

    sum->fitTo(*dataset,Range(kslow,ksup));
    dataset->plotOn(xframe,Binning(nBins));
    sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
    sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
    sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
    sum->plotOn(xframe);
    xframe->Draw();
    TPaveText *pt = new TPaveText(0.65,0.50,0.95,0.95,"BRNDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(4000);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    sprintf(name,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
    pt->AddText(name);
    sprintf(name,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
    pt->AddText(name);
    sprintf(name,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
    pt->AddText(name);
    sprintf(name,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
    pt->AddText(name);
    sprintf(name,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
    pt->AddText(name);
    sprintf(name,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
    pt->AddText(name);
    sprintf(name,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
    pt->AddText(name);
    double factor4err=TMath::Sqrt(TMath::Power(mean.getError()/facfit->GetParameter(1),2) + TMath::Power(facfit->GetParError(0),2));
    sprintf(name,"factor = %.6f #pm %.6f",factor,factor4err);
    pt->AddText(name);

    pt->Draw();
    c1->Update();
    sprintf(name,"part%d_final_fit",partj);
    c1->SetName(name);
    c1->Write();

    partid[loopj] = partj;//pcut[partj]+(pcut[partj+1]-pcut[partj])/2;
	//partjd[loopj] = coscut[cosj]+(coscut[cosj+1]-coscut[cosj])/2;
    ofpar<<"p="<<partid[loopj]<<"\tpeak estimate : "<<peakest[partj]<<std::endl;
    parter[loopj] = 0;
	//partje[loopj] = 0;
    corfac[loopj] = factor;
    corerr[loopj] = factor4err;
   	
	double s1 = sigma.getVal();
	double s2 = sigma2.getVal();
	double n1 = signal.getVal();
	double n2 = signal2.getVal();
	comsgm[loopj][1] = sqrt((s1*s1*n1+s2*s2*n2)/(n1+n2));
    
	delete sum;
    delete dataset;
    delete xframe;
 
  }// n part end

  realsize = partmap.size();
  TGraphErrors *graph = new TGraphErrors(realsize,partid,corfac,parter,corerr);
  graph->SetTitle("compare factors");
  graph->Draw("AP");
  graph->SetMarkerStyle(5);
  gStyle->SetOptFit(1111);
  c1->SetName("compare_factors");
  c1->Write();

  for (int i=0;i<realsize;i++){
    ofpar<<"p="<<partid[i]<<"\tfactor: "<<corfac[i]<<"\t +/- \t"<< corerr[i]<<"\t ";
	ofpar<<'\t'<<"old_sigma: "<< comsgm[i][0]<<"\t cor_sigma: "<< comsgm[i][1]<<std::endl;
    ofparpur<<partid[i]<<"\t"<<corfac[i]<<"\t"<< corerr[i]<<std::endl;
  }

}

#ifdef Ks0Alg_cxx
Ks0Alg::Ks0Alg(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ksto2pi_583_140226_0035862.root");
      if (!f) {
         f = new TFile("Ksto2pi_583_140226_0035862.root");
      }
      tree = (TTree*)gDirectory->Get("Ks_info");

   }
   Init(tree);
}

Ks0Alg::~Ks0Alg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Ks0Alg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ks0Alg::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ks0Alg::Init(TTree *tree)
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

   fChain->SetBranchAddress("rec_truth_Mks", &rec_truth_Mks, &b_rec_truth_Mks);
   fChain->SetBranchAddress("rec_truth_Pks", &rec_truth_Pks, &b_rec_truth_Pks);
   fChain->SetBranchAddress("rec_truth_Eks", &rec_truth_Eks, &b_rec_truth_Eks);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("nKsCan", &nKsCan, &b_nKsCan);
   fChain->SetBranchAddress("m_chisq0", m_chisq0, &b_m_chisq0);
   fChain->SetBranchAddress("m_chisqvtxnd", m_chisqvtxnd, &b_m_chisqvtxnd);
   fChain->SetBranchAddress("m_decayL", m_decayL, &b_m_decayL);
   fChain->SetBranchAddress("m_decayLerr", m_decayLerr, &b_m_decayLerr);
   fChain->SetBranchAddress("m_ctau", m_ctau, &b_m_ctau);
   fChain->SetBranchAddress("Ks_id", Ks_id, &b_Ks_id);
   fChain->SetBranchAddress("Mpippim", Mpippim, &b_Mpippim);
   fChain->SetBranchAddress("Ppippim", Ppippim, &b_Ppippim);
   fChain->SetBranchAddress("Epippim", Epippim, &b_Epippim);
   fChain->SetBranchAddress("Thepippim", Thepippim, &b_Thepippim);
   fChain->SetBranchAddress("Phipippim", Phipippim, &b_Phipippim);
   fChain->SetBranchAddress("Ppip", Ppip, &b_Ppip);
   fChain->SetBranchAddress("Ppim", Ppim, &b_Ppim);
   fChain->SetBranchAddress("Thepip", Thepip, &b_Thepip);
   fChain->SetBranchAddress("Thepim", Thepim, &b_Thepim);
   fChain->SetBranchAddress("Phipip", Phipip, &b_Phipip);
   fChain->SetBranchAddress("Phipim", Phipim, &b_Phipim);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   Notify();
}

Bool_t Ks0Alg::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ks0Alg::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ks0Alg::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Ks0Alg_cxx
