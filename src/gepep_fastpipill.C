#define gepep_fastpipill_cxx
#include "gepep_fastpipill.h"
#include <TH1D.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
//#include <TStyle.h>
#include <TF1.h>
#include "function.h"
#include "bes3plotstyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"
//#include <stdlib.h>
#include <fstream>
//#include <string>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
extern std::string outputdir;
using RooFit::Title;
using RooFit::Components;
using RooFit::LineStyle;
using RooFit::LineColor;
using RooFit::Range;

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
   if (fChain == 0) return false;
   Long64_t nentries = fChain->GetEntriesFast();

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=12;
   double NBratio=((double)nentries)/nBins;
   double psilow=3.0;
   double psiup=3.2;
   double psiplow=3.6;
   double psipup=3.8;
   //TF1 *brewig = new TF1("BreWig",BreitWigner,3.05,3.15,3);
   //######## for psi(2S)
   TF1 *fit = new TF1("fit",fitfun,3.0,3.2,5);
   fit->SetLineColor(kMagenta);
   fit->SetParNames("Const","E0","sigma","intercept","slope");
   // set limits for better fit of jpsi -> mu- + mu+
   fit->SetParLimits(0,0.05,1000);
   fit->SetParLimits(2,0.03,0.04);
   fit->SetParLimits(4,-5.0,5.0);
   TF1 *sig = new TF1("signal",BreitWigner,3.6,3.8,3);
   sig->SetLineColor(kBlue);
   TF1 *bck = new TF1("background",line,3.6,3.8,2);
   bck->SetLineColor(kRed);
   // for factor fit
   TF1 *facfit = new TF1("facfit",line2,3.0,3.2,2);

   // try to use roofit
   RooRealVar x("x","energy",3.0,3.2,"GeV");
   RooRealVar mean("mean","mean of gaussian",3.0,3.8);
   RooRealVar sigma("sigma","width of gaussian",0.03,0.04);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooRealVar co1("co1","coefficient #1",0);
   RooRealVar co2("co2","coefficient #2",0);
   RooRealVar co3("co3","coefficient #3",0);
   RooRealVar co4("co4","coefficient #4",0);
   RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",3.1,3.0,3.8);
   RooRealVar background("background"," ",0);
   RooAddPdf sum("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));

   double par1[6];// for e
   double parerr1[6];
   double par2[6];// for mu
   double parerr2[6];
   double par3[6];// for pi
   double parerr3[6];
   ofstream ofpar;
   ofpar.open("par.txt",std::ios::app);
   ofstream ofpardetail;
   ofpardetail.open("detail.txt",std::ios::app);

   TCanvas *c1=new TCanvas("","",800,600);
   TH1D *h1   = new TH1D("h1","2 electron invariant mass",nBins,3,3.2);
   TH1D *h1_1 = new TH1D("h1_1","2 electron invariant mass",nBins,3,3.2);
   TH1D *h2   = new TH1D("h2","2 muon invariant mass",nBins,3,3.2);
   TH1D *h2_1 = new TH1D("h2_1","2 muon invariant mass",nBins,3,3.2);
   TH1D *h3   = new TH1D("h3","total invariant mass",nBins,3.6,3.8);
   TH1D *h3_1 = new TH1D("h3_1","total invariant mass",nBins,3.6,3.8);

   // Roofit part
   RooDataHist data_e("data_e","data_e",x,h1);
   RooDataHist data_mu("data_mu","data_mu",x,h2);
   //RooDataHist data_pi("data_pi","data_pi",x,h3);
   RooPlot *xframe=x.frame(Title("xxxxxxx"));

   //Long64_t nentries = fChain->GetEntriesFast();
   //std::cout<<"Toral entry is "<<nentries<<std::endl;
   // for initial spectrum
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
	  if(cos(angle4)>0.90) continue; // cut bhabha 
	  if(decay_ee==1)
	    h1_1->Fill(llm4);
	  else
	    h2_1->Fill(llm4);
	  h3_1->Fill(psipm4);
	  // 2 pi invariant mass
      // if (Cut(ientry) < 0) continue;
   }  

   // try to correct the spectrum
   // iniialize the fit function
   double me=0.000511;
   double mmu=0.105658;
   double mpi=0.13957;
   double factor1,factor1err;// for e+,e-
   double factor2,factor2err;// for mu+,mu-
   double factor3,factor3err;// for pi

   // psi 3.097
   //fit->SetParameters(1.08e-3*NBratio,3.097,0.04,-6.0e-2*NBratio-150,2.2e-2*NBratio+42.4);
   //fit->SetParameters(1.08e-3*NBratio,3.097,0.04,0,0);
   fit->SetParLimits(1,3.0,3.2);
   const int pointNo=20;
   double factor=0.99;
   double factorstep=(1.-factor)*2/pointNo;
   double peakvalue=3.09690;
   double deltapeak=0.002;
   double fittimes=0;
   //const int pointNo=100;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   std::string fitepsname;
   std::string fiteps_start;
   std::string fiteps_stop;
   std::string tmpstr;
   TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

//~~~~~~~~~electron part~~~~~~~~~~
   // for saving the fit result
   fitepsname  =outputdir+"/fitee.eps";
   fiteps_start=fitepsname+"[";
   fiteps_stop =fitepsname+"]";
   c1->Print(fiteps_start.c_str());
   //
   for (int i=0;i<pointNo;i++){
      h1->Reset();
	  //h1->SetAxisRange(psilow+3.1*(factor-1.0),psiup+3.1*(factor-1.0));
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
		 //if(ngam>0) continue;
	     double mass;
	     double totpx,totpy,totpz,tote;
		 double lee[2];
	     // total invariant mass
	     totpx=factor*(lepx4[0]+lepx4[1]);
	     totpy=factor*(lepy4[0]+lepy4[1]);
	     totpz=factor*(lepz4[0]+lepz4[1]);
		 if(cos(angle4)>0.90) continue; // cut bhabha 
		 if(decay_ee==0){
		   continue;
		 }
		 else if(decay_ee==1){
		   lee[0]=TMath::Sqrt(me*me + 
		          factor*factor*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
		   lee[1]=TMath::Sqrt(me*me +
		          factor*factor*(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
		   //lee[0]=TMath::Sqrt(mmu*mmu+
		   //       factor*factor*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
		 }
		 else{
		   std::cout<<"can not identify lepton. "<<std::endl;
		   return false;
		 }

	     tote=lee[0]+lee[1];
	     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	     h1->Fill(mass);
         // if (Cut(ientry) < 0) continue;
      }

      /*  
      fit->SetParameters(2.0e-2*NBratio,3.097+3.1*(factor-1.0),0.04,0,0);
	  h1->Fit(fit,"Q","",psilow,psiup);
      fit->GetParameters(par1);
	  double *tmperr;
      tmperr = fit->GetParErrors();
	  parerr1[1]=tmperr[1];
	  deltapeak = par1[1] - peakvalue;
	  factors[i]=factor;
	  factorserr[i]=0.;//factorstep*0.34;
	  deltapeaks[i]=deltapeak;
	  deltapeakserr[i]=parerr1[1];
	  ofpardetail<<factor<<"\t"<<factorserr[i]<<std::endl;
	  ofpardetail<<deltapeak<<"\t"<<parerr1[1]<<std::endl;
	  
      // save the histogram
      c1->Clear();
      h1->Draw("E");// or "scatter"
      fit->SetParameters(par1);
      fit->SetRange(psilow,psiup);
      //fit->GetParErrors(parerr);
      fit->Draw("same");
      sig->SetParameters(par1);
      sig->SetRange(psilow,psiup);
      bck->SetParameters(&par1[3]);
      bck->SetRange(psilow,psiup);
      sig->Draw("same");
      bck->Draw("same");
      legend->Clear();
      legend->AddEntry(h1,"data");
      legend->AddEntry(fit,"signal+backgorund");
      legend->AddEntry(sig,"signal");
      legend->AddEntry(bck,"backgorund");
      legend->Draw();
	  */
      sum.fitTo(data_e,Range(psilow,psiup));
	  data_e.plotOn(xframe);
	  xframe->Draw();
	  sum.plotOn(xframe);
	  sum.plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	  sum.plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
      c1->Print(fitepsname.c_str());
	  xframe->Clear();
	  fittimes++;
	  factor += factorstep;
   }
   std::cout<<"entry is "<<nentries<<std::endl;
   c1->Print(fiteps_stop.c_str());
   c1->Clear();
   /* 
   TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph1->SetTitle("delta peak");
   graph1->Draw("AP");
   graph1->SetMarkerStyle(5);
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   factor1=facfit->GetParameter(0);
   factor1err=facfit->GetParError(0);
   ofpar<<factor1<<"\t"<<factor1err<<std::endl;
   std::cout<<"fit factor: "<<factor1<<", error is "<<factor1err<<std::endl;
   std::string tmpstr=outputdir+"/factore.eps";
   c1->Print(tmpstr.c_str());*/
//~~~~~~~~~~~~~electron part end~~~~~~~~~~

//~~~~~~~~~~~~~muon part start~~~~~~~~~
   factor = 0.99;
   fittimes =0;
   deltapeak=0.002;
   // for saving the fit result
   fitepsname  =outputdir+"/fitmu.eps";
   fiteps_start=fitepsname+"[";
   fiteps_stop =fitepsname+"]";
   c1->Print(fiteps_start.c_str());
   //
   for (int i=0;i<pointNo;i++){
      h2->Reset();
      std::cout<<"factor is "<<factor<<std::endl;

      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
		 //if(ngam>0) continue;
	     double mass;
	     double totpx,totpy,totpz,tote;
		 double lee[2];
	     totpx=factor*(lepx4[0]+lepx4[1]);
	     totpy=factor*(lepy4[0]+lepy4[1]);
	     totpz=factor*(lepz4[0]+lepz4[1]);
		 if(cos(angle4)>0.90) continue; // cut bhabha 
		 if(decay_ee==0){
		   lee[0]=TMath::Sqrt(mmu*mmu + 
		          factor*factor*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
		   lee[1]=TMath::Sqrt(mmu*mmu+
		          factor*factor*(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
		 }
		 else if(decay_ee==1){
		   continue;
		 }
		 else{
		   std::cout<<"can not identify lepton. "<<std::endl;
		   return false;
		 }

	     //tote=lee4[0]+lee4[1]+pie4[0]+pie4[1];
	     tote=lee[0]+lee[1];
	     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	     h2->Fill(mass);
         // if (Cut(ientry) < 0) continue;
      }
      fit->SetParameters(2.0e-2*NBratio,3.097+3.1*(factor-1.0),0.04,0,0);
	  h2->Fit(fit,"Q","",psilow,psiup);
      fit->GetParameters(par2);
	  double *tmperr;
      tmperr = fit->GetParErrors();
	  parerr2[1]=tmperr[1];
	  deltapeak = par2[1] - peakvalue;
	  factors[i]=factor;
	  factorserr[i]=0.;//factorstep*0.34;
	  deltapeaks[i]=deltapeak;
	  deltapeakserr[i]=parerr2[1];
	  ofpardetail<<factor<<"\t"<<factorserr[i]<<std::endl;
	  ofpardetail<<deltapeak<<"\t"<<parerr2[1]<<std::endl;
	  
	  fittimes ++;
	  factor += factorstep;
      // save the histogram
      c1->Clear();
      h2->Draw("E");// or "scatter"
      fit->SetParameters(par2);
      fit->SetRange(psilow,psiup);
      fit->Draw("same");
      sig->SetParameters(par2);
      sig->SetRange(psilow,psiup);
      bck->SetParameters(&par2[3]);
      bck->SetRange(psilow,psiup);
      sig->Draw("same");
      bck->Draw("same");
      legend->Clear();
      legend->AddEntry(h2,"data");
      legend->AddEntry(fit,"signal+backgorund");
      legend->AddEntry(sig,"signal");
      legend->AddEntry(bck,"backgorund");
      legend->Draw();
      c1->Print(fitepsname.c_str());

      //std::cout<<"factor is "<<factor-1<<std::endl;
   }
   std::cout<<"entry is "<<nentries<<std::endl;
   c1->Print(fiteps_stop.c_str());
   c1->Clear();
   
   TGraphErrors *graph2 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph2->SetTitle("delta peak");
   graph2->Draw("AP");
   graph2->SetMarkerStyle(5);
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph2->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   factor2=facfit->GetParameter(0);
   factor2err=facfit->GetParError(0);
   ofpar<<factor2<<"\t"<<factor2err<<std::endl;
   std::cout<<"fit factor: "<<factor2<<", error is "<<factor2err<<std::endl;
   tmpstr=outputdir+"/factormu.eps";
   c1->Print(tmpstr.c_str());
//~~~~~~~~~~~muon part end~~~~~~~~~

//~~~~~~~~~~~pion part start~~~~~~~~
   // iniialize the fit function
   // fit->SetParameters(2.08e-3*NBratio,3.68,0.04,-NBratio*0.19-26,-7.4e-2*NBratio+9.5);
   //fit->SetParameters(5.0e-3*NBratio,3.68+3.0*(factor-1),0.04,0,0);
   fit->SetParLimits(1,3.6,3.8);
   // psi(2S) 3.686
   factor=0.99;
   //factorstep=(1.0-factor)*2/pointNo;
   peakvalue=3.686109;
   //peakerror=0.0001;
   deltapeak=0.002;
   fittimes=0;

   // for saving the fit result
   fitepsname=outputdir+"/fitpi.eps";
   fiteps_start=fitepsname+"[";
   fiteps_stop =fitepsname+"]";
   c1->Print(fiteps_start.c_str());
   //
   for (int i=0;i<pointNo;i++){
      h3->Reset();
	  std::cout<<"factor is "<<factor<<std::endl;

      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;

         //if(ngam>0) continue;
	     double mass;
	     double totpx,totpy,totpz,tote;
		 double lee[2],pie[2];
	     // total invariant mass
	     totpx=factor1*(lepx4[0]+lepx4[1])+factor*(pipx4[0]+pipx4[1]);
	     totpy=factor1*(lepy4[0]+lepy4[1])+factor*(pipy4[0]+pipy4[1]);
	     totpz=factor1*(lepz4[0]+lepz4[1])+factor*(pipz4[0]+pipz4[1]);
		 if(cos(angle4)>0.90) continue; // cut bhabha 
		 if (decay_ee == 0){
		   lee[0]=TMath::Sqrt(mmu*mmu + 
		          factor2*factor2*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
		   lee[1]=TMath::Sqrt(mmu*mmu +
		          factor2*factor2*(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
		 }
		 else if(decay_ee==1){
     	   lee[0]=TMath::Sqrt(me*me + 
		          factor1*factor1*(lepx4[0]*lepx4[0]+lepy4[0]*lepy4[0]+lepz4[0]*lepz4[0]) );
		   lee[1]=TMath::Sqrt(me*me +
		          factor1*factor1*(lepx4[1]*lepx4[1]+lepy4[1]*lepy4[1]+lepz4[1]*lepz4[1]) );
		 }
		 else {
		   std::cout<<"can not identify lepton. "<<std::endl;
		   return false;
		 }
		 pie[0]=TMath::Sqrt(mpi*mpi+
		        factor*factor*(pipx4[0]*pipx4[0]+pipy4[0]*pipy4[0]+pipz4[0]*pipz4[0]) );
		 pie[1]=TMath::Sqrt(mpi*mpi+
		        factor*factor*(pipx4[1]*pipx4[1]+pipy4[1]*pipy4[1]+pipz4[1]*pipz4[1]) );

	     //tote=lee4[0]+lee4[1]+pie4[0]+pie4[1];
	     tote=lee[0]+lee[1]+pie[0]+pie[1];
	     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	     h3->Fill(mass);

         // if (Cut(ientry) < 0) continue;
      } 
      fit->SetParameters(2.0e-2*NBratio,3.68+0.4*(factor-1),0.04,0,0);
	  h3->Fit(fit,"Q","",psiplow,psipup);
      fit->GetParameters(par3);
      double *tmperr;
      tmperr = fit->GetParErrors();
	  parerr3[1]=tmperr[1];
	  deltapeak = par3[1] - peakvalue;
	  factors[i]=factor;
	  factorserr[i]=0.;//factorstep*0.34;
	  deltapeaks[i]=deltapeak;
	  deltapeakserr[i]=parerr3[1];
	  
	  fittimes ++;
	  //factor =1.0;
	  factor += factorstep;
	  
	  // save the histogram
      c1->Clear();
      h3->Draw("E");// or "scatter"
      fit->SetParameters(par3);
      fit->SetRange(psiplow,psipup);
      //fit->GetParErrors(parerr);
      fit->Draw("same");
      sig->SetParameters(par3);
      sig->SetRange(psiplow,psipup);
      bck->SetParameters(&par3[3]);
      bck->SetRange(psiplow,psipup);
      sig->Draw("same");
      bck->Draw("same");
      legend->Clear();
      legend->AddEntry(h1,"data");
      legend->AddEntry(fit,"signal+backgorund");
      legend->AddEntry(sig,"signal");
      legend->AddEntry(bck,"backgorund");
      legend->Draw();
      c1->Print(fitepsname.c_str());

   }
   c1->Print(fiteps_stop.c_str());

   c1->Clear();
   TGraphErrors *graph3 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph3->SetTitle("delta peak");
   graph3->Draw("AP");
   graph3->SetMarkerStyle(5);
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph3->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   factor3=facfit->GetParameter(0);
   factor3err=facfit->GetParError(0);
   ofpar<<factor3<<"\t"<<factor3err<<std::endl;
   std::cout<<"fit factor: "<<factor3<<", error is "<<factor3err<<std::endl;
   tmpstr=outputdir+"/factorpi.eps";
   c1->Print(tmpstr.c_str());
//~~~~~~~~~~~pion part end~~~~~~~~


   // save the histogram
   std::string epsname=outputdir+"/fpipill.eps";
   std::string eps_start=epsname+"[";
   std::string eps_stop =epsname+"]";
   c1->Print(eps_start.c_str());
   //gStyle->SetOptFit(1111);
   gStyle->SetOptStat(0);
   //SetStyle();

   
   c1->Clear();
   h1->Draw("E");// or "scatter"
   fit->SetParameters(par1);
   fit->SetRange(psilow,psiup);
   //fit->GetParErrors(parerr);
   fit->Draw("same");
   sig->SetParameters(par1);
   sig->SetRange(3.0,3.2);
   bck->SetParameters(&par1[3]);
   bck->SetRange(3.0,3.2);
   sig->Draw("same");
   bck->Draw("same");
   legend->Clear();
   legend->AddEntry(h1,"data");
   legend->AddEntry(fit,"signal+backgorund");
   legend->AddEntry(sig,"signal");
   legend->AddEntry(bck,"backgorund");
   legend->Draw();
   c1->Print(epsname.c_str());
   c1->Clear();
   h1_1->Draw("E");
   //std::cout<<"before fitting"<<std::endl;
   fit->SetParameters(1.08e-3*NBratio,3.097,0.04,0,0);
   fit->SetParLimits(1,3.0,3.2);
   h1_1->Fit(fit,"","",3.0,3.2);
   c1->Print(epsname.c_str());
   
   c1->Clear();
   h2->Draw("E");
   fit->SetParameters(par2);
   fit->SetRange(3.0,3.2);
   fit->Draw("same");
   sig->SetParameters(par2);
   sig->SetRange(3.02,3.18);
   bck->SetParameters(&par2[3]);
   bck->SetRange(3.0,3.2);
   sig->Draw("same");
   bck->Draw("same");
   c1->Print(epsname.c_str());
   c1->Clear();
   h2_1->Draw("E");
   fit->SetParameters(1.08e-3*NBratio,3.097,0.04,0,0);
   fit->SetParLimits(1,3.0,3.2);
   h2_1->Fit(fit,"","",3.0,3.2);
   c1->Print(epsname.c_str());
   
   c1->Clear();
   h3->Draw("E");
   fit->SetParameters(par3);
   fit->SetRange(3.6,3.8);
   fit->Draw("same");
   sig->SetParameters(par3);
   sig->SetRange(3.6,3.8);
   bck->SetParameters(&par3[3]);
   bck->SetRange(3.6,3.8);
   sig->Draw("same");
   bck->Draw("same");
   c1->Print(epsname.c_str());
   c1->Clear();
   h3_1->Draw("E");
   fit->SetParameters(1.08e-3*NBratio,3.687,0.04,0,0);
   fit->SetParLimits(1,3.6,3.8);
   h3_1->Fit(fit,"","",3.6,3.8);
   c1->Print(epsname.c_str());

   ofpar.close();
   ofpardetail.close();
   c1->Print(eps_stop.c_str());
   delete h1;
   delete h1_1;
   delete h2;
   delete h2_1;
   delete legend;
   //delete h3;
   delete c1;
   delete fit;
   delete sig;
   delete bck;
   return true;
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
