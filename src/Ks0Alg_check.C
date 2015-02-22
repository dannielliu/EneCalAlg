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
namespace KS{
  void FitSpectrum(TTree *&dataraw, char* namesfx);
}

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

   std::cout<<"loop"<<std::endl;

   char fname[100];
   sprintf(fname,"%s/plot_ks.root",outputdir.c_str());
   TFile *f=new TFile(fname,"RECREATE");
   double mpi=0.13957;
   double kslow=0.45;
   double ksup =0.55;
   
   TTree *vars = new TTree("vars","vars");
   double mass;
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
   
   TH1D *hmass  = new TH1D("hmass" ,"Ks mass",200,0.45,0.55);
   TH1D *hmassU = new TH1D("hmassU","Ks mass",200,0.45,0.55);
   TH1D *hmassc = new TH1D("hmassc","Ks mass",200,0.45,0.55);
   //TH2D *h2pos = new TH2D("h2pos","pos",10, 0,5, 10,0,5);
   TH2D *h2p = new TH2D("h2p","P",200, 0,2, 200,0,2);

   const int Npart=1;
   double start=0.1;
   double stop =1.0;
   double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
//  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range
   pcut[0] = start;
   pcut[1] = stop;
 //pcut[0] =0.0 ;  
 //pcut[1] =0.05;  
 //pcut[2] =0.10;  
 //pcut[3] =0.15;  
 //pcut[4] =0.20;  
 //pcut[5] =0.25;  
 //pcut[6] =0.30;  
 //pcut[7] =0.35;  
 //pcut[8] =0.40;  
 //pcut[9] =0.45;  
 //pcut[10]=0.50;  
 //pcut[11]=0.60;  
 //pcut[12]=0.70;  
 //pcut[13]=0.80;  
 //pcut[14]=0.90;  
 //pcut[15]=1.00;  
 //pcut[16]=1.20;  
 //pcut[17]=1.40;  
 //pcut[18]=1.60;       
 //pcut[19]=1.80;       
 //pcut[20]=2.00;  //facmap[20]=2.00;

   Event evt;
   std::vector<Event> evts[Npart][Npart];
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
      int ipip=0,ipim=0;
      for (int ipip=0; ipip<npip;ipip++)
      for (int ipim=0; ipim<npim;ipim++){
        mass = CalInvMass(mpi,pippx[ipip],pippy[ipip],pippz[ipip],mpi,pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        double chi2 = m_chisq0[ipip*npim+ipim]+m_chisqvtxnd[ipip*npim+ipim];
        if (chi2 > 100) continue;
        double ratio = m_decayL[ipip*npim+ipim]/m_decayLerr[ipip*npim+ipim];
        if (ratio<2) continue;
        if (Ks_id[ipip*npim+ipim]==0){
          hmassU->Fill(mass);
          continue;
        }
        hmass->Fill(Mpippim[ipip*npim+ipim]);
        hmassc->Fill(mass);
        p1=CalMom(pippx[ipip],pippy[ipip],pippz[ipip]);
        p2=CalMom(pimpx[ipim],pimpy[ipim],pimpz[ipim]);
        if (p1<0.15 || p1>1.0) continue;
        if (p2<0.15 || p2>1.0) continue;
        costheta1 = pippz[ipip]/p1;
        costheta2 = pimpz[ipim]/p2;
        if (pippy[ipip]>0)  phi1 = acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
        else  phi1 = 2*TMath::Pi()-acos(pippx[ipip]/CalMom(pippx[ipip],pippy[ipip]));
        if (pimpy[ipim]>0)  phi2 = acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
        else  phi2 = 2*TMath::Pi()-acos(pimpx[ipim]/CalMom(pimpx[ipim],pimpy[ipim]));
        int parti=-1,partj=-1;
		for (int i=0;i<Npart;i++){
          if (p1>=pcut[i]&&p1<pcut[i+1]){
            parti=i;
            break;
          }
        }
        for (int i=0;i<Npart;i++){
          if (p2>=pcut[i]&&p2<pcut[i+1]){
            partj = i;
            break;
          }
        }
        if(parti==-1 || partj==-1) continue;
        if (mass>kslow && mass<ksup){
          vars->Fill();
		  evt.SetVal(pippx[ipip],pippy[ipip],pippz[ipip],
		             pimpx[ipim],pimpy[ipim],pimpz[ipim]);
		  evts[parti][partj].push_back(evt);
        }
        // allocate factor according to p, average factor

      }
      
   }// loop data end
   vars->Write();

   for (int parti=0; parti<Npart;parti++)
   for (int partj=0; partj<Npart;partj++){
     sprintf(fname,"part%02d%02d",parti,partj);
	 if(evts[parti][partj].size()<200) continue;
	 FitSpe(evts[parti][partj],fname);
   }

   return;
}

void Ks0Alg::FitSpe(std::vector<Event> &evts, const char* namesfx)
 {
   double mpi=0.13957;
   double kslow=0.45;
   double ksup =0.55;
   //TF1 *facfit = new TF1("facfit",line2,kslow,ksup,2);
   TTree *datarawo= new TTree("datarawo" ,"dataraw");
   TTree *dataraw = new TTree("dataraw" ,"dataraw");
   TTree *datarawl= new TTree("datarawl" ,"dataraw");
   TTree *datarawu= new TTree("datarawu" ,"dataraw");
   double mass;
   datarawo->Branch("x",&mass,"x/D");
   dataraw->Branch("x",&mass,"x/D");
   datarawl->Branch("x",&mass,"x/D");
   datarawu->Branch("x",&mass,"x/D");
   double p1,p2;


   const int Npart=20;
   double start=0.1;
   double stop =1.0;
   double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   //pcut[0]=0.1;
   //pcut[1]=0.4;
   //pcut[2]=0.9;
   double facmap[Npart]; // for pi+
   double facemap[Npart];
   double facmapm[Npart]; // for pi-
   double facemapm[Npart];

//  set normal factor in (0.1, 0.4) to 1.00082, get factor in different range
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00263 ;  facemap[2] =0.000106524;
   pcut[3] =0.15;  facmap[3] =1.00126 ;  facemap[3] =5.65503e-05;
   pcut[4] =0.20;  facmap[4] =1.00087 ;  facemap[4] =5.36283e-05;
   pcut[5] =0.25;  facmap[5] =1.00060 ;  facemap[5] =4.96854e-05;
   pcut[6] =0.30;  facmap[6] =1.00063 ;  facemap[6] =5.4565e-05 ;
   pcut[7] =0.35;  facmap[7] =1.00034 ;  facemap[7] =6.70578e-05;
   pcut[8] =0.40;  facmap[8] =0.999913;  facemap[8] =8.50691e-05;
   pcut[9] =0.45;  facmap[9] =0.999925;  facemap[9] =0.000103853;
   pcut[10]=0.50;  facmap[10]=0.999614;  facemap[10]=0.000101297;
   pcut[11]=0.60;  facmap[11]=0.999918;  facemap[11]=9.75854e-05;
   pcut[12]=0.70;  facmap[12]=1.000260;  facemap[12]=0.000228099;
   pcut[13]=0.80;  facmap[13]=1.000340;  facemap[13]=0.000356495;
   pcut[14]=0.90;  facmap[14]=1.00056 ;  facemap[14]=0.000554875;
   pcut[15]=1.00;  facmap[15]=1.00185 ;  facemap[15]=0.000650392;
   pcut[16]=1.20;  facmap[16]=0.994443;  facemap[16]=0.00288715 ;
   pcut[17]=1.40;  facmap[17]=0.994232;  facemap[17]=0.00301238 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
//  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00263 ;  facemap[2] =0.000104711;
   pcut[3] =0.15;  facmap[3] =1.00127 ;  facemap[3] =5.66523e-05;
   pcut[4] =0.20;  facmap[4] =1.00089 ;  facemap[4] =5.40327e-05;
   pcut[5] =0.25;  facmap[5] =1.00061 ;  facemap[5] =4.96881e-05;
   pcut[6] =0.30;  facmap[6] =1.00063 ;  facemap[6] =5.45699e-05;
   pcut[7] =0.35;  facmap[7] =1.00034 ;  facemap[7] =6.70568e-05;
   pcut[8] =0.40;  facmap[8] =0.999916;  facemap[8] =8.51944e-05;
   pcut[9] =0.45;  facmap[9] =0.999916;  facemap[9] =0.00010567 ;
   pcut[10]=0.50;  facmap[10]=0.99955 ;  facemap[10]=8.6275e-05 ;
   pcut[11]=0.60;  facmap[11]=0.999922;  facemap[11]=0.000112065;
   pcut[12]=0.70;  facmap[12]=1.00024 ;  facemap[12]=0.000229724;
   pcut[13]=0.80;  facmap[13]=1.00037 ;  facemap[13]=0.000353896;
   pcut[14]=0.90;  facmap[14]=1.00057 ;  facemap[14]=0.000551636;
   pcut[15]=1.00;  facmap[15]=1.00208 ;  facemap[15]=0.000718673;
   pcut[16]=1.20;  facmap[16]=0.994435;  facemap[16]=0.00286605 ;
   pcut[17]=1.40;  facmap[17]=0.98542 ;  facemap[17]=0.00410906 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
//  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range for pi-
/* 
 //pcut[0] =0.0 ;//facmapm[0] =1.0;       facemapm[0] =1.0;     
 //pcut[1] =0.05;//facmapm[1] =1.0;       facemapm[1] =1.0;     
 //pcut[2] =0.10;//facmapm[2] =1.0;       facemapm[2] =1.0;
 //pcut[3] =0.15;//facmapm[3] =1.00196 ;  facemapm[3] =6.16454e-05;
 //pcut[4] =0.20;//facmapm[4] =1.00103 ;  facemapm[4] =6.03674e-05;
 //pcut[5] =0.25;//facmapm[5] =1.00037 ;  facemapm[5] =7.83072e-05;
 //pcut[6] =0.30;//facmapm[6] =1.00029 ;  facemapm[6] =0.000100928;
 //pcut[7] =0.35;//facmapm[7] =1.00006 ;  facemapm[7] =0.000131676;
 //pcut[8] =0.40;//facmapm[8] =0.999687;  facemapm[8] =0.000174975;
 //pcut[9] =0.45;//facmapm[9] =0.999197;  facemapm[9] =0.000210539;
 //pcut[10]=0.50;//facmapm[10]=0.999238;  facemapm[10]=0.000194209;
 //pcut[11]=0.60;//facmapm[11]=0.999154;  facemapm[11]=0.00028029 ;
 //pcut[12]=0.70;//facmapm[12]=0.999249;  facemapm[12]=0.00044126 ;
 //pcut[13]=0.80;//facmapm[13]=0.999974;  facemapm[13]=0.000672389;
 //pcut[14]=0.90;//facmapm[14]=1.00274 ;  facemapm[14]=0.00108137 ;
 //pcut[15]=1.00;//facmapm[15]=0.999595;  facemapm[15]=0.00119913 ;
 //pcut[16]=1.20;//facmapm[16]=1.00089 ;  facemapm[16]=0.00206829 ;
 //pcut[17]=1.40;//facmapm[17]=0.988005;  facemapm[17]=0.00665777 ;
 //pcut[18]=1.60;//facmapm[18]=1.0;       facemapm[18]=1.0;     
 //pcut[19]=1.80;//facmapm[19]=1.0;       facemapm[19]=1.0;     
 //pcut[20]=2.00;////facmap[20]=2.00;

                 //facmapp[0] =1.0;       facemapp[0] =1.0;     
                 //facmapp[1] =1.0;       facemapp[1] =1.0;     
                 //facmapp[2] =1.0;       facemapp[2] =1.0;
                 //facmapp[3] =1.00187 ;  facemapp[3] =6.54054e-05;
                 //facmapp[4] =1.00108 ;  facemapp[4] =6.1776e-05 ;
                 //facmapp[5] =1.00033 ;  facemapp[5] =8.0688e-05 ;
                 //facmapp[6] =1.00034 ;  facemapp[6] =8.32413e-05;
                 //facmapp[7] =0.999868;  facemapp[7] =0.000105336;
                 //facmapp[8] =0.999483;  facemapp[8] =0.000156073;
                 //facmapp[9] =0.999932;  facemapp[9] =0.000191548;
                 //facmapp[10]=0.999581;  facemapp[10]=0.000181117;
                 //facmapp[11]=1.00036 ;  facemapp[11]=0.000271994;
                 //facmapp[12]=1.00085 ;  facemapp[12]=0.000409303;
                 //facmapp[13]=1.00128 ;  facemapp[13]=0.000622945;
                 //facmapp[14]=1.00197 ;  facemapp[14]=0.000878579;
                 //facmapp[15]=1.00456 ;  facemapp[15]=0.00127855 ;
                 //facmapp[16]=0.999691;  facemapp[16]=0.000841315;
                 //facmapp[17]=1.01137 ;  facemapp[17]=0.0145154  ;
                 //facmapp[18]=1.0;       facemapp[18]=1.0;     
                 //facmapp[19]=1.0;       facemapp[19]=1.0;     
*/

//  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
//  only for r value data, limit pi+ split pi-
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00283 ;  facemap[2] =0.000250076;
   pcut[3] =0.15;  facmap[3] =1.00132 ;  facemap[3] =0.000131411;
   pcut[4] =0.20;  facmap[4] =1.00073 ;  facemap[4] =0.000129736;
   pcut[5] =0.25;  facmap[5] =1.00052 ;  facemap[5] =0.000116304;
   pcut[6] =0.30;  facmap[6] =1.00051 ;  facemap[6] =0.00012692 ;
   pcut[7] =0.35;  facmap[7] =1.00015 ;  facemap[7] =0.000155205;
   pcut[8] =0.40;  facmap[8] =0.999884;  facemap[8] =0.000200586;
   pcut[9] =0.45;  facmap[9] =0.999119;  facemap[9] =0.000245204;
   pcut[10]=0.50;  facmap[10]=0.99943 ;  facemap[10]=0.000135434;
   pcut[11]=0.60;  facmap[11]=0.999845;  facemap[11]=0.000352511;
   pcut[12]=0.70;  facmap[12]=0.999938;  facemap[12]=0.000551319;
   pcut[13]=0.80;  facmap[13]=0.999888;  facemap[13]=0.000923838;
   pcut[14]=0.90;  facmap[14]=1.00093 ;  facemap[14]=0.00142603 ;
   pcut[15]=1.00;  facmap[15]=0.999717;  facemap[15]=0.00181855 ;
   pcut[16]=1.20;  facmap[16]=0.993496;  facemap[16]=0.00501287 ;
   pcut[17]=1.40;  facmap[17]=0.981287;  facemap[17]=0.0131427  ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
//  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
//  only for r value data, limit pi- split pi+
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00172 ;  facemap[2] =0.000197672;
   pcut[3] =0.15;  facmap[3] =1.00101 ;  facemap[3] =0.000120255;
   pcut[4] =0.20;  facmap[4] =1.00081 ;  facemap[4] =0.000122614;
   pcut[5] =0.25;  facmap[5] =1.00048 ;  facemap[5] =0.000122483;
   pcut[6] =0.30;  facmap[6] =1.00074 ;  facemap[6] =0.000145777;
   pcut[7] =0.35;  facmap[7] =1.0003  ;  facemap[7] =0.000182639;
   pcut[8] =0.40;  facmap[8] =0.999589;  facemap[8] =0.000226131;
   pcut[9] =0.45;  facmap[9] =0.999426;  facemap[9] =0.000266261;
   pcut[10]=0.50;  facmap[10]=0.998955;  facemap[10]=0.000260378;
   pcut[11]=0.60;  facmap[11]=0.999262;  facemap[11]=0.000221771;
   pcut[12]=0.70;  facmap[12]=0.999362;  facemap[12]=0.000590433;
   pcut[13]=0.80;  facmap[13]=0.999897;  facemap[13]=0.000976447;
   pcut[14]=0.90;  facmap[14]=1.00053 ;  facemap[14]=0.00145681 ;
   pcut[15]=1.00;  facmap[15]=0.998737;  facemap[15]=0.00191489 ;
   pcut[16]=1.20;  facmap[16]=1.01641 ;  facemap[16]=0.0126093  ;
   pcut[17]=1.40;  facmap[17]=1.00226 ;  facemap[17]=0.00331346 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/   
   //  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
//  only for r value data, combine both part
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00215 ;  facemap[2] =0.000155076;
   pcut[3] =0.15;  facmap[3] =1.00115 ;  facemap[3] =8.87154e-05;
   pcut[4] =0.20;  facmap[4] =1.00077 ;  facemap[4] =8.91127e-05;
   pcut[5] =0.25;  facmap[5] =1.0005  ;  facemap[5] =8.43392e-05;
   pcut[6] =0.30;  facmap[6] =1.00061 ;  facemap[6] =9.57233e-05;
   pcut[7] =0.35;  facmap[7] =1.00021 ;  facemap[7] =0.000118269;
   pcut[8] =0.40;  facmap[8] =0.999754;  facemap[8] =0.000150058;
   pcut[9] =0.45;  facmap[9] =0.99926 ;  facemap[9] =0.000180371;
   pcut[10]=0.50;  facmap[10]=0.999329;  facemap[10]=0.000120152;
   pcut[11]=0.60;  facmap[11]=0.999427;  facemap[11]=0.000187713;
   pcut[12]=0.70;  facmap[12]=0.99967 ;  facemap[12]=0.00040296 ;
   pcut[13]=0.80;  facmap[13]=0.999892;  facemap[13]=0.00067108 ;
   pcut[14]=0.90;  facmap[14]=1.00073 ;  facemap[14]=0.00101906 ;
   pcut[15]=1.00;  facmap[15]=0.999252;  facemap[15]=0.00131865 ;
   pcut[16]=1.20;  facmap[16]=0.996623;  facemap[16]=0.00465825 ;
   pcut[17]=1.40;  facmap[17]=1.00101 ;  facemap[17]=0.00321292 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
   //  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range
//  only for r value data, combine both part

   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.0;       facemap[2] =1.0;
   pcut[3] =0.15;  facmap[3] =1.00193 ;  facemap[3] =0.00014646 ;
   pcut[4] =0.20;  facmap[4] =1.00092 ;  facemap[4] =0.00014185 ;
   pcut[5] =0.25;  facmap[5] =1.00027 ;  facemap[5] =0.00017408 ;
   pcut[6] =0.30;  facmap[6] =1.00028 ;  facemap[6] =0.00025041 ;
   pcut[7] =0.35;  facmap[7] =0.999665;  facemap[7] =0.00028309 ;
   pcut[8] =0.40;  facmap[8] =0.999421;  facemap[8] =0.00040939 ;
   pcut[9] =0.45;  facmap[9] =0.998626;  facemap[9] =0.00048584 ;
   pcut[10]=0.50;  facmap[10]=0.999499;  facemap[10]=0.00043914 ;
   pcut[11]=0.60;  facmap[11]=0.998844;  facemap[11]=0.00070548 ;
   pcut[12]=0.70;  facmap[12]=0.998429;  facemap[12]=0.00139006 ;
   pcut[13]=0.80;  facmap[13]=1.00149 ;  facemap[13]=0.00189529 ;
   pcut[14]=0.90;  facmap[14]=1.00292 ;  facemap[14]=0.00278649 ;
   pcut[15]=1.00;  facmap[15]=1.00281 ;  facemap[15]=0.00383812 ;
   pcut[16]=1.20;  facmap[16]=1.02283 ;  facemap[16]=0.0085055  ;
   pcut[17]=1.40;  facmap[17]=1.0;       facemap[17]=1.0;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
 
   char name[200];
   double factori=1;
   double factorj=1;
   // ~~~~~~~~ draw nxn histogram, m distribution in different range
 
   // loop data
   for (Long64_t jentry=0; jentry<evts.size();jentry++) {
      double px1=evts.at(jentry).px1;
      double py1=evts.at(jentry).py1;
      double pz1=evts.at(jentry).pz1;
      double px2=evts.at(jentry).px2;
      double py2=evts.at(jentry).py2;
      double pz2=evts.at(jentry).pz2;
     
      mass = CalInvMass(mpi,px1,py1,pz1,mpi,px2,py2,pz2);
      p1=CalMom(px1,py1,pz1);
      p2=CalMom(px2,py2,pz2);
      if (p1<0.1  || p1>1.0) continue;
      if (p2<0.1  || p2>1.0) continue;
      
      if (mass>kslow && mass<ksup){
        datarawo->Fill();
      }
	  
	  // allocate factor according to p, average factor
      for (int i=0;i<Npart;i++){
        if (p1>=pcut[i]&&p1<pcut[i+1]){
          factori = facmap[i];
          break;
        }
      }
      for (int i=0;i<Npart;i++){
        if (p2>=pcut[i]&&p2<pcut[i+1]){
          factorj = facmap[i];
          break;
        }
      }
      mass = CalInvMass(mpi,factori*px1,factori*py1,factori*pz1,
	  	               mpi,factorj*px2,factorj*py2,factorj*pz2);
      if (mass>kslow && mass<ksup){
        dataraw->Fill();
      }
 
      // allocate factor according to p, factor at low edge
      for (int i=0;i<Npart;i++){
        if (p1>=pcut[i]&&p1<pcut[i+1]){
          factori = facmap[i]-facemap[i];
          break;
        }
      }
      for (int i=0;i<Npart;i++){
        if (p2>=pcut[i]&&p2<pcut[i+1]){
          factorj = facmap[i]-facemap[i];
          break;
        }
      }
      mass = CalInvMass(mpi,factori*px1,factori*py1,factori*pz1,
	  	               mpi,factorj*px2,factorj*py2,factorj*pz2);
      if (mass>kslow && mass<ksup){
        datarawl->Fill();
      }

      // allocate factor according to p, factor at up edge
      for (int i=0;i<Npart;i++){
        if (p1>=pcut[i]&&p1<pcut[i+1]){
          factori = facmap[i]+facemap[i];
          break;
        }
      }
      for (int i=0;i<Npart;i++){
        if (p2>=pcut[i]&&p2<pcut[i+1]){
          factorj = facmap[i]+facemap[i];
          break;
        }
      }
      mass = CalInvMass(mpi,factori*px1,factori*py1,factori*pz1,
	  	               mpi,factorj*px2,factorj*py2,factorj*pz2);
      if (mass>kslow && mass<ksup){
        datarawu->Fill();
      }
      
   }// loop data end
   std::cout<<"fffffffffffffffff"<<std::endl;
   
   // no correction
   sprintf(name,"raw_%s",namesfx);
   KS::FitSpectrum(datarawo,name);

   // average factor part
   sprintf(name,"%s",namesfx);
   KS::FitSpectrum(dataraw,name);

   // low edge part
   sprintf(name,"low_%s",namesfx);
   KS::FitSpectrum(datarawl,name);
   // low part end

   //  up part 
   sprintf(name,"up_%s",namesfx);
   KS::FitSpectrum(datarawu,name);
   // up part end
   return;
}

void KS::FitSpectrum(TTree *&dataraw, char* namesfx)
{  
   TCanvas *c1 = new TCanvas();
   int Npar;
   char name[200];
   int nBins=100;
   double mparticle=0.497614;
   double kslow=0.45;
   double ksup =0.55;

   // roofit variables and functions
   RooRealVar x("x","energy",mparticle,kslow,ksup,"GeV");
   RooRealVar mean("mean","mean of gaussian",mparticle,kslow,ksup);
   RooRealVar sigma("sigma","width of gaussian",0.0030,0.0010,0.0050);
   RooRealVar sigma2("sigma2","width of gaussian",0.009,0.005,0.02);
   RooRealVar brewid("brewid","width of breit wigner",0.0023,0.0010,0.05);
   mean.setError(0.000014);
   sigma.setError(0.000032);
   sigma2.setError(0.000078);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
 
   RooRealVar a0("a0","coefficient #0",1000,-1000000,1000000);
   RooRealVar a1("a1","coefficient #1",-1000,-1000000,1000000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
 
   RooRealVar signal("signal"," ",100000,0,10000000);//event number
   RooRealVar signal2("signal2"," ",100000,0,10000000);//event number
   RooRealVar background("background"," ",300000,0,10000000);
   signal.setError(sqrt(signal.getVal()));
   signal2.setError(sqrt(signal2.getVal()));
   background.setError(sqrt(background.getVal()));

   RooAddPdf *sum;
   RooDataHist *datahist;
   RooDataSet *dataset;
   RooPlot *xframe;  
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   
   if (dataraw->GetEntries("x>0.49 && x<0.51")<50) return;
   xframe = x.frame(Title("fit Ks"));
   // fit data
   sprintf(name,"data_Ks_%s",namesfx);
   dataset = new RooDataSet(name,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
   Npar=8;
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

   pt->Draw();
   c1->Update();
   sprintf(name,"final_fit_%s",namesfx);
   c1->SetName(name);
   c1->Write();

   delete sum;
   delete dataset;
   delete xframe;
   delete pt;
   delete c1;
   //average part end

   return;
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
