#define gepep_fast6pi_cxx
#include "gepep_fast6pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "function.h"
#include <fstream>

void gepep_fast6pi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fast6pi.C
//      Root > gepep_fast6pi t
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
   //TH1D *h = new TH1D("h","invariant mass",100,1.,4.);
   TH1D *h = new TH1D("h","invariant mass",100,3.,4.);
   double factore,factoreerr;
   double factorpi,factorpierr;
   double me=0.000511;
   double mmu=0.105658;
   double mpi=0.13957;
   
   ifstream f;
   f.open("par.txt");
   f >> factore >> factoreerr;
   f >> factorpi>> factorpierr;
   f.close();
   //factorpi =1.0;
   std::cout<<"factor is "<<factorpi<<std::endl;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;	  
	  double totpx,totpy,totpz,tote;
	  double pie[6];
	  double mass;
	  totpx = pippx[0]+pippx[1]+pippx[2]+pimpx[0]+pimpx[1]+pimpx[2];
	  totpy = pippy[0]+pippy[1]+pippy[2]+pimpy[0]+pimpy[1]+pimpy[2];
	  totpz = pippz[0]+pippz[1]+pippz[2]+pimpz[0]+pimpz[1]+pimpz[2];
	  //tote  = pipe[0] +pipe[1] +pipe[2] +pime[0] +pime[1] +pime[2];
	  pie[0]=TMath::Sqrt(mpi*mpi+
	         factorpi*factorpi*(pippx[0]*pippx[0]+pippy[0]*pippy[0]+pippz[0]*pippz[0]) );
	  pie[1]=TMath::Sqrt(mpi*mpi+
	         factorpi*factorpi*(pippx[1]*pippx[1]+pippy[1]*pippy[1]+pippz[1]*pippz[1]) );	
	  pie[2]=TMath::Sqrt(mpi*mpi+
	         factorpi*factorpi*(pippx[2]*pippx[2]+pippy[2]*pippy[2]+pippz[2]*pippz[2]) );
	  pie[3]=TMath::Sqrt(mpi*mpi+
	         factorpi*factorpi*(pimpx[0]*pimpx[0]+pimpy[0]*pimpy[0]+pimpz[0]*pimpz[0]) );
	  pie[4]=TMath::Sqrt(mpi*mpi+
	         factorpi*factorpi*(pimpx[1]*pimpx[1]+pimpy[1]*pimpy[1]+pimpz[1]*pimpz[1]) );
	  pie[5]=TMath::Sqrt(mpi*mpi+
	         factorpi*factorpi*(pimpx[2]*pimpx[2]+pimpy[2]*pimpy[2]+pimpz[2]*pimpz[2]) );
	  tote  = pie[0] +pie[1] +pie[2] +pie[3] +pie[4] +pie[5];
  mass = TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	  h->Fill(mass);
   }

   TCanvas *c1=new TCanvas("","",800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);
   TF1 *fit=new TF1("fit",GausLineBack,3.75,4.0,5);
   fit->SetParameters(1,3.8,0.02,10,1);
   h->Fit(fit,"","",3.75,4.0);
   h->Draw();
   //h->Fit("gaus","","",3.75,4.0);
   c1->Print("f6pi.eps");
}

#ifdef gepep_fast6pi_cxx
gepep_fast6pi::gepep_fast6pi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_Rvalue_f6pi_e3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_Rvalue_f6pi_e3850.root");
      }
      f->GetObject("gepep_fast6pi",tree);

   }
   Init(tree);
}

gepep_fast6pi::~gepep_fast6pi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fast6pi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fast6pi::LoadTree(Long64_t entry)
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

void gepep_fast6pi::Init(TTree *tree)
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
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_fast6pi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fast6pi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fast6pi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif 
