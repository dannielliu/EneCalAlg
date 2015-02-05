#define gepep_3pi_cxx
//#include "gepep_3pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void gepep_3pi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_3pi.C
//      Root > gepep_3pi t
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
   
   TCanvas *c1 = new TCanvas();
   TH2D *hntrack = new TH2D("hntrack","",7,0,7,7,0,7);
   TH1D *hmass = new TH1D("hmass","invariant mass",1000,0.4,0.6);
   TH2D *h2p = new TH2D("h2p","h2p",100,0.0,2.0,100,0.0,2.0);
   TH1D *hp = new TH1D("hp","hp",100,0.0,2.0);

   double mass;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	  hntrack->Fill(npip,npim);
	  if (npip==1 && npim==1){
	    double px,py,pz;
		double p1,p2,p3;
		double e;
		px = pippx[0]+pimpx[0];
		py = pippy[0]+pimpy[0];
		pz = pippz[0]+pimpz[0];
		e  = pipe[0] +pime[0] ;
		p1 = sqrt(pippx[0]*pippx[0]+pippy[0]*pippy[0]+pippz[0]*pippz[0]);
		p2 = sqrt(pimpx[0]*pimpx[0]+pimpy[0]*pimpy[0]+pimpz[0]*pimpz[0]);
		mass = sqrt(e*e-px*px-py*py-pz*pz);
		//if(p1>1.8||p2>1.8)
		  hmass->Fill(mass);
		if (mass > 0.49 && mass < 0.505){
		  h2p->Fill(p1,p2);
		  hp->Fill(p1);
		  hp->Fill(p2);
		}
	  }
   }
   //hntrack->Draw();
   hmass->Draw();
   TCanvas *c2 = new TCanvas();
   h2p->Draw();
   TCanvas *c3 = new TCanvas();
   hp->Draw();
}

#ifdef gepep_3pi_cxx
gepep_3pi::gepep_3pi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RValue_Dto3pi_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RValue_Dto3pi_3850.root");
      }
      f->GetObject("gepep_3pi",tree);

   }
   Init(tree);
}
gepep_3pi::gepep_3pi(string file) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
 
   //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RValue_Dto3pi_3850.root");
   TFile *f = new TFile(file.c_str());
   if (!f || !f->IsOpen()) {
      f = new TFile("RValue_Dto3pi_3850.root");
   }
   TTree *tree;
   f->GetObject("gepep_3pi",tree);

   Init(tree);
}

gepep_3pi::~gepep_3pi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_3pi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_3pi::LoadTree(Long64_t entry)
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

void gepep_3pi::Init(TTree *tree)
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
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_3pi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_3pi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_3pi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_3pi_cxx
