#define gepep_4k_cxx
#include "gepep_4k.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

void gepep_4k::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_4k.C
//      Root > gepep_4k t
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

   TH1D *hmass = new TH1D("hmass","mass",200,2.0,5.0);
   double mass;
   TFile tf("kkkk.root","RECREATE");
   TTree tree("tree","tree");
   tree.Branch("mass",&mass,"mass/D");

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
	  mass=CalMass();
	  hmass->Fill(mass);
	  tree.Fill();
   }
   tf.Write();
   TCanvas *c1 = new TCanvas();
   hmass->Draw();
   c1->Print("kkkkmass.pdf");
}

#ifdef gepep_4k_cxx
gepep_4k::gepep_4k(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RValue_4k_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RValue_4k_3850.root");
      }
      f->GetObject("gepep_4k",tree);

   }
   Init(tree);
}

gepep_4k::~gepep_4k()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_4k::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_4k::LoadTree(Long64_t entry)
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

void gepep_4k::Init(TTree *tree)
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
   fChain->SetBranchAddress("kappx", kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", kappz, &b_kappz);
   fChain->SetBranchAddress("kape", kape, &b_kape);
   fChain->SetBranchAddress("kampx", kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", kampz, &b_kampz);
   fChain->SetBranchAddress("kame", kame, &b_kame);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_4k::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_4k::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_4k::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
double gepep_4k::CalMass()
{
   double mass;
   double px,py,pz,p,e;
   double mk=0.493;
   px=kappx[0]+kappx[1]+kampx[0]+kampx[1];
   py=kappy[0]+kappy[1]+kampy[0]+kampy[1];
   pz=kappz[0]+kappz[1]+kampz[0]+kampz[1];
   p=TMath::Sqrt(px*px+py*py+pz*pz);
   e=kape[0]+kape[1]+kame[0]+kame[1];
   mass=TMath::Sqrt(e*e-p*p);
   return mass;
}
#endif // #ifdef gepep_4k_cxx
