//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May  3 23:08:51 2015 by ROOT version 5.34/28
// from TTree gepep_pipipp/ks N-Tuple example
// found on file: RValue_pipipp_3850.root
//////////////////////////////////////////////////////////

#ifndef gepep_pipipp_h
#define gepep_pipipp_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gepep_pipipp {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           rec;
   Int_t           evttag;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Int_t           ngch;
   Int_t           ncharg;
   Int_t           nneu;
   Double_t        pippx;
   Double_t        pippy;
   Double_t        pippz;
   Double_t        pipe;
   Double_t        pimpx;
   Double_t        pimpy;
   Double_t        pimpz;
   Double_t        pime;
   Double_t        protonpx;
   Double_t        protonpy;
   Double_t        protonpz;
   Double_t        protone;
   Double_t        pbarpx;
   Double_t        pbarpy;
   Double_t        pbarpz;
   Double_t        pbare;
   Double_t        totalm;
   Double_t        pipippm;
   Double_t        pipippbarm;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_rec;   //!
   TBranch        *b_evttag;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_ngch;   //!
   TBranch        *b_ncharg;   //!
   TBranch        *b_nneu;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!
   TBranch        *b_protonpx;   //!
   TBranch        *b_protonpy;   //!
   TBranch        *b_protonpz;   //!
   TBranch        *b_protone;   //!
   TBranch        *b_pbarpx;   //!
   TBranch        *b_pbarpy;   //!
   TBranch        *b_pbarpz;   //!
   TBranch        *b_pbare;   //!
   TBranch        *b_totalm;   //!
   TBranch        *b_pipippm;   //!
   TBranch        *b_pipippbarm;   //!

   gepep_pipipp(TTree *tree=0);
   virtual ~gepep_pipipp();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gepep_pipipp_cxx
gepep_pipipp::gepep_pipipp(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RValue_pipipp_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RValue_pipipp_3850.root");
      }
      f->GetObject("gepep_pipipp",tree);

   }
   Init(tree);
}

gepep_pipipp::~gepep_pipipp()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_pipipp::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_pipipp::LoadTree(Long64_t entry)
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

void gepep_pipipp::Init(TTree *tree)
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
   fChain->SetBranchAddress("pippx", &pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", &pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", &pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", &pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", &pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", &pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", &pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", &pime, &b_pime);
   fChain->SetBranchAddress("protonpx", &protonpx, &b_protonpx);
   fChain->SetBranchAddress("protonpy", &protonpy, &b_protonpy);
   fChain->SetBranchAddress("protonpz", &protonpz, &b_protonpz);
   fChain->SetBranchAddress("protone", &protone, &b_protone);
   fChain->SetBranchAddress("pbarpx", &pbarpx, &b_pbarpx);
   fChain->SetBranchAddress("pbarpy", &pbarpy, &b_pbarpy);
   fChain->SetBranchAddress("pbarpz", &pbarpz, &b_pbarpz);
   fChain->SetBranchAddress("pbare", &pbare, &b_pbare);
   fChain->SetBranchAddress("totalm", &totalm, &b_totalm);
   fChain->SetBranchAddress("pipippm", &pipippm, &b_pipippm);
   fChain->SetBranchAddress("pipippbarm", &pipippbarm, &b_pipippbarm);
   Notify();
}

Bool_t gepep_pipipp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_pipipp::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_pipipp::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_pipipp_cxx
