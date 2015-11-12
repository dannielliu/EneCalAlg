//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 11 11:10:22 2015 by ROOT version 5.34/28
// from TTree mctruth/ks N-Tuple example
// found on file: mcKpi.root
//////////////////////////////////////////////////////////

#ifndef mctruth_h
#define mctruth_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class mctruth {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        truthpkap[4];
   Double_t        truthpkam[4];
   Int_t           mcnpip;
   Double_t        mcpippx[10];   //[mcnpip]
   Double_t        mcpippy[10];   //[mcnpip]
   Double_t        mcpippz[10];   //[mcnpip]
   Double_t        mcpippe[10];   //[mcnpip]
   Int_t           mcnpim;
   Double_t        mcpimpx[10];   //[mcnpim]
   Double_t        mcpimpy[10];   //[mcnpim]
   Double_t        mcpimpz[10];   //[mcnpim]
   Double_t        mcpimpe[10];   //[mcnpim]
   Double_t        mcgpx[3];
   Double_t        mcgpy[3];
   Double_t        mcgpz[3];
   Double_t        mcge[3];
   Double_t        dang_min;
   Double_t        dang_e;

   // List of branches
   TBranch        *b_truthpkap;   //!
   TBranch        *b_truthpkam;   //!
   TBranch        *b_mcnpip;   //!
   TBranch        *b_mcpippx;   //!
   TBranch        *b_mcpippy;   //!
   TBranch        *b_mcpippz;   //!
   TBranch        *b_mcpippe;   //!
   TBranch        *b_mcnpim;   //!
   TBranch        *b_mcpimpx;   //!
   TBranch        *b_mcpimpy;   //!
   TBranch        *b_mcpimpz;   //!
   TBranch        *b_mcpimpe;   //!
   TBranch        *b_mcgpx;   //!
   TBranch        *b_mcgpy;   //!
   TBranch        *b_mcgpz;   //!
   TBranch        *b_mcge;   //!
   TBranch        *b_dang_min;   //!
   TBranch        *b_dang_e;   //!

   mctruth(TTree *tree=0);
   virtual ~mctruth();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mctruth_cxx
mctruth::mctruth(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcKpi.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("mcKpi.root");
      }
      f->GetObject("mctruth",tree);

   }
   Init(tree);
}

mctruth::~mctruth()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mctruth::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mctruth::LoadTree(Long64_t entry)
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

void mctruth::Init(TTree *tree)
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

   fChain->SetBranchAddress("truthpkap", truthpkap, &b_truthpkap);
   fChain->SetBranchAddress("truthpkam", truthpkam, &b_truthpkam);
   fChain->SetBranchAddress("mcnpip", &mcnpip, &b_mcnpip);
   fChain->SetBranchAddress("mcpippx", mcpippx, &b_mcpippx);
   fChain->SetBranchAddress("mcpippy", mcpippy, &b_mcpippy);
   fChain->SetBranchAddress("mcpippz", mcpippz, &b_mcpippz);
   fChain->SetBranchAddress("mcpippe", mcpippe, &b_mcpippe);
   fChain->SetBranchAddress("mcnpim", &mcnpim, &b_mcnpim);
   fChain->SetBranchAddress("mcpimpx", mcpimpx, &b_mcpimpx);
   fChain->SetBranchAddress("mcpimpy", mcpimpy, &b_mcpimpy);
   fChain->SetBranchAddress("mcpimpz", mcpimpz, &b_mcpimpz);
   fChain->SetBranchAddress("mcpimpe", mcpimpe, &b_mcpimpe);
   fChain->SetBranchAddress("mcgpx", mcgpx, &b_mcgpx);
   fChain->SetBranchAddress("mcgpy", mcgpy, &b_mcgpy);
   fChain->SetBranchAddress("mcgpz", mcgpz, &b_mcgpz);
   fChain->SetBranchAddress("mcge", mcge, &b_mcge);
   fChain->SetBranchAddress("dang_min", &dang_min, &b_dang_min);
   fChain->SetBranchAddress("dang_e", &dang_e, &b_dang_e);
   Notify();
}

Bool_t mctruth::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mctruth::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mctruth::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mctruth_cxx
