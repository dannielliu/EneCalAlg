//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jun 21 09:32:17 2015 by ROOT version 5.34/28
// from TTree gepep_lambdac/ks N-Tuple example
// found on file: xyz_lambdac_4575_v2.root
//////////////////////////////////////////////////////////

#ifndef gepep_lambdac_h
#define gepep_lambdac_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gepep_lambdac {
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
   Int_t           nkap;
   Double_t        kappx[5];   //[nkap]
   Double_t        kappy[5];   //[nkap]
   Double_t        kappz[5];   //[nkap]
   Double_t        kape[5];   //[nkap]
   Int_t           nkam;
   Double_t        kampx[5];   //[nkam]
   Double_t        kampy[5];   //[nkam]
   Double_t        kampz[5];   //[nkam]
   Double_t        kame[5];   //[nkam]
   Int_t           npip;
   Double_t        pippx[8];   //[npip]
   Double_t        pippy[8];   //[npip]
   Double_t        pippz[8];   //[npip]
   Double_t        pipe[8];   //[npip]
   Int_t           npim;
   Double_t        pimpx[7];   //[npim]
   Double_t        pimpy[7];   //[npim]
   Double_t        pimpz[7];   //[npim]
   Double_t        pime[7];   //[npim]
   Int_t           npp;
   Double_t        protonpx[6];   //[npp]
   Double_t        protonpy[6];   //[npp]
   Double_t        protonpz[6];   //[npp]
   Double_t        protone[6];   //[npp]
   Int_t           npm;
   Double_t        pbarpx[5];   //[npm]
   Double_t        pbarpy[5];   //[npm]
   Double_t        pbarpz[5];   //[npm]
   Double_t        pbare[5];   //[npm]
   Double_t        totalm;
   Double_t        lambdacm;
   Double_t        lambdacbarm;

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
   TBranch        *b_nkap;   //!
   TBranch        *b_kappx;   //!
   TBranch        *b_kappy;   //!
   TBranch        *b_kappz;   //!
   TBranch        *b_kape;   //!
   TBranch        *b_nkam;   //!
   TBranch        *b_kampx;   //!
   TBranch        *b_kampy;   //!
   TBranch        *b_kampz;   //!
   TBranch        *b_kame;   //!
   TBranch        *b_npip;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!
   TBranch        *b_npp;   //!
   TBranch        *b_protonpx;   //!
   TBranch        *b_protonpy;   //!
   TBranch        *b_protonpz;   //!
   TBranch        *b_protone;   //!
   TBranch        *b_npm;   //!
   TBranch        *b_pbarpx;   //!
   TBranch        *b_pbarpy;   //!
   TBranch        *b_pbarpz;   //!
   TBranch        *b_pbare;   //!
   TBranch        *b_totalm;   //!
   TBranch        *b_lambdacm;   //!
   TBranch        *b_lambdacbarm;   //!

   gepep_lambdac(TTree *tree=0);
   virtual ~gepep_lambdac();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gepep_lambdac_cxx
gepep_lambdac::gepep_lambdac(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("xyz_lambdac_4575_v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("xyz_lambdac_4575_v2.root");
      }
      f->GetObject("gepep_lambdac",tree);

   }
   Init(tree);
}

gepep_lambdac::~gepep_lambdac()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_lambdac::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_lambdac::LoadTree(Long64_t entry)
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

void gepep_lambdac::Init(TTree *tree)
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
   fChain->SetBranchAddress("nkap", &nkap, &b_nkap);
   fChain->SetBranchAddress("kappx", kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", kappz, &b_kappz);
   fChain->SetBranchAddress("kape", kape, &b_kape);
   fChain->SetBranchAddress("nkam", &nkam, &b_nkam);
   fChain->SetBranchAddress("kampx", kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", kampz, &b_kampz);
   fChain->SetBranchAddress("kame", kame, &b_kame);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   fChain->SetBranchAddress("npp", &npp, &b_npp);
   fChain->SetBranchAddress("protonpx", protonpx, &b_protonpx);
   fChain->SetBranchAddress("protonpy", protonpy, &b_protonpy);
   fChain->SetBranchAddress("protonpz", protonpz, &b_protonpz);
   fChain->SetBranchAddress("protone", protone, &b_protone);
   fChain->SetBranchAddress("npm", &npm, &b_npm);
   fChain->SetBranchAddress("pbarpx", pbarpx, &b_pbarpx);
   fChain->SetBranchAddress("pbarpy", pbarpy, &b_pbarpy);
   fChain->SetBranchAddress("pbarpz", pbarpz, &b_pbarpz);
   fChain->SetBranchAddress("pbare", pbare, &b_pbare);
   fChain->SetBranchAddress("totalm", &totalm, &b_totalm);
   fChain->SetBranchAddress("lambdacm", &lambdacm, &b_lambdacm);
   fChain->SetBranchAddress("lambdacbarm", &lambdacbarm, &b_lambdacbarm);
   Notify();
}

Bool_t gepep_lambdac::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_lambdac::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_lambdac::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_lambdac_cxx
