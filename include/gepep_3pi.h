//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  2 16:31:50 2015 by ROOT version 5.34/24
// from TTree gepep_3pi/ks N-Tuple example
// found on file: RValue_Dto3pi_3850.root
//////////////////////////////////////////////////////////

#ifndef gepep_3pi_h
#define gepep_3pi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gepep_3pi {
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
   Double_t        pippx[3];
   Double_t        pippy[3];
   Double_t        pippz[3];
   Double_t        pipe[3];
   Double_t        pimpx[3];
   Double_t        pimpy[3];
   Double_t        pimpz[3];
   Double_t        pime[3];
   Int_t           npip;
   Int_t           npim;
   Double_t        kkm4;

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
   TBranch        *b_npip;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_kkm4;   //!

   gepep_3pi(TTree *tree=0);
   gepep_3pi(string file);
   virtual ~gepep_3pi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

