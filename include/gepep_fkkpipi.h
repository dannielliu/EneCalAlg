//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov  3 16:35:11 2014 by ROOT version 5.34/19
// from TTree gepep_fastkkpipi/ks N-Tuple example
// found on file: data/RValue_fkkpipi_3850.root
//////////////////////////////////////////////////////////

#ifndef gepep_fkkpipi_h
#define gepep_fkkpipi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gepep_fkkpipi {
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
   Double_t        kappx;
   Double_t        kappy;
   Double_t        kappz;
   Double_t        kape;
   Double_t        kampx;
   Double_t        kampy;
   Double_t        kampz;
   Double_t        kame;
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
   TBranch        *b_kappx;   //!
   TBranch        *b_kappy;   //!
   TBranch        *b_kappz;   //!
   TBranch        *b_kape;   //!
   TBranch        *b_kampx;   //!
   TBranch        *b_kampy;   //!
   TBranch        *b_kampz;   //!
   TBranch        *b_kame;   //!
   TBranch        *b_kkm4;   //!

   gepep_fkkpipi(TTree *tree=0);
   virtual ~gepep_fkkpipi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

