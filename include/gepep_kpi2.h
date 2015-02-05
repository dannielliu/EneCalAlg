//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov  5 12:52:41 2014 by ROOT version 5.34/19
// from TTree gepep_kpi/ks N-Tuple example
// found on file: RValue_kpi2_3850.root
//////////////////////////////////////////////////////////

#ifndef gepep_kpi2_h
#define gepep_kpi2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gepep_kpi2 {
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
   Double_t        pippx[5];   //[npip]
   Double_t        pippy[5];   //[npip]
   Double_t        pippz[5];   //[npip]
   Double_t        pipe[5];   //[npip]
   Double_t        delthe[5];   //[npip]
   Int_t           npim;
   Double_t        pimpx[5];   //[npim]
   Double_t        pimpy[5];   //[npim]
   Double_t        pimpz[5];   //[npim]
   Double_t        pime[5];   //[npim]
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
   TBranch        *b_delthe;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!
   TBranch        *b_kkm4;   //!

   gepep_kpi2(TTree *tree=0);
   virtual ~gepep_kpi2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

