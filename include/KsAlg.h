//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  4 14:44:48 2015 by ROOT version 5.34/24
// from TTree gepep_2pi/ks N-Tuple example
// found on file: multipi_001_131210_0034011.root
//////////////////////////////////////////////////////////

#ifndef KsAlg_h
#define KsAlg_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class KsAlg {
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
   Int_t           npip;
   Int_t           npim;
   Int_t           nKsCan;
   Double_t        pippx[5];   //[npip]
   Double_t        pippy[5];   //[npip]
   Double_t        pippz[5];   //[npip]
   Double_t        pipe[5];   //[npip]
   Double_t        pimpx[5];   //[npim]
   Double_t        pimpy[5];   //[npim]
   Double_t        pimpz[5];   //[npim]
   Double_t        pime[5];   //[npim]
   Double_t        Ks_mass[25];   //[nKsCan]
   Double_t        Ks_ratio[25];   //[nKsCan]
   Double_t        Ks_decayL[25];   //[nKsCan]
   Double_t        Ks_chi2[25];   //[nKsCan]
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
   TBranch        *b_npip;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_nKsCan;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!
   TBranch        *b_Ks_mass;   //!
   TBranch        *b_Ks_ratio;   //!
   TBranch        *b_Ks_decayL;   //!
   TBranch        *b_Ks_chi2;   //!
   TBranch        *b_kkm4;   //!

   KsAlg(TTree *tree=0);
   KsAlg(std::string filename);
   virtual ~KsAlg();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(double cutd=2);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

