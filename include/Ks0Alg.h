//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb  6 09:56:41 2015 by ROOT version 5.24/00b
// from TTree Ks_info/track N-Tuple example
// found on file: Ksto2pi_583_140226_0035862.root
//////////////////////////////////////////////////////////

#ifndef Ks0Alg_h
#define Ks0Alg_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

struct Event;

class Ks0Alg {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        rec_truth_Mks;
   Double_t        rec_truth_Pks;
   Double_t        rec_truth_Eks;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Int_t           npip;
   Int_t           npim;
   Int_t           nKsCan;
   Double_t        m_chisq0[100];   //[nKsCan]
   Double_t        m_chisqvtxnd[100];   //[nKsCan]
   Double_t        m_decayL[100];   //[nKsCan]
   Double_t        m_decayLerr[100];   //[nKsCan]
   Double_t        m_ctau[100];   //[nKsCan]
   Int_t           Ks_id[100];   //[nKsCan]
   Double_t        Mpippim[100];   //[nKsCan]
   Double_t        Ppippim[100];   //[nKsCan]
   Double_t        Epippim[100];   //[nKsCan]
   Double_t        Thepippim[100];   //[nKsCan]
   Double_t        Phipippim[100];   //[nKsCan]
   Double_t        Ppip[100];   //[nKsCan]
   Double_t        Ppim[100];   //[nKsCan]
   Double_t        Thepip[100];   //[nKsCan]
   Double_t        Thepim[100];   //[nKsCan]
   Double_t        Phipip[100];   //[nKsCan]
   Double_t        Phipim[100];   //[nKsCan]
   Double_t        pippx[10];   //[npip]
   Double_t        pippy[10];   //[npip]
   Double_t        pippz[10];   //[npip]
   Double_t        pipe[10];   //[npip]
   Double_t        pimpx[10];   //[npim]
   Double_t        pimpy[10];   //[npim]
   Double_t        pimpz[10];   //[npim]
   Double_t        pime[10];   //[npim]

   // List of branches
   TBranch        *b_rec_truth_Mks;   //!
   TBranch        *b_rec_truth_Pks;   //!
   TBranch        *b_rec_truth_Eks;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_npip;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_nKsCan;   //!
   TBranch        *b_m_chisq0;   //!
   TBranch        *b_m_chisqvtxnd;   //!
   TBranch        *b_m_decayL;   //!
   TBranch        *b_m_decayLerr;   //!
   TBranch        *b_m_ctau;   //!
   TBranch        *b_Ks_id;   //!
   TBranch        *b_Mpippim;   //!
   TBranch        *b_Ppippim;   //!
   TBranch        *b_Epippim;   //!
   TBranch        *b_Thepippim;   //!
   TBranch        *b_Phipippim;   //!
   TBranch        *b_Ppip;   //!
   TBranch        *b_Ppim;   //!
   TBranch        *b_Thepip;   //!
   TBranch        *b_Thepim;   //!
   TBranch        *b_Phipip;   //!
   TBranch        *b_Phipim;   //!
   TBranch        *b_pippx;   //!
   TBranch        *b_pippy;   //!
   TBranch        *b_pippz;   //!
   TBranch        *b_pipe;   //!
   TBranch        *b_pimpx;   //!
   TBranch        *b_pimpy;   //!
   TBranch        *b_pimpz;   //!
   TBranch        *b_pime;   //!

   Ks0Alg(TTree *tree=0);
   virtual ~Ks0Alg();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   void FitSpe(std::vector<Event> &evts, const char* namesfx);
};

#endif

